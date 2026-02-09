#include "AlignmentReaderFactory.h"
#include "SeqLibAlignmentReader.h"
#include <fstream>
#include <algorithm>
#include <cctype>
#include <cstring>
#include <stdexcept>

namespace freebayes {

// Helper function to convert string to lowercase
static std::string toLower(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return result;
}

// Helper function to check if string ends with suffix
static bool endsWith(const std::string& str, const std::string& suffix) {
    if (suffix.length() > str.length()) {
        return false;
    }
    return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
}

// Helper function to check if string starts with prefix
static bool startsWith(const std::string& str, const std::string& prefix) {
    if (prefix.length() > str.length()) {
        return false;
    }
    return str.compare(0, prefix.length(), prefix) == 0;
}

std::unique_ptr<IAlignmentReader> AlignmentReaderFactory::create(
    AlignmentFormat format,
    const std::map<std::string, std::string>& options)
{
    switch (format) {
        case AlignmentFormat::BAM:
        case AlignmentFormat::CRAM:
        case AlignmentFormat::SAM: {
            auto reader = std::make_unique<SeqLibAlignmentReader>();

            // Apply options
            auto refIt = options.find("reference");
            if (refIt != options.end()) {
                reader->setReferenceFile(refIt->second);
            }

            return reader;
        }

        case AlignmentFormat::PAF:
            throw std::runtime_error("PAF format not yet implemented");

        case AlignmentFormat::AUTO_DETECT:
            throw std::runtime_error(
                "Cannot create reader with AUTO_DETECT format. "
                "Use detectFormat() first to determine actual format."
            );

        case AlignmentFormat::UNKNOWN:
        default:
            throw std::runtime_error("Unknown or unsupported alignment format");
    }
}

AlignmentFormat AlignmentReaderFactory::detectFormat(const std::string& path) {
    // First try extension-based detection
    AlignmentFormat format = detectFormatFromExtension(path);

    if (format != AlignmentFormat::UNKNOWN) {
        return format;
    }

    // Fall back to magic byte detection if extension didn't help
    format = detectFormatFromMagic(path);

    return format;
}

AlignmentFormat AlignmentReaderFactory::detectFormatFromExtension(const std::string& path) {
    std::string lowerPath = toLower(path);

    // Handle special cases
    if (path == "-" || path == "stdin") {
        // Stdin could be any format, default to BAM
        return AlignmentFormat::BAM;
    }

    // Check for compressed extensions first (.bam.gz, .sam.gz, etc.)
    // Note: .bam files are already gzip-compressed, so .bam.gz is redundant but may exist

    // Remove trailing .gz if present
    std::string checkPath = lowerPath;
    if (endsWith(checkPath, ".gz")) {
        checkPath = checkPath.substr(0, checkPath.length() - 3);
    }

    // Check extensions
    if (endsWith(checkPath, ".bam")) {
        return AlignmentFormat::BAM;
    }
    if (endsWith(checkPath, ".cram")) {
        return AlignmentFormat::CRAM;
    }
    if (endsWith(checkPath, ".sam")) {
        return AlignmentFormat::SAM;
    }
    if (endsWith(checkPath, ".paf")) {
        return AlignmentFormat::PAF;
    }

    return AlignmentFormat::UNKNOWN;
}

AlignmentFormat AlignmentReaderFactory::detectFormatFromMagic(const std::string& path) {
    // Skip magic detection for stdin/special files
    if (path == "-" || path == "stdin") {
        return AlignmentFormat::UNKNOWN;
    }

    std::ifstream file(path, std::ios::binary);
    if (!file.good()) {
        return AlignmentFormat::UNKNOWN;
    }

    // Read first 4 bytes
    char magic[4];
    file.read(magic, 4);

    if (!file.good()) {
        return AlignmentFormat::UNKNOWN;
    }

    // Check magic bytes
    // BAM: 0x1f 0x8b 0x08 (gzip header, but we need more context)
    // Actually BAM after decompression starts with "BAM\1"
    // CRAM: "CRAM" (0x43 0x52 0x41 0x4D)
    // SAM: "@HD" or "@SQ" or "@RG" (text header)
    // PAF: Text format, typically starts with read name (no fixed magic)

    // Check for CRAM
    if (memcmp(magic, "CRAM", 4) == 0) {
        return AlignmentFormat::CRAM;
    }

    // Check for SAM (text starting with @)
    if (magic[0] == '@' && (
        (magic[1] == 'H' && magic[2] == 'D') ||
        (magic[1] == 'S' && magic[2] == 'Q') ||
        (magic[1] == 'R' && magic[2] == 'G') ||
        (magic[1] == 'P' && magic[2] == 'G') ||
        (magic[1] == 'C' && magic[2] == 'O')
    )) {
        return AlignmentFormat::SAM;
    }

    // Check for gzip header (likely BAM)
    if ((unsigned char)magic[0] == 0x1f &&
        (unsigned char)magic[1] == 0x8b) {
        // This is gzip-compressed, likely BAM
        return AlignmentFormat::BAM;
    }

    // Could be PAF (text format) if it's plain text
    // PAF lines have tab-separated fields
    // Read a bit more to check for tab characters
    file.seekg(0);
    char buffer[100];
    file.read(buffer, sizeof(buffer) - 1);
    int bytesRead = file.gcount();
    buffer[bytesRead] = '\0';

    // Check if it looks like tab-delimited text
    bool hasTab = false;
    bool hasOnlyPrintable = true;
    for (int i = 0; i < bytesRead; ++i) {
        if (buffer[i] == '\t') {
            hasTab = true;
        }
        if (buffer[i] != '\t' && buffer[i] != '\n' && buffer[i] != '\r' &&
            !std::isprint(static_cast<unsigned char>(buffer[i]))) {
            hasOnlyPrintable = false;
            break;
        }
    }

    if (hasTab && hasOnlyPrintable) {
        // Likely PAF or similar text format
        return AlignmentFormat::PAF;
    }

    return AlignmentFormat::UNKNOWN;
}

std::string AlignmentReaderFactory::formatToString(AlignmentFormat format) {
    switch (format) {
        case AlignmentFormat::BAM:         return "BAM";
        case AlignmentFormat::CRAM:        return "CRAM";
        case AlignmentFormat::SAM:         return "SAM";
        case AlignmentFormat::PAF:         return "PAF";
        case AlignmentFormat::AUTO_DETECT: return "AUTO_DETECT";
        case AlignmentFormat::UNKNOWN:     return "UNKNOWN";
        default:                           return "INVALID";
    }
}

AlignmentFormat AlignmentReaderFactory::stringToFormat(const std::string& name) {
    std::string lowerName = toLower(name);

    if (lowerName == "bam")            return AlignmentFormat::BAM;
    if (lowerName == "cram")           return AlignmentFormat::CRAM;
    if (lowerName == "sam")            return AlignmentFormat::SAM;
    if (lowerName == "paf")            return AlignmentFormat::PAF;
    if (lowerName == "auto" ||
        lowerName == "auto_detect")    return AlignmentFormat::AUTO_DETECT;

    return AlignmentFormat::UNKNOWN;
}

} // namespace freebayes
