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

        case AlignmentFormat::GBAM:
            throw std::runtime_error("GBAM format not yet implemented");

        default:
            throw std::runtime_error("Unknown alignment format");
    }
}

AlignmentFormat AlignmentReaderFactory::detectFormat(const std::string& path) {
    // First try extension-based detection
    AlignmentFormat format = detectFormatFromExtension(path);

    if (format != AlignmentFormat::BAM) {
        return format;
    }

    // Fall back to magic byte detection if extension gave us default BAM
    AlignmentFormat magicFormat = detectFormatFromMagic(path);

    // If magic detection found something specific, use it
    if (magicFormat != AlignmentFormat::BAM ||
        path == "-" || path == "stdin") {
        return magicFormat;
    }

    return format;
}

AlignmentFormat AlignmentReaderFactory::detectFormatFromExtension(const std::string& path) {
    std::string lowerPath = toLower(path);

    // Handle special cases
    if (path == "-" || path == "stdin") {
        // Stdin could be any format, default to BAM
        return AlignmentFormat::BAM;
    }

    // Remove trailing .gz if present
    std::string checkPath = lowerPath;
    if (endsWith(checkPath, ".gz")) {
        checkPath = checkPath.substr(0, checkPath.length() - 3);
    }

    // Check extensions
    if (endsWith(checkPath, ".cram")) {
        return AlignmentFormat::CRAM;
    }
    if (endsWith(checkPath, ".sam")) {
        return AlignmentFormat::SAM;
    }
    if (endsWith(checkPath, ".gbam") || endsWith(checkPath, ".gam")) {
        return AlignmentFormat::GBAM;
    }
    if (endsWith(checkPath, ".bam")) {
        return AlignmentFormat::BAM;
    }

    // Default to BAM
    return AlignmentFormat::BAM;
}

AlignmentFormat AlignmentReaderFactory::detectFormatFromMagic(const std::string& path) {
    // Skip magic detection for stdin/special files
    if (path == "-" || path == "stdin") {
        return AlignmentFormat::BAM;
    }

    std::ifstream file(path, std::ios::binary);
    if (!file.good()) {
        return AlignmentFormat::BAM;
    }

    // Read first 4 bytes
    char magic[4];
    file.read(magic, 4);

    if (!file.good()) {
        return AlignmentFormat::BAM;
    }

    // Check magic bytes
    // CRAM: "CRAM" (0x43 0x52 0x41 0x4D)
    if (memcmp(magic, "CRAM", 4) == 0) {
        return AlignmentFormat::CRAM;
    }

    // SAM: text starting with @ header
    if (magic[0] == '@' && (
        (magic[1] == 'H' && magic[2] == 'D') ||
        (magic[1] == 'S' && magic[2] == 'Q') ||
        (magic[1] == 'R' && magic[2] == 'G') ||
        (magic[1] == 'P' && magic[2] == 'G') ||
        (magic[1] == 'C' && magic[2] == 'O')
    )) {
        return AlignmentFormat::SAM;
    }

    // BAM: gzip header (0x1f 0x8b)
    if ((unsigned char)magic[0] == 0x1f &&
        (unsigned char)magic[1] == 0x8b) {
        return AlignmentFormat::BAM;
    }

    // Default to BAM
    return AlignmentFormat::BAM;
}

std::string AlignmentReaderFactory::formatToString(AlignmentFormat format) {
    switch (format) {
        case AlignmentFormat::BAM:  return "BAM";
        case AlignmentFormat::CRAM: return "CRAM";
        case AlignmentFormat::SAM:  return "SAM";
        case AlignmentFormat::GBAM: return "GBAM";
        default:                    return "UNKNOWN";
    }
}

AlignmentFormat AlignmentReaderFactory::stringToFormat(const std::string& name) {
    std::string lowerName = toLower(name);

    if (lowerName == "cram")           return AlignmentFormat::CRAM;
    if (lowerName == "sam")            return AlignmentFormat::SAM;
    if (lowerName == "gbam" ||
        lowerName == "gam")            return AlignmentFormat::GBAM;
    if (lowerName == "bam")            return AlignmentFormat::BAM;

    // Default to BAM
    return AlignmentFormat::BAM;
}

} // namespace freebayes
