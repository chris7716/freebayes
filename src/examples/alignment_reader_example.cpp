/**
 * Example demonstrating the use of the IAlignmentReader interface.
 *
 * This example shows how to:
 * 1. Auto-detect alignment format
 * 2. Create a reader using the factory
 * 3. Open alignment files
 * 4. Iterate through alignments
 * 5. Access alignment properties through the interface
 *
 * Usage:
 *   alignment_reader_example <alignment.bam> [region]
 *
 * Examples:
 *   alignment_reader_example alignments.bam
 *   alignment_reader_example alignments.cram chr1:1000-2000
 */

#include "../AlignmentReaderFactory.h"
#include "../IAlignmentReader.h"
#include "../IAlignment.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace freebayes;

void printUsage(const char* progName) {
    std::cerr << "Usage: " << progName << " <alignment_file> [chr:start-end]\n";
    std::cerr << "\n";
    std::cerr << "Examples:\n";
    std::cerr << "  " << progName << " alignments.bam\n";
    std::cerr << "  " << progName << " alignments.cram chr1:1000-2000\n";
    std::cerr << "  " << progName << " alignments.sam chr2:5000-10000\n";
}

bool parseRegion(const std::string& regionStr, std::string& chr, long& start, long& end) {
    size_t colonPos = regionStr.find(':');
    if (colonPos == std::string::npos) {
        return false;
    }

    chr = regionStr.substr(0, colonPos);

    size_t dashPos = regionStr.find('-', colonPos);
    if (dashPos == std::string::npos) {
        return false;
    }

    try {
        start = std::stol(regionStr.substr(colonPos + 1, dashPos - colonPos - 1));
        end = std::stol(regionStr.substr(dashPos + 1));
    } catch (...) {
        return false;
    }

    return true;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    std::string alignmentFile = argv[1];
    std::string region = (argc > 2) ? argv[2] : "";

    // Step 1: Detect format
    std::cout << "Detecting format for: " << alignmentFile << "\n";
    AlignmentFormat format = AlignmentReaderFactory::detectFormat(alignmentFile);
    std::cout << "Detected format: " << AlignmentReaderFactory::formatToString(format) << "\n\n";

    if (format == AlignmentFormat::UNKNOWN) {
        std::cerr << "Error: Could not detect alignment format\n";
        return 1;
    }

    // Step 2: Create reader using factory
    std::unique_ptr<IAlignmentReader> reader;
    try {
        reader = AlignmentReaderFactory::create(format);
    } catch (const std::exception& e) {
        std::cerr << "Error creating reader: " << e.what() << "\n";
        return 1;
    }

    // Step 3: Open alignment file
    std::vector<std::string> files = {alignmentFile};
    if (!reader->open(files)) {
        std::cerr << "Error opening file: " << reader->getErrorString() << "\n";
        return 1;
    }

    std::cout << "Successfully opened alignment file\n";

    // Print header information
    std::cout << "\n=== Reference Sequences ===\n";
    auto refSeqs = reader->getReferenceSequences();
    for (size_t i = 0; i < refSeqs.size() && i < 10; ++i) {
        std::cout << i << ": " << refSeqs[i].name
                  << " (length: " << refSeqs[i].length << ")\n";
    }
    if (refSeqs.size() > 10) {
        std::cout << "... and " << (refSeqs.size() - 10) << " more\n";
    }

    // Step 4: Set region if specified
    if (!region.empty()) {
        std::string chr;
        long start, end;
        if (!parseRegion(region, chr, start, end)) {
            std::cerr << "Error: Invalid region format. Use chr:start-end\n";
            return 1;
        }

        int refId = reader->getReferenceId(chr);
        if (refId < 0) {
            std::cerr << "Error: Reference sequence '" << chr << "' not found\n";
            return 1;
        }

        std::cout << "\nSetting region: " << chr << ":" << start << "-" << end << "\n";
        if (!reader->setRegion(refId, start, end)) {
            std::cerr << "Error setting region: " << reader->getErrorString() << "\n";
            return 1;
        }
    }

    // Step 5: Iterate through alignments
    std::cout << "\n=== Alignments ===\n";
    std::shared_ptr<IAlignment> alignment;
    int count = 0;
    int mapped = 0;
    int unmapped = 0;

    while (reader->getNextAlignment(alignment)) {
        count++;

        // Print first 10 alignments in detail
        if (count <= 10) {
            std::cout << "\nAlignment " << count << ":\n";
            std::cout << "  Name: " << alignment->queryName() << "\n";
            std::cout << "  Position: " << alignment->position() << "\n";
            std::cout << "  End: " << alignment->endPosition() << "\n";
            std::cout << "  MapQ: " << alignment->mappingQuality() << "\n";
            std::cout << "  Mapped: " << (alignment->isMapped() ? "yes" : "no") << "\n";
            std::cout << "  Paired: " << (alignment->isPaired() ? "yes" : "no") << "\n";
            std::cout << "  Reverse: " << (alignment->isReverseStrand() ? "yes" : "no") << "\n";
            std::cout << "  Duplicate: " << (alignment->isDuplicate() ? "yes" : "no") << "\n";

            std::string rg;
            if (alignment->getTag("RG", rg)) {
                std::cout << "  Read Group: " << rg << "\n";
            }

            // Print CIGAR string
            auto cigar = alignment->getCigar();
            std::cout << "  CIGAR: ";
            for (const auto& op : cigar) {
                std::cout << op.length << op.type;
            }
            std::cout << "\n";

            // Print first 50bp of sequence
            std::string seq = alignment->sequence();
            if (seq.length() > 50) {
                seq = seq.substr(0, 50) + "...";
            }
            std::cout << "  Sequence: " << seq << "\n";
        }

        if (alignment->isMapped()) {
            mapped++;
        } else {
            unmapped++;
        }

        // Stop after 1000 for this example
        if (count >= 1000) {
            std::cout << "\n(stopping after 1000 alignments)\n";
            break;
        }
    }

    // Print summary
    std::cout << "\n=== Summary ===\n";
    std::cout << "Total alignments processed: " << count << "\n";
    std::cout << "Mapped: " << mapped << "\n";
    std::cout << "Unmapped: " << unmapped << "\n";

    reader->close();

    return 0;
}
