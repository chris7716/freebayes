/**
 * Minimal test to reproduce the interface crash
 */

#include "../src/SeqLibAlignmentReader.h"
#include "../src/AlignmentReaderFactory.h"
#include <iostream>
#include <cassert>

using namespace freebayes;

int main(int argc, char** argv) {
    std::cout << "=== Testing Interface-based Reader ===\n";

    // Test 1: Create reader directly
    std::cout << "\n1. Creating SeqLibAlignmentReader...\n";
    auto reader = std::make_unique<SeqLibAlignmentReader>();
    std::cout << "   ✓ Created\n";

    // Test 2: Open BAM file
    std::cout << "\n2. Opening BAM file...\n";
    std::vector<std::string> files = {"test/data/test.bam"};
    if (!reader->open(files)) {
        std::cerr << "   ✗ Failed to open: " << reader->getErrorString() << "\n";
        return 1;
    }
    std::cout << "   ✓ Opened\n";

    // Test 3: Get header
    std::cout << "\n3. Getting header...\n";
    std::string header = reader->getHeaderText();
    std::cout << "   ✓ Header length: " << header.size() << " bytes\n";

    // Test 4: Get reference sequences
    std::cout << "\n4. Getting reference sequences...\n";
    auto refSeqs = reader->getReferenceSequences();
    std::cout << "   ✓ Found " << refSeqs.size() << " reference sequences\n";
    if (!refSeqs.empty()) {
        std::cout << "     First: " << refSeqs[0].name << " (" << refSeqs[0].length << " bp)\n";
    }

    // Test 5: Read first alignment
    std::cout << "\n5. Reading first alignment...\n";
    std::shared_ptr<IAlignment> alignment;
    if (!reader->getNextAlignment(alignment)) {
        std::cerr << "   ✗ Failed to get alignment\n";
        return 1;
    }
    std::cout << "   ✓ Got alignment\n";

    // Test 6: Check alignment is valid
    std::cout << "\n6. Checking alignment validity...\n";
    if (!alignment) {
        std::cerr << "   ✗ Alignment is null!\n";
        return 1;
    }
    std::cout << "   ✓ Alignment is not null\n";

    // Test 7: Access alignment properties
    std::cout << "\n7. Accessing alignment properties...\n";
    std::cout << "   Query name: " << alignment->queryName() << "\n";
    std::cout << "   Position: " << alignment->position() << "\n";
    std::cout << "   Sequence length: " << alignment->sequenceLength() << "\n";
    std::cout << "   ✓ Properties accessible\n";

    // Test 8: Cast to SeqLibAlignment
    std::cout << "\n8. Casting to SeqLibAlignment...\n";
    auto* seqLibAlign = dynamic_cast<SeqLibAlignment*>(alignment.get());
    if (!seqLibAlign) {
        std::cerr << "   ✗ Cast failed!\n";
        return 1;
    }
    std::cout << "   ✓ Cast succeeded\n";

    // Test 9: Get underlying BamRecord
    std::cout << "\n9. Getting underlying BamRecord...\n";
    const SeqLib::BamRecord& record = seqLibAlign->getRecord();
    std::cout << "   ✓ Got BamRecord reference\n";

    // Test 10: Copy BamRecord
    std::cout << "\n10. Copying BamRecord...\n";
    SeqLib::BamRecord recordCopy = record;
    std::cout << "   ✓ Copied BamRecord\n";
    std::cout << "   Copy position: " << recordCopy.Position() << "\n";

    std::cout << "\n=== All tests passed! ===\n";
    return 0;
}
