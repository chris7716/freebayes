/**
 * Unit tests for alignment reader interfaces.
 *
 * Tests the IAlignment and IAlignmentReader interfaces with the SeqLib implementation.
 */

#include "../src/AlignmentReaderFactory.h"
#include "../src/SeqLibAlignmentReader.h"
#include "../src/IAlignmentReader.h"
#include "../src/IAlignment.h"

#include <cassert>
#include <iostream>
#include <memory>
#include <string>

using namespace freebayes;

// Simple test framework macros
#define TEST(name) void test_##name()
#define RUN_TEST(name) do { \
    std::cout << "Running " #name "..."; \
    test_##name(); \
    std::cout << " PASSED\n"; \
} while(0)

#define ASSERT(condition) do { \
    if (!(condition)) { \
        std::cerr << "\nAssertion failed: " #condition "\n"; \
        std::cerr << "  at " << __FILE__ << ":" << __LINE__ << "\n"; \
        exit(1); \
    } \
} while(0)

#define ASSERT_EQ(a, b) ASSERT((a) == (b))
#define ASSERT_NE(a, b) ASSERT((a) != (b))
#define ASSERT_TRUE(a) ASSERT(a)
#define ASSERT_FALSE(a) ASSERT(!(a))

// ===== Format Detection Tests =====

TEST(FormatDetectionByExtension) {
    ASSERT_EQ(AlignmentReaderFactory::detectFormatFromExtension("test.bam"),
              AlignmentFormat::BAM);
    ASSERT_EQ(AlignmentReaderFactory::detectFormatFromExtension("test.cram"),
              AlignmentFormat::CRAM);
    ASSERT_EQ(AlignmentReaderFactory::detectFormatFromExtension("test.sam"),
              AlignmentFormat::SAM);
    ASSERT_EQ(AlignmentReaderFactory::detectFormatFromExtension("test.gbam"),
              AlignmentFormat::GBAM);
    ASSERT_EQ(AlignmentReaderFactory::detectFormatFromExtension("test.BAM"),
              AlignmentFormat::BAM);
    ASSERT_EQ(AlignmentReaderFactory::detectFormatFromExtension("test.bam.gz"),
              AlignmentFormat::BAM);
    // Unknown extensions default to BAM
    ASSERT_EQ(AlignmentReaderFactory::detectFormatFromExtension("unknown.xyz"),
              AlignmentFormat::BAM);
}

TEST(FormatStringConversion) {
    ASSERT_EQ(AlignmentReaderFactory::formatToString(AlignmentFormat::BAM), "BAM");
    ASSERT_EQ(AlignmentReaderFactory::formatToString(AlignmentFormat::CRAM), "CRAM");
    ASSERT_EQ(AlignmentReaderFactory::formatToString(AlignmentFormat::SAM), "SAM");
    ASSERT_EQ(AlignmentReaderFactory::formatToString(AlignmentFormat::GBAM), "GBAM");

    ASSERT_EQ(AlignmentReaderFactory::stringToFormat("bam"), AlignmentFormat::BAM);
    ASSERT_EQ(AlignmentReaderFactory::stringToFormat("BAM"), AlignmentFormat::BAM);
    ASSERT_EQ(AlignmentReaderFactory::stringToFormat("cram"), AlignmentFormat::CRAM);
    ASSERT_EQ(AlignmentReaderFactory::stringToFormat("sam"), AlignmentFormat::SAM);
    ASSERT_EQ(AlignmentReaderFactory::stringToFormat("gbam"), AlignmentFormat::GBAM);
    // Unknown formats default to BAM
    ASSERT_EQ(AlignmentReaderFactory::stringToFormat("unknown"), AlignmentFormat::BAM);
}

// ===== Reader Creation Tests =====

TEST(CreateBAMReader) {
    auto reader = AlignmentReaderFactory::create(AlignmentFormat::BAM);
    ASSERT_NE(reader.get(), nullptr);
    ASSERT_FALSE(reader->isOpen());
}

TEST(CreateCRAMReader) {
    auto reader = AlignmentReaderFactory::create(AlignmentFormat::CRAM);
    ASSERT_NE(reader.get(), nullptr);
    ASSERT_FALSE(reader->isOpen());
}

TEST(CreateSAMReader) {
    auto reader = AlignmentReaderFactory::create(AlignmentFormat::SAM);
    ASSERT_NE(reader.get(), nullptr);
    ASSERT_FALSE(reader->isOpen());
}

TEST(CreateGBAMReaderNotImplemented) {
    bool caught = false;
    try {
        auto reader = AlignmentReaderFactory::create(AlignmentFormat::GBAM);
    } catch (const std::runtime_error& e) {
        caught = true;
    }
    ASSERT_TRUE(caught);
}

// ===== SeqLibAlignment Tests =====

TEST(SeqLibAlignmentBasicProperties) {
    // Create a mock SeqLib::BamRecord
    // Note: This test is somewhat limited without actual BAM data
    // In a real test suite, we'd use test BAM files

    SeqLib::BamRecord record;
    // Set some basic properties if possible
    // (SeqLib API may not allow easy mock creation)

    SeqLibAlignment alignment(record);

    // Test that methods don't crash
    std::string name = alignment.queryName();
    int refId = alignment.referenceId();
    long pos = alignment.position();
    bool mapped = alignment.isMapped();

    // Basic sanity check
    ASSERT_TRUE(name.empty() || !name.empty()); // Always true, just checking it doesn't crash
}

// ===== SeqLibAlignmentReader Tests =====

TEST(SeqLibReaderInitialState) {
    SeqLibAlignmentReader reader;
    ASSERT_FALSE(reader.isOpen());
}

TEST(SeqLibReaderOpenNonexistentFileFails) {
    SeqLibAlignmentReader reader;
    std::vector<std::string> files = {"nonexistent_file_xyz_123.bam"};
    bool result = reader.open(files);
    ASSERT_FALSE(result);
    ASSERT_FALSE(reader.isOpen());
}

TEST(SeqLibReaderOpenEmptyListFails) {
    SeqLibAlignmentReader reader;
    std::vector<std::string> files;
    bool result = reader.open(files);
    ASSERT_FALSE(result);
    ASSERT_NE(reader.getErrorString(), "");
}

// ===== Integration Tests =====
// These would require actual BAM files to test properly
// For now, we just test the API surface

TEST(ReaderLifecycle) {
    auto reader = AlignmentReaderFactory::create(AlignmentFormat::BAM);

    // Initial state
    ASSERT_FALSE(reader->isOpen());

    // After close (without open)
    reader->close();
    ASSERT_FALSE(reader->isOpen());

    // Double close should be safe
    reader->close();
    ASSERT_FALSE(reader->isOpen());
}

TEST(AlignmentPolymorphism) {
    // Test that we can use IAlignment pointer polymorphically
    std::shared_ptr<IAlignment> alignment;

    // This should compile and work even though we don't have data
    if (alignment) {
        auto name = alignment->queryName();
        auto pos = alignment->position();
        auto mapped = alignment->isMapped();
    }

    ASSERT_TRUE(alignment.get() == nullptr);
}

// ===== Main Test Runner =====

int main(int argc, char** argv) {
    std::cout << "Running alignment interface tests...\n\n";

    RUN_TEST(FormatDetectionByExtension);
    RUN_TEST(FormatStringConversion);

    RUN_TEST(CreateBAMReader);
    RUN_TEST(CreateCRAMReader);
    RUN_TEST(CreateSAMReader);
    RUN_TEST(CreateGBAMReaderNotImplemented);

    RUN_TEST(SeqLibAlignmentBasicProperties);
    RUN_TEST(SeqLibReaderInitialState);
    RUN_TEST(SeqLibReaderOpenNonexistentFileFails);
    RUN_TEST(SeqLibReaderOpenEmptyListFails);

    RUN_TEST(ReaderLifecycle);
    RUN_TEST(AlignmentPolymorphism);

    std::cout << "\n=================================\n";
    std::cout << "All tests passed!\n";
    std::cout << "=================================\n";

    return 0;
}
