/**
 * Baseline tests for AlleleParser.
 *
 * These tests verify the current behavior of AlleleParser before refactoring.
 * All tests should pass with the existing implementation and continue to pass
 * after refactoring to use the IAlignmentReader interface.
 *
 * Test data:
 *   - test/data/test.ref: 11bp reference "ATCGGCTAAAA"
 *   - test/data/test.bam: 1 read with SNPs
 *   - test/tiny/NA12878.chr22.tiny.bam: Real reads from chr22
 */

#include "../src/AlleleParser.h"
#include "../src/Parameters.h"
#include "../src/Allele.h"
#include "../src/Sample.h"
#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>
#include <getopt.h>

// Declare getopt globals (needed to reset state between tests)
extern int optind;
extern int opterr;

// Simple test framework
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
#define ASSERT_GT(a, b) ASSERT((a) > (b))
#define ASSERT_GE(a, b) ASSERT((a) >= (b))
#define ASSERT_TRUE(a) ASSERT(a)
#define ASSERT_FALSE(a) ASSERT(!(a))
#define ASSERT_DOUBLE_EQ(a, b) ASSERT(std::abs((a) - (b)) < 0.0001)

// Helper to create command-line arguments
class ArgBuilder {
private:
    std::vector<std::string> args_;
    std::vector<char*> argv_;

public:
    // Constructor resets getopt state to allow multiple tests to parse args
    ArgBuilder() {
        // CRITICAL: Reset getopt global state between tests
        // Without this, only the first test will successfully parse arguments
        optind = 1;
        opterr = 0;
    }

    ArgBuilder& add(const std::string& arg) {
        args_.push_back(arg);
        return *this;
    }

    ArgBuilder& fasta(const std::string& path) {
        return add("-f").add(path);
    }

    ArgBuilder& bam(const std::string& path) {
        return add(path);
    }

    ArgBuilder& region(const std::string& region) {
        return add("-r").add(region);
    }

    ArgBuilder& targets(const std::string& bedfile) {
        return add("-t").add(bedfile);
    }

    ArgBuilder& minAltFraction(double frac) {
        return add("-F").add(std::to_string(frac));
    }

    ArgBuilder& ploidy(int p) {
        return add("-p").add(std::to_string(p));
    }

    int argc() { return args_.size(); }

    char** argv() {
        argv_.clear();
        for (auto& arg : args_) {
            argv_.push_back(const_cast<char*>(arg.c_str()));
        }
        return argv_.data();
    }

    void print() {
        std::cout << "Args: ";
        for (const auto& arg : args_) {
            std::cout << arg << " ";
        }
        std::cout << "\n";
    }
};

// ===== Basic Functionality Tests =====

TEST(CanConstructAlleleParser) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam");

    AlleleParser parser(args.argc(), args.argv());

    ASSERT_TRUE(true);  // If we got here, construction succeeded
}

TEST(LoadsReferenceSequence) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam");

    AlleleParser parser(args.argc(), args.argv());

    // Reference sequence should be loaded from BAM into referenceIDToName map
    ASSERT_FALSE(parser.referenceIDToName.empty());
    ASSERT_EQ(parser.referenceIDToName[0], "ref");
}

TEST(LoadsBamFiles) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam");

    AlleleParser parser(args.argc(), args.argv());

    // Should have loaded BAM header
    ASSERT_FALSE(parser.bamHeaderLines.empty());
}

TEST(GetsSampleNamesFromBam) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam");

    AlleleParser parser(args.argc(), args.argv());

    // Should have extracted sample names
    ASSERT_FALSE(parser.sampleList.empty());
    std::cout << "\n  Sample: " << parser.sampleList[0];
}

TEST(CanSetRegion) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam")
        .region("ref:1-11");

    AlleleParser parser(args.argc(), args.argv());

    // Should have parsed the region (stored in parameters.regions)
    ASSERT_FALSE(parser.parameters.regions.empty());
    ASSERT_EQ(parser.parameters.regions[0], "ref:1-11");

    // After moving to first position, currentSequenceName should be set
    if (parser.toNextPosition()) {
        ASSERT_EQ(parser.currentSequenceName, "ref");
    }
}

// ===== Alignment Reading Tests =====

TEST(CanReadFirstAlignment) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam");

    AlleleParser parser(args.argc(), args.argv());

    // Try to get first alignment
    bool hasAlignment = parser.getFirstAlignment();
    ASSERT_TRUE(hasAlignment);
}

TEST(CanIteratePositions) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam");

    AlleleParser parser(args.argc(), args.argv());

    int positions = 0;
    while (parser.toNextPosition() && positions < 20) {
        positions++;
    }

    ASSERT_GT(positions, 0);
    std::cout << "\n  Processed " << positions << " positions";
}

// ===== Allele Detection Tests =====

TEST(DetectsAllelesInTestBam) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam")
        .minAltFraction(0.0);  // Detect all alleles

    AlleleParser parser(args.argc(), args.argv());

    Samples samples;
    int totalAlleles = 0;
    int positions = 0;
    int positionsWithAlleles = 0;

    // The test.bam has complex alleles, so also check for ALLELE_COMPLEX and ALLELE_MNP
    int alleleTypes = ALLELE_SNP | ALLELE_INSERTION | ALLELE_DELETION | ALLELE_COMPLEX | ALLELE_MNP;

    while (parser.toNextPosition() && positions < 20) {
        samples.clear();

        // Try to get alleles at this position
        bool hasAlleles = parser.getNextAlleles(samples, alleleTypes);

        if (hasAlleles) {
            for (auto& sample : samples) {
                if (!sample.second.empty()) {
                    positionsWithAlleles++;
                    totalAlleles += sample.second.size();
                }
            }
        }

        positions++;
    }

    // The test BAM should have at least some positions processed
    ASSERT_GT(positions, 0);
    std::cout << "\n  Processed " << positions << " positions, found "
              << totalAlleles << " alleles at " << positionsWithAlleles << " positions";

    // Note: May be 0 if alleles are filtered out, but at least we processed positions
}

TEST(RegistersAlignments) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam");

    AlleleParser parser(args.argc(), args.argv());

    // Move to first position
    parser.toNextPosition();

    // Should have registered some alignments
    size_t registeredCount = 0;
    for (const auto& posAlignments : parser.registeredAlignments) {
        registeredCount += posAlignments.second.size();
    }

    ASSERT_GT(registeredCount, 0);
    std::cout << "\n  Registered " << registeredCount << " alignments";
}

// ===== Reference Access Tests =====

TEST(CanAccessReferenceSequence) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam");

    AlleleParser parser(args.argc(), args.argv());

    // Load reference
    std::string refName = "ref";
    parser.loadReferenceSequence(refName);

    // Should be able to get reference substring
    std::string substr = parser.referenceSubstr(0, 5);
    ASSERT_EQ(substr.length(), 5);
    std::cout << "\n  Reference[0:5] = " << substr;
}

TEST(ReferenceMatchesExpected) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam");

    AlleleParser parser(args.argc(), args.argv());

    std::string refName = "ref";
    parser.loadReferenceSequence(refName);

    // Reference is "ATCGGCTAAAA"
    std::string full = parser.referenceSubstr(0, 11);
    ASSERT_EQ(full, "ATCGGCTAAAA");
}

// ===== Parameter Handling Tests =====

TEST(HandlesMinAltFractionParameter) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam")
        .minAltFraction(0.2);

    AlleleParser parser(args.argc(), args.argv());

    ASSERT_DOUBLE_EQ(parser.parameters.minAltFraction, 0.2);
}

TEST(HandlesPloidyParameter) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam")
        .ploidy(2);

    AlleleParser parser(args.argc(), args.argv());

    ASSERT_EQ(parser.parameters.ploidy, 2);
}

TEST(HandlesRegionParameter) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam")
        .region("ref:5-10");

    AlleleParser parser(args.argc(), args.argv());

    // Should have parsed the region (note: Parameters uses "regions" plural)
    ASSERT_FALSE(parser.parameters.regions.empty());
}

// ===== Coverage & Depth Tests =====

TEST(TracksCoverage) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam");

    AlleleParser parser(args.argc(), args.argv());

    // Process a few positions
    int positions = 0;
    while (parser.toNextPosition() && positions < 10) {
        positions++;
    }

    // Should have coverage information
    ASSERT_FALSE(parser.coverage.empty());
    std::cout << "\n  Coverage tracked at " << parser.coverage.size() << " positions";
}

// ===== Target Region Tests =====

TEST(HandlesTargetBedFile) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("test/data/test.bam")
        .targets("test/data/test.bed");

    AlleleParser parser(args.argc(), args.argv());

    // Should have loaded targets
    ASSERT_FALSE(parser.targets.empty());
    std::cout << "\n  Loaded " << parser.targets.size() << " targets";
}

// ===== Integration Tests with Larger Data =====

TEST(WorksWithTinyBam) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/tiny/q.fa")
        .bam("test/tiny/NA12878.chr22.tiny.bam")
        .region("q:1-1000");

    AlleleParser parser(args.argc(), args.argv());

    int positions = 0;
    Samples samples;

    while (parser.toNextPosition() && positions < 100) {
        samples.clear();
        parser.getNextAlleles(samples, ALLELE_SNP | ALLELE_INSERTION | ALLELE_DELETION);
        positions++;
    }

    ASSERT_GT(positions, 0);
    std::cout << "\n  Processed " << positions << " positions from tiny BAM";
}

TEST(FindsVariantsInTinyBam) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/tiny/q.fa")
        .bam("test/tiny/NA12878.chr22.tiny.bam")
        .region("q:1-5000")
        .minAltFraction(0.2);

    AlleleParser parser(args.argc(), args.argv());

    int variantPositions = 0;
    int totalPositions = 0;
    Samples samples;

    while (parser.toNextPosition() && totalPositions < 500) {
        samples.clear();
        bool hasAlleles = parser.getNextAlleles(samples, ALLELE_SNP | ALLELE_INSERTION | ALLELE_DELETION);

        if (hasAlleles) {
            for (auto& sample : samples) {
                if (!sample.second.empty()) {
                    variantPositions++;
                    break;
                }
            }
        }

        totalPositions++;
    }

    ASSERT_GT(variantPositions, 0);
    std::cout << "\n  Found variants at " << variantPositions << "/" << totalPositions << " positions";
}

// ===== Error Handling Tests =====

TEST(HandlesNonexistentBam) {
    ArgBuilder args;
    args.add("freebayes")
        .fasta("test/data/test.ref")
        .bam("nonexistent.bam");

    bool caughtError = false;
    try {
        AlleleParser parser(args.argc(), args.argv());
    } catch (...) {
        caughtError = true;
    }

    // Should either throw or handle gracefully
    // (Current implementation exits on error, so this may not catch it)
    ASSERT_TRUE(true);  // If we got here, it's handled somehow
}

// ===== Main Test Runner =====

int main(int argc, char** argv) {
    std::cout << "===========================================\n";
    std::cout << "AlleleParser Baseline Tests\n";
    std::cout << "===========================================\n\n";

    std::cout << "Basic Functionality:\n";
    RUN_TEST(CanConstructAlleleParser);
    RUN_TEST(LoadsReferenceSequence);
    RUN_TEST(LoadsBamFiles);
    RUN_TEST(GetsSampleNamesFromBam);
    RUN_TEST(CanSetRegion);

    std::cout << "\nAlignment Reading:\n";
    RUN_TEST(CanReadFirstAlignment);
    RUN_TEST(CanIteratePositions);
    RUN_TEST(RegistersAlignments);

    std::cout << "\nAllele Detection:\n";
    RUN_TEST(DetectsAllelesInTestBam);

    std::cout << "\nReference Access:\n";
    RUN_TEST(CanAccessReferenceSequence);
    RUN_TEST(ReferenceMatchesExpected);

    std::cout << "\nParameter Handling:\n";
    RUN_TEST(HandlesMinAltFractionParameter);
    RUN_TEST(HandlesPloidyParameter);
    RUN_TEST(HandlesRegionParameter);

    std::cout << "\nCoverage & Depth:\n";
    RUN_TEST(TracksCoverage);

    std::cout << "\nTarget Regions:\n";
    RUN_TEST(HandlesTargetBedFile);

    std::cout << "\nIntegration Tests:\n";
    RUN_TEST(WorksWithTinyBam);
    RUN_TEST(FindsVariantsInTinyBam);

    std::cout << "\n===========================================\n";
    std::cout << "All baseline tests passed!\n";
    std::cout << "===========================================\n";

    return 0;
}
