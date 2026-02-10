/**
 * Simple syntax check for the alignment interfaces.
 * This file verifies that the interfaces compile without linking.
 * Compile with: c++ -std=c++17 -fsyntax-only -I../src test_interfaces_syntax.cpp
 */

#include "../src/IAlignment.h"
#include "../src/IAlignmentReader.h"
#include <memory>
#include <vector>
#include <string>

using namespace freebayes;

// Mock implementation to verify interface is complete
class TestAlignment : public IAlignment {
public:
    std::string queryName() const override { return "test"; }
    std::string readGroup() const override { return "RG1"; }
    int referenceId() const override { return 0; }
    long position() const override { return 100; }
    long endPosition() const override { return 200; }
    int mappingQuality() const override { return 60; }
    std::string sequence() const override { return "ACGT"; }
    std::string qualities() const override { return "IIII"; }
    int alignmentLength() const override { return 100; }
    int sequenceLength() const override { return 4; }

    bool isMapped() const override { return true; }
    bool isPaired() const override { return true; }
    bool isMateMapped() const override { return true; }
    bool isProperPair() const override { return true; }
    bool isReverseStrand() const override { return false; }
    bool isDuplicate() const override { return false; }

    std::vector<CigarOp> getCigar() const override {
        return {CigarOp('M', 4)};
    }

    bool getTag(const std::string& tag, std::string& value) const override {
        if (tag == "RG") { value = "RG1"; return true; }
        return false;
    }

    bool getTag(const std::string& tag, int& value) const override {
        if (tag == "NM") { value = 0; return true; }
        return false;
    }

    bool getTag(const std::string& tag, double& value) const override {
        return false;
    }
};

// Mock reader to verify interface is complete
class TestReader : public IAlignmentReader {
private:
    bool open_ = false;

public:
    bool open(const std::vector<std::string>& paths) override {
        open_ = true;
        return true;
    }

    void close() override { open_ = false; }

    bool isOpen() const override { return open_; }

    bool setRegion(int refId, long start, long end) override {
        return true;
    }

    bool getNextAlignment(std::shared_ptr<IAlignment>& alignment) override {
        alignment = std::make_shared<TestAlignment>();
        return false; // EOF
    }

    std::string getHeaderText() const override {
        return "@HD\tVN:1.6\n";
    }

    std::vector<ReferenceSequence> getReferenceSequences() const override {
        return {ReferenceSequence("chr1", 1000000)};
    }

    int getReferenceId(const std::string& name) const override {
        if (name == "chr1") return 0;
        return -1;
    }

    std::string getErrorString() const override {
        return "";
    }

    bool locateIndexes() override { return true; }

    void setReferenceFile(const std::string& fastaPath) override {}
};

// Test function to verify polymorphism
void processAlignments(IAlignmentReader& reader) {
    std::shared_ptr<IAlignment> aln;
    while (reader.getNextAlignment(aln)) {
        auto pos = aln->position();
        auto name = aln->queryName();
        auto cigar = aln->getCigar();
    }
}

int main() {
    // Test interface instantiation
    TestReader reader;

    // Test polymorphic usage
    IAlignmentReader* readerPtr = &reader;
    readerPtr->open({"test.bam"});

    // Test alignment access
    std::shared_ptr<IAlignment> aln = std::make_shared<TestAlignment>();
    auto pos = aln->position();
    auto mapped = aln->isMapped();

    // Test function accepting interface
    processAlignments(reader);

    return 0;
}
