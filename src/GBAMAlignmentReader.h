#ifndef FREEBAYES_GBAMALIGNMENTREADER_H
#define FREEBAYES_GBAMALIGNMENTREADER_H

#include "IAlignment.h"
#include "IAlignmentReader.h"
#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include <htslib/sam.h>

// Forward declaration to avoid pulling in gbam headers everywhere.
class Reader;

namespace freebayes {

/**
 * IAlignment implementation backed by a htslib bam1_t decoded from a GBAM file.
 * The object owns its bam1_t and destroys it on destruction.
 */
class GBAMAlignment : public IAlignment {
public:
    // Takes ownership of a bam1_t allocated with bam_init1() / bam_dup1().
    explicit GBAMAlignment(bam1_t* record);
    ~GBAMAlignment() override;

    // Non-copyable (owns a raw pointer).
    GBAMAlignment(const GBAMAlignment&) = delete;
    GBAMAlignment& operator=(const GBAMAlignment&) = delete;

    std::string queryName() const override;
    std::string readGroup() const override;
    int referenceId() const override;
    long position() const override;
    long endPosition() const override;
    int mappingQuality() const override;
    std::string sequence() const override;
    std::string qualities() const override;
    int alignmentLength() const override;
    int sequenceLength() const override;

    bool isMapped() const override;
    bool isPaired() const override;
    bool isMateMapped() const override;
    bool isProperPair() const override;
    bool isReverseStrand() const override;
    bool isDuplicate() const override;

    std::vector<CigarOp> getCigar() const override;

    bool getTag(const std::string& tag, std::string& value) const override;
    bool getTag(const std::string& tag, int& value) const override;
    bool getTag(const std::string& tag, double& value) const override;

private:
    bam1_t* record_;
};

/**
 * IAlignmentReader implementation for GBAM files.
 *
 * Opens a single .gbam file via mmap and wraps the gbam::Reader.
 * Region seeking is implemented via binary search over the refID/pos columns.
 */
class GBAMAlignmentReader : public IAlignmentReader {
public:
    GBAMAlignmentReader();
    ~GBAMAlignmentReader() override;

    GBAMAlignmentReader(const GBAMAlignmentReader&) = delete;
    GBAMAlignmentReader& operator=(const GBAMAlignmentReader&) = delete;

    bool open(const std::vector<std::string>& paths) override;
    void close() override;
    bool isOpen() const override;

    bool setRegion(int refId, long start, long end) override;
    bool getNextAlignment(std::shared_ptr<IAlignment>& alignment) override;

    std::string getHeaderText() const override;
    std::vector<ReferenceSequence> getReferenceSequences() const override;
    int getReferenceId(const std::string& name) const override;

    std::string getErrorString() const override;

    bool locateIndexes() override;
    void setReferenceFile(const std::string& fastaPath) override;

private:
    // mmap state
    void*  mmap_addr_  = nullptr;
    size_t mmap_size_  = 0;

    std::unique_ptr<Reader> reader_;
    bam1_t* buf_ = nullptr;   // reusable decode buffer

    // Region state
    int64_t current_rec_    = 0;
    int64_t region_end_rec_ = 0;  // exclusive upper bound

    bool   isOpen_    = false;
    std::string lastError_;
};

} // namespace freebayes

#endif // FREEBAYES_GBAMALIGNMENTREADER_H
