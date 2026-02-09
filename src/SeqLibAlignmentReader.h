#ifndef FREEBAYES_SEQLIBALIGNMENTREADER_H
#define FREEBAYES_SEQLIBALIGNMENTREADER_H

#include "IAlignment.h"
#include "IAlignmentReader.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include <memory>
#include <string>
#include <vector>

namespace freebayes {

/**
 * Adapter that wraps a SeqLib::BamRecord to implement IAlignment interface.
 * This allows SeqLib BAM/CRAM records to be used polymorphically.
 */
class SeqLibAlignment : public IAlignment {
public:
    /**
     * Construct from a SeqLib BamRecord.
     * The record is copied, so the original can be reused.
     */
    explicit SeqLibAlignment(const SeqLib::BamRecord& record);

    // Copy constructor and assignment
    SeqLibAlignment(const SeqLibAlignment& other) = default;
    SeqLibAlignment& operator=(const SeqLibAlignment& other) = default;

    // IAlignment interface implementation
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

    /**
     * Get the underlying SeqLib::BamRecord (for performance-critical code)
     * @return const reference to the wrapped record
     */
    const SeqLib::BamRecord& getRecord() const { return record_; }

private:
    SeqLib::BamRecord record_;
    mutable std::string cachedReadGroup_;  // Cache RG tag lookup
    mutable bool readGroupCached_;
};

/**
 * IAlignmentReader implementation using SeqLib::BamReader.
 * Supports BAM, CRAM, and SAM formats through SeqLib.
 *
 * This is a wrapper around the existing SeqLib reader used by freebayes,
 * providing a format-agnostic interface.
 */
class SeqLibAlignmentReader : public IAlignmentReader {
public:
    SeqLibAlignmentReader();
    ~SeqLibAlignmentReader() override;

    // Prevent copying (reader manages file handles)
    SeqLibAlignmentReader(const SeqLibAlignmentReader&) = delete;
    SeqLibAlignmentReader& operator=(const SeqLibAlignmentReader&) = delete;

    // Allow moving
    SeqLibAlignmentReader(SeqLibAlignmentReader&&) = default;
    SeqLibAlignmentReader& operator=(SeqLibAlignmentReader&&) = default;

    // IAlignmentReader interface implementation
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

    /**
     * Get the underlying SeqLib::BamReader (for advanced use cases)
     * @return pointer to the wrapped reader
     */
    SeqLib::BamReader* getReader() { return reader_.get(); }

private:
    std::unique_ptr<SeqLib::BamReader> reader_;
    std::string lastError_;
    bool isOpen_;

    // Reusable record for efficiency
    SeqLib::BamRecord currentRecord_;
};

} // namespace freebayes

#endif // FREEBAYES_SEQLIBALIGNMENTREADER_H
