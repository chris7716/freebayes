#include "SeqLibAlignmentReader.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamRecord.h"
#include <stdexcept>

namespace freebayes {

// ===== SeqLibAlignment Implementation =====

SeqLibAlignment::SeqLibAlignment(const SeqLib::BamRecord& record)
    : record_(record)
    , cachedReadGroup_()
    , readGroupCached_(false)
{
}

std::string SeqLibAlignment::queryName() const {
    return record_.Qname();
}

std::string SeqLibAlignment::readGroup() const {
    if (!readGroupCached_) {
        record_.GetZTag("RG", cachedReadGroup_);
        readGroupCached_ = true;
    }
    return cachedReadGroup_;
}

int SeqLibAlignment::referenceId() const {
    return record_.ChrID();
}

long SeqLibAlignment::position() const {
    return record_.Position();
}

long SeqLibAlignment::endPosition() const {
    return record_.PositionEnd();
}

int SeqLibAlignment::mappingQuality() const {
    return record_.MapQuality();
}

std::string SeqLibAlignment::sequence() const {
    return record_.Sequence();
}

std::string SeqLibAlignment::qualities() const {
    return record_.Qualities();
}

int SeqLibAlignment::alignmentLength() const {
    return record_.Length();
}

int SeqLibAlignment::sequenceLength() const {
    return record_.Length();
}

bool SeqLibAlignment::isMapped() const {
    return record_.MappedFlag();
}

bool SeqLibAlignment::isPaired() const {
    return record_.PairedFlag();
}

bool SeqLibAlignment::isMateMapped() const {
    return record_.MateMappedFlag();
}

bool SeqLibAlignment::isProperPair() const {
    return record_.ProperPair();
}

bool SeqLibAlignment::isReverseStrand() const {
    return record_.ReverseFlag();
}

bool SeqLibAlignment::isDuplicate() const {
    return record_.DuplicateFlag();
}

std::vector<CigarOp> SeqLibAlignment::getCigar() const {
    std::vector<CigarOp> result;
    const SeqLib::Cigar& cigar = record_.GetCigar();

    for (const auto& field : cigar) {
        result.emplace_back(field.Type(), field.Length());
    }

    return result;
}

bool SeqLibAlignment::getTag(const std::string& tag, std::string& value) const {
    return record_.GetZTag(tag, value);
}

bool SeqLibAlignment::getTag(const std::string& tag, int& value) const {
    return record_.GetIntTag(tag, value);
}

bool SeqLibAlignment::getTag(const std::string& tag, double& value) const {
    // SeqLib doesn't have GetDoubleTag, use GetSmartDoubleTag or parse from string
    std::string strValue;
    if (record_.GetZTag(tag, strValue)) {
        try {
            value = std::stod(strValue);
            return true;
        } catch (...) {
            return false;
        }
    }
    return false;
}

// ===== SeqLibAlignmentReader Implementation =====

SeqLibAlignmentReader::SeqLibAlignmentReader()
    : reader_(std::make_unique<SeqLib::BamReader>())
    , lastError_()
    , isOpen_(false)
{
}

SeqLibAlignmentReader::~SeqLibAlignmentReader() {
    close();
}

bool SeqLibAlignmentReader::open(const std::vector<std::string>& paths) {
    if (paths.empty()) {
        lastError_ = "No input files provided";
        return false;
    }

    // Handle stdin
    if (paths.size() == 1 && paths[0] == "-") {
        if (!reader_->Open("-")) {
            lastError_ = "Failed to open stdin for reading";
            return false;
        }
        isOpen_ = true;
        return true;
    }

    // For multiple files or single file, use Open with vector
    if (!reader_->Open(paths)) {
        lastError_ = "Failed to open BAM/CRAM files";
        return false;
    }

    isOpen_ = true;
    return true;
}

void SeqLibAlignmentReader::close() {
    if (isOpen_) {
        // SeqLib::BamReader destructor handles cleanup
        // We can recreate the reader if needed
        reader_ = std::make_unique<SeqLib::BamReader>();
        isOpen_ = false;
    }
}

bool SeqLibAlignmentReader::isOpen() const {
    return isOpen_;
}

bool SeqLibAlignmentReader::setRegion(int refId, long start, long end) {
    if (!isOpen_) {
        lastError_ = "Cannot set region: reader not open";
        return false;
    }

    SeqLib::GenomicRegion region(refId, start, end);
    if (!reader_->SetRegion(region)) {
        lastError_ = "Failed to set region";
        return false;
    }

    return true;
}

bool SeqLibAlignmentReader::getNextAlignment(std::shared_ptr<IAlignment>& alignment) {
    if (!isOpen_) {
        lastError_ = "Cannot read alignment: reader not open";
        return false;
    }

    if (!reader_->GetNextRecord(currentRecord_)) {
        return false;
    }

    // Create new alignment wrapper
    alignment = std::make_shared<SeqLibAlignment>(currentRecord_);
    return true;
}

std::string SeqLibAlignmentReader::getHeaderText() const {
    if (!isOpen_) {
        return "";
    }
    return reader_->HeaderConcat();
}

std::vector<ReferenceSequence> SeqLibAlignmentReader::getReferenceSequences() const {
    std::vector<ReferenceSequence> result;

    if (!isOpen_) {
        return result;
    }

    const auto& headerSeqs = reader_->Header().GetHeaderSequenceVector();

    for (const auto& seq : headerSeqs) {
        result.emplace_back(seq.Name, seq.Length, true);
    }

    return result;
}

int SeqLibAlignmentReader::getReferenceId(const std::string& name) const {
    if (!isOpen_) {
        return -1;
    }

    int id = reader_->Header().Name2ID(name);
    return id;
}

std::string SeqLibAlignmentReader::getErrorString() const {
    return lastError_;
}

bool SeqLibAlignmentReader::locateIndexes() {
    // SeqLib automatically handles index loading when SetRegion is called
    // We can try to load explicitly if needed
    // For now, return true as SeqLib handles this internally
    return true;
}

void SeqLibAlignmentReader::setReferenceFile(const std::string& fastaPath) {
    if (!fastaPath.empty()) {
        reader_->SetCramReference(fastaPath);
    }
}

} // namespace freebayes
