#include "GBAMAlignmentReader.h"
#include "../gbam/reader.hpp"
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdexcept>
#include <cstring>

namespace freebayes {

// ── GBAMAlignment ─────────────────────────────────────────────────────────────

GBAMAlignment::GBAMAlignment(bam1_t* record) : record_(record) {}

GBAMAlignment::~GBAMAlignment() {
    if (record_) bam_destroy1(record_);
}

std::string GBAMAlignment::queryName() const {
    return bam_get_qname(record_);
}

std::string GBAMAlignment::readGroup() const {
    const uint8_t* rg = bam_aux_get(record_, "RG");
    if (!rg) return "";
    const char* val = bam_aux2Z(rg);
    return val ? val : "";
}

int GBAMAlignment::referenceId() const {
    return record_->core.tid;
}

long GBAMAlignment::position() const {
    return record_->core.pos;
}

long GBAMAlignment::endPosition() const {
    return bam_endpos(record_);
}

int GBAMAlignment::mappingQuality() const {
    return record_->core.qual;
}

std::string GBAMAlignment::sequence() const {
    static const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
    int l = record_->core.l_qseq;
    const uint8_t* seq = bam_get_seq(record_);
    std::string result(l, 'N');
    for (int i = 0; i < l; ++i)
        result[i] = seq_nt16_str[bam_seqi(seq, i)];
    return result;
}

std::string GBAMAlignment::qualities() const {
    int l = record_->core.l_qseq;
    const uint8_t* qual = bam_get_qual(record_);
    // 0xff sentinel means qualities absent
    if (l > 0 && qual[0] == 0xff)
        return std::string(l, 'I');
    std::string result(l, '\0');
    for (int i = 0; i < l; ++i)
        result[i] = static_cast<char>(qual[i] + 33);
    return result;
}

int GBAMAlignment::alignmentLength() const {
    return bam_cigar2rlen(record_->core.n_cigar, bam_get_cigar(record_));
}

int GBAMAlignment::sequenceLength() const {
    return record_->core.l_qseq;
}

bool GBAMAlignment::isMapped() const {
    return !(record_->core.flag & BAM_FUNMAP);
}

bool GBAMAlignment::isPaired() const {
    return record_->core.flag & BAM_FPAIRED;
}

bool GBAMAlignment::isMateMapped() const {
    return !(record_->core.flag & BAM_FMUNMAP);
}

bool GBAMAlignment::isProperPair() const {
    return record_->core.flag & BAM_FPROPER_PAIR;
}

bool GBAMAlignment::isReverseStrand() const {
    return record_->core.flag & BAM_FREVERSE;
}

bool GBAMAlignment::isDuplicate() const {
    return record_->core.flag & BAM_FDUP;
}

std::vector<CigarOp> GBAMAlignment::getCigar() const {
    static const char cigar_ops[] = "MIDNSHP=X";
    std::vector<CigarOp> result;
    const uint32_t* cigar = bam_get_cigar(record_);
    int n = record_->core.n_cigar;
    result.reserve(n);
    for (int i = 0; i < n; ++i) {
        int op_idx = bam_cigar_op(cigar[i]);
        char op = (op_idx < 9) ? cigar_ops[op_idx] : '?';
        result.emplace_back(op, bam_cigar_oplen(cigar[i]));
    }
    return result;
}

bool GBAMAlignment::getTag(const std::string& tag, std::string& value) const {
    if (tag.size() != 2) return false;
    const uint8_t* p = bam_aux_get(record_, tag.c_str());
    if (!p) return false;
    char type = *p;
    if (type == 'Z' || type == 'H') {
        const char* s = bam_aux2Z(p);
        if (!s) return false;
        value = s;
        return true;
    }
    // For numeric types, convert to string
    if (type == 'i' || type == 'I' || type == 's' || type == 'S' ||
        type == 'c' || type == 'C') {
        value = std::to_string(bam_aux2i(p));
        return true;
    }
    if (type == 'f' || type == 'd') {
        value = std::to_string(bam_aux2f(p));
        return true;
    }
    return false;
}

bool GBAMAlignment::getTag(const std::string& tag, int& value) const {
    if (tag.size() != 2) return false;
    const uint8_t* p = bam_aux_get(record_, tag.c_str());
    if (!p) return false;
    char type = *p;
    if (type == 'i' || type == 'I' || type == 's' || type == 'S' ||
        type == 'c' || type == 'C') {
        value = bam_aux2i(p);
        return true;
    }
    return false;
}

bool GBAMAlignment::getTag(const std::string& tag, double& value) const {
    if (tag.size() != 2) return false;
    const uint8_t* p = bam_aux_get(record_, tag.c_str());
    if (!p) return false;
    char type = *p;
    if (type == 'f' || type == 'd') {
        value = bam_aux2f(p);
        return true;
    }
    if (type == 'i' || type == 'I' || type == 's' || type == 'S' ||
        type == 'c' || type == 'C') {
        value = static_cast<double>(bam_aux2i(p));
        return true;
    }
    return false;
}

// ── GBAMAlignmentReader ───────────────────────────────────────────────────────

GBAMAlignmentReader::GBAMAlignmentReader() = default;

GBAMAlignmentReader::~GBAMAlignmentReader() {
    close();
}

bool GBAMAlignmentReader::open(const std::vector<std::string>& paths) {
    if (paths.size() != 1) {
        lastError_ = "GBAM reader supports exactly one file at a time";
        return false;
    }

    const std::string& path = paths[0];
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        lastError_ = "Cannot open file: " + path;
        return false;
    }

    struct stat st;
    if (::fstat(fd, &st) < 0) {
        ::close(fd);
        lastError_ = "Cannot stat file: " + path;
        return false;
    }
    mmap_size_ = static_cast<size_t>(st.st_size);

    mmap_addr_ = ::mmap(nullptr, mmap_size_, PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);
    if (mmap_addr_ == MAP_FAILED) {
        mmap_addr_ = nullptr;
        lastError_ = "mmap failed for: " + path;
        return false;
    }

    try {
        reader_ = std::make_unique<Reader>(static_cast<const char*>(mmap_addr_));
    } catch (const std::exception& e) {
        ::munmap(mmap_addr_, mmap_size_);
        mmap_addr_ = nullptr;
        lastError_ = std::string("Failed to parse GBAM file: ") + e.what();
        return false;
    }

    buf_ = bam_init1();
    int64_t n = reader_->totalRecords();
    current_rec_    = 0;
    region_end_rec_ = n;
    isOpen_ = true;
    return true;
}

void GBAMAlignmentReader::close() {
    reader_.reset();
    if (buf_) { bam_destroy1(buf_); buf_ = nullptr; }
    if (mmap_addr_) {
        ::munmap(mmap_addr_, mmap_size_);
        mmap_addr_ = nullptr;
        mmap_size_ = 0;
    }
    isOpen_ = false;
}

bool GBAMAlignmentReader::isOpen() const {
    return isOpen_;
}

bool GBAMAlignmentReader::setRegion(int refId, long start, long end) {
    if (!isOpen_) {
        lastError_ = "Cannot set region: reader not open";
        return false;
    }

    int64_t n = reader_->totalRecords();

    // Binary search: first record where refID >= refId
    int64_t lo = 0, hi = n;
    while (lo < hi) {
        int64_t mid = (lo + hi) / 2;
        if (reader_->getRefId(mid) < refId) lo = mid + 1;
        else hi = mid;
    }
    if (lo == n || reader_->getRefId(lo) != refId) {
        // No records for this reference
        current_rec_    = n;
        region_end_rec_ = n;
        return true;
    }
    int64_t ref_start = lo;

    // Binary search: first record where refID > refId (upper bound of refId range)
    lo = ref_start; hi = n;
    while (lo < hi) {
        int64_t mid = (lo + hi) / 2;
        if (reader_->getRefId(mid) == refId) lo = mid + 1;
        else hi = mid;
    }
    int64_t ref_end = lo;  // exclusive

    // Binary search within [ref_start, ref_end): first pos >= start
    lo = ref_start; hi = ref_end;
    while (lo < hi) {
        int64_t mid = (lo + hi) / 2;
        if (reader_->getPos(mid) < start) lo = mid + 1;
        else hi = mid;
    }
    current_rec_ = lo;

    // Binary search: first record with pos >= end (still within refId range)
    lo = current_rec_; hi = ref_end;
    while (lo < hi) {
        int64_t mid = (lo + hi) / 2;
        if (reader_->getPos(mid) < end) lo = mid + 1;
        else hi = mid;
    }
    region_end_rec_ = lo;

    return true;
}

bool GBAMAlignmentReader::getNextAlignment(std::shared_ptr<IAlignment>& alignment) {
    if (!isOpen_) {
        lastError_ = "Cannot read alignment: reader not open";
        return false;
    }
    if (current_rec_ >= region_end_rec_) {
        return false;
    }

    try {
        reader_->read_record(current_rec_, buf_);
    } catch (const std::exception& e) {
        lastError_ = std::string("read_record failed: ") + e.what();
        return false;
    }
    ++current_rec_;

    // Deep-copy the decoded record and hand ownership to GBAMAlignment.
    bam1_t* copy = bam_dup1(buf_);
    alignment = std::make_shared<GBAMAlignment>(copy);
    return true;
}

std::string GBAMAlignmentReader::getHeaderText() const {
    if (!isOpen_ || !reader_->header()) return "";
    char* text = sam_hdr_str(reader_->header());
    return text ? std::string(text) : "";
}

std::vector<ReferenceSequence> GBAMAlignmentReader::getReferenceSequences() const {
    std::vector<ReferenceSequence> result;
    if (!isOpen_ || !reader_->header()) return result;
    sam_hdr_t* hdr = reader_->header();
    int n = sam_hdr_nref(hdr);
    for (int i = 0; i < n; ++i) {
        const char* name = sam_hdr_tid2name(hdr, i);
        hts_pos_t   len  = sam_hdr_tid2len(hdr, i);
        result.emplace_back(name ? name : "", static_cast<int>(len), true);
    }
    return result;
}

int GBAMAlignmentReader::getReferenceId(const std::string& name) const {
    if (!isOpen_ || !reader_->header()) return -1;
    return sam_hdr_name2tid(reader_->header(), name.c_str());
}

std::string GBAMAlignmentReader::getErrorString() const {
    return lastError_;
}

bool GBAMAlignmentReader::locateIndexes() {
    // GBAM uses in-file metadata; no external index needed.
    return true;
}

void GBAMAlignmentReader::setReferenceFile(const std::string& /*fastaPath*/) {
    // GBAM stores full sequences; no external reference required.
}

} // namespace freebayes
