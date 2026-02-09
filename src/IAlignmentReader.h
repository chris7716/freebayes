#ifndef FREEBAYES_IALIGNMENTREADER_H
#define FREEBAYES_IALIGNMENTREADER_H

#include <string>
#include <vector>
#include <memory>
#include "IAlignment.h"

namespace freebayes {

/**
 * Structure representing a reference sequence in the header
 */
struct ReferenceSequence {
    std::string name;
    int length;
    bool hasAlignments;

    ReferenceSequence(const std::string& n, int len, bool has = true)
        : name(n), length(len), hasAlignments(has) {}
};

/**
 * Interface for reading alignments from various formats.
 * Implementations can support BAM, CRAM, SAM, PAF, or custom formats.
 *
 * Usage pattern:
 *   auto reader = AlignmentReaderFactory::create(format);
 *   reader->open(files);
 *   reader->setRegion(refId, start, end);
 *
 *   std::shared_ptr<IAlignment> aln;
 *   while (reader->getNextAlignment(aln)) {
 *       // process alignment
 *   }
 *
 *   reader->close();
 */
class IAlignmentReader {
public:
    virtual ~IAlignmentReader() = default;

    // ===== Lifecycle Management =====

    /**
     * Open alignment file(s) for reading.
     *
     * @param paths Vector of file paths to open. For multi-file formats like BAM,
     *              all files should be listed. For single-file formats, vector
     *              should contain one element. Special paths like "-" for stdin
     *              may be supported by specific implementations.
     * @return true on success, false on failure (check getErrorString())
     */
    virtual bool open(const std::vector<std::string>& paths) = 0;

    /**
     * Close the alignment reader and release resources.
     */
    virtual void close() = 0;

    /**
     * Check if reader is currently open
     * @return true if open and ready to read
     */
    virtual bool isOpen() const = 0;

    // ===== Navigation =====

    /**
     * Set region to read from (enables random access).
     * For formats that support indexing (BAM, CRAM), this jumps to the region.
     * For non-indexed formats (SAM, PAF), this may require sequential scan.
     *
     * @param refId Reference sequence ID (index into getReferenceSequences())
     * @param start 0-based start position (inclusive)
     * @param end 0-based end position (exclusive)
     * @return true on success, false if region setting not supported or failed
     */
    virtual bool setRegion(int refId, long start, long end) = 0;

    /**
     * Read the next alignment from the current position/region.
     *
     * @param alignment Shared pointer that will be set to the next alignment.
     *                  The pointer is reused/reset on each call for efficiency.
     * @return true if an alignment was read, false if EOF or error
     */
    virtual bool getNextAlignment(std::shared_ptr<IAlignment>& alignment) = 0;

    // ===== Metadata Access =====

    /**
     * Get the complete header text from the alignment file.
     * For SAM/BAM/CRAM, this is the SAM header.
     * For other formats, may return format-specific metadata.
     *
     * @return Header text as a string
     */
    virtual std::string getHeaderText() const = 0;

    /**
     * Get list of reference sequences from the file header.
     *
     * @return Vector of reference sequences with names and lengths
     */
    virtual std::vector<ReferenceSequence> getReferenceSequences() const = 0;

    /**
     * Get reference sequence ID from name.
     *
     * @param name Reference sequence name
     * @return Reference ID (non-negative) if found, -1 if not found
     */
    virtual int getReferenceId(const std::string& name) const = 0;

    // ===== Error Handling =====

    /**
     * Get description of last error that occurred.
     *
     * @return Human-readable error message
     */
    virtual std::string getErrorString() const = 0;

    // ===== Optional: Index Management =====

    /**
     * Attempt to locate/load index files for random access.
     * For BAM, looks for .bai or .bam.bai files.
     * For CRAM, looks for .crai files.
     *
     * @return true if index found/loaded, false otherwise
     */
    virtual bool locateIndexes() = 0;

    /**
     * Set reference FASTA file (required for CRAM format).
     *
     * @param fastaPath Path to reference FASTA file
     */
    virtual void setReferenceFile(const std::string& fastaPath) = 0;
};

} // namespace freebayes

#endif // FREEBAYES_IALIGNMENTREADER_H
