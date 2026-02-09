#ifndef FREEBAYES_IALIGNMENT_H
#define FREEBAYES_IALIGNMENT_H

#include <string>
#include <vector>
#include <memory>

namespace freebayes {

// Forward declaration for CIGAR operations
struct CigarOp {
    char type;      // 'M', 'I', 'D', 'S', 'H', 'N', 'P', '=', 'X'
    uint32_t length;

    CigarOp(char t, uint32_t len) : type(t), length(len) {}
};

/**
 * Interface representing a single alignment/read.
 * Abstracts away the underlying format (BAM, CRAM, SAM, PAF, etc.)
 * and provides a uniform API for accessing alignment properties.
 */
class IAlignment {
public:
    virtual ~IAlignment() = default;

    // ===== Core Properties =====

    /**
     * @return Query (read) name
     */
    virtual std::string queryName() const = 0;

    /**
     * @return Read group identifier (from RG tag)
     */
    virtual std::string readGroup() const = 0;

    /**
     * @return Reference sequence ID (index into reference list)
     */
    virtual int referenceId() const = 0;

    /**
     * @return 0-based leftmost mapping position
     */
    virtual long position() const = 0;

    /**
     * @return 0-based rightmost mapping position (exclusive)
     */
    virtual long endPosition() const = 0;

    /**
     * @return Mapping quality (0-255)
     */
    virtual int mappingQuality() const = 0;

    /**
     * @return Read sequence (A/C/G/T/N)
     */
    virtual std::string sequence() const = 0;

    /**
     * @return Base quality scores (Phred+33 encoded string)
     */
    virtual std::string qualities() const = 0;

    /**
     * @return Length of the alignment on the reference
     */
    virtual int alignmentLength() const = 0;

    /**
     * @return Length of the query sequence
     */
    virtual int sequenceLength() const = 0;

    // ===== SAM Flags =====

    /**
     * @return true if read is mapped to reference
     */
    virtual bool isMapped() const = 0;

    /**
     * @return true if read is paired in sequencing
     */
    virtual bool isPaired() const = 0;

    /**
     * @return true if mate is mapped
     */
    virtual bool isMateMapped() const = 0;

    /**
     * @return true if read is proper pair (both mates mapped correctly)
     */
    virtual bool isProperPair() const = 0;

    /**
     * @return true if read maps to reverse strand
     */
    virtual bool isReverseStrand() const = 0;

    /**
     * @return true if read is marked as duplicate
     */
    virtual bool isDuplicate() const = 0;

    // ===== CIGAR Operations =====

    /**
     * @return CIGAR operations describing alignment
     */
    virtual std::vector<CigarOp> getCigar() const = 0;

    // ===== Tags (Optional Fields) =====

    /**
     * Get a string tag value
     * @param tag Two-character tag name (e.g., "RG", "NM")
     * @param value Output parameter to store tag value
     * @return true if tag exists, false otherwise
     */
    virtual bool getTag(const std::string& tag, std::string& value) const = 0;

    /**
     * Get an integer tag value
     * @param tag Two-character tag name
     * @param value Output parameter to store tag value
     * @return true if tag exists, false otherwise
     */
    virtual bool getTag(const std::string& tag, int& value) const = 0;

    /**
     * Get a float tag value
     * @param tag Two-character tag name
     * @param value Output parameter to store tag value
     * @return true if tag exists, false otherwise
     */
    virtual bool getTag(const std::string& tag, double& value) const = 0;
};

} // namespace freebayes

#endif // FREEBAYES_IALIGNMENT_H
