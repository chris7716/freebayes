#ifndef FREEBAYES_ALIGNMENTREADERFACTORY_H
#define FREEBAYES_ALIGNMENTREADERFACTORY_H

#include "IAlignmentReader.h"
#include <memory>
#include <string>
#include <map>

namespace freebayes {

/**
 * Supported alignment file formats
 */
enum class AlignmentFormat {
    BAM,            // Binary Alignment Map
    CRAM,           // Reference-based compressed alignment
    SAM,            // Sequence Alignment/Map (text)
    GBAM          // GABM format
};

/**
 * Factory for creating alignment readers based on format.
 *
 * This factory provides:
 * - Format auto-detection from file paths
 * - Reader instantiation for different formats
 * - Configuration option passing to readers
 *
 * Example usage:
 *   auto format = AlignmentReaderFactory::detectFormat("alignments.bam");
 *   auto reader = AlignmentReaderFactory::create(format);
 *   reader->open({"alignments.bam"});
 */
class AlignmentReaderFactory {
public:
    /**
     * Create an alignment reader for the specified format.
     *
     * @param format The alignment format to read
     * @param options Optional configuration parameters (format-specific)
     *                Common options:
     *                - "reference": Path to reference FASTA (for CRAM)
     *                - "threads": Number of threads for decompression
     * @return Unique pointer to reader implementation
     * @throws std::runtime_error if format is not supported
     */
    static std::unique_ptr<IAlignmentReader> create(
        AlignmentFormat format,
        const std::map<std::string, std::string>& options = {}
    );

    /**
     * Detect alignment format from file path.
     * Checks file extension and, if available, magic bytes.
     *
     * @param path Path to alignment file
     * @return Detected format (defaults to BAM if cannot determine)
     */
    static AlignmentFormat detectFormat(const std::string& path);

    /**
     * Detect format from file extension only.
     *
     * @param path File path
     * @return Format based on extension (defaults to BAM if unknown)
     */
    static AlignmentFormat detectFormatFromExtension(const std::string& path);

    /**
     * Detect format from file magic bytes.
     * Requires file to exist and be readable.
     *
     * @param path File path
     * @return Format based on magic bytes (defaults to BAM if cannot read)
     */
    static AlignmentFormat detectFormatFromMagic(const std::string& path);

    /**
     * Convert format enum to string name.
     *
     * @param format Format enum value
     * @return Human-readable format name
     */
    static std::string formatToString(AlignmentFormat format);

    /**
     * Convert string name to format enum.
     *
     * @param name Format name (case-insensitive: "bam", "cram", "sam", "gbam")
     * @return Format enum value (defaults to BAM if not recognized)
     */
    static AlignmentFormat stringToFormat(const std::string& name);
};

} // namespace freebayes

#endif // FREEBAYES_ALIGNMENTREADERFACTORY_H
