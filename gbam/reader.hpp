#pragma once

#include "defs.hpp"
#include "codec.hpp"
#include <array>
#include <cstdint>
#include <vector>
#include <htslib/sam.h>

class Reader {
public:
    explicit Reader(const char* mmaped_file);
    ~Reader();

    // Read record rec_num into a pre-allocated bam1_t.
    void read_record(int64_t rec_num, bam1_t* aln);

    // Total number of records in the file.
    int64_t totalRecords() const { return rec_num_; }

    // Fast single-field accessors for region seeking (reads only the needed column).
    int32_t getRefId(int64_t rec_num);
    int32_t getPos(int64_t rec_num);

    // Access the embedded SAM header.
    sam_hdr_t* header() const { return header_; }

private:
    void fetch_field(int64_t rec_num, int col);

    const char* mmaped_file_;
    std::array<std::vector<ColumnChunkMeta>, COLUMNTYPE_SIZE> metadatas_;
    sam_hdr_t* header_ = nullptr;

    std::array<std::vector<int64_t>, COLUMNTYPE_SIZE> record_counts_per_chunk_;
    int64_t rec_num_ = 0;

    std::array<int64_t, COLUMNTYPE_SIZE> currently_loaded_chunk_;
    std::array<int64_t, COLUMNTYPE_SIZE> loaded_up_to_rec_num_;
    std::array<int64_t, COLUMNTYPE_SIZE> loaded_since_rec_num_;
    std::array<int64_t, COLUMNTYPE_SIZE> m_chunk_memory_;
    std::array<Column, COLUMNTYPE_SIZE>  columns_;
};
