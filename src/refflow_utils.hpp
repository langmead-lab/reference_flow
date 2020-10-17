#ifndef REFFLOW_UTILS_H__
#define REFFLOW_UTILS_H__
#include <htslib/sam.h>
#include <vector>

static int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };


static std::string make_cmd(int argc, char** argv) {
    std::string cmd("");
    for (auto i = 0; i < argc; ++i) {
        cmd += std::string(argv[i]) + " ";
    }
    return cmd;
}

static std::string get_read(const bam1_t *rec);

void write_fq_from_bam(bam1_t* aln, std::ofstream& out_fq);

int sam_read1_selective(samFile* sam_fp, bam_hdr_t* hdr, bam1_t* aln,
                        const std::vector<int>& exclude_flag);

#endif /* REFFLOW_UTILS_H__ */
