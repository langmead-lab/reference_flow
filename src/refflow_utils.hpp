#ifndef REFFLOW_UTILS_H__
#define REFFLOW_UTILS_H__
static int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };


const int MERGE_PE_SUM = 0;
const int MERGE_PE_MAX = 1;
// const int MERGE_PE_PAIR_TLEN_THRESHOLD = 2000;

std::string make_cmd(int argc, char** argv);

static std::string get_read(const bam1_t *rec);

void write_fq_from_bam(bam1_t* aln, std::ofstream& out_fq);

int sam_read1_selective(samFile* sam_fp, bam_hdr_t* hdr, bam1_t* aln,
                        const std::vector<int>& exclude_flag);

struct merge_sam_opts{
    bool paired_end = false;
    std::string cmd = "";
    int decoy_threshold = 0;
    int rand_seed = 0;
    std::string sam_list = "";
    std::string id_list = "";
    std::string decoy_list = "";
    std::string output_prefix = "";
};

#endif /* REFFLOW_UTILS_H__ */
