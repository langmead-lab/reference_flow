
const int SPLIT_OPT = 0;
const int SPLIT_PES = 1;
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

struct split_sam_opts{
    std::string cmd = "";
    bool write_hq_to_stdout = false;
    int split_strategy = SPLIT_OPT;
    int mapq_threshold = 0;
    std::string sam_fn = "";
    std::string output_prefix = "";
};

const int MERGE_PE_SUM = 0;
const int MERGE_PE_MAX = 1;
// const int MERGE_PE_PAIR_TLEN_THRESHOLD = 2000;

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
