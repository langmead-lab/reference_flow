#ifndef SPLIT_SAM_H__
#define SPLIT_SAM_H__

const int SPLIT_OPT = 0;
const int SPLIT_PES = 1;


struct split_sam_opts{
    std::string cmd = "";
    bool write_hq_to_stdout = false;
    int split_strategy = SPLIT_PES;
    int mapq_threshold = 10;
    std::string sam_fn = "";
    std::string output_prefix = "";
    std::string output_ext = "sam";
};

static void split_sam_help(){
    std::cerr << "\n";
    std::cerr << "Usage: refflow_utils split [options] -s <input> -o <out>\n";
    std::cerr << "\n";
    std::cerr << "  <input>     SAM or BAM file to split\n";
    std::cerr << "  <out>       Prefix of the output files\n";
    std::cerr << "\n";
    std::cerr << "Options (defaults in parentheses):\n";
    std::cerr << "  -q <int>    MAPQ threshold of splitting [10]\n";
    std::cerr << "  -O <text>   Output alignment format, can be sam or bam [sam]\n";
    std::cerr << "  -p <int>    Split strategy: set to 0 for the optimistic mode\n";
    std::cerr << "                              set to 1 for the pessimistic mode (1)\n";
    std::cerr << "              The optmistic mode treats a pair as high-quality if either alignments has >= threshold MAPQ\n";
    std::cerr << "              The pessimistic mode requires both segments with >= threshold MAPQ to be considered high-quality\n";
    std::cerr << "  -d          Set to write high-quality alignments to stdout\n";
    std::cerr << "\n";
}

void split_sam(split_sam_opts args);

void split_sam_main(int argc, char** argv);


#endif /* SPLIT_SAM_H__ */
