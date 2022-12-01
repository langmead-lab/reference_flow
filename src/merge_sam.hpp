#ifndef MERGE_SAM_H__
#define MERGE_SAM_H__
#include <tuple>
#include <vector>

#include <htslib/sam.h>


const int MERGE_PE_SUM = 0;
const int MERGE_PE_MAX = 1;
// const int MERGE_PE_PAIR_TLEN_THRESHOLD = 2000;

struct merge_sam_opts{
    bool paired_end = false;
    std::string cmd = "";
    int decoy_threshold = 0;
    int rand_seed = 0;
    // std::string sam_list = "";
    // std::string id_list = "";
    std::string decoy_list = "";
    std::string output_prefix = "";
    std::vector<std::string> inputs;
};

static void merge_sam_help(){
    std::cerr << "\n";
    std::cerr << "Usage: refflow_utils merge [options] -s <input_list> -i <id_list> -o <out>\n";
    std::cerr << "\n";
    std::cerr << "  <input_list>Path to a list of SAM/BAM files to merge. One path in each line.\n";
    std::cerr << "  <id_list>   Path to a list of labels corresponding to lines in <input_list>\n";
    std::cerr << "  <out>       Prefix of the output files\n";
    std::cerr << "\n";
    std::cerr << "Options (defaults in parentheses):\n";
    std::cerr << "  -m          Set to perform merging in pairs [off]\n";
    std::cerr << "  -r <int>    Random seed used by the program [0]\n";
    // TODO
    // std::cerr << "  -O <text>   Output alignment format, can be sam or bam [sam]\n";
    // decoy_list
    // decoy_threshold
    std::cerr << "\n";
}

template <typename A, typename B, typename C, typename D>
std::vector<std::tuple<A, B, C, D>> zip(const std::vector<A> &a,
                                        const std::vector<B> &b,
                                        const std::vector<C> &c,
                                        const std::vector<D> &d);

template <typename A, typename B, typename C, typename D>
void unzip(const std::vector<std::tuple<A, B, C, D>> &zipped,
           std::vector<A> &a,
           std::vector<B> &b,
           std::vector<C> &c,
           std::vector<D> &d);

/* Random generator function
 *
 * From http://www.cplusplus.com/reference/algorithm/random_shuffle/
 */
static int myrandom (int i) { return std::rand() % i;}

std::vector<std::string> read_file_as_vector(std::string list_fn);

int select_best_aln(const std::vector<bool>& pair_indicators,
                    const std::vector<int>& scores,
                    const std::vector<int>& mapqs,
                    int& num_tied_best);

int select_best_aln_paired_end(const std::vector<bam1_t*>& aln1s,
                               const std::vector<bam1_t*>& aln2s,
                               const int merge_pe_mode);

void merge_sam(merge_sam_opts args);

void merge_sam_main(int argc, char** argv);

#endif /* MERGE_SAM_H__ */
