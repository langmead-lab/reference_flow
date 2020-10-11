#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

#include <getopt.h>
#include <htslib/sam.h>
#include <limits.h>
#include <string.h>

#include <refflow_utils.hpp>

std::string make_cmd(int argc, char** argv) {
    std::string cmd("");
    for (auto i = 0; i < argc; ++i) {
        cmd += std::string(argv[i]) + " ";
    }
    return cmd;
}


/* Return the read, reverse complemented if necessary
   Adapted from: https://github.com/samtools/samtools/blob/develop/bam_fastq.c 
*/
static std::string get_read(const bam1_t *rec){
    int len = rec->core.l_qseq + 1;
    char *seq = (char *)bam_get_seq(rec);
    std::string read = "";

    for (int n = 0; n < rec->core.l_qseq; n++) {
        if (rec->core.flag & BAM_FREVERSE)
            read.append(1, seq_nt16_str[seq_comp_table[bam_seqi(seq, n)]]);
        else
            read.append(1, seq_nt16_str[bam_seqi(seq, n)]);
    }
    if (rec->core.flag & BAM_FREVERSE)
        std::reverse(read.begin(), read.end());
    return read;
}


/* Write a bam1_t object to a FASTQ record.
 */
void write_fq_from_bam(bam1_t* aln, std::ofstream& out_fq){
    out_fq << "@" << bam_get_qname(aln) << "\n";
    out_fq << get_read(aln) << "\n+\n";
    std::string qual_seq("");
    uint8_t* qual = bam_get_qual(aln);
    if (qual[0] == 255) qual_seq = "*";
    else {
        for (auto i = 0; i < aln->core.l_qseq; ++i) {
            qual_seq += (char) (qual[i] + 33);
        }
    }
    if (aln->core.flag & BAM_FREVERSE)
        std::reverse(qual_seq.begin(), qual_seq.end());
    out_fq << qual_seq << "\n";
}


/* Fectch alignments that are unapped or mapped with low MAPQ.
 * 
 */
void split_sam(split_sam_opts args){
    // Read raw SAM file.
    samFile* sam_fp = (args.sam_fn == "")?
        sam_open("-", "r") :
        sam_open(args.sam_fn.data(), "r");
    bam_hdr_t* hdr = sam_hdr_read(sam_fp);

    // High-quality alignments (SAM).
    std::string out_sam_hq_fn = args.output_prefix + "-high_qual.sam";
    samFile* out_sam_hq_fp = (args.write_hq_to_stdout)?
        sam_open("-", "w") :
        sam_open(out_sam_hq_fn.data(), "w");
    int write_hdr = sam_hdr_write(out_sam_hq_fp, hdr);
    // Low-quality alignments (SAM).
    std::string out_sam_lq_fn = args.output_prefix + "-low_qual.sam";
    samFile* out_sam_lq_fp = sam_open(out_sam_lq_fn.data(), "w");
    write_hdr = sam_hdr_write(out_sam_lq_fp, hdr);
    // Low-quality reads (paired-end FQ).
    std::ofstream out_fq1, out_fq2;
    out_fq1.open(args.output_prefix + "-R1.fq");
    out_fq2.open(args.output_prefix + "-R2.fq");

    // Read two reads in one iteration to for the paired-end mode.
    bam1_t* aln1 = bam_init1(), * aln2 = bam_init1();
    while(true){
        if (sam_read1(sam_fp, hdr, aln1) < 0 || sam_read1(sam_fp, hdr, aln2) < 0)
            break;
        // Check read names: they should be identical.
        if (strcmp(bam_get_qname(aln1), bam_get_qname(aln2)) != 0){
            std::cerr << "[Error] Input SAM file should be sorted by read name.\n";
            exit(1);
        }
        bam1_core_t c_aln1 = aln1->core, c_aln2 = aln2->core;

        bool is_low_quality = false;
        if (c_aln1.tid != c_aln1.mtid)
            is_low_quality = true;
        else if (args.split_strategy == SPLIT_OPT)
            is_low_quality = (c_aln1.qual < args.mapq_threshold && c_aln2.qual < args.mapq_threshold);
        else if (args.split_strategy == SPLIT_PES)
            is_low_quality = (c_aln1.qual < args.mapq_threshold || c_aln2.qual < args.mapq_threshold);

        if (is_low_quality){
            if (aln1->core.flag & BAM_FREAD1){
                write_fq_from_bam(aln1, out_fq1);
                write_fq_from_bam(aln2, out_fq2);
            }
            else{
                write_fq_from_bam(aln1, out_fq2);
                write_fq_from_bam(aln2, out_fq1);
            }
            if (sam_write1(out_sam_lq_fp, hdr, aln1) < 0 || sam_write1(out_sam_lq_fp, hdr, aln2) < 0){
                std::cerr << "[Error] Failure when writing low-qualiy alignments to SAM.\n";
                exit(1);
            }
        }
        else{
            if (sam_write1(out_sam_hq_fp, hdr, aln1) < 0 || sam_write1(out_sam_hq_fp, hdr, aln2) < 0){
                std::cerr << "[Error] Failure when writing high-qualiy alignments to SAM.\n";
                exit(1);
            }
        }
    }
    sam_close(sam_fp);
    sam_close(out_sam_hq_fp);
    sam_close(out_sam_lq_fp);
    out_fq1.close();
    out_fq2.close();
}


void split_sam_main(int argc, char** argv){
    int c;
    split_sam_opts args;
    args.cmd = make_cmd(argc, argv);
    static struct option long_options[]{
        {"write_hq_to_stdout", no_argument, 0, 'd'},
        {"sam_fn", required_argument, 0, 's'},
        {"output_prefix", required_argument, 0, 'o'},
        {"split_strategy", required_argument, 0, 'p'},
        {"mapq_threshold", required_argument, 0, 'q'}
    };
    int long_idx = 0;
    while((c = getopt_long(argc, argv, "dho:p:q:s:", long_options, &long_idx)) != -1){
        switch (c){
            case 'd':
                args.write_hq_to_stdout = true;
                break;
            case 'o':
                args.output_prefix = optarg;
                break;
            case 'p':
                args.split_strategy = atoi(optarg);
                break;
            case 'q':
                args.mapq_threshold = atoi(optarg);
                break;
            case 's':
                args.sam_fn = optarg;
                break;
            default:
                std::cerr << "Ignoring option " << c << " \n";
                exit(1);
        }
    }
    split_sam(args);
}


/* Permutation used to sort a vector.
 *
 * From Timothy Shields,
 * https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
 */
template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
    const std::vector<T>& vec,
    Compare compare)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
        [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
    return p;
}


/* Build a new std::vector<T> that is reordered according to the permutation.
 *
 * Example:
 *     vector<MyObject> vectorA;
 *     vector<int> vectorB;
 *     auto p = sort_permutation(vectorA,
 *         [](T const& a, T const& b){ / *some comparison* / });
 *     vectorA = apply_permutation(vectorA, p);
 *     vectorB = apply_permutation(vectorB, p);
 *
 * From Timothy Shields,
 * https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
 */
template <typename T>
std::vector<T> apply_permutation(
    const std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
        [&](std::size_t i){ return vec[i]; });
    return sorted_vec;
}


// Fill the zipped vector with pairs consisting of the
// corresponding elements of a and b. (This assumes 
// that the vectors have equal length)
// https://stackoverflow.com/questions/37368787/c-sort-one-vector-based-on-another-one/46370189
template <typename A, typename B, typename C, typename D>
std::vector<std::tuple<A, B, C, D>> zip(const std::vector<A> &a,
                                        const std::vector<B> &b,
                                        const std::vector<C> &c,
                                        const std::vector<D> &d){
    std::vector<std::tuple<A, B, C, D>> zipped;
    for(int i = 0; i < a.size(); ++i){
        zipped.push_back(std::make_tuple(a[i], b[i], c[i], d[i]));
    }
    return (zipped);
}


// Write the first and second element of the pairs in 
// the given zipped vector into a and b. (This assumes 
// that the vectors have equal length)
// https://stackoverflow.com/questions/37368787/c-sort-one-vector-based-on-another-one/46370189
template <typename A, typename B, typename C, typename D>
void unzip(const std::vector<std::tuple<A, B, C, D>> &zipped,
           std::vector<A> &a,
           std::vector<B> &b,
           std::vector<C> &c,
           std::vector<D> &d){
    for(int i = 0; i < a.size(); i++){
        a[i] = std::get<0>(zipped[i]);
        b[i] = std::get<1>(zipped[i]);
        c[i] = std::get<2>(zipped[i]);
        d[i] = std::get<3>(zipped[i]);
    }
}


/* Random generator function
 *
 * From http://www.cplusplus.com/reference/algorithm/random_shuffle/
 */
int myrandom (int i) { return std::rand() % i;}


std::vector<std::string> read_file_as_vector(std::string list_fn){
    std::vector<std::string> list;
    std::ifstream read_list(list_fn.c_str());
    if (!read_list){
        std::cerr << "[Error] Cannot open file " << list_fn << "\n";
        exit(1);
    }
    std::string line;
    while (std::getline(read_list, line)){
        if (line.size() > 0){
            list.push_back(line);
        }
    }
    read_list.close();
    return (list);
}

/* Rank a set of alignments.
 * Order by: is_proper_pair > score > MAPQ
 */
int rank_alns(const std::vector<bool>& pair_indicators,
              const std::vector<int>& scores,
              const std::vector<int>& mapqs){
    int vec_size = pair_indicators.size();
    std::vector<int> ranks(vec_size);
    std::iota(ranks.begin(), ranks.end(), 0);

    std::vector<std::tuple<int, int, bool, int>> zipped = zip(mapqs, scores, pair_indicators, ranks);
    // for(int i = 0; i < vec_size; i++){
    //     std::cerr << std::get<1>(zipped[i]) << " ";
    // }
    // std::cerr << "\n";
    // Sort the vector of tuples.
    std::sort(zipped.begin(), zipped.end(),
        [&](const auto& a, const auto& b)
        {
            return (std::get<2>(a) != std::get<2>(b))? std::get<2>(a) :
                   (std::get<1>(a) != std::get<1>(b))? (std::get<1>(a) > std::get<1>(b)) :
                   (std::get<0>(a) > std::get<0>(b));
        });
    int num_ties = vec_size;
    for (int i = 0; i < vec_size - 1; i++){
        // Compare elements in tuples (ranks should be excluded so cannot simply compare tuples).
        if (std::get<2>(zipped[i]) != std::get<2>(zipped[i+1]) ||
            std::get<1>(zipped[i]) != std::get<1>(zipped[i+1]) ||
            std::get<0>(zipped[i]) != std::get<0>(zipped[i+1])){
            num_ties = i + 1;
            break;
        }
    }
    ranks.clear();

    for(int i = 0; i < vec_size; i++){
        ranks.push_back(std::get<3>(zipped[i]));
    }
    // std::cerr << "\n";
    // if (num_ties != vec_size){
    //     std::cerr << "num_ties " << num_ties << "\n";
    //     for (auto rs : ranks) { std::cerr << rs << ' ';};
    //     std::cerr << "\n";
    // }
    std::random_shuffle(ranks.begin(), ranks.begin() + num_ties, myrandom);
    // if (num_ties != vec_size){
    //     for (auto rs : ranks) { std::cerr << rs << ' ';};
    //     std::cerr << "\n";
    //     std::cerr << "best: " << ranks[0] << "\n";

    //     std::cerr << "Enter 0 to exit\n";
    //     int cin_a;
    //     std::cin >> cin_a;
    //     if (cin_a == 0) exit(0);
    // }
    return (ranks[0]);
}

/* Rank a set of alignments.
 * Order by: is_proper_pair > score > MAPQ
 */
//int rank_alns(const std::vector<int>& scores,
//              const std::vector<int>& mapqs,
//              const std::vector<bool>& pair_indicators){
//    std::vector<int> rank(scores.size());
//    std::iota(rank.begin(), rank.end(), 0);
//
//    auto sort_mapqs = sort_permutation(mapqs, std::greater<int>());
//    auto rank_sorted_mapq = apply_permutation(rank, sort_mapqs);
//    std:: cerr << "sort by mapq\n";
//    for (auto rs : rank_sorted_mapq) { std::cerr << rs << ' ';};
//    std:: cerr << "\n";
//
//    auto sort_scores = sort_permutation(scores, std::greater<int>());
//    auto rank_sorted_mapq_score = apply_permutation(rank_sorted_mapq, sort_scores);
//    std:: cerr << "sort by score\n";
//    for (auto rs : rank_sorted_mapq_score) { std::cerr << rs << ' ';};
//    std:: cerr << "\n";
//
//    auto sort_pind = sort_permutation(pair_indicators, std::greater<int>());
//    auto rank_sorted_all = apply_permutation(rank_sorted_mapq_score, sort_pind);
//    std:: cerr << "sort by pairness\n";
//    for (auto rs : rank_sorted_all) { std::cerr << rs << ' ';};
//    std:: cerr << "\n";
//
//    int num_ties = rank_sorted_all.size();
//    for (int i = 0; i < rank_sorted_all.size() - 1; i++){
//        int actual_idx = rank_sorted_all[i];
//        int next_idx = rank_sorted_all[i + 1];
//        if (pair_indicators[actual_idx] != pair_indicators[next_idx] ||
//            scores[actual_idx] != scores[next_idx] ||
//            mapqs[actual_idx] != mapqs[next_idx]){
//            num_ties = i;
//            break;
//        }
//    }
//    std::cerr << "num_ties = " << num_ties << "\n";
//    for (auto r : rank) { std::cerr << r << ' ';};
//    std:: cerr << "\n";
//    std::random_shuffle(rank_sorted_all.begin(), rank_sorted_all.begin() + num_ties, myrandom);
//    for (auto rs : rank_sorted_all) { std::cerr << rs << ' ';};
//    std:: cerr << "\n";
//    return (rank_sorted_all[0]);
//}


int select_best_aln_paired_end(const std::vector<bam1_t*>& aln1s,
                               const std::vector<bam1_t*>& aln2s,
                               const int merge_pe_mode = MERGE_PE_SUM){
    std::vector<int> mapqs, scores;
    // Indicator 
    //      1 if a pair is properly paired (TLEN != 0). We apply a loose threshold here.
    //      0 otherwise
    std::vector<bool> pair_indicators;
    for (int i = 0; i < aln1s.size(); i++){
        bam1_core_t c_aln1 = aln1s[i]->core, c_aln2 = aln2s[i]->core;
        pair_indicators.push_back(c_aln1.isize != 0);
        // if (c_aln1.size == 0){
        //     // Set MAPQ and AS to 0 for an unaligned read.
        //     maps.push_back(0);
        //     scores.push_back(0);
        // } else 
        if (merge_pe_mode == MERGE_PE_SUM){
            // MERGE_PE_SUM mode sums MAPQ and AS. AS is set to 0 for an unaligned read.
            mapqs.push_back(c_aln1.qual + c_aln2.qual);
            int score = 0;
            if (!(c_aln1.flag & BAM_FUNMAP))
                score += bam_aux2i(bam_aux_get(aln1s[i], "AS"));
            if (!(c_aln2.flag & BAM_FUNMAP))
                score += bam_aux2i(bam_aux_get(aln2s[i], "AS"));
            scores.push_back(score);
        } else if (merge_pe_mode == MERGE_PE_MAX){
            // MERGE_PE_MAX mode takes max MAPQ and AS.
            if (c_aln1.qual > c_aln2.qual)
                mapqs.push_back(c_aln1.qual);
            else
                mapqs.push_back(c_aln2.qual);
            int score = INT_MIN;
            if (!(c_aln1.flag & BAM_FUNMAP))
                score = bam_aux2i(bam_aux_get(aln1s[i], "AS"));
            if (!(c_aln2.flag & BAM_FUNMAP))
                score = (score > bam_aux2i(bam_aux_get(aln2s[i], "AS")))?
                    score : bam_aux2i(bam_aux_get(aln2s[i], "AS"));
            scores.push_back(score);
        } else{
            std::cerr << "[Error] Invalid merging mode for paired-end alignments "
                      << merge_pe_mode << "\n";
        }
        std::cerr << bam_get_qname(aln1s[i]) << " " << scores[i] << " " << mapqs[i] << " " << pair_indicators[i] << "\n";
    }
    int best_index = rank_alns(pair_indicators=pair_indicators, scores=scores, mapqs=mapqs);
    return (best_index);
}


void merge_sam(merge_sam_opts args){
    std::srand(args.rand_seed);
    std::vector<std::string> sam_fns = read_file_as_vector(args.sam_list);
    std::vector<std::string> ids = read_file_as_vector(args.id_list);

    // Read each SAM file listed in `--sam_list`.
    std::vector<samFile*> sam_fps;
    std::vector<bam_hdr_t*> hdrs;
    std::vector<bam1_t*> aln1s, aln2s;
    for (int i = 0; i < sam_fns.size(); i++){
        sam_fps.push_back(sam_open(sam_fns[i].data(), "r"));
        hdrs.push_back(sam_hdr_read(sam_fps[i]));
        aln1s.push_back(bam_init1());
        aln2s.push_back(bam_init1());
    }
    while(true){
        // Paired-end mode.
        // Read two reads from each of the SAM files in each iteration.
        for (int i = 0; i < sam_fns.size(); i++){
            if (sam_read1(sam_fps[i], hdrs[i], aln1s[i]) < 0 ||
                sam_read1(sam_fps[i], hdrs[i], aln2s[i]) < 0)
                break;
            // Check read names: they should be identical.
            if (strcmp(bam_get_qname(aln1s[i]), bam_get_qname(aln2s[i])) != 0){
                std::cerr << "[Error] SAM file should be sorted by read name.\n";
                std::cerr << "Mismatched reads: " << bam_get_qname(aln1s[i]) <<
                             " and " << bam_get_qname(aln2s[i]) << "\n";
                exit(1);
            }
            // Check if read names from all files are identical.
            if (i > 0)
                if (strcmp(bam_get_qname(aln1s[0]), bam_get_qname(aln1s[i])) != 0){
                    std::cerr << "[Error] Reads mismatch across SAM files.\n";
                    std::cerr << "Mismatched reads: " << bam_get_qname(aln1s[0]) <<
                                 " and " << bam_get_qname(aln1s[i]) << "\n";
                    exit(1);
                }
        }
        select_best_aln_paired_end(aln1s, aln2s, MERGE_PE_SUM);
    }
}


void merge_sam_main(int argc, char** argv){
    int c;
    merge_sam_opts args;
    args.cmd = make_cmd(argc, argv);
    static struct option long_options[]{
        {"sam_list", required_argument, 0, 's'},
        {"id_list", required_argument, 0, 'i'},
        {"output_prefix", required_argument, 0, 'o'},
        {"rand_seed", required_argument, 0, 'r'},
        {"decoy_list", required_argument, 0, 'd'},
        {"decoy_threshold", required_argument, 0, 'p'},
    };
    int long_idx = 0;
    while((c = getopt_long(argc, argv, "d:i:o:p:r:s:", long_options, &long_idx)) != -1){
        switch (c){
            case 'd':
                args.decoy_list = optarg;
                break;
            case 'i':
                args.id_list = optarg;
                break;
            case 'o':
                args.output_prefix = optarg;
                break;
            case 'p':
                args.decoy_threshold = atoi(optarg);
                break;
            case 'r':
                args.rand_seed = atoi(optarg);
                break;
            case 's':
                args.sam_list = optarg;
                break;
            default:
                std::cerr << "Ignoring option " << c << " \n";
                exit(1);
        }
    }
    merge_sam(args);
}

int main(int argc, char** argv) {
    if (!strcmp(argv[optind], "split_sam")){
        std::cerr << "[split_sam] Split a SAM file into high- and low-quality sub SAM files and " <<
            "generate FASTQ files for low-quality reads.\n";
        split_sam_main(argc, argv);
        std::cerr << "[split_sam] Completed.\n";
    } else if (!strcmp(argv[optind], "merge")){
        // TODO
        merge_sam_main(argc, argv);
    }
    return 0;
}
