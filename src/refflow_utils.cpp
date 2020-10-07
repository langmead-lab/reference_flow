#include <algorithm>
#include <fstream>
#include <getopt.h>
#include <htslib/sam.h>
#include <iostream>
#include <string.h>


const int SPLIT_OPT = 0;
const int SPLIT_PES = 1;
int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

struct refflow_utils_opts{
    std::string cmd = "";
    std::string sam_fn = "";
    std::string output_prefix = "";
    int split_strategy = SPLIT_OPT;
    int mapq_threshold = 0;
};


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
void fetch_low_mapq(refflow_utils_opts args){
    samFile* sam_fp = (args.sam_fn == "")?
        sam_open("-", "r") :
        sam_open(args.sam_fn.data(), "r");

    std::cerr << args.sam_fn.data() << "\n";

    bam_hdr_t* hdr = sam_hdr_read(sam_fp);
    // Read two reads in one iteration to for the paired-end mode.
    bam1_t* aln1 = bam_init1(), * aln2 = bam_init1();
    std::ofstream out_fq1, out_fq2;
    out_fq1.open(args.output_prefix + "-R1.fq");
    out_fq2.open(args.output_prefix + "-R2.fq");
    while(true){
        if (sam_read1(sam_fp, hdr, aln1) < 0 || sam_read1(sam_fp, hdr, aln2) < 0)
            break;
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
        }
    }
    sam_close(sam_fp);
    out_fq1.close();
    out_fq2.close();
}


std::string make_cmd(int argc, char** argv) {
    std::string cmd("");
    for (auto i = 0; i < argc; ++i) {
        cmd += std::string(argv[i]) + " ";
    }
    return cmd;
}


int main(int argc, char** argv) {
    int c;
    refflow_utils_opts args;
    args.cmd = make_cmd(argc, argv);
    static struct option long_options[]{
        {"sam_fn", required_argument, 0, 's'},
        {"output_prefix", required_argument, 0, 'o'}
    };
    int long_idx = 0;
    while((c = getopt_long(argc, argv, "ho:p:q:s:", long_options, &long_idx)) != -1){
        switch (c){
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
    if (!strcmp(argv[optind], "fetch_low_mapq")){
        std::cerr << "fetch_low_mapq\n";
        fetch_low_mapq(args);
        // TODO
    } else if (!strcmp(argv[optind], "merge")){
        // TODO
        std::cerr << "Merge is not yet supported.\n";
    }
    return 0;
}
