#include <iostream>

#include <getopt.h>
#include <htslib/sam.h>

#include <add_aux.hpp>
#include <refflow_utils.hpp>


/* Parse a raw AUX tag into right htslib format.
 *
 * Adapted from Vincenzo Pii
 * https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
 */
void parse_aux_tag(const std::string &aux, char (&aux_tag)[2], char &aux_type,
                   size_t &aux_len, std::string &aux_content) {
    std::string token, s = aux, delimiter = ":";
    size_t pos = 0, tag_idx = 0;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        if (tag_idx == 0) {
            if (token.length() != 2) {
                std::cerr << "[Error] Aux tag name should have exactly two characters.\n";
                std::cerr << "Given tag " << token << " doesn't satisfy\n";
                exit(1);
            }
            aux_tag[0] = token[0];
            aux_tag[1] = token[1];
        } else if (tag_idx == 1) {
            if (token.length() != 1) {
                std::cerr << "[Error] Aux tag type should have exactly one character.\n";
                std::cerr << "Given aux type " << token << " doesn't satisfy\n";
                exit(1);
            }
            aux_type = token[0];
        } else {
            std::cerr << "[Error] Aux tag should have exactly three fields.\n";
            std::cerr << "Given tag " << aux_tag << " doesn't satisfy\n";
            exit(1);
        }
        s.erase(0, pos + delimiter.length());
        tag_idx ++;
    }
    aux_content = s;
    aux_len = aux_content.length() + 1;
}

void add_aux(add_aux_opts args) {
    // Read raw SAM file.
    samFile* sam_fp = sam_open(args.sam_fn.data(), "r");
    bam_hdr_t* hdr = sam_hdr_read(sam_fp);

    char const *out_mode = (args.output_ext == "bam")? "wb" : "w";
    std::string out_sam_fn = (args.output_prefix == "-")?
        "-" : args.output_prefix + "." + args.output_ext;
    samFile* out_sam_fp = sam_open(out_sam_fn.data(), out_mode);
    int write_hdr = sam_hdr_write(out_sam_fp, hdr);

    char aux_tag[2];
    char aux_type;
    std::string aux_content;
    size_t aux_len;
    parse_aux_tag(args.added_tag, aux_tag, aux_type, aux_len, aux_content);

    bam1_t* aln = bam_init1();
    while(true){
        if (sam_read1(sam_fp, hdr, aln) < 0) break;

        uint8_t* aux_data = (uint8_t*) aux_content.data();
        // If AUX type is 'i' (int)
        if (aux_type == 'i') {
            aux_len = 4;
            int32_t int_aux_content = std::stoi(aux_content);
            aux_data = (uint8_t*) &int_aux_content;
        }
        if (bam_aux_append(aln, aux_tag, aux_type, aux_len, aux_data) < 0) {
            std::cerr << "[Error] Failed to add aux tag for " << bam_get_qname(aln) << "\n";
            exit(1);
        }

        if (sam_write1(out_sam_fp, hdr, aln) < 0) {
            std::cerr << "[Error] Failure when writing an alignment.\n";
            std::cerr << "Read name: " << bam_get_qname(aln) << "\n";
            exit(1);
        }
    }
    bam_destroy1(aln);
    sam_close(sam_fp);
    sam_close(out_sam_fp);
}


void add_aux_main(int argc, char** argv) {
    int c;
    add_aux_opts args;
    args.cmd = make_cmd(argc, argv);
    static struct option long_options[]{
        {"sam_fn", required_argument, 0, 's'},
        {"output_prefix", required_argument, 0, 'o'},
        {"output_ext", required_argument, 0, 'O'},
        {"added_tag", required_argument, 0, 'g'},
    };
    int long_idx = 0;
    while ((c = getopt_long(argc, argv, "hg:o:O:s:", long_options, &long_idx)) != -1){
        switch (c){
            case 'o':
                args.output_prefix = optarg;
                break;
            case 'O':
                args.output_ext = optarg;
                break;
            case 'g':
                args.added_tag = optarg;
                break;
            case 's':
                args.sam_fn = optarg;
                break;
            case 'h':
                add_aux_help();
                exit(1);
            default:
                std::cerr << "Ignoring option " << c << " \n";
                add_aux_help();
                exit(1);
        }
    }
    if (args.output_ext != "sam" && args.output_ext != "bam") {
        std::cerr << "[Error] Unsupported output extension " << args.output_ext << "\n";
        add_aux_help();
        exit(1);
    }
    if (args.added_tag == "") {
        std::cerr << "[Error] Empty added tag. Set it using the `-g` option\n";
        add_aux_help();
        exit(1);
    }
    add_aux(args);
}

