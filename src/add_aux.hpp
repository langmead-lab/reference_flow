#ifndef ADD_AUX_H__
#define ADD_AUX_H__

struct add_aux_opts{
    std::string cmd = "";
    std::string added_tag = "";
    std::string sam_fn = "-";
    std::string output_ext = "sam";
    std::string output_prefix = "";
};

static void add_aux_help(){
    std::cerr << "\n";
    std::cerr << "Usage: refflow_utils add_aux [options] -s <input> -o <out> -g <aux>\n";
    std::cerr << "\n";
    std::cerr << "  -s <input>  SAM or BAM file to split\n";
    std::cerr << "  -o <out>    Prefix of the output files\n";
    std::cerr << "  -g <aux>    AUX tag to be added ('i' and 'Z' formats are supported)\n";
    std::cerr << "\n";
    std::cerr << "Options (defaults in parentheses):\n";
    std::cerr << "  -O <text>   Output alignment format, can be sam or bam [sam]\n";
    std::cerr << "\n";
}

void parse_aux_tag(const std::string &aux, char (&aux_tag)[2], char &aux_type,
                   int &aux_len, std::string &aux_content);

void add_aux(add_aux_opts args);

void add_aux_main(int argc, char** argv);


#endif /* ADD_AUX_H__ */
