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


#endif /* REFFLOW_UTILS_H__ */
