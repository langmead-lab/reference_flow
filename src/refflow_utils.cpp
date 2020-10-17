#include <chrono>
#include <ctime>
#include <iostream>

#include <getopt.h>

#include <merge_sam.hpp>
#include <refflow_utils.hpp>
#include <split_sam.hpp>


int main(int argc, char** argv) {
    /* Time measuring is borrowed from carlduke
     * https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
     */
    double start_cputime = std::clock();
    auto start_walltime = std::chrono::system_clock::now();
    if (!strcmp(argv[optind], "split")){
        std::cerr <<
            "[split] Split a SAM/BAM file into high- and low-quality sub SAM/BAM files and " <<
            "generate FASTQ files for low-quality reads.\n";
        split_sam_main(argc, argv);
        std::cerr << "[split] Completed.\n";
    } else if (!strcmp(argv[optind], "merge")){
        std::cerr << "[merge] Merge a list of SAM/BAM files that contain the same set of reads " <<
            "in the same order.\n";
        std::cerr << "        Best-ranked alignments will be selected in separate SAM/BAM " <<
            "files, labelled according to a list of IDs.\n";
        merge_sam_main(argc, argv);
        std::cerr << "[merge] Completed.\n";
    } else {
        std::cerr << "Subcommand " << argv[optind] << " not found.\n";
        std::cerr << "\n";
        std::cerr << "Currently supported subcommands:\n";
        std::cerr << " - split: split alignments by MAPQ\n";
        std::cerr << " - merge: select best alignments from a set of plural alignments\n";
        std::cerr << "\n";
        std::cerr << "Please try `refflow_utils <subcommand> -h` to see usages\n";
        std::cerr << "\n";
        exit(0);
    }
    double cpu_duration = (std::clock() - start_cputime) / (double)CLOCKS_PER_SEC;
    std::chrono::duration<double> wall_duration = (std::chrono::system_clock::now() - start_walltime);
    std::cerr << "\n";
    std::cerr << "Finished in " << cpu_duration << " CPU seconds, or " << 
                                   wall_duration.count() << " wall clock seconds\n";

    return 0;
}
