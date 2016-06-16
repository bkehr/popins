#include <iostream>
#include <ctime>
#include <stdio.h>
#include <string.h>

#include "command_line_parsing.h"
#include "assemble/popins_assemble.h"
#include "merge/popins_merge.h"
#include "contigmap/popins_contigmap.h"
#include "place/popins_place.h"
#include "genotype/popins_genotype.h"

// ==========================================================================

void printHelp(char const * name)
{
    std::cerr << "PopIns - population-scale detection of novel sequence insertions" << std::endl;
    std::cerr << "================================================================" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mSYNOPSIS\033[0m" << std::endl;
    std::cerr << "    \033[1m" << name << " COMMAND\033[0m [\033[4mOPTIONS\033[0m]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mCOMMAND\033[0m" << std::endl;
    std::cerr << "    \033[1massemble\033[0m   Crop unmapped reads from a bam file and assemble them." << std::endl;
    std::cerr << "    \033[1mmerge\033[0m      Merge contigs from assemblies of unmapped reads into supercontigs." << std::endl;
    std::cerr << "    \033[1mcontigmap\033[0m  Map unmapped reads to (super-)contigs." << std::endl;
    std::cerr << "    \033[1mplace\033[0m      Find position of (super-)contigs in the reference genome." << std::endl;
    std::cerr << "    \033[1mgenotype\033[0m   Genotype insertions for an individual." << std::endl;
    std::cerr << std::endl;
    std::cerr << "\033[1mVERSION\033[0m" << std::endl;
    std::cerr << "    " << VERSION << ", Date: " << VERSION_DATE << std::endl;
    std::cerr << std::endl;
    std::cerr << "Try `" << name << " COMMAND --help' for more information on each command." << std::endl;
    std::cerr << std::endl;
}

// ==========================================================================
// Function main()
// ==========================================================================

int main(int argc, char const ** argv)
{
    std::time_t start_time = std::time(0);

    int ret = 0;
    const char * prog_name = argv[0];
    if (argc < 2)
    {
        printHelp(prog_name);
        return 1;
    }

    const char * command = argv[1];
    if (strcmp(command,"assemble") == 0) ret = popins_assemble(argc, argv);
    else if (strcmp(command,"merge") == 0) ret = popins_merge(argc, argv);
    else if (strcmp(command,"contigmap") == 0) ret = popins_contigmap(argc, argv);
    else if (strcmp(command,"place") == 0) ret = popins_place(argc, argv);
    else if (strcmp(command,"genotype") == 0) ret = popins_genotype(argc, argv);
    else if (strcmp(command, "--help") == 0 || strcmp(command, "-h") == 0)
    {
        printHelp(prog_name);
        return 1;
    }
    else
    {
        std::cerr << "ERROR: Unknown command: " << command << std::endl;
        printHelp(prog_name);
        return 1;
    }

    if (ret == 0)
    {
        std::ostringstream msg;
        msg << "popins " << command << " finished in " << (std::time(0) - start_time) << " seconds.";
        printStatus(msg);
    }
    return ret;
}

