//
// panfreebayes — throwaway CLI wrapper around the extracted core (Milestone 1)
//
// This exists only to exercise panfreebayes::callVariants from the command line
// while the core is being validated. It is NOT the Milestone 3 CLI.
//
//   panfreebayes --ref <ref.fasta> --bam <aln.bam> [--pooled-continuous]
//                [--min-alternate-count N] [--min-alternate-fraction X]
//                [--limit-coverage N] [-- <extra freebayes args...>] > out.vcf
//
// Or pass a full freebayes command line straight through:
//
//   panfreebayes --passthrough -f ref.fa --pooled-continuous aln.bam > out.vcf
//

#include "panfreebayes_core.h"

#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <csignal>

#include "SegfaultHandler.h"

static void usage(const char* prog) {
    std::cerr
        << "usage: " << prog << " --ref <ref.fasta> --bam <aln.bam> [options] > out.vcf\n"
        << "\n"
        << "  --ref FILE                  local reference FASTA (required)\n"
        << "  --bam FILE                  BAM aligned to --ref (required)\n"
        << "  --pooled-continuous\n"
        << "  --min-alternate-count N\n"
        << "  --min-alternate-fraction X\n"
        << "  --limit-coverage N\n"
        << "  --                          pass all following tokens to freebayes verbatim\n"
        << "  --passthrough <args...>      treat the entire remaining command line as a\n"
        << "                              freebayes argv (argv[0] synthesised)\n";
}

int main(int argc, char** argv) {
    signal(SIGSEGV, segfaultHandler);

    if (argc < 2) { usage(argv[0]); return 1; }

    // passthrough mode: everything after --passthrough is a freebayes argv
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--passthrough") == 0) {
            std::vector<std::string> a;
            a.push_back("freebayes");
            for (int j = i + 1; j < argc; ++j) a.push_back(argv[j]);
            return panfreebayes::callVariantsArgv(a, std::cout);
        }
    }

    panfreebayes::Options opt;
    bool sawExtra = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto need = [&](const char* what) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "error: " << what << " requires an argument\n";
                std::exit(1);
            }
            return argv[++i];
        };
        if (sawExtra) {
            opt.extraArgs.push_back(a);
        } else if (a == "--") {
            sawExtra = true;
        } else if (a == "--ref" || a == "-f") {
            opt.fasta = need("--ref");
        } else if (a == "--bam" || a == "-b") {
            opt.bam = need("--bam");
        } else if (a == "--pooled-continuous") {
            opt.pooledContinuous = true;
        } else if (a == "--min-alternate-count" || a == "-C") {
            opt.minAlternateCount = std::stoi(need("--min-alternate-count"));
        } else if (a == "--min-alternate-fraction" || a == "-F") {
            opt.minAlternateFraction = std::stod(need("--min-alternate-fraction"));
        } else if (a == "--limit-coverage") {
            opt.limitCoverage = std::stoi(need("--limit-coverage"));
        } else if (a == "-h" || a == "--help") {
            usage(argv[0]);
            return 0;
        } else {
            std::cerr << "error: unrecognised argument '" << a << "'\n";
            usage(argv[0]);
            return 1;
        }
    }

    if (opt.fasta.empty() || opt.bam.empty()) {
        std::cerr << "error: --ref and --bam are both required\n";
        usage(argv[0]);
        return 1;
    }

    return panfreebayes::callVariants(opt, std::cout);
}
