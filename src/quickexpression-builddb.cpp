#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "config.h"
#include "dbbuilder.hpp"
#include "quickcommon.hpp"

#include <unistd.h>
#include <stdbool.h>

static const int DEFAULT_NGRAM = 32;

int main(int argc, char *argv[])
{
    bool showHelp = false;
    bool onMemory = false;
    int exitStatus = 0;
    int ngram = DEFAULT_NGRAM;
    int ch;
    const char *coveragepath = NULL;
    const char *program = argv[0];

    while ((ch = getopt(argc, argv, "mhn:c:f:")) != -1) {
        switch (ch) {
        case 'h':
            showHelp = true;
            break;
        case 'm':
            onMemory = true;
            break;
        case 'c':
            coveragepath = optarg;
            break;
        case 'n': {
            if (*optarg == '\0') {
                fprintf(stderr, "a number is required for '-n'\n");
                showHelp = true;
                exitStatus = 1;
                break;
            }
            
            char *endptr;
            long tmp = strtol(optarg, &endptr, 10);
            if (*endptr != '\0') {
                fprintf(stderr, "a number is required for '-n'\n");
                showHelp = true;
                exitStatus = 1;
                break;
            }

            ngram = tmp;
            break;
        }
        default:
            showHelp = true;
            exitStatus = 1;
            break;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc != 3) {
        fprintf(stderr, "arguments are missing\n");
        showHelp = true;
        exitStatus = 1;
    }

    if (showHelp) {
        FILE *out = exitStatus ? stderr : stdout;
        fprintf(out, "%s [-h] GENE2REFSEQ REFSEQ_FASTA DBFILE\n", program);
        fprintf(out,
                "Options:\n"
                "             -h : show this help\n"
                "             -m : do not use file memory mapping\n"
                "     -c COVERAGE: output coverage file\n"
                "        -n NGRAM: set n-gram to build index (default: %d)\n",
                DEFAULT_NGRAM);
        return exitStatus;
    }

    bool success = buildDatabase(argv[2], coveragepath, argv[1], argv[0], ngram, onMemory);
    if (success)
        return 0;
    fprintf(stderr, "Failed to build index\n");
    return 1;
}
