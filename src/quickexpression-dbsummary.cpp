#include "dbsummary.hpp"
#include "filewriter.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "dbcommon.hpp"

int main(int argc, char *argv[])
{
    bool showHelp = false;
    bool showDetail = false;
    int exitStatus = 0;
    int ch;
    const char *program = argv[0];
    const char *outputpath = NULL;

    while ((ch = getopt(argc, argv, "ho:d")) != -1) {
        switch (ch) {
        case 'h':
            showHelp = true;
            break;
        case 'd':
            showDetail = true;
            break;
        case 'o':
            outputpath = optarg;
            break;
        default:
            showHelp = true;
            exitStatus = 1;
            break;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc != 1) {
        fprintf(stderr, "arguments are missing\n");
        showHelp = true;
        exitStatus = 1;
    }

    if (showHelp) {
        FILE *out = exitStatus ? stderr : stdout;
        fprintf(out, "%s [options] INDEXDB\n", program);
        fprintf(out,
                "Options:\n"
                "        -h : show this help\n"
                "        -d : show detail\n"
                "  -o OUTPUT: summary result output\n");
        return exitStatus;
    }

    std::map<std::string, int> summary;
    int ngram;
    int fragment;

    GenomeHashWithMetadata map;
    FILE *file = fopen(argv[0], "r");
    if (file == NULL) {
        perror("Cannot open database");
        return 1;
    }
    
    bool success = map.load(file);
    if (!success) {
        fprintf(stderr, "Cannot open database\n");
        return 1;
    }
    
    success = dbsummary(argv[0], &summary, &fragment, &ngram);
    if (!success) {
        fprintf(stderr, "Cannot summarize database\n");
        return 1;
    }
    
    cppfasta::BasicWriter *output = NULL;
    if (output) {
        output = cppfasta::openFileForWrite(outputpath);
    } else {
        output = new cppfasta::FileWriter(stdout);
    }

    if (!output->printf("        N-Gram: %d\n", ngram))
        goto onioerror;
    if (!output->printf("# of fragments: %d\n", map.numberOfEntries()))
        goto onioerror;
    if (!output->printf("    # of genes: %d\n", summary.size()))
        goto onioerror;

    if (showDetail) {
        if (!output->printf("\n"))
            goto onioerror;
        for (std::map<std::string, int>::const_iterator it = summary.begin();
             it != summary.end(); ++it) {
            if (!output->printf("%s\t%d\n", it->first.c_str(), it->second))
                goto onioerror;
        }
    }
    delete output;
        
    return 0;
    
onioerror:
    delete output;
    fprintf(stderr, "Cannot write expression result\n");
    return 1;
}

