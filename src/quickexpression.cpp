#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include "mapping.hpp"
#include "filewriter.hpp"
#include "quickcommon.hpp"
#include "pipedopen.h"

int main(int argc, char *argv[])
{
    bool showHelp = false;
    bool checkAll = false;
    bool unifiedPairFile = false;
    int exitStatus = 0;
    int ch;
    long stepSize = 1;
    long acceptcount = 1;
    const char *program = argv[0];
    const char *outputpath = NULL;
    const char *resultpath = NULL;

    while ((ch = getopt(argc, argv, "ho:r:s:a:cU")) != -1) {
        switch (ch) {
        case 'h':
            showHelp = true;
            break;
        case 'c':
            checkAll = true;
            break;
        case 'U':
            unifiedPairFile = true;
            break;
        case 'o':
            outputpath = optarg;
            break;
        case 'r':
            resultpath = optarg;
            break;
        case 's': {
            bool ok;
            stepSize = string2long(optarg, &ok);
            if (!ok) {
                exitStatus = 1;
                showHelp = true;
            }
            break;
        }
        case 'a': {
            bool ok;
            acceptcount = string2long(optarg, &ok);
            if (!ok) {
                exitStatus = 1;
                showHelp = true;
            }
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

    if (argc != 2 && argc != 3) {
        fprintf(stderr, "arguments are missing\n");
        showHelp = true;
        exitStatus = 1;
    }

    if (stepSize < 1) {
        fprintf(stderr, "Step Size should be 1 or larger\n");
        showHelp = true;
        exitStatus = 1;
    }

    if (showHelp) {
        FILE *out = exitStatus ? stderr : stdout;
        fprintf(out, "%s [options] INDEXDB FASTQ [FASTQ]\n", program);
        fprintf(out,
                "Options:\n"
                "        -h : show this help\n"
                "        -m : do not use file mapped memory\n"
                "        -U : mate pair in one file\n"
                "-s STEPSIZE: Step size\n"
                "   -a COUNT: Accept count\n"
                "  -o OUTPUT: Expression result output\n"
                "  -r RESULT: Mapping result for each fragment\n");
        return exitStatus;
    }

    cppfasta::BasicWriter *result = NULL;
    if (resultpath) {
        result = cppfasta::openFileForWrite(resultpath);
    }

    std::map<std::string, int> expression;
    std::map<std::string, double> fpkm;
    uint64_t processedSequecnes;

    cppfasta::FastqReader2 *reader;
    
    PFILE *file = pipedopen(argv[1], false);
    if (file == NULL) {
        perror("cannot open FASTQ\n");
        return false;
    }

    bool pairedend = false;
    PFILE *file2 = NULL;
    if (argc == 3) {
        file2 = pipedopen(argv[2], false);
        if (file2 == NULL) {
            perror("cannot open FASTQ\n");
            return false;
        }
        reader = new cppfasta::FastqReader2(file->file, file2->file);
        pairedend = true;
    } else {
        reader = new cppfasta::FastqReader2(file->file);
        pairedend = unifiedPairFile;
    }
    
    
    bool success = domapping(argv[0], reader, &expression, &fpkm, result, &processedSequecnes, stepSize, acceptcount, checkAll, pairedend);

    pipedclose(file);
    if (file2)
        pipedclose(file2);
    delete result;
    delete reader;
    
    if (success) {
        cppfasta::BasicWriter *output = NULL;
        if (outputpath) {
            output = cppfasta::openFileForWrite(outputpath);
        } else {
            output = new cppfasta::FileWriter(stdout);
        }
        if (output == NULL) {
            perror("Cannot open result file");
            return 1;
        }
        
        output->write("GeneID\tMapped\tFPKM\n");

        for (std::map<std::string, int>::const_iterator it = expression.begin();
            it != expression.end(); ++it) {
            if (!output->printf("%s\t%d\t%lf\n", it->first.c_str(), it->second, fpkm[it->first])) {
                fprintf(stderr, "Cannot write expression result\n");
                return 1;
            }
        }
        delete output;

        fprintf(stderr, "Mapping rate: %.2f%%\n", 100*(1-static_cast<double>(expression["UNMAPPED"])/processedSequecnes));
        
        return 0;
    }
    
    fprintf(stderr, "Failed to mapping\n");
    return 1;
}
