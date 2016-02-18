#include <stdio.h>
#include <cppfasta.hpp>
#include <unistd.h>
#include <stdbool.h>

#include "quickcommon.hpp"
#include "gene2refseq.hpp"

int main(int argc, char *argv[])
{
    bool showHelp = false;
    int exitStatus = 0;
    int ch;
    const char *program = argv[0];

    const char *output_fasta = NULL;
    const char *taxid = NULL;

    while ((ch = getopt(argc, argv, "o:t:h")) != -1) {
        switch (ch) {
        case 'h':
            showHelp = true;
            break;
        case 'o':
            output_fasta = optarg;
            break;
        case 't':
            taxid = optarg;
            break;
        default:
            showHelp = true;
            exitStatus = 1;
            break;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc < 2) {
        fprintf(stderr, "GENE2REFSEQ or/and REFSEQ_FATA are missing\n");
        showHelp = true;
        exitStatus = 1;
    }

    if (output_fasta == NULL) {
        fprintf(stderr, "OUTPUT_FASTA is missing\n");
        showHelp = true;
        exitStatus = 1;
    }


    if (taxid == NULL) {
        fprintf(stderr, "TAXID is missing\n");
        showHelp = true;
        exitStatus = 1;
    }

    if (showHelp) {
        FILE *out = exitStatus ? stderr : stdout;
        fprintf(out, "%s [-h] -o OUTPUT_FASTA -t TAXID GENE2REFSEQ REFSEQ_FASTA ...\n", program);
        fprintf(out,
                "Options:\n"
                "              -h : show this help\n"
                " -o OUTPUT_FASTA : filtered fasta\n"
                "        -t TAXID : Taxonomy ID\n"
                );
        return exitStatus;
    }

    Gene2Refseq gene2refseq;
    bool ok = gene2refseq.load(argv[0]);
    if (!ok) {fprintf(stderr, "Failed to load gene2refseq\n"); return 1;}

    cppfasta::FastaWriter writer;
    ok = writer.open(output_fasta);
    if (!ok) {fprintf(stderr, "Failed to open %s\n", output_fasta); return 1;}

    for (int i = 1; i < argc; i++) {
        fprintf(stderr, "Loading %s\n", argv[i]);
        cppfasta::FastaReader reader;
        ok = reader.open(argv[i]);
        if (!ok) {fprintf(stderr, "Failed to open %s\n", argv[i]); return 1;}
        cppfasta::SequenceRecord *record;
        while ((record = reader.nextRecord()) != NULL) {
            std::string refseqid = stringsplit(record->name(), "|")[3];
            if (taxid == gene2refseq.refseq2tax(refseqid)) {
                writer.write(*record);
            }
        }
    }

    return 0;
}
