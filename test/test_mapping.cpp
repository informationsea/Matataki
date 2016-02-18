#include "gtest/gtest.h"
#include "quickcommon.hpp"
#include "dbbuilder.hpp"
#include "mapping.hpp"
#include "tabtable.hpp"
#include "filewriter.hpp"
#include "dbcommon.hpp"
#include "dbsummary.hpp"
#include <stdio.h>
#include <stdlib.h>
#include "pipedopen.h"

TEST(SKIPARRAY, SKIPARRAY) {
    uint8_t steplen[CPPFASTA_MAXIMUM_SEQUENCE_LENGTH];

    bzero(steplen, sizeof(steplen));
    QuickMapping::createSkipArray(steplen, "ATCGNGGA", 8, 4);
    uint8_t expected1[] = {0, 4, 3, 2, 1, 0, 0, 0};
    for (size_t i = 0; i < sizeof(expected1); i++) fprintf(stderr, "%u ", steplen[i]);
    fprintf(stderr, "\n");
    ASSERT_TRUE(memcmp(expected1, steplen, sizeof(expected1)) == 0);

    
    bzero(steplen, sizeof(steplen));
    QuickMapping::createSkipArray(steplen, "ATCGNGNA", 8, 4);
    for (size_t i = 0; i < sizeof(expected1); i++) fprintf(stderr, "%u ", steplen[i]);
    fprintf(stderr, "\n");
    uint8_t expected2[] = {0, 6, 5, 4, 3, 2, 1, 0};
    ASSERT_TRUE(memcmp(expected2, steplen, sizeof(expected2)) == 0);


    bzero(steplen, sizeof(steplen));
    QuickMapping::createSkipArray(steplen, "NNNNNNNN", 8, 4);
    for (size_t i = 0; i < sizeof(expected1); i++) fprintf(stderr, "%u ", steplen[i]);
    fprintf(stderr, "\n");
    uint8_t expected3[] = {8, 7, 6, 5, 4, 3, 2, 1};
    ASSERT_TRUE(memcmp(expected3, steplen, sizeof(expected3)) == 0);
    
}

TEST(MAPPING, MAPPING) {
    mkdir("tmp", 0700);
    TemporaryFile tmpfile("dbbuilder", ".dat", "./tmp");

    FILE *filew = fopen(tmpfile.path(), "w");
    GenomeHash genomeHash;
    genomeHash.open(4, 4);
    genomeHash.insert("ATCG", 1);
    genomeHash.insert("GCTA", 2);
    genomeHash.save(filew);
    fprintf(filew, "1\tX\t1.0\n");
    fprintf(filew, "2\tY\t2.0\n");

    ASSERT_EQ(1U, genomeHash.get("ATCG"));
    ASSERT_EQ(2U, genomeHash.get("GCTA"));
    ASSERT_EQ(1U, genomeHash.get("ATCGGA"));
    
    fclose(filew);

    QuickMapping mapping;
    ASSERT_TRUE(mapping.loadIndex(tmpfile.path()));
    ASSERT_EQ(1U, mapping.addSequence(cppfasta::SequenceRecord2("A", "ATCGGA")));
    ASSERT_EQ(1U, mapping.addSequence(cppfasta::SequenceRecord2("B", "GATCGG")));
    ASSERT_EQ(2U, mapping.addSequence(cppfasta::SequenceRecord2("C", "GCTAAG")));
    ASSERT_EQ(2U, mapping.addSequence(cppfasta::SequenceRecord2("D", "AGGCTA")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("E", "AGGCTT")));
    ASSERT_EQ(2U, mapping.addSequence(cppfasta::SequenceRecord2("F", "GNGCTA")));
    ASSERT_EQ(2U, mapping.addSequence(cppfasta::SequenceRecord2("G", "NNGCTA")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("H", "NTCGGA")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("I", "ANTCGG")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("J", "AXTCGC")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("K", "GCTXTT")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("L", "GCTNAA")));
    ASSERT_EQ(2U, mapping.addSequence(cppfasta::SequenceRecord2("M", "GCTANN")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("N", "GCNNAT")));

    std::map<std::string, int> result = mapping.result();

    for (std::map<std::string, int>::const_iterator it = result.begin();
         it != result.end(); ++it) {
        DEBUG("%s %d", it->first.c_str(), it->second);
    }

    ASSERT_EQ(2, result["X"]);
    ASSERT_EQ(5, result["Y"]);
    ASSERT_EQ(7, result["UNMAPPED"]);
    int total = result["X"] + result["Y"];

    std::map<std::string, double> fpkm = mapping.fpkm();
    ASSERT_EQ(1U, fpkm.count("X"));
    ASSERT_DOUBLE_EQ(((double)result["X"])/total*1000000000, fpkm["X"]);
    ASSERT_EQ(1U, fpkm.count("Y"));
    ASSERT_DOUBLE_EQ(((double)result["Y"])/2/total*1000000000, fpkm["Y"]);
    rmdir("tmp");
}

TEST(MAPPING, MAPPING2) {
    mkdir("tmp", 0700);
    TemporaryFile tmpfile("dbbuilder", ".dat", "./tmp");

    FILE *filew = fopen(tmpfile.path(), "w");
    GenomeHash genomeHash;
    genomeHash.open(4, 4);
    genomeHash.insert("ATCG", 1);
    genomeHash.insert("GCTA", 2);
    genomeHash.save(filew);
    fprintf(filew, "1\tX\t1.0\n");
    fprintf(filew, "2\tY\t2.0\n");
    fclose(filew);

    QuickMapping mapping(4, 2);
    ASSERT_TRUE(mapping.loadIndex(tmpfile.path()));
    ASSERT_EQ(1U, mapping.addSequence(cppfasta::SequenceRecord2("A", "ATCGATCGTTTT")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("B", "ATCGATCATTTT")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("C", "TATCGATCGTTT")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("D", "TTATCGATCGTT")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("E", "TTTATCGATCGT")));
    ASSERT_EQ(1U, mapping.addSequence(cppfasta::SequenceRecord2("F", "TTTTATCGATCG")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("G", "GCTNTTTTGCTA")));
    ASSERT_EQ(2U, mapping.addSequence(cppfasta::SequenceRecord2("H", "GCTATTTTGCTA")));
    ASSERT_EQ(2U, mapping.addSequence(cppfasta::SequenceRecord2("I", "GCTAGCTAGCTA")));
    ASSERT_EQ(0U, mapping.addSequence(cppfasta::SequenceRecord2("J", "NGCTAGCTAGCT")));
    ASSERT_EQ(2U, mapping.addSequence(cppfasta::SequenceRecord2("K", "NNNGGCTAGCTA")));
    rmdir("tmp");
}

TEST(MAPPING, MAPPING3) {
    mkdir("tmp", 0700);
    TemporaryFile file("mapping2", ".dat", "./tmp");
    ASSERT_TRUE(buildDatabase(file.path(), NULL, "./data/refseq.rna.bz2", "./data/gene2refseq_test.txt.bz2", 10, false));

    //frozenhashmap::FrozenMap map;
    //ASSERT_TRUE(map.open(file.path()));

    //mapping
    fprintf(stderr, "Mapping...\n");
    cppfasta::FileWriter fileWriter;
    ASSERT_TRUE(fileWriter.open("tmp/mapping-result.txt"));
    std::map<std::string, int> result;
    std::map<std::string, double> fpkm;
    uint64_t num;

    PFILE *pfile = pipedopen("./data/simulation.fastq.bz2", false);
    ASSERT_TRUE(pfile);
    cppfasta::FastqReader2 reader(pfile->file);
    
    ASSERT_TRUE(domapping(file.path(), &reader, &result, &fpkm, &fileWriter, &num, 1, 1, false, false));
    ASSERT_EQ(493150UL, num);

    pipedclose(pfile);

    for (std::map<std::string, int>::iterator it = result.begin();
         it != result.end(); ++it) {
        fprintf(stderr, "%10s : %d\n", it->first.c_str(), it->second);
    }
    fprintf(stderr, "OK...\n");

    TabTableReader answer_table;
    std::map<std::string, int> answer;
    ASSERT_TRUE(answer_table.open("./data/python-mapping-result.txt"));
    answer_table.next(); // skip header
    std::vector<std::string> row;
    while ((row = answer_table.next()).size() >= 2) {
        char *endptr;
        long num = strtol(row[1].c_str(), &endptr, 10);
        ASSERT_EQ('\0', *endptr);
        ASSERT_EQ(num, result[row[0]]) << row[0];
        //DEBUG("%s %ld %ld", row[0].c_str(), num, result[row[0]]);

        if (row[0] != "UNMAPPED" || row[0] != "MULTIMAPPED") {
            char *endptr;
            double fpkm_expected = strtod(row[2].c_str(), &endptr);
            ASSERT_EQ('\0', *endptr);
            //fprintf(stderr, "FPKM Expected[%s]: %lf  Result: %lf\n", row[0].c_str(), fpkm_expected, fpkm[row[0]]);
            ASSERT_NEAR(fpkm_expected, fpkm[row[0]], 1e-19) << "GeneID " << row[0];
        }
    }
}


TEST(MAPPING, MAPPING4) {
    mkdir("tmp", 0700);
    TemporaryFile tmpfile("dbbuilder4", ".dat", "./tmp");

    FILE *filew = fopen(tmpfile.path(), "w");
    GenomeHash genomeHash;
    genomeHash.open(4, 4);
    genomeHash.insert("ATCG", 1);
    genomeHash.insert("GCTA", 2);
    genomeHash.save(filew);
    fprintf(filew, "1\tX\t1.0\n");
    fprintf(filew, "2\tY\t2.0\n");
    fclose(filew);

    {
        QuickMapping mapping(1, 1);
        ASSERT_TRUE(mapping.loadIndex(tmpfile.path()));
        ASSERT_EQ(1U, mapping.mapSequence(cppfasta::SequenceRecord2("A", "ATCGATCGTTTT")));
        ASSERT_EQ(1U, mapping.mapSequence(cppfasta::SequenceRecord2("A", "AAAAATCGTTTT")));
        ASSERT_EQ(1U, mapping.addSequencePair(cppfasta::SequenceRecord2("A", "ATCGATCGTTTT"), cppfasta::SequenceRecord2("B", "AAAAATCGTTTT")));
        ASSERT_EQ(1U, mapping.addSequencePair(cppfasta::SequenceRecord2("A", "ATCGATCGTTTT"), cppfasta::SequenceRecord2("B", "AAAAAAAATTTT")));
        ASSERT_EQ(1U, mapping.addSequencePair(cppfasta::SequenceRecord2("A", "TTTTTTTTTTTT"), cppfasta::SequenceRecord2("B", "AAAAATCGTTTT")));
        ASSERT_EQ(0U, mapping.addSequencePair(cppfasta::SequenceRecord2("A", "TTTTTTTTTTTT"), cppfasta::SequenceRecord2("B", "AAAAAAAATTTT")));
        ASSERT_EQ(GENOMEHASH_MULTIMAP_ID, mapping.addSequencePair(cppfasta::SequenceRecord2("A", "ATCGTTTTTTTT"), cppfasta::SequenceRecord2("B", "AAAAAAAAGCTA")));
        ASSERT_EQ(3, mapping.result()["X"]);
    }

    {
        QuickMapping mapping(4, 2);
        ASSERT_TRUE(mapping.loadIndex(tmpfile.path()));
        ASSERT_EQ(GENOMEHASH_MULTIMAP_ID, mapping.mapSequence(cppfasta::SequenceRecord2("A", "ATCGGCTATTTT")));
        ASSERT_EQ(GENOMEHASH_MULTIMAP_ID, mapping.addSequencePair(cppfasta::SequenceRecord2("A", "ATCGGCTATTTT"), cppfasta::SequenceRecord2("A", "ATCGATCGTTTT")));
    }

    
    rmdir("tmp");
}

TEST(MAPPING, MAPPING5) {
    mkdir("tmp", 0700);
    TemporaryFile tmpfile("dbbuilder5", ".dat", "./tmp");

    FILE *filew = fopen(tmpfile.path(), "w");
    GenomeHash genomeHash;
    genomeHash.open(4, 4);
    genomeHash.insert("ATCG", 1);
    genomeHash.insert("GCTA", 2);
    genomeHash.save(filew);
    fprintf(filew, "1\tX\t1.0\n");
    fprintf(filew, "2\tY\t2.0\n");
    fclose(filew);

    {
        PFILE *pfile = pipedopen("./data/paired1.fastq", false);
        ASSERT_TRUE(pfile);
        cppfasta::FastqReader2 reader(pfile->file);

        std::map<std::string, int> result;
        std::map<std::string, double> fpkm;
        uint64_t num;
    
        ASSERT_TRUE(domapping(tmpfile.path(), &reader, &result, &fpkm, NULL, &num, 1, 1, false, false));
        ASSERT_EQ(4U, num);
        ASSERT_EQ(2, result["X"]);
        ASSERT_EQ(2, result["UNMAPPED"]);
        pipedclose(pfile);
    }

    {
        PFILE *pfile1 = pipedopen("./data/paired1.fastq", false);
        ASSERT_TRUE(pfile1);
        PFILE *pfile2 = pipedopen("./data/paired2.fastq", false);
        ASSERT_TRUE(pfile2);
        cppfasta::FastqReader2 reader(pfile1->file, pfile2->file);

        std::map<std::string, int> result;
        std::map<std::string, double> fpkm;
        uint64_t num;
    
        ASSERT_TRUE(domapping(tmpfile.path(), &reader, &result, &fpkm, NULL, &num, 1, 1, false, true));
        ASSERT_EQ(4U, num);
        ASSERT_EQ(3, result["X"]);
        ASSERT_EQ(1, result["UNMAPPED"]);

        pipedclose(pfile1);
        pipedclose(pfile2);
    }

    {
        PFILE *pfile = pipedopen("./data/paired-unified.fastq", false);
        ASSERT_TRUE(pfile);
        cppfasta::FastqReader2 reader(pfile->file);

        std::map<std::string, int> result;
        std::map<std::string, double> fpkm;
        uint64_t num;
    
        ASSERT_TRUE(domapping(tmpfile.path(), &reader, &result, &fpkm, NULL, &num, 1, 1, false, true));
        ASSERT_EQ(4U, num);
        ASSERT_EQ(3, result["X"]);
        ASSERT_EQ(1, result["UNMAPPED"]);

        pipedclose(pfile);
    }

    rmdir("tmp");
}






