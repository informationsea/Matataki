#include "gtest/gtest.h"
#include "dbbuilder.hpp"
#include "tabtable.hpp"
#include "gene2refseq.hpp"
#include "quickcommon.hpp"
#include "dbcommon.hpp"
#include <sys/stat.h>
#include "cppfasta.hpp"

TEST(QUICK, BUILDER2) {
    const char *path = "./tmp/dbbuilder2.dat";
    mkdir("tmp", 0700);

    DBBuilder builder(4, true);
    ASSERT_TRUE(builder.isReady());

    // warning: leaking memory
    ASSERT_TRUE(builder.addSequence("A", new cppfasta::SequenceRecord("A1", "ATCGG"))); // CCGAT
    ASSERT_TRUE(builder.addSequence("A", new cppfasta::SequenceRecord("A2", "GATCG"))); // CGATC

    ASSERT_TRUE(builder.addSequence("B", new cppfasta::SequenceRecord("B1", "CGGAT"))); // ATCCG
    ASSERT_TRUE(builder.addSequence("B", new cppfasta::SequenceRecord("B2", "ATCCG"))); // CGGAT

    ASSERT_TRUE(builder.addSequence("C", new cppfasta::SequenceRecord("C1", "GCAGT"))); // ACTGC
    ASSERT_TRUE(builder.addSequence("C", new cppfasta::SequenceRecord("C2", "AGTTA"))); // TAACT

    ASSERT_TRUE(builder.addSequence("D", new cppfasta::SequenceRecord("D1", "TCCGA"))); // TCGGA

    ASSERT_TRUE(builder.build(path));

    FILE *f = fopen(path, "r");
    ASSERT_TRUE(f);
    GenomeHashWithMetadata map;
    ASSERT_TRUE(map.load(f));
    fclose(f);

    ASSERT_EQ(std::string("A"), map.geneIdForIndex(map.get("ATCG")));
    ASSERT_EQ(std::string("A"), map.geneIdForIndex(map.get("CGAT")));
    
    ASSERT_EQ(std::string("B"), map.geneIdForIndex(map.get("GGAT")));
    ASSERT_EQ(std::string("B"), map.geneIdForIndex(map.get("ATCC")));

    ASSERT_EQ(4UL, map.numberOfEntries());

    ASSERT_EQ(1.0, map.foundFragmentsForIndex(map.indexForGeneID("A")));
    ASSERT_EQ(1.0, map.foundFragmentsForIndex(map.indexForGeneID("B")));
    ASSERT_EQ(0.0, map.foundFragmentsForIndex(map.indexForGeneID("C")));
    ASSERT_EQ(0.0, map.foundFragmentsForIndex(map.indexForGeneID("D")));
}

#if 1
TEST(QUICK, BUILDER) {
    const char *path = "./tmp/dbbuilder1.dat";
    mkdir("tmp", 0700);
    //TemporaryFile file("dbbuilder", ".dat", "./tmp");
    ASSERT_TRUE(buildDatabase(path, "./tmp/dbbuilder1-coverage.txt", "./data/refseq.rna.bz2", "./data/gene2refseq_test.txt.bz2", 10, true));

    FILE *f = fopen(path, "r");
    ASSERT_TRUE(f);
    GenomeHashWithMetadata map;
    ASSERT_TRUE(map.load(f));
    fclose(f);
    
    //frozenhashmap::FrozenMap map;
    //ASSERT_TRUE(map.open(path));
    //fprintf(stderr, "%llu\n", map.count());
    
   
    TabTableReader expected;
    ASSERT_TRUE(expected.open("./data/unique-10gram.txt.bz2"));

    size_t numberOfCount = 0; // # of fragment
    do {
        std::vector<std::string> row = expected.next();
        if (row.size() == 0) break;
        ASSERT_EQ(row[1], map.geneIdForIndex(map.get(row[0].c_str()))) << "key is " << row[0];
        numberOfCount += 1;
    } while(1);

    ASSERT_TRUE(numberOfCount == map.numberOfEntries()) << "expected: " << numberOfCount << "  actual: " << map.numberOfEntries();

    TabTableReader coverageExpected;
    ASSERT_TRUE(coverageExpected.open("data/coverage.txt"));
    coverageExpected.next(); // skip header

    std::map<std::string, int> totalSelectedFragments;
    std::map<std::string, int> numberOfTranscript;
    std::vector<std::string> row;
    while ((row = coverageExpected.next()).size() != 0) {
        char *endptr;
        long value = strtol(row[5].c_str(), &endptr, 10);
        ASSERT_EQ('\0', *endptr);
        totalSelectedFragments[row[0]] += value;
        numberOfTranscript[row[0]] += 1;
    }

    for (std::map<std::string, int>::const_iterator it = totalSelectedFragments.begin();
         it != totalSelectedFragments.end(); ++it) {
        double averageNumberOfSelected = (double)(it->second) / numberOfTranscript[it->first];
        ASSERT_EQ(averageNumberOfSelected, map.foundFragmentsForIndex(map.indexForGeneID(it->first)))
            << "GeneID " << it->first << " Internal ID " << map.indexForGeneID(it->first);
    }

    ASSERT_EQ(10U, map.sequenceLength());
    
    rmdir("tmp");
}
#endif

TEST(GENE2REFSEQ, LOAD) {
    Gene2Refseq gene2refseq;
    ASSERT_TRUE(gene2refseq.load("data/gene2refseq_test.txt.bz2"));
    ASSERT_EQ("10090", gene2refseq.gene2tax("12295"));
    ASSERT_EQ("10090", gene2refseq.gene2tax("12297"));
    ASSERT_EQ("10090", gene2refseq.gene2tax("12298"));

    std::set<std::string> set1;
    set1.insert("XM_006526759.1");
    set1.insert("XM_006526760.1");
    set1.insert("XM_006526761.1");
    set1.insert("XR_386446.1");
    set1.insert("NM_001163144.1");
    set1.insert("NM_001190483.1");
    ASSERT_EQ(set1, gene2refseq.gene2refseq("18552"));

    //std::set<std::string> set2;
    //set2.insert("NM_001013853.1");
    //ASSERT_EQ(set2, gene2refseq.gene2refseq("287167"));

    ASSERT_EQ("12297", gene2refseq.refseq2gene("XM_006520361.1"));
    ASSERT_EQ("10090", gene2refseq.refseq2tax("NR_045533.1"));
}

TEST(GENE2REFSEQ, LOAD2) {
    Gene2Refseq gene2refseq;
    ASSERT_TRUE(gene2refseq.load("data/gene2refseq_test.txt.bz2", "9606"));
    ASSERT_EQ("", gene2refseq.gene2tax("24440"));

    std::set<std::string> set1;
    ASSERT_EQ(set1, gene2refseq.gene2refseq("361619"));

    std::set<std::string> set2;
    ASSERT_EQ(set2, gene2refseq.gene2refseq("287167"));

    ASSERT_EQ("", gene2refseq.refseq2gene("NM_033234.1"));
    ASSERT_EQ("", gene2refseq.refseq2tax("NM_001013853.1"));
}

TEST(SEQUENCE, REVERSECOMPLEMENT) {
    char buf[256];
    
    bzero(buf, sizeof(buf));
    ASSERT_TRUE(reverseComplement(buf, sizeof(buf), "ATCG", 4));
    ASSERT_STREQ("CGAT", buf);

    ASSERT_TRUE(reverseComplement(buf, sizeof(buf), "CGNT", 4));
    ASSERT_STREQ("ANCG", buf);
}


TEST(SEQUENCE, REVERSECOMPLEMENT2) {
    ASSERT_EQ(std::string("CGAT"), reverseComplement("ATCG"));
    ASSERT_EQ(std::string("ANCG"), reverseComplement("CGNT"));
}
