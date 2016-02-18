#include "gtest/gtest.h"
#include "quickcommon.hpp"
#include "genomehash.hpp"
#include "tabtable.hpp"
#include <stdlib.h>
#include <map>
#include <string>

TEST(GENOMEHASH, OPEN) {
    GenomeHash genomeHash;
    ASSERT_TRUE(genomeHash.open(16, 4));
    ASSERT_TRUE(genomeHash.ready());
    ASSERT_FALSE(genomeHash.open(16, 4));

    GenomeHash genomeHash2;
    ASSERT_FALSE(genomeHash2.open(10, 4));

    GenomeHash genomeHash3;
    ASSERT_FALSE(genomeHash3.open(16, 400));
}

TEST(GENOMEHASH, TEST1) {
    GenomeHash genomeHash;
    ASSERT_TRUE(genomeHash.open(16, 4));
    ASSERT_TRUE(genomeHash.insert("ATCG", 1));
    ASSERT_TRUE(genomeHash.insert("TCGA", 2));

    ASSERT_EQ(1U, genomeHash.get("ATCG"));
    ASSERT_EQ(2U, genomeHash.get("TCGA"));
    ASSERT_EQ(GENOMEHASH_UNMAPPED_ID, genomeHash.get("CGAA"));

    ASSERT_EQ(0, system("mkdir -p ./tmp"));
    FILE *f = fopen("./tmp/genomehash1.dat", "w");
    ASSERT_TRUE(genomeHash.save(f));
    fclose(f);
    genomeHash.printInformation(stdout);
}

TEST(GENOMEHASH, TEST2) {
    GenomeHash genomeHash;
    ASSERT_TRUE(genomeHash.open(4, 4));
    ASSERT_TRUE(genomeHash.insert("ATCG", 1));
    ASSERT_TRUE(genomeHash.insert("TCGA", 2));
    ASSERT_TRUE(genomeHash.insert("ACGA", 3));
    ASSERT_TRUE(genomeHash.insert("GCGA", 4));


    ASSERT_EQ(1U, genomeHash.get("ATCG"));
    ASSERT_EQ(2U, genomeHash.get("TCGA"));
    ASSERT_EQ(3U, genomeHash.get("ACGA"));
    ASSERT_EQ(4U, genomeHash.get("GCGA"));

    ASSERT_EQ(GENOMEHASH_FAILED_ID, genomeHash.get("GCGT"));
    
    ASSERT_FALSE(genomeHash.insert("GCGT", 5));
}

TEST(GENOMEHASH, LARGE) {
    GenomeHash genomeHash;
    genomeHash.open(0x10000, 10);

    {
        TabTableReader expected;
        ASSERT_TRUE(expected.open("./data/unique-10gram.txt.bz2"));

        size_t numberOfCount = 0; // # of fragment
        do {
            std::vector<std::string> row = expected.next();
            if (row.size() == 0) break;
            bool ok;
            long geneid = string2long(row[1], &ok);
            ASSERT_TRUE(ok);
            genomeHash.insert(row[0].c_str(), geneid);
        
            numberOfCount += 1;
        } while(1);
        
        ASSERT_EQ(numberOfCount, genomeHash.numberOfEntries());
        fprintf(stderr, "Hash size: %"UINT64F"u\n", genomeHash.hashtableSize());
        fprintf(stderr, "# of Fragment: %zu\n", numberOfCount);
    }

    {
        TabTableReader expected;
        ASSERT_TRUE(expected.open("./data/unique-10gram.txt.bz2"));

        do {
            std::vector<std::string> row = expected.next();
            if (row.size() == 0) break;
            bool ok;
            long geneid = string2long(row[1], &ok);
            ASSERT_TRUE(ok);
            ASSERT_EQ((uint32_t)geneid, genomeHash.get(row[0].c_str()));
        } while(1);
    }

    {
        ASSERT_EQ(0, system("mkdir -p ./tmp"));
        FILE *f = fopen("./tmp/genomehash2.dat", "w");
        ASSERT_TRUE(genomeHash.save(f));
        fclose(f);
    }
    
    {
        GenomeHash genomeHash2;
        FILE *f = fopen("./tmp/genomehash2.dat", "r");
        ASSERT_TRUE(genomeHash2.load(f));
        fclose(f);

                TabTableReader expected;
        ASSERT_TRUE(expected.open("./data/unique-10gram.txt.bz2"));

        do {
            std::vector<std::string> row = expected.next();
            if (row.size() == 0) break;
            bool ok;
            long geneid = string2long(row[1], &ok);
            ASSERT_TRUE(ok);
            ASSERT_EQ((uint32_t)geneid, genomeHash.get(row[0].c_str()));

        } while(1);

        ASSERT_FALSE(genomeHash2.open(16, 10));
    }
}

TEST(GENOMEHASH, BESTSIZE) {
    ASSERT_EQ(32U, bestHashSizeForExpectedEntries(16));
    ASSERT_EQ(32U, bestHashSizeForExpectedEntries(15));
    ASSERT_EQ(32U, bestHashSizeForExpectedEntries(14));
    ASSERT_EQ(32U, bestHashSizeForExpectedEntries(13));
    ASSERT_EQ(32U, bestHashSizeForExpectedEntries(12));
    ASSERT_EQ(32U, bestHashSizeForExpectedEntries(11));
    ASSERT_EQ(32U, bestHashSizeForExpectedEntries(10));
    ASSERT_EQ(32U, bestHashSizeForExpectedEntries(9));
    ASSERT_EQ(16U, bestHashSizeForExpectedEntries(8));
    ASSERT_EQ(2048U, bestHashSizeForExpectedEntries(1000));
}

TEST(GENOMEHASH, CURSOR) {
    GenomeHash genomeHash;
    ASSERT_TRUE(genomeHash.open(16, 4));
    ASSERT_TRUE(genomeHash.insert("ATCG", 1));
    ASSERT_TRUE(genomeHash.insert("TCGA", 2));
    ASSERT_TRUE(genomeHash.insert("ACGA", 3));
    ASSERT_TRUE(genomeHash.insert("GCGA", 4));

    GenomeHashCursor cursor(&genomeHash);
    std::map<std::string, uint32_t> cursorResult;

    char seq[10];
    uint32_t id;
    while ((id = cursor.next(seq, sizeof(seq))) != GENOMEHASH_FAILED_ID) {
        cursorResult[seq] = id;
    }

    ASSERT_EQ(4U, cursorResult.size());
    ASSERT_EQ(1U, cursorResult["ATCG"]);
    ASSERT_EQ(2U, cursorResult["TCGA"]);
    ASSERT_EQ(3U, cursorResult["ACGA"]);
    ASSERT_EQ(4U, cursorResult["GCGA"]);
}
