#include "gtest/gtest.h"
#include "quickcommon.hpp"
#include "dbbuilder.hpp"
#include "gene2refseq.hpp"
#include "dbcoverage.hpp"
#include "tabtable.hpp"

#include <stdio.h>

TEST(COVERAGE, Coverage1) {
    const char *path = "./tmp/dbcoverage1.dat";
    const int NGRAM = 10;
    
    mkdir("tmp", 0700);
    //TemporaryFile file("dbbuilder", ".dat", "./tmp");
    ASSERT_TRUE(buildDatabase(path, "./tmp/dbcoverage1-coverage.txt", "./data/refseq.rna.bz2", "./data/gene2refseq_test.txt.bz2", NGRAM, true));

    std::map<std::string, TranscriptCoverage> coverages;
    //ASSERT_TRUE(dbcoverage(path, "./data/refseq.rna.bz2", "./data/gene2refseq_test.txt.bz2", &coverages));

    TabTableReader resultReader;
    ASSERT_TRUE(resultReader.open("./tmp/dbcoverage1-coverage.txt"));
    resultReader.next(); // skip header
    std::vector<std::string> row;
    
    while ((row = resultReader.next()).size() > 0) {
        bool ok;
        long fragment = string2long(row[4], &ok);
        ASSERT_TRUE(ok);
        //fprintf(stderr, "%s .. %s\n%s\n", row[1].c_str(), stringsplit(row[1], "|")[3].c_str(), row[7].c_str());
        ASSERT_EQ(8UL, row.size());
        coverages[stringsplit(row[1], "|")[3]] = TranscriptCoverage(row[1], row[0], std::string(row[7].size(), 'A'),
                                                                    row[7], fragment, NGRAM);
    }

    /*
    TabTableReader reader;
    ASSERT_TRUE(reader.open("./data/index-coverage.txt"));
    reader.next(); // skip header

    while ((row = reader.next()).size() > 0) {
        std::string refseqid = row[1];
        TranscriptCoverage record = coverages[refseqid];
        ASSERT_EQ(row[4], record.cover());
        ASSERT_EQ(row[0], record.geneId());

        bool ok;
        long value;

        value = string2long(row[3], &ok);
        ASSERT_TRUE(ok);
        ASSERT_EQ(value, record.selectedFragments());
    }
    */

    TabTableReader reader2;
    ASSERT_TRUE(reader2.open("./data/coverage.txt"));
    reader2.next(); // skip header
    
    while ((row = reader2.next()).size() > 0) {
        std::string refseqid = row[1];
        TranscriptCoverage record = coverages[refseqid];

        char *endptr;
        double coverage = strtod(row[6].c_str(), &endptr);
        ASSERT_EQ('\0', *endptr);
        ASSERT_NEAR(coverage, record.coverage(), 0.0000000001);

        ASSERT_EQ(row[0].c_str(), record.geneId());

        size_t length = strtod(row[2].c_str(), &endptr);
        ASSERT_EQ('\0', *endptr);
        ASSERT_EQ(length, record.cover().size());

        int coveredLength = strtod(row[3].c_str(), &endptr);
        ASSERT_EQ('\0', *endptr);
        ASSERT_EQ(coveredLength, record.coveredLength());

        int fragment = strtod(row[4].c_str(), &endptr);
        ASSERT_EQ('\0', *endptr);
        ASSERT_EQ(fragment, record.fragment());

        int foundFragment = strtod(row[5].c_str(), &endptr);
        ASSERT_EQ('\0', *endptr);
        ASSERT_EQ(foundFragment, record.foundFragments());

        ASSERT_EQ(row[7], record.cover());
    }
}
