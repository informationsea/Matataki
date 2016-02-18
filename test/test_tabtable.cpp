#include "gtest/gtest.h"
#include "tabtable.hpp"
#include "filereader.hpp"
#include "quickcommon.hpp"
#include <limits.h>

char EXPECTED_LINES[][16][40] =
    {
//        {"10116", "24440", "PROVISIONAL", "NM_033234.1", "17985948", "NP_150237.1", "17985949", "AC_000069.1", "109644432", "156261275", "156262686", "-", "Alternate Rn_Celera", "-", "-", "Hbb"},
//        {"10116", "24440", "PROVISIONAL", "NM_033234.1", "17985948", "NP_150237.1", "17985949", "NC_005100.4", "666184579", "168971268", "168972679", "-", "Reference Rnor_6.0 Primary Assembly", "-", "-", "Hbb"},
//        {"10116", "24440", "PROVISIONAL", "NM_033234.1", "17985948", "NP_150237.1", "17985949", "NW_001084773.1", "109463005", "18469486", "18470897", "-", "Alternate Rn_Celera", "-", "-", "Hbb"},

{"10090", "12295", "MODEL", "XM_006532091.1", "568971380", "XP_006532154.1", "568971381", "NC_000077.6", "372099099", "98001507", "98022886", "-", "Reference GRCm38.p2 C57BL/6J", "-", "-", "Cacnb1"},
{"10090", "12295", "MODEL", "XM_006532091.1", "568971380", "XP_006532154.1", "568971381", "NT_096135.6", "372098990", "63319887", "63341266", "-", "Reference GRCm38.p2 C57BL/6J", "-", "-", "Cacnb1"},
{"10090", "12295", "MODEL", "XM_006532092.1", "568971382", "XP_006532155.1", "568971383", "NC_000077.6", "372099099", "98001507", "98022886", "-", "Reference GRCm38.p2 C57BL/6J", "-", "-", "Cacnb1"}


    };

TEST(TABTABLE, READER) {
    TabTableReader reader;
    ASSERT_TRUE(reader.open("data/gene2refseq_test.txt.bz2"));

    for (size_t i = 0; i < LENGTH(EXPECTED_LINES); i++) {
        std::vector<std::string> row = reader.next();

        for (size_t j = 0; j < LENGTH(EXPECTED_LINES[0]); j++) {
            ASSERT_EQ(EXPECTED_LINES[i][j], row[j]);
        }
    }

    for (size_t i = 0; i < 171 - LENGTH(EXPECTED_LINES); i++) {
        std::vector<std::string> row = reader.next();
        ASSERT_TRUE(16 == row.size());
    }

    std::vector<std::string> row = reader.next();
    ASSERT_TRUE(0 == row.size());
}

TEST(TABTABLE, WRITER) {
    ASSERT_EQ(0, system("mkdir -p tmp"));

    char tempfilepath[PATH_MAX];
    sprintf(tempfilepath, "tmp/normal-XXXXXX");
    ASSERT_NE(-1, mkstemp(tempfilepath));

    TabTableWriter *writer = new TabTableWriter;
    ASSERT_TRUE(writer->open(tempfilepath));
    std::vector<std::string> row;
    row.push_back("Hello");
    row.push_back("world");
    writer->write(row);
    row.erase(row.begin(), row.end());
    row.push_back("!");
    row.push_back("?");
    writer->write(row);
    delete writer;

    cppfasta::FileReader reader;
    ASSERT_TRUE(reader.open(tempfilepath));
    ASSERT_EQ("Hello\tworld\n", reader.readLine());
    ASSERT_EQ("!\t?\n", reader.readLine());
    ASSERT_EQ("", reader.readLine());

    unlink(tempfilepath);
    rmdir("tmp");
}
