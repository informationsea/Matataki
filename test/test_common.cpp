#include "gtest/gtest.h"
#include "quickcommon.hpp"

TEST(COMMON, SPLIT) {
    std::vector<std::string> expected;
    expected.push_back("A");
    expected.push_back("B");
    expected.push_back("C");
    ASSERT_EQ(expected, stringsplit("A,B,C", ","));
    ASSERT_EQ(expected, stringsplit("A_!B_!C", "_!"));
}

TEST(COMMON, CONVERT) {
    bool ok;
    long value;

    value = string2long("12", &ok);
    ASSERT_EQ(12, value);
    ASSERT_TRUE(ok);

    value = string2long("2147483647", &ok);
    ASSERT_EQ(2147483647, value);
    ASSERT_TRUE(ok);

    value = string2long("9223372036854775808", &ok);
    ASSERT_EQ(9223372036854775807L, value);
    ASSERT_FALSE(ok);

    value = string2long("-10", &ok);
    ASSERT_EQ(-10, value);
    ASSERT_TRUE(ok);

    value = string2long("-10k", &ok);
    ASSERT_EQ(0, value);
    ASSERT_FALSE(ok);

    value = string2long("abc", &ok);
    ASSERT_EQ(0, value);
    ASSERT_FALSE(ok);

    value = string2long("0x10", &ok);
    ASSERT_EQ(16, value);
    ASSERT_TRUE(ok);
}

TEST(COMMON, ATCG2BIT) {
    uint8_t buf[3];
    ASSERT_TRUE(atcg2bit(buf, sizeof(buf), "A", 1));
    ASSERT_EQ(SEQ2BIT_A << 6, buf[0]);

    ASSERT_TRUE(atcg2bit(buf, sizeof(buf), "T", 1));
    ASSERT_EQ(SEQ2BIT_T << 6, buf[0]);

    ASSERT_TRUE(atcg2bit(buf, sizeof(buf), "C", 1));
    ASSERT_EQ(SEQ2BIT_C << 6, buf[0]);

    ASSERT_TRUE(atcg2bit(buf, sizeof(buf), "G", 1));
    ASSERT_EQ(SEQ2BIT_G << 6, buf[0]);

    ASSERT_TRUE(atcg2bit(buf, sizeof(buf), "ATCG", 4));
    ASSERT_EQ((SEQ2BIT_A << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_C << 2) | (SEQ2BIT_G << 0), buf[0]);

    ASSERT_TRUE(atcg2bit(buf, sizeof(buf), "ATCGAT", 4));
    ASSERT_EQ((SEQ2BIT_A << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_C << 2) | (SEQ2BIT_G << 0), buf[0]);
    
    ASSERT_TRUE(atcg2bit(buf, sizeof(buf), "CTGA", 4));
    ASSERT_EQ((SEQ2BIT_C << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_G << 2) | (SEQ2BIT_A << 0), buf[0]);

    ASSERT_TRUE(atcg2bit(buf, sizeof(buf), "CTGAATCG", 8));
    ASSERT_EQ((SEQ2BIT_C << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_G << 2) | (SEQ2BIT_A << 0), buf[0]);
    ASSERT_EQ((SEQ2BIT_A << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_C << 2) | (SEQ2BIT_G << 0), buf[1]);

    ASSERT_TRUE(atcg2bit(buf, sizeof(buf), "CTGAAG", 6));
    ASSERT_EQ((SEQ2BIT_C << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_G << 2) | (SEQ2BIT_A << 0), buf[0]);
    ASSERT_EQ((SEQ2BIT_A << 6) | (SEQ2BIT_G << 4) , buf[1]);

    ASSERT_TRUE(atcg2bit(buf, sizeof(buf), "CTGAAGTAG", 9));
    ASSERT_EQ((SEQ2BIT_C << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_G << 2) | (SEQ2BIT_A << 0), buf[0]);
    ASSERT_EQ((SEQ2BIT_A << 6) | (SEQ2BIT_G << 4) | (SEQ2BIT_T << 2) | (SEQ2BIT_A << 0), buf[1]);
    ASSERT_EQ((SEQ2BIT_G << 6), buf[2]);

    ASSERT_FALSE(atcg2bit(buf, sizeof(buf), "CTGAAGK", 7));
    ASSERT_FALSE(atcg2bit(buf, sizeof(buf), "ATCGATCGATCGATCG", 16));
}

TEST(COMMON, ATCG2BIT2) {
    uint8_t buf[3];
    ASSERT_TRUE(atcg2bit_nocheck(buf, sizeof(buf), "A", 1));
    ASSERT_EQ(SEQ2BIT_A << 6, buf[0]);

    ASSERT_TRUE(atcg2bit_nocheck(buf, sizeof(buf), "T", 1));
    ASSERT_EQ(SEQ2BIT_T << 6, buf[0]);

    ASSERT_TRUE(atcg2bit_nocheck(buf, sizeof(buf), "C", 1));
    ASSERT_EQ(SEQ2BIT_C << 6, buf[0]);

    ASSERT_TRUE(atcg2bit_nocheck(buf, sizeof(buf), "G", 1));
    ASSERT_EQ(SEQ2BIT_G << 6, buf[0]);

    ASSERT_TRUE(atcg2bit_nocheck(buf, sizeof(buf), "ATCG", 4));
    ASSERT_EQ((SEQ2BIT_A << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_C << 2) | (SEQ2BIT_G << 0), buf[0]);

    ASSERT_TRUE(atcg2bit_nocheck(buf, sizeof(buf), "ATCGAT", 4));
    ASSERT_EQ((SEQ2BIT_A << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_C << 2) | (SEQ2BIT_G << 0), buf[0]);
    
    ASSERT_TRUE(atcg2bit_nocheck(buf, sizeof(buf), "CTGA", 4));
    ASSERT_EQ((SEQ2BIT_C << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_G << 2) | (SEQ2BIT_A << 0), buf[0]);

    ASSERT_TRUE(atcg2bit_nocheck(buf, sizeof(buf), "CTGAATCG", 8));
    ASSERT_EQ((SEQ2BIT_C << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_G << 2) | (SEQ2BIT_A << 0), buf[0]);
    ASSERT_EQ((SEQ2BIT_A << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_C << 2) | (SEQ2BIT_G << 0), buf[1]);

    ASSERT_TRUE(atcg2bit_nocheck(buf, sizeof(buf), "CTGA.G", 6));
    ASSERT_EQ((SEQ2BIT_C << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_G << 2) | (SEQ2BIT_A << 0), buf[0]);
    ASSERT_EQ((SEQ2BIT_A << 6) | (SEQ2BIT_G << 4) , buf[1]);

    ASSERT_TRUE(atcg2bit_nocheck(buf, sizeof(buf), "CTGAAGTXG", 9));
    ASSERT_EQ((SEQ2BIT_C << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_G << 2) | (SEQ2BIT_A << 0), buf[0]);
    ASSERT_EQ((SEQ2BIT_A << 6) | (SEQ2BIT_G << 4) | (SEQ2BIT_T << 2) | (SEQ2BIT_A << 0), buf[1]);
    ASSERT_EQ((SEQ2BIT_G << 6), buf[2]);
}

TEST(COMMON, BIT2ATCG) {
    char buf[10];
    uint8_t seq_a[] = {SEQ2BIT_A << 6,};
    ASSERT_TRUE(bit2atcg(buf, sizeof(buf), seq_a, 1));
    ASSERT_STREQ("A", buf);

    uint8_t seq_t[] = {SEQ2BIT_T << 6,};
    ASSERT_TRUE(bit2atcg(buf, sizeof(buf), seq_t, 1));
    ASSERT_STREQ("T", buf);

    uint8_t seq_c[] = {SEQ2BIT_C << 6,};
    ASSERT_TRUE(bit2atcg(buf, sizeof(buf), seq_c, 1));
    ASSERT_STREQ("C", buf);

    uint8_t seq_g[] = {SEQ2BIT_G << 6,};
    ASSERT_TRUE(bit2atcg(buf, sizeof(buf), seq_g, 1));
    ASSERT_STREQ("G", buf);

    uint8_t seq2[] = {(SEQ2BIT_C << 6) | (SEQ2BIT_T << 4) | (SEQ2BIT_G << 2) | (SEQ2BIT_A << 0),
                      (SEQ2BIT_A << 6) | (SEQ2BIT_G << 4) | (SEQ2BIT_T << 2) | (SEQ2BIT_A << 0),
                      (SEQ2BIT_G << 6)};
    ASSERT_TRUE(bit2atcg(buf, sizeof(buf), seq2, 9));
    ASSERT_STREQ("CTGAAGTAG", buf);

    char smallbuf[9];
    ASSERT_FALSE(bit2atcg(smallbuf, sizeof(smallbuf), seq2, 9));
}
