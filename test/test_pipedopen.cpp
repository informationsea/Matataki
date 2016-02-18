#include "gtest/gtest.h"
#include "pipedopen.h"
#include <errno.h>
#include <string.h>

void common_test(const char *filename) {
    PFILE *file = pipedopen(filename, false);
    ASSERT_TRUE(file);
    ASSERT_TRUE(file->file) << strerror(errno);

    char buf[300];
    ASSERT_EQ(buf, fgets(buf, sizeof(buf), file->file));
    ASSERT_STREQ("Hello, world!\n", buf);
    ASSERT_EQ(buf, fgets(buf, sizeof(buf), file->file));
    ASSERT_STREQ("This is sample text file\n", buf);
    ASSERT_EQ(NULL, fgets(buf, sizeof(buf), file->file));
    ASSERT_TRUE(feof(file->file));
    
    pipedclose(file);
}

TEST(PIPEDOPEN, READ_GZIP) {
    common_test("./data/sample.txt.gz");
}

TEST(PIPEDOPEN, READ_BZIP2) {
    common_test("./data/sample.txt.bz2");
}

TEST(PIPEDOPEN, READ_XZ) {
    common_test("./data/sample.txt.xz");
}

TEST(PIPEDOPEN, READ_PLAIN) {
    common_test("./data/sample.txt");
}

TEST(PIPEDOPEN, READ_PLAIN_FAIL) {
    ASSERT_FALSE(pipedopen("./this/file/is/not/exist", false));
}

TEST(PIPEDOPEN, READ_XZ_FAIL) {
    ASSERT_FALSE(pipedopen("./this/file/is/not/exist.xz", false));
}

TEST(PIPEDOPEN, READ_GZIP_FAIL) {
    ASSERT_FALSE(pipedopen("./this/file/is/not/exist.gz;", false));
}

TEST(PIPEDOPEN, WRITE_GZIP) {
    PFILE *file = pipedopen("./tmp/sample.txt.gz", true);
    fprintf(stderr, "writing..\n");
    fprintf(file->file, "Hello, world!\n" "This is sample text file\n");
    fprintf(stderr, "closing..\n");
    pipedclose(file);
}
