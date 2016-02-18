#ifndef QUICKCOMMON_H
#define QUICKCOMMON_H

#include "config.h"
#include <string>
#include <vector>
#include <string.h>
#include <limits.h>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <stdint.h>

#define _CONCAT(x, y) x ## y
#define CONCAT(x, y) _CONCAT(x, y)

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define LENGTH(x) (sizeof(x)/sizeof(x[0]))

std::vector<std::string> stringsplit(const std::string &str, const std::string &delimiter);
long string2long(const std::string &value, bool *ok);
double string2double(const std::string &value, bool *ok);

#define SEQ2BIT_A 0
#define SEQ2BIT_T 1
#define SEQ2BIT_C 2
#define SEQ2BIT_G 3

#ifdef __APPLE__
#define UINT64F "ll"
#else
#define UINT64F "l"
#endif

bool atcg2bit(uint8_t *dest, size_t destbufsize, const char *atcg, size_t len);
bool atcg2bit_nocheck(uint8_t *dest, size_t destbufsize, const char *atcg, size_t len); // treat invalid character as A
bool bit2atcg(char *dest, size_t destbufsize, const uint8_t *bit, size_t len);

#define DEBUG(fmt,...) (fprintf(stderr, __FILE__ ": %3d: " fmt "\n" ,__LINE__, ## __VA_ARGS__))

class TemporaryFile
{
public:
    TemporaryFile(const std::string &prefix, const std::string &suffix = "", std::string tempdir = "") : autodelete(true) {
        if (tempdir.size() == 0) {
            const char* tmpdir_parent = std::getenv("TMPDIR");
            if (tmpdir_parent == NULL)
                tempdir = "/tmp";
            else
                tempdir = tmpdir_parent;
        }
        
        std::snprintf(tempfilepath, sizeof(tempfilepath)-1, "%s/%s-XXXXXX%s", tempdir.c_str(), prefix.c_str(), suffix.c_str());
        if (mkstemps(tempfilepath, suffix.size()) == -1) {
            bzero(tempfilepath, sizeof(tempfilepath));
        }
    }

    virtual ~TemporaryFile() {
        if (autodelete)
            unlink(tempfilepath);
    }

    bool ready() const {return tempfilepath[0] != '\0';}
    const char *path() const {return tempfilepath;}
    void setAutoDelete(bool enable) {autodelete = enable;}
    
private:
    char tempfilepath[PATH_MAX+1];
    bool autodelete;
};

#endif /* QUICKCOMMON_H */
