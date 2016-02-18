#include "quickcommon.hpp"
#include <stdlib.h>
#include <errno.h>

std::vector<std::string> stringsplit(const std::string &str, const std::string &delimiter)
{
    std::vector<std::string> result;
    std::string::size_type pos = - delimiter.size();
    do {
        std::string::size_type newpos = str.find(delimiter, pos + delimiter.size());
        if (newpos == std::string::npos) {
            result.push_back(str.substr(pos + delimiter.size()));
            break;
        }
        result.push_back(str.substr(pos + delimiter.size(), newpos - pos - delimiter.size()));
        pos = newpos;
    } while(1);
    
    return result;
}

long string2long(const std::string &value, bool *ok)
{
    *ok = false;
    if (value.size() == 0) return 0;

    errno = 0;
    char *endptr;
    long ret = strtol(value.c_str(), &endptr, 0);
    if (*endptr != '\0')
        return 0;
    if (errno != 0)
        return ret;
    *ok = true;
    return ret;
}

double string2double(const std::string &value, bool *ok)
{
    *ok = false;
    if (value.size() == 0) return 0;

    errno = 0;
    char *endptr;
    double ret = strtod(value.c_str(), &endptr);
    if (*endptr != '\0')
        return 0;
    if (errno != 0)
        return ret;
    *ok = true;
    return ret;
}

static uint8_t atcg2bit_converttable[4][256];

static uint8_t atcg2bit_initailize() {
    for (int i = 0; i < 4; i++) {
        atcg2bit_converttable[i][(int)'T'] = SEQ2BIT_T << ((3-i)*2);
        atcg2bit_converttable[i][(int)'C'] = SEQ2BIT_C << ((3-i)*2);
        atcg2bit_converttable[i][(int)'G'] = SEQ2BIT_G << ((3-i)*2);
    }
    return 0;
}

static uint8_t atcg2bit_dummy = atcg2bit_initailize();

bool atcg2bit_nocheck(uint8_t *dest, size_t destbufsize, const char *atcg, size_t len)
{
    for (size_t i = 0; i < len/4; ++i) {
        uint8_t ch = 0;
        ch |= atcg2bit_converttable[0][(int)atcg[i*4]];
        ch |= atcg2bit_converttable[1][(int)atcg[i*4 + 1]];
        ch |= atcg2bit_converttable[2][(int)atcg[i*4 + 2]];
        ch |= atcg2bit_converttable[3][(int)atcg[i*4 + 3]];
        dest[i] = ch;
    }
    
    uint8_t ch = 0;
    for (size_t i = 0; i < len % 4; ++i) {
        ch |= atcg2bit_converttable[i][(int)atcg[len/4*4 + i]];
    }

    dest[len/4] = ch;
    return true;
}

bool atcg2bit(uint8_t *dest, size_t destbufsize, const char *atcg, size_t len)
{
    if ((len + 3)/4 > destbufsize) // buffer overflow
        return false;

    uint8_t ch = 0;
    for (size_t i = 0; i < len; ++i) {
        if ((i % 4) == 0)
            ch = 0;

        int shift = 2 * (3 - (i % 4));
        switch (atcg[i]) {
        case 'A':
            ch |= SEQ2BIT_A << shift;
            break;
        case 'T':
            ch |= SEQ2BIT_T << shift;
            break;
        case 'C':
            ch |= SEQ2BIT_C << shift;
            break;
        case 'G':
            ch |= SEQ2BIT_G << shift;
            break;
        default:
            return false;
        }

        if ((i % 4) == 3)
            dest[i/4] = ch;
    }

    dest[len/4] = ch;
    return true;    
}

bool bit2atcg(char *dest, size_t destbufsize, const uint8_t *bit, size_t len)
{
    if (destbufsize - 1 < len) // buffer overflow with null termination
        return false;

    for (size_t i = 0; i < len; ++i) {
        switch ((bit[i/4] >> (2*(3-(i%4)))) & 0x3) {
        case SEQ2BIT_A:
            dest[i] = 'A';
            break;
        case SEQ2BIT_T:
            dest[i] = 'T';
            break;
        case SEQ2BIT_C:
            dest[i] = 'C';
            break;
        case SEQ2BIT_G:
            dest[i] = 'G';
            break;
        }
    }

    dest[len] = 0;
    
    return true;
}
