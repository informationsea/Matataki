#ifndef GENOMEHASH_H
#define GENOMEHASH_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#ifndef UINT32_MAX
#define UINT32_MAX         (4294967295U)
#endif

#define GENOMEHASH_UNMAPPED_ID 0U
#define GENOMEHASH_FAILED_ID (UINT32_MAX)

uint64_t bestHashSizeForExpectedEntries(uint64_t expected);

struct GenomeHashHeader;
class GenomeHashCursor;

class GenomeHash
{
    friend class GenomeHashCursor;
public:
    GenomeHash();
    virtual ~GenomeHash();

    bool open(size_t hashtableSize, size_t seqlen); // create new empty hash
    virtual bool load(FILE *file); // load from file
    virtual bool save(FILE *file) const; // write to file

    bool insert(const char *seq, uint32_t id);
    bool insert2bit(const uint8_t *seq, uint32_t id);

    uint32_t get(const char *seq) const;
    uint32_t get2bit(const uint8_t *seq) const;

    bool ready() const {return m_ready;}
    void printInformation(FILE *output) const;

    uint32_t sequenceLength() const;
    uint64_t hashtableSize() const;
    uint64_t numberOfEntries() const;

private:
    struct GenomeHashHeader *m_header;
    bool m_mmap;
    uint8_t *m_buffer;
    bool m_ready;
};

class GenomeHashCursor
{
public:
    GenomeHashCursor(const GenomeHash *genomeHash);
    virtual ~GenomeHashCursor() {};

    uint32_t next(char *seqbuf, size_t seqbuflen);
    uint32_t next2bit(uint8_t *buf, size_t buflen);
private:
    const GenomeHash *m_genomeHash;
    size_t m_pos;
};

#endif /* GENOMEHASH_H */
