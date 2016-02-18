#include "genomehash.hpp"
#include "quickcommon.hpp"
#include <sys/mman.h>
#include <string.h>

#define GENOMEHASH_HASHTABLE_OFFSET (4*1024)
#define GENOMEHASH_MAGIC "GENOMEHASH"
#define GENOMEHASH_VERSION "0.0.1"
#define GENOMEHASH_IDLEN (sizeof(uint32_t))
#define GENOMEHASH_ALIGNMENT GENOMEHASH_IDLEN
#define GENOMEHASH_MAXIMUM_SEQLEN 128
#define GENOMEHASH_MURMURSEED 3590585660
#define SEQLENGTH_IN_BYTE 4

extern "C" {
    void MurmurHash3_x64_128 ( const void * key, int len, uint32_t seed, void * out );
}

static inline size_t alignmentedSize(size_t size)
{
    return ((size + GENOMEHASH_ALIGNMENT - 1)/GENOMEHASH_ALIGNMENT)*GENOMEHASH_ALIGNMENT;
}

static inline size_t bitlen2bytelen(size_t size)
{
    return ((size + SEQLENGTH_IN_BYTE - 1)/SEQLENGTH_IN_BYTE);
}

uint64_t bestHashSizeForExpectedEntries(uint64_t expected)
{
    uint64_t minimumHashTableSize = expected * 2;
    for (int i = 0; i < 64; ++i) {
        if (minimumHashTableSize <= (1 << i)) {
            return 1 << i;
        }
    }
    return 0; // Failed to suggest
}

struct GenomeHashHeader {
    char magicbytes[32];
    char version[32];
    uint64_t hashtableSize;
    uint64_t hashMask;
    uint32_t seqlen;
    uint32_t seqbytelen;
    uint32_t entrylen;
    uint64_t hashtableOffset;
    uint64_t hashtableBufferSize;
    uint64_t numberOfEntries;

    GenomeHashHeader(){}
    
    GenomeHashHeader(size_t hashtableSize, size_t seqlen) :
        hashtableSize(hashtableSize),
        seqlen(seqlen), hashtableOffset(GENOMEHASH_HASHTABLE_OFFSET), numberOfEntries(0) {

        strncpy(magicbytes, GENOMEHASH_MAGIC, sizeof(magicbytes));
        strncpy(version, GENOMEHASH_VERSION, sizeof(version));

        seqbytelen = bitlen2bytelen(seqlen);
        entrylen = alignmentedSize(seqbytelen + GENOMEHASH_IDLEN);
        hashtableBufferSize = entrylen * hashtableSize;
    }
};

GenomeHash::GenomeHash() : m_header(0), m_mmap(false), m_buffer(0), m_ready(false)
{}

GenomeHash::~GenomeHash()
{
    if (m_mmap) {
        munmap(m_buffer, m_header->hashtableBufferSize);
    } else {
        delete m_buffer;
    }
    delete m_header;
}

bool GenomeHash::open(size_t hashtableSize, size_t seqlen)
{
    if (m_ready) {
        fprintf(stderr, "Already opened\n");
        return false; // cannot open twice
    }
    
    uint64_t hashmask = 0;
    // check whether bufferSize is power of 2 or not
    for (size_t i = 0; i < sizeof(size_t)*8; ++i) {
        if (hashtableSize & (1 << i)) {
            if ((hashtableSize ^ (1 << i)) == 0) {
                hashmask = (1 << i) - 1;
                break;
            } else {
                fprintf(stderr, "Invalid buffer size %zu\n", hashtableSize);
                return false;
            }
        }
    }

    // check maximum sequence length
    if (seqlen > GENOMEHASH_MAXIMUM_SEQLEN) {
        fprintf(stderr, "Too long sequence length\n");
        return false;
    }

    m_header = new GenomeHashHeader(hashtableSize, seqlen);
    m_header->hashMask = hashmask;
    m_buffer = new uint8_t[m_header->hashtableBufferSize];
    bzero(m_buffer, m_header->hashtableBufferSize);
    m_mmap = false;
    m_ready = true;
    return true;
}

bool GenomeHash::load(FILE *file)
{
    if (m_ready) {
        fprintf(stderr, "Already opened\n");
        return false; // cannot open twice
    }

    m_header = new GenomeHashHeader;
    if (fread(m_header, sizeof(struct GenomeHashHeader), 1, file) != 1) {
        perror("Cannot load genome hash header");
        return false;
    }

    if (strcmp(GENOMEHASH_MAGIC, m_header->magicbytes) != 0) {
        fprintf(stderr, "Invalid signature. May be this file is not index file.\n");
        return false;
    }

    if (strcmp(GENOMEHASH_VERSION, m_header->version) != 0) {
        fprintf(stderr, "Invalid version %s\n", m_header->version);
        return false;
    }

    size_t remainBytes = m_header->hashtableOffset - sizeof(struct GenomeHashHeader);
    char empty[256];
    while (remainBytes != 0) {
        size_t readBytes = fread(empty, 1, MIN(remainBytes, sizeof(empty)), file);
        if (readBytes != MIN(remainBytes, sizeof(empty)))
            return false;
        remainBytes -= readBytes;
    }

    m_buffer = new uint8_t[m_header->hashtableBufferSize];
    if (fread(m_buffer, m_header->hashtableBufferSize, 1, file) != 1) {
        perror("Cannot load genome hash data");
        return false;
    }

    m_ready = true;
    return true;
}

bool GenomeHash::save(FILE *file) const
{
    // write header
    if (fwrite(m_header, sizeof(struct GenomeHashHeader), 1, file) != 1)
        return false;

    // fill padding
    size_t remainBytes = m_header->hashtableOffset - sizeof(struct GenomeHashHeader);
    char empty[256];
    bzero(empty, sizeof(empty));
    while (remainBytes != 0) {
        size_t wroteBytes = fwrite(empty, 1, MIN(remainBytes, sizeof(empty)), file);
        if (wroteBytes != MIN(remainBytes, sizeof(empty)))
            return false;
        remainBytes -= wroteBytes;
    }

    // write data buffer
    if (fwrite(m_buffer, 1, m_header->hashtableBufferSize, file) != m_header->hashtableBufferSize)
        return false;
    
    return true;
}


bool GenomeHash::insert(const char *seq, uint32_t id)
{
    uint8_t buf[GENOMEHASH_MAXIMUM_SEQLEN/4];
    if (atcg2bit(buf, sizeof(buf), seq, m_header->seqlen)) {
        return insert2bit(buf, id);
    }
    return false;
}    

bool GenomeHash::insert2bit(const uint8_t *seq, uint32_t id)
{
    if (id == GENOMEHASH_UNMAPPED_ID || id == GENOMEHASH_FAILED_ID) {
        DEBUG("Invalid ID");
        return false;
    }

    uint64_t hashvalue[2];
    MurmurHash3_x64_128(seq, m_header->seqbytelen, GENOMEHASH_MURMURSEED, hashvalue);
    uint64_t hashpos = hashvalue[0] & m_header->hashMask;
    uint64_t originalpos = hashpos;

    do {
        uint8_t *entry = m_buffer + m_header->entrylen * hashpos;
        if (*((uint32_t *)entry) == GENOMEHASH_UNMAPPED_ID) {
            *((uint32_t *)entry) = id;
            memcpy(entry + GENOMEHASH_IDLEN, seq, m_header->seqbytelen);
            m_header->numberOfEntries += 1;
            return true;
        }
        hashpos += 1;
        if (originalpos == hashpos) {
            fprintf(stderr, "No empty location in hash table\n");
            return false;
        }
        
        if (hashpos >= m_header->hashtableSize)
            hashpos = 0;
    } while(1);
}

uint32_t GenomeHash::get(const char *seq) const
{
    uint8_t buf[GENOMEHASH_MAXIMUM_SEQLEN/4];
    if (atcg2bit(buf, sizeof(buf), seq, m_header->seqlen)) {
        return get2bit(buf);
    }
    return GENOMEHASH_FAILED_ID;
}

uint32_t GenomeHash::get2bit(const uint8_t *seq) const
{
    uint64_t hashvalue[2];
    MurmurHash3_x64_128(seq, m_header->seqbytelen, GENOMEHASH_MURMURSEED, hashvalue);
    uint64_t hashpos = hashvalue[0] & m_header->hashMask;
    uint64_t originalpos = hashpos;
    do {
        uint8_t *entry = m_buffer + m_header->entrylen * hashpos;

        if (*((uint32_t *)entry) == GENOMEHASH_UNMAPPED_ID)
            return GENOMEHASH_UNMAPPED_ID;

        if (memcmp(seq, entry +  GENOMEHASH_IDLEN, m_header->seqbytelen) == 0) {
            return  *((uint32_t *)entry);
        }
        
        hashpos += 1;
        if (hashpos == originalpos)
            return GENOMEHASH_FAILED_ID;
        
        if (hashpos >= m_header->hashtableSize)
            hashpos = 0;
    } while(1);
}

void GenomeHash::printInformation(FILE *output) const
{
    fprintf(output, "               Version: %s\n", m_header->version);
    fprintf(output, "       Hash Table Size: %"UINT64F"u\n", m_header->hashtableSize);
    fprintf(output, "          # of entries: %"UINT64F"u\n", m_header->numberOfEntries);
    fprintf(output, "             Hash Mask: 0x%"UINT64F"x\n", m_header->hashMask);
    fprintf(output, "       Sequence Length: %u\n", m_header->seqlen);
    fprintf(output, "   2Bit Sequence Bytes: %u\n", m_header->seqbytelen);
    fprintf(output, "          Entry Length: %u\n", m_header->entrylen);
    fprintf(output, "     Hash Table Offset: 0x%"UINT64F"x\n", m_header->hashtableOffset);
    fprintf(output, "Hash Table Buffer Size: 0x%"UINT64F"x\n", m_header->hashtableBufferSize);
}

uint32_t GenomeHash::sequenceLength() const
{
    return m_header->seqlen;
}

uint64_t GenomeHash::hashtableSize() const
{
    return m_header->hashtableSize;
}

uint64_t GenomeHash::numberOfEntries() const
{
    return m_header->numberOfEntries;
}

GenomeHashCursor::GenomeHashCursor(const GenomeHash *genomeHash) : m_genomeHash(genomeHash), m_pos(0) {}

uint32_t GenomeHashCursor::next(char *seqbuf, size_t seqbuflen)
{
    uint8_t buf[GENOMEHASH_MAXIMUM_SEQLEN/4];
    uint32_t id = next2bit(buf, sizeof(buf));
    bool ok = bit2atcg(seqbuf, seqbuflen, buf, m_genomeHash->m_header->seqlen);
    if (!ok)
        return GENOMEHASH_FAILED_ID;
    return id;
}

uint32_t GenomeHashCursor::next2bit(uint8_t *buf, size_t buflen)
{
    if (buflen < m_genomeHash->m_header->seqbytelen) // too small buffer size
        return GENOMEHASH_FAILED_ID;
    
    while (m_pos < m_genomeHash->m_header->hashtableSize) {
        uint8_t *entry = m_genomeHash->m_buffer + m_genomeHash->m_header->entrylen * m_pos;
        
        if (*((uint32_t *)entry) == GENOMEHASH_UNMAPPED_ID) {
            m_pos += 1;
            continue;
        } else {
            uint32_t value = *((uint32_t *)entry);
            memcpy(buf, entry +  GENOMEHASH_IDLEN, m_genomeHash->m_header->seqbytelen);
            m_pos += 1;
            return value;
        }
    }
    
    return GENOMEHASH_FAILED_ID; // end of file
}
