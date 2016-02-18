#include "mapping.hpp"
#include "quickcommon.hpp"
#include "dbcommon.hpp"
#include "filewriter.hpp"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "pipedopen.h"

using namespace std;

bool QuickMapping::loadIndex(const std::string &path)
{
    FILE *file = fopen(path.c_str(), "r");
    if (file == NULL) {
        perror("Cannot open index file");
        return false;
    }

    if (!m_index.load(file)) {
        fprintf(stderr, "Failed to load index file");
        return false;
    }

    m_ngram = m_index.sequenceLength();

    m_mapped.clear();
    size_t genenum = m_index.gene2index().size();
    for (size_t i = 0; i < genenum + 1; ++i) {
        m_mapped.push_back(0);
    }
    
    fclose(file);
    return true;
}

void QuickMapping::createSkipArray(uint8_t *skipArray, const char* seq, size_t len, long ngram)
{
    for (int i = len - 1; i >= 0; i--) {
        switch (seq[i]) {
        case 'A':
        case 'T':
        case 'C':
        case 'G':
            break;
        default:
            for (int j = 0; j < MIN(ngram, i+1); j++) {
                skipArray[i - j] = skipArray[i - j + 1] + 1;
            }
            break;
        }
    }
}

uint32_t QuickMapping::mapSequence(const cppfasta::SequenceRecord2 &record) const
{
    //DEBUG("Add Sequence %s", record.sequence);
    size_t len = strlen(record.sequence);
    const char *seq = record.sequence;
    uint32_t candidate = GENOMEHASH_UNMAPPED_ID;
    long foundcount = 0;

    if (len < (size_t)m_ngram) {
        static bool small_warned = false;

        if (!small_warned) {
            fprintf(stderr, "Supplied sequence is smaller than N-gram\n");
            small_warned = true;
        }
        return GENOMEHASH_FAILED_ID;
    }

#if CPPFASTA_MAXIMUM_SEQUENCE_LENGTH > 255
#error Too long CPPFASTA_MAXIMUM_SEQUENCE_LENGTH
#endif
    
    uint8_t steplen[CPPFASTA_MAXIMUM_SEQUENCE_LENGTH+1];
    bzero(steplen, sizeof(steplen));
    createSkipArray(steplen, seq, len, m_ngram);
    
    uint8_t bitseq[4][CPPFASTA_MAXIMUM_SEQUENCE_LENGTH/4];
    if ((m_stepsize % 4) == 0) {
        int i = 0;
        atcg2bit_nocheck(bitseq[i], CPPFASTA_MAXIMUM_SEQUENCE_LENGTH, seq + i, len - i);
    } else {
        for (int i = 0; i < 4; i++) {
            atcg2bit_nocheck(bitseq[i], CPPFASTA_MAXIMUM_SEQUENCE_LENGTH, seq + i, len - i);
        }
    }

    uint8_t mappingresult[CPPFASTA_MAXIMUM_SEQUENCE_LENGTH];
    memset(mappingresult, ' ', sizeof(mappingresult));
    mappingresult[len] = '\0';

    uint8_t tmp[CPPFASTA_MAXIMUM_SEQUENCE_LENGTH/4];
    uint8_t lastByteMask = 0;
    switch (m_index.sequenceLength() % 4) {
    case 0:
        lastByteMask = 0x00;
        break;
    case 1:
        lastByteMask = 0xc0;
        break;
    case 2:
        lastByteMask = 0xf0;
        break;
    case 3:
        lastByteMask = 0xfc;
        break;
    }

    //DEBUG("Step Size: %ld", m_stepsize);
    for (size_t i = 0; i < len - m_ngram + 1; i += m_stepsize) {
        uint32_t value;

        i += ((steplen[i] + m_stepsize - 1)/m_stepsize) * m_stepsize;
        if (i >= len - m_ngram + 1) break;

        //DEBUG("  PROCESSING %zu", i);

        char tmp2[11];
        bzero(tmp2, sizeof(tmp2));
        memcpy(tmp, bitseq[i%4] + i/4, m_index.sequenceLength()/4);
        tmp[m_index.sequenceLength()/4] = bitseq[i%4][i/4 + m_index.sequenceLength()/4] & lastByteMask;

        value = m_index.get2bit(tmp);

        if (value == GENOMEHASH_FAILED_ID) {
            fprintf(stderr, "Failed to map %s\n", seq);
            return GENOMEHASH_FAILED_ID;
        }
        
        if (value == GENOMEHASH_UNMAPPED_ID) {
            mappingresult[i] = '.';
            continue;
        }

        foundcount += 1;
        
        if (candidate != GENOMEHASH_UNMAPPED_ID) {
            if (candidate != value) {
                mappingresult[i] = 'X';
                if (m_resultwriter) m_resultwriter->printf("%s\t%s\t%ld\t%s\t%s - %s\n",
                                                           record.name, "MULTIMAPPED", foundcount, mappingresult,
                                                           m_index.geneIdForIndex(candidate).c_str(), m_index.geneIdForIndex(value).c_str());
                //m_multimapped += 1;
                return GENOMEHASH_MULTIMAP_ID;
            }
            mappingresult[i] = 'O';
        } else {
            mappingresult[i] = 'O';
            candidate = value;
            if (!m_checkall && foundcount >= m_acceptcount) // Accept
                break;
        }
    }

    if (candidate == GENOMEHASH_UNMAPPED_ID || foundcount < m_acceptcount) {
        if (m_resultwriter) m_resultwriter->printf("%s\t%s\t%ld\t%s\n", record.name, "UNMAPPED", foundcount, mappingresult);
        //m_unmapped += 1;
        return GENOMEHASH_UNMAPPED_ID;
    }

    //DEBUG("MAPPED %s", m_index.geneIdForIndex(candidate).c_str());
    if (m_resultwriter) m_resultwriter->printf("%s\t%s\t%ld\t%s\n", record.name, m_index.geneIdForIndex(candidate).c_str(), foundcount, mappingresult);
    //m_nread += 1;
    //m_mapped[candidate] += 1;
    //DEBUG("MAPPED %zu %zu %zu", candidate, m_mapped[candidate], m_mapped.size());
    return candidate;
}

uint32_t QuickMapping::addSequence(const cppfasta::SequenceRecord2 &record)
{
    uint32_t gene = mapSequence(record);
    if (gene == GENOMEHASH_UNMAPPED_ID) {
        m_unmapped += 1;
    } else if (gene == GENOMEHASH_FAILED_ID) {
        // ignore
    } else {
        m_nread += 1;
        m_mapped[gene] += 1;
    }
    return gene;
}

uint32_t QuickMapping::addSequencePair(const cppfasta::SequenceRecord2 &record1, const cppfasta::SequenceRecord2 &record2)
{
    uint32_t gene1 = mapSequence(record1);
    uint32_t gene2 = mapSequence(record2);

    if (gene1 == GENOMEHASH_FAILED_ID && gene2 == GENOMEHASH_FAILED_ID)
        return GENOMEHASH_FAILED_ID;
    if (gene1 == GENOMEHASH_FAILED_ID)
        gene1 = GENOMEHASH_UNMAPPED_ID;
    if (gene2 == GENOMEHASH_FAILED_ID)
        gene2 = GENOMEHASH_UNMAPPED_ID;

    uint32_t accepted_gene;
    if (gene1 == gene2) {
        accepted_gene = gene1;
    } else if (gene1 == GENOMEHASH_MULTIMAP_ID || gene2 == GENOMEHASH_MULTIMAP_ID) {
        accepted_gene = GENOMEHASH_MULTIMAP_ID;
    } else if (gene1 == GENOMEHASH_UNMAPPED_ID) {
        accepted_gene = gene2;
    } else if (gene2 == GENOMEHASH_UNMAPPED_ID) {
        accepted_gene = gene1;
    } else { // conflict
        accepted_gene = GENOMEHASH_MULTIMAP_ID;
    }

    if (accepted_gene == GENOMEHASH_UNMAPPED_ID) {
        m_unmapped += 1;
    } else if (accepted_gene == GENOMEHASH_MULTIMAP_ID) {
        m_multimapped += 1;
    } else if (accepted_gene == GENOMEHASH_FAILED_ID) {
        // ignore
    } else {
        m_nread += 1;
        m_mapped[accepted_gene] += 1;
    }

    
    return accepted_gene;
}

std::map<std::string, int> QuickMapping::result() const
{
    std::map<std::string, int> result;

    for (std::map<std::string, size_t>::const_iterator it = m_index.gene2index().begin();
         it != m_index.gene2index().end(); ++it) {
        //DEBUG("Result for %zu %zu", it->second, m_mapped.at(it->second));
        result[it->first] = m_mapped.at(it->second);
    }

    result["UNMAPPED"] = m_unmapped;
    result["MULTIMAPPED"] = m_multimapped;

    return result;
}

std::map<std::string, double> QuickMapping::fpkm() const
{
    std::map<std::string, double> fpkm_map;

    for (std::map<std::string, size_t>::const_iterator it = m_index.gene2index().begin();
         it != m_index.gene2index().end(); ++it) {
        if (m_nread == 0 || m_index.foundFragmentsForIndex(it->second) == 0)
            fpkm_map[it->first] = nan("");
        else
            fpkm_map[it->first] = m_mapped.at(it->second)/m_index.foundFragmentsForIndex(it->second)*1000/m_nread*1000000;
    }

    return fpkm_map;
}

bool domapping(const char *dbpath, cppfasta::FastqReader2 *reader, std::map<std::string, int> *result, std::map<std::string, double> *fpkm_result, cppfasta::BasicWriter* resultWriter, uint64_t *processedSequences, long stepsize, long acceptcount, bool checkall, bool pairedend)
{
    QuickMapping mapping(stepsize, acceptcount, checkall);
    *processedSequences = 0;
    mapping.setResultWriter(resultWriter);

    fprintf(stderr, "Loading index...\n");
    if (!mapping.loadIndex(dbpath)) {
        fprintf(stderr, "cannot open index\n");
        return false;
    }

    cppfasta::SequenceRecord2 record;
    cppfasta::SequenceRecord2 record2;

    fprintf(stderr, "Start mapping...\n");
    while (reader->next(&record)) {
        if (pairedend) {
            if (!reader->next(&record2)) {
                break;
            }
            mapping.addSequencePair(record, record2);
        } else {
            mapping.addSequence(record);
        }
        *processedSequences += 1;
    }

    fprintf(stderr, "Writing results...\n");

    *result = mapping.result();
    fprintf(stderr, "FPKM...\n");
    *fpkm_result = mapping.fpkm();
    return true;
}
