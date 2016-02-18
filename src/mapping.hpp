#ifndef MAPPING_H
#define MAPPING_H

#include <string>
#include <map>
#include "dbcommon.hpp"
#include "cppfasta2.hpp"
#include "filewriter.hpp"
#include <stdint.h>

#define GENOMEHASH_MULTIMAP_ID (UINT32_MAX-1)

class QuickMapping
{
public:
    QuickMapping(long stepsize = 1, long acceptcount = 1, bool checkall = false) :
        m_resultwriter(0), m_stepsize(stepsize), m_acceptcount(acceptcount),
        m_unmapped(0), m_multimapped(0), m_nread(0), m_ngram(0), m_checkall(checkall) {};
    virtual ~QuickMapping(){};

    bool loadIndex(const std::string &path);
    uint32_t mapSequence(const cppfasta::SequenceRecord2 &record) const;
    uint32_t addSequence(const cppfasta::SequenceRecord2 &record);
    uint32_t addSequencePair(const cppfasta::SequenceRecord2 &record1, const cppfasta::SequenceRecord2 &record2);
    std::map<std::string, int> result() const;
    std::map<std::string, double> fpkm() const;
    void setResultWriter(cppfasta::BasicWriter *writer) {m_resultwriter = writer;}

    static void createSkipArray(uint8_t *skipArray, const char* seq, size_t len, long ngram);
    
private:
    cppfasta::BasicWriter *m_resultwriter;
    GenomeHashWithMetadata m_index;

    long m_stepsize;
    long m_acceptcount;
    long m_unmapped;
    long m_multimapped;
    std::vector<long> m_mapped;
    
    uint64_t m_nread;
    int m_ngram;
    bool m_checkall;
};

/**
 * Run mapping
 * @param dbpath Index database path
 * @param reader Input FASTQ reader
 * @param result Map of GeneID and count of reads
 * @param fpkm_result Map of GeneID and estimated FPKM
 * @param resultWriter write out mapping result. set null to disable
 * @param processedSequecnes The number of FASTQ reads
 * @param stepsize Skip size of check hash
 * @param acceptcount The nubmer of found count to accept
 * @param checkall continue checking N-gram after accepting the read
 * @param pairedend set true if input data is paired end FASTQ
 */
bool domapping(const char *dbpath, cppfasta::FastqReader2 *reader, std::map<std::string, int> *result, std::map<std::string, double> *fpkm_result,
               cppfasta::BasicWriter* resultWriter, uint64_t *processedSequences, long stepsize, long acceptcount, 
               bool checkall, bool pairedend);

#endif /* MAPPING_H */
