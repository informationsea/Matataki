#ifndef DBBUILDER_H
#define DBBUILDER_H

#include <limits.h>
#include <map>
#include <string>
#include <set>
#include <vector>
#include "dbcoverage.hpp"

namespace cppfasta {
    class SequenceRecord;
}

namespace kyotocabinet {
    class PolyDB;
}

class TemporaryFile;

class DBBuilder
{
public:
    DBBuilder(int ngram, bool onmemory);
    virtual ~DBBuilder();

    bool isReady() const {return m_ready;}
    bool addSequence(const char *geneid, const cppfasta::SequenceRecord *sequence);
    bool build(const char *filepath);
    const std::map<std::string, TranscriptCoverage>& transcriptomeCoverage() const {
        return m_transcriptome_coverage;
    }
    
private:
    int m_ngram;
    bool m_ready;
    bool m_onmemory;

    TemporaryFile *m_uniqueMapPath;
    TemporaryFile *m_unique2isoformPath;
    long long m_next_sequence_id;
    kyotocabinet::PolyDB *m_uniquemap;
    kyotocabinet::PolyDB *m_unique2isoform;
    std::set<std::string> m_geneset;
    std::map<std::string, std::vector<std::string> > avaliable_gene2refseq;
    std::map<std::string, TranscriptCoverage> m_transcriptome_coverage;
};

bool buildDatabase(const char *dbpath, const char *coveragepath, const char *refseqpath,
                   const char *gene2refseqpath, int ngram, bool onmemory);

bool reverseComplement(char *dest, size_t destsize, const char *src, size_t srcsize);
std::string reverseComplement(std::string seq);

#endif /* DBBUILDER_H */
