#ifndef DBCOVERAGE_H
#define DBCOVERAGE_H

#include <string>
#include <map>

const static char UNIQUE_START = 'X';
const static char UNIQUE_COVERED = '*';
const static char UNIQUE_UNCOVERED = '.';
const static char UNIQUE_NO_COMMON_START = '_';
const static char UNIQUE_NO_COMMON_COVERED = ',';


class TranscriptCoverage
{
public:
    TranscriptCoverage();
    TranscriptCoverage(const std::string &refseq_id, const std::string &gene_id, const std::string &sequence);
    TranscriptCoverage(const std::string &refseq_id, const std::string &gene_id, const std::string &sequence, const std::string &cover, int fragment, int ngram);
    virtual ~TranscriptCoverage();

    const std::string &sequence() const {return m_seq;}
    const std::string &geneId() const {return m_gene_id;}
    const std::string &refseqId() const {return m_refseq_id;}
    const std::string &cover() const {return m_cover;}
    
    int coveredLength() const;
    int foundFragments() const;
    int fragment() const {return m_fragment;}
    double coverage() const;

    bool addNoCommonFragment(const std::string fragment);
    bool addFragment(const std::string fragment);
    
private:
    std::string m_refseq_id;
    std::string m_gene_id;
    std::string m_seq;
    std::string m_cover;
    int m_fragment;
    int m_ngram;
};

typedef void (*callbackTranscriptCoverage)(const TranscriptCoverage &coverage, void *userobj);

//bool dbcoverage(const char *dbpath, const char *refseq_path, const char* gene2refseq_path, std::map<std::string, TranscriptCoverage> *coverages);
bool exportCoverage(const char *export_path, const std::map<std::string, TranscriptCoverage> *coverages);

#endif /* DBCOVERAGE_H */
