#include "dbcoverage.hpp"
#include <map>

#include "gene2refseq.hpp"
#include "cppfasta.hpp"
#include "dbcommon.hpp"
#include "quickcommon.hpp"
#include "filewriter.hpp"

#define DEBUG(fmt,...) (fprintf(stderr, __FILE__ ": %3d: " fmt "\n" ,__LINE__, ## __VA_ARGS__))

TranscriptCoverage::TranscriptCoverage() :
    m_refseq_id(), m_gene_id(), m_seq(), m_cover(), m_fragment(0), m_ngram(0){}

TranscriptCoverage::TranscriptCoverage(const std::string &refseq_id, const std::string &gene_id, const std::string &sequence) :
    m_refseq_id(refseq_id), m_gene_id(gene_id), m_seq(sequence), m_cover(sequence.size(), UNIQUE_UNCOVERED), m_fragment(0), m_ngram(0)
{
    
}

TranscriptCoverage::TranscriptCoverage(const std::string &refseq_id, const std::string &gene_id, const std::string &sequence,
                   const std::string &cover, int fragment, int ngram) :
    m_refseq_id(refseq_id), m_gene_id(gene_id), m_seq(sequence), m_cover(cover), m_fragment(fragment), m_ngram(ngram)
{
    
}

TranscriptCoverage::~TranscriptCoverage()
{

}

bool TranscriptCoverage::addNoCommonFragment(const std::string fragment)
{
    std::string::size_type pos = 0;
    while ((pos = m_seq.find(fragment, pos)) != std::string::npos) {
        //DEBUG("Found %s at %lu %c", fragment.c_str(), pos, m_cover[pos]);
        size_t ngram = fragment.size();

        if (m_cover[pos] == UNIQUE_UNCOVERED || m_cover[pos] == UNIQUE_NO_COMMON_COVERED)
            m_cover[pos] = UNIQUE_NO_COMMON_START;
        for (size_t i = 1; i < ngram; ++i) {
            if (m_cover[pos + i] == UNIQUE_UNCOVERED)
                m_cover[pos + i] = UNIQUE_NO_COMMON_COVERED;
        }
        pos += 1;
    }
    return true;
}

bool TranscriptCoverage::addFragment(const std::string fragment)
{
    //DEBUG("Add Fragment %s", fragment.c_str());
    std::string::size_type pos = 0;
    while ((pos = m_seq.find(fragment, pos)) != std::string::npos) {
        //DEBUG("Found %s at %lu %c", fragment.c_str(), pos, m_cover[pos]);
        size_t ngram = fragment.size();
        
        m_cover[pos] = UNIQUE_START;
        for (size_t i = 1; i < ngram; ++i) {
            if (m_cover[pos + i] != UNIQUE_START)
                m_cover[pos + i] = UNIQUE_COVERED;
        }
        pos += 1;
    }

    m_fragment += 1;
    m_ngram = (int)fragment.size();
    
    return true;
}

int TranscriptCoverage::foundFragments() const
{
    int covered = 0;
    for (std::string::const_iterator it = m_cover.begin(); it != m_cover.end(); ++it) {
        if (*it == UNIQUE_START) covered += 1;
    }
    return covered;
}

int TranscriptCoverage::coveredLength() const
{
    int covered = 0;
    for (std::string::const_iterator it = m_cover.begin(); it != m_cover.end(); ++it) {
        if (*it == UNIQUE_START || *it == UNIQUE_COVERED) covered += 1;
    }
    return covered;
}

double TranscriptCoverage::coverage() const
{
    int covered = foundFragments();
    return ((double) covered)/(m_seq.size() - m_ngram + 1);
}

bool exportCoverage(const char *export_path, const std::map<std::string, TranscriptCoverage> *coverages)
{
    cppfasta::BasicWriter *writer;
    if ((writer = cppfasta::openFileForWrite(export_path)) == NULL) {
        fprintf(stderr, "Failed to write data\n");
        return false;
    }

    writer->write("GeneID\tRefSeqID\tLength\t# of covered base\t# of unique fragment\t# of found fragment\tcoverage\tcover\n");

    for (std::map<std::string, TranscriptCoverage>::const_iterator it = coverages->begin();
         it != coverages->end(); ++it) {

        writer->printf("%s\t%s\t%zu\t%d\t%d\t%d\t%.12lf\t",
                       it->second.geneId().c_str(),
                       it->second.refseqId().c_str(),
                       it->second.cover().size(),
                       it->second.coveredLength(),
                       it->second.fragment(),
                       it->second.foundFragments(),
                       it->second.coverage()
            );
        writer->write(it->second.cover());
        writer->write("\n");
    }

    delete writer;
    return true;
}
