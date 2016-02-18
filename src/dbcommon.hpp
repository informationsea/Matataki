#ifndef DBCOMMON_H
#define DBCOMMON_H

#include "genomehash.hpp"
#include <vector>
#include <string>
#include <map>

class GenomeHashWithMetadata : public GenomeHash
{
public:
    GenomeHashWithMetadata();
    virtual ~GenomeHashWithMetadata() {}

    virtual bool load(FILE *file); // load from file
    virtual bool save(FILE *file) const; // write to file

    std::string geneIdForIndex(size_t i) const {return m_geneidlist[i];}
    double foundFragmentsForIndex(size_t i) const {return m_foundFragments[i];}
    size_t indexForGeneID(std::string id) const {return m_gene2index.at(id);}

    const std::map<std::string, size_t>& gene2index() const {return m_gene2index;}

private:
    std::map<std::string, size_t> m_gene2index;
    std::vector<std::string> m_geneidlist;
    std::vector<double> m_foundFragments;
};

#endif /* DBCOMMON_H */
