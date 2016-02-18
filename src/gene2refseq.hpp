#ifndef GENE2REFSEQ_H
#define GENE2REFSEQ_H

#include <map>
#include <set>
#include <string>

class Gene2Refseq
{
public:
    Gene2Refseq();
    virtual ~Gene2Refseq();

    bool load(const char *filepath, const char *taxid = 0);
    const std::set<std::string> &gene2refseq(const std::string &geneid) {return m_gene2refseq[geneid];};
    const std::string &refseq2gene(const std::string &refseqid) {return m_refseq2gene[refseqid];}
    const std::string &refseq2tax(const std::string &refseqid) {return m_gene2tax[m_refseq2gene[refseqid]];};
    const std::string &gene2tax(const std::string &geneid) {return m_gene2tax[geneid];};

private:
    std::map<std::string, std::set<std::string> > m_gene2refseq;
    std::map<std::string, std::string> m_refseq2gene;
    std::map<std::string, std::string> m_gene2tax;
    
};


#endif /* GENE2REFSEQ_H */
