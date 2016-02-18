#include "gene2refseq.hpp"
#include "tabtable.hpp"

Gene2Refseq::Gene2Refseq() {}
Gene2Refseq::~Gene2Refseq() {}

bool Gene2Refseq::load(const char *filepath, const char *taxid)
{
    TabTableReader reader;
    if (!reader.open(filepath))
        return false;
    
    std::vector<std::string> row;
    do {
        row = reader.next();
        if (row.size() == 0) break;
        if (row.size() < 4) continue;

        if (taxid && (taxid != row[0])) continue;

        m_gene2tax[row[1]] = row[0];
        m_refseq2gene[row[3]] = row[1];
        m_gene2refseq[row[1]].insert(row[3]);
        
    } while(1);
    return true;
}
