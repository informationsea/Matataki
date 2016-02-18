#include "dbcommon.hpp"
#include "quickcommon.hpp"
#include "tabtable.hpp"
#include <cstdio>

using namespace std;

GenomeHashWithMetadata::GenomeHashWithMetadata()
{
}

bool GenomeHashWithMetadata::save(FILE *file) const
{
    return false; // not implemented
}

bool GenomeHashWithMetadata::load(FILE *file)
{
    if (!GenomeHash::load(file)) return false;

    //DEBUG("Ready for loading");

    TabTableReader tabTableReader;
    if (!tabTableReader.open(file)) return false;

    m_geneidlist.push_back("DUMMY"); // Insert dummy data first
    m_foundFragments.push_back(0);

    long expectedIndex = 1;
    std::vector<std::string> row;
    while ((row = tabTableReader.next()).size() > 0) {
        //DEBUG("Loading list %s %s %s", row[0].c_str(), row[1].c_str(), row[2].c_str());
        bool ok;
        long index = string2long(row[0], &ok);
        if (!ok) {
            fprintf(stderr, "Failed to convert gene index number\n");
            return false;
        }
        
        if (index != expectedIndex) {
            fprintf(stderr, "Invalid Index in gene list\n");
            return false;
        }

        m_geneidlist.push_back(row[1]);
        m_gene2index[row[1]] = index;

        double foundFragments = string2double(row[2], &ok);
        if (!ok) {
            fprintf(stderr, "Cannot read # of found fragments in gene list\n");
            return false;
        }

        m_foundFragments.push_back(foundFragments);

        expectedIndex += 1;
    }

    //DEBUG("Finish loading");
    
    return true;
}
