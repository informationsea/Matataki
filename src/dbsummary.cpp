#include "dbsummary.hpp"
#include "dbcommon.hpp"
#include "quickcommon.hpp"
#include <cstdio>

using namespace std;

bool dbsummary(const char *dbpath, std::map<std::string, int> *gene2num, int *fragments, int *ngram)
{
    FILE *db = fopen(dbpath, "r");
    if (db == NULL) {
        perror("Cannot open database");
        return false;
    }
    
    GenomeHashWithMetadata genomeHash;
    if (!genomeHash.load(db)) {
        fprintf(stderr, "Cannot open database");
        return false;
    }

    *ngram = genomeHash.sequenceLength();
    *fragments = 0;

    GenomeHashCursor cursor(&genomeHash);
    char buf[128];
    uint32_t value;
    while ((value = cursor.next(buf, sizeof(buf))) != GENOMEHASH_FAILED_ID) {
        (*gene2num)[genomeHash.geneIdForIndex(value)] += 1;
        *fragments += 1;
    }
    
    return true;
}
