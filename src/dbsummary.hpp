#ifndef DBSUMMARY_H
#define DBSUMMARY_H

#include <string>
#include <map>

bool dbsummary(const char *dbpath, std::map<std::string, int> *gene2num, int *fragments, int *ngram);

#endif /* DBSUMMARY_H */
