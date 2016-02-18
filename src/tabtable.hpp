#ifndef TABTABLE_H
#define TABTABLE_H

#include <stdio.h>
#include <string>
#include <vector>

namespace cppfasta{
    class BasicReader;
    class BasicWriter;
}

class TabTableReader
{
public:
    TabTableReader();
    virtual ~TabTableReader();

    virtual bool open(const char *filepath);
    virtual bool open(FILE *file);
    virtual std::vector<std::string> next();
    
private:
    cppfasta::BasicReader *reader;
};

class TabTableWriter
{
public:
    TabTableWriter();
    virtual ~TabTableWriter();

    virtual bool open(const char *filepath);
    virtual bool write(const std::vector<std::string> &row);
private:
    cppfasta::BasicWriter *writer;
};



#endif /* TABTABLE_H */
