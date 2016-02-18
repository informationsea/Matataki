#include "tabtable.hpp"
#include "filereader.hpp"
#include "filewriter.hpp"
#include "quickcommon.hpp"

TabTableReader::TabTableReader() : reader(0) {}
TabTableReader::~TabTableReader()
{
    delete reader;
}

bool TabTableReader::open(const char *filepath)
{
    reader = cppfasta::openFile(filepath);
    return reader != NULL;
}

bool TabTableReader::open(FILE *file)
{
    cppfasta::FileReader *fileReader = new cppfasta::FileReader;
    if (!fileReader->open(file)) {
        delete fileReader;
        return false;
    }
    reader = fileReader;
    return true;
}

std::vector<std::string> TabTableReader::next()
{
    std::string line = reader->readLine();
    if (line.size() == 0)
        return std::vector<std::string>();
    
    if (line.rfind('\n') == line.size()-1)
        line = line.substr(0, line.size()-1);
    if (line.rfind('\r') == line.size()-1)
        line = line.substr(0, line.size()-1);

    std::vector<std::string> row;
    std::string::size_type pos = -1;
    while (1) {
        std::string::size_type newpos = line.find('\t', pos+1);
        if (newpos == std::string::npos) {
            row.push_back(line.substr(pos+1));
            break;
        }

        row.push_back(line.substr(pos+1, newpos - pos - 1));
        pos = newpos;
    }
    return row;
}

TabTableWriter::TabTableWriter() : writer(0) {}
TabTableWriter::~TabTableWriter()
{
    delete writer;
}

bool TabTableWriter::open(const char *filepath)
{
    writer = cppfasta::openFileForWrite(filepath);
    return writer != NULL;
}

bool TabTableWriter::write(const std::vector<std::string> &row)
{
    for (std::vector<std::string>::const_iterator it = row.begin();
         it < row.end(); ++it) {
        if (it != row.begin())
            if (!writer->write("\t")) return false;
        if (!writer->write(*it)) return false;
    }
    if (!writer->write("\n")) return false;
    return true;
}
