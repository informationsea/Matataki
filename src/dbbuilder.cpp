#include "dbbuilder.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <cppfasta.hpp>

#include "gene2refseq.hpp"
#include "cppfasta.hpp"
#include "quickcommon.hpp"
#include "genomehash.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#include <kcpolydb.h>
#pragma clang diagnostic pop
#pragma GCC diagnostic pop

#define DEBUG(fmt,...) (fprintf(stderr, __FILE__ ": %3d: " fmt "\n" ,__LINE__, ## __VA_ARGS__))
#define NOT_UNIQUE_ID "@"

DBBuilder::DBBuilder(int ngram, bool onmemory) :
    m_ngram(ngram), m_ready(false), m_onmemory(onmemory),
    m_uniqueMapPath(0), m_unique2isoformPath(0), m_next_sequence_id(0),
    m_uniquemap(new kyotocabinet::PolyDB), m_unique2isoform(new kyotocabinet::PolyDB)
{

    if (onmemory) {
        if (!m_uniquemap->open())
            return;
        if (!m_unique2isoform->open())
            return;
    } else {
        m_uniqueMapPath = new TemporaryFile("QE_DBBUILDER_UM", ".kch");
        if (!m_uniqueMapPath->ready()) {
            DEBUG("Failed to create temporary file");
            return;
        }
        
        if (!m_uniquemap->open(m_uniqueMapPath->path())) {
            DEBUG("Failed to open kyoto cabinet");
            return;
        }

        m_unique2isoformPath = new TemporaryFile("QE_DBBUILDER_U2S", ".kch");
        if (!m_unique2isoformPath->ready()) {
            DEBUG("Failed to create temporary file");
            return;
        }
    
        if (!m_unique2isoform->open(m_unique2isoformPath->path())) {
            DEBUG("Failed to open kyoto cabinet");
            return;
        }
    }

    m_ready = true;
}

DBBuilder::~DBBuilder()
{
    if (m_uniquemap) {
        m_uniquemap->close();
        delete m_uniquemap;
    }
    delete m_uniqueMapPath;

    if (m_unique2isoform) {
        m_unique2isoform->close();
        delete m_unique2isoform;
    }
    delete m_unique2isoformPath;
}

bool DBBuilder::addSequence(const char *geneid, const cppfasta::SequenceRecord *sequence)
{
    if (!m_ready) return false;
    //DEBUG("ADD SEQUENCE  : %s %s", geneid, sequence->name().c_str());
    size_t sequenceLength = sequence->sequence().size();

    if ((int)sequenceLength < m_ngram) {
        fprintf(stderr, "Length of sequence is smaller than N-gram. %s\n", sequence->name().c_str());
        return true;
    }
    
    char sequence_id[100];
    snprintf(sequence_id, sizeof(sequence_id) - 1, "%llu", m_next_sequence_id);
    avaliable_gene2refseq[geneid].push_back(sequence_id);
    m_geneset.insert(geneid);
    m_transcriptome_coverage[sequence_id] = TranscriptCoverage(sequence->name(), geneid, sequence->sequence());
    m_next_sequence_id += 1;

    //DEBUG("LENGTH : %zu  NGRAM : %d", sequenceLength, m_ngram);
    for (size_t i = 0; i < sequenceLength  - m_ngram + 1; i++) {
        std::string keys[2];
        keys[0] = sequence->sequence().substr(i, m_ngram);
        keys[1] = reverseComplement(keys[0]);
        //DEBUG("TESTING<%s -> %s> %s %s %d %d", geneid, sequence_id, keys[0].c_str(), keys[1].c_str(), m_uniquemap->check(keys[0]), m_uniquemap->check(keys[1]));

        for (size_t j = 0; j < 2; ++j) {
            if (m_uniquemap->check(keys[j]) > 0) {
                std::string value;
                if (!m_uniquemap->get(keys[j], &value)) return false;
                if (value == NOT_UNIQUE_ID)
                    continue;
                if (value != geneid) {
                    m_uniquemap->replace(keys[j], NOT_UNIQUE_ID);
                    if (!m_unique2isoform->remove(keys[j])) {
                        DEBUG("Failed on deleting data from kyoto cabinet");
                        return false;
                    }
                } else {
                    if (!m_unique2isoform->append(keys[j], std::string(sequence_id) + ",")) {
                        DEBUG("Failed on adding data from kyoto cabinet");
                        return false;
                    }
                }
            } else {
                if (!m_uniquemap->add(keys[j], geneid))
                    return false;
                if (!m_unique2isoform->add(keys[j], std::string(sequence_id) + ","))
                    return false;
            }
        }
    }
    return true;
}

bool DBBuilder::build(const char *filepath)
{
    if (!m_ready) return false;
    kyotocabinet::DB::Cursor *cursor = m_uniquemap->cursor();
    //frozenhashmap::FrozenMapBuilder builder;
    //if (!builder.open()) return false;
    if (!cursor->jump()) return false;

    // list up genes
    std::vector<std::string> genelist;
    genelist.push_back("UNMAPPED");
    genelist.insert(genelist.end(), m_geneset.begin(), m_geneset.end());
    std::map<std::string, int> geneindex;
    {
        int i = 0;
        for (std::vector<std::string>::iterator it = genelist.begin();
            it != genelist.end(); ++it, ++i) {
            geneindex[*it] = i;
        }
    }
    

    // check whether a N-gram available commonly among isoforms
    fprintf(stderr, "checking whether a N-gram available commonly among isoforms...\n");
    size_t numberOfNgrams = 0;
    while (1) {
        std::string key, value;
        if (!cursor->get(&key, &value)) break;
        //DEBUG("%s : %s", key.c_str(), value.c_str());
        if (value == NOT_UNIQUE_ID) {
            cursor->remove();
            continue;
        }

        std::string fragments;
        if (!m_unique2isoform->get(key, &fragments)) {
            return false;
        }

        std::vector<std::string> fragmentvector = stringsplit(fragments, ",");
        std::set<std::string> fragmentset = std::set<std::string>(fragmentvector.begin(), fragmentvector.end()-1); // ignore last element
        std::set<std::string> available_ids = std::set<std::string>(avaliable_gene2refseq[value].begin(), avaliable_gene2refseq[value].end());
        
        if (fragmentset != available_ids) {
            for (std::vector<std::string>::const_iterator it = avaliable_gene2refseq[value].begin();
                 it != avaliable_gene2refseq[value].end(); ++it) {
                m_transcriptome_coverage[*it].addNoCommonFragment(key);
            }

            cursor->remove();
            continue;
        }

        for (std::vector<std::string>::const_iterator it = avaliable_gene2refseq[value].begin();
             it != avaliable_gene2refseq[value].end(); ++it) {
            m_transcriptome_coverage[*it].addFragment(key);
        }

        //builder.put(key, value);
        
        numberOfNgrams += 1;

        if (!cursor->step())
            break;
    }

    //DEBUG("preparing genome hash %zu -> %llu", m_uniquemap->count(), bestHashSizeForExpectedEntries(m_uniquemap->count()));
  
    GenomeHash genomeHash;
    if (!genomeHash.open(bestHashSizeForExpectedEntries(m_uniquemap->count()), m_ngram)) return false;

    if (!cursor->jump()) return false;
    while (1) {
        std::string key, value;
        if (!cursor->get(&key, &value)) break;
        genomeHash.insert(key.c_str(), geneindex[value]);
        if (!cursor->step()) break;
    }

    FILE *file = fopen(filepath, "w");
    if (file == NULL) {
        perror("Failed to open file");
    }
        
    genomeHash.save(file);

    // write out coverage
    fprintf(stderr, "Writing coverage\n");
    {
        size_t i = 1;
        for (std::vector<std::string>::const_iterator it = genelist.begin() + 1; // skip dummy
             it != genelist.end(); ++it, ++i) {

            std::vector<std::string> &refseqlist = avaliable_gene2refseq[*it];
            
            double foundFragment = 0;
            for (std::vector<std::string>::const_iterator it2 = refseqlist.begin();
                 it2 != refseqlist.end(); ++it2) {
                foundFragment += m_transcriptome_coverage[*it2].foundFragments();
            }

            fprintf(file, "%zu\t%s\t%.22lf\n", i, it->c_str(), foundFragment/refseqlist.size());
        }
    }

    fclose(file);

    return true;
}

bool buildDatabase(const char *dbpath, const char *coveragepath, const char *refseqpath,
                   const char *gene2refseqpath, int ngram, bool onmemory)
{
    DBBuilder builder(ngram, onmemory);
    Gene2Refseq gene2refseq;
    cppfasta::FastaReader fasta;

    if (!gene2refseq.load(gene2refseqpath)) return false;
    if (!fasta.open(refseqpath)) return false;

    fprintf(stderr, "Loading FASTA...\n");
    cppfasta::SequenceRecord *record;
    while ((record = fasta.nextRecord()) != NULL) {
        std::string refseqid = stringsplit(record->name(), "|")[3];
        const std::string &geneid = gene2refseq.refseq2gene(refseqid);
        bool ok = true;
        if (geneid.size() == 0) {
            fprintf(stderr, "WARNING!: Unknown RefSeq ID: %s\n", refseqid.c_str());
        } else {
            ok = builder.addSequence(geneid.c_str(), record);
        }
        delete record;
        if (!ok) return false;
    }

    if (!builder.build(dbpath))
        return false;
    
    if (coveragepath)
        exportCoverage(coveragepath, &builder.transcriptomeCoverage());
    return true;
}

static char reverseCompelemntArray[256];

static bool initializeReverseCompelemntArray()
{
    memset(reverseCompelemntArray, 'N', sizeof(reverseCompelemntArray));
    reverseCompelemntArray[(int)'A'] = 'T';
    reverseCompelemntArray[(int)'T'] = 'A';
    reverseCompelemntArray[(int)'C'] = 'G';
    reverseCompelemntArray[(int)'G'] = 'C';
    return true;
}

static bool dummy1 = initializeReverseCompelemntArray();

bool reverseComplement(char *dest, size_t destsize, const char *src, size_t srcsize)
{
    if (destsize < srcsize) return false; // small destination buffer
    
    for (size_t i = 0; i < srcsize; i++) {
        dest[i] = reverseCompelemntArray[(size_t)src[srcsize - i - 1]];
    }
    return true;
}

std::string reverseComplement(std::string seq)
{
    std::string complement(seq.size(), 'N');

    size_t index = 0;
    for (std::string::const_reverse_iterator it = seq.rbegin();
         it != seq.rend(); ++it) {
        complement[index] = reverseCompelemntArray[(size_t) *it];
        index += 1;
    }
    return complement;
}
