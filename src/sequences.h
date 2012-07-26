/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2012

  Common sequence functions

=============================================================================*/


#ifndef ARGHMM_SEQUENCES_H
#define ARGHMM_SEQUENCES_H

// c++ includes
#include <string>
#include <vector>


namespace arghmm {

using namespace std;


// The alignment of sequences
class Sequences
{
public:
    Sequences() : owned(false) {}

    Sequences(char **_seqs, int nseqs, int seqlen) :
        seqlen(seqlen), owned(false)
    {
        extend(_seqs, nseqs);
    }

    // initialize from a subset of another Sequences alignment
    Sequences(Sequences *sequences, int nseqs=-1, int _seqlen=-1, 
              int offset=0) :
        seqlen(_seqlen), owned(false)
    {
        // use same nseqs and/or seqlen by default
        if (nseqs == -1)
            nseqs = sequences->get_nseqs();
        if (seqlen == -1)
            seqlen = sequences->length();
        
        for (int i=0; i<nseqs; i++)
            append(sequences->names[i], &sequences->seqs[i][offset]);
    }

    ~Sequences()
    {
        clear();
    }
    
    inline int get_nseqs() const
    {
        return seqs.size();
    }

    inline int length() const
    {
        return seqlen;
    }

    inline int set_length()
    {
        seqlen = -1;
        for (unsigned int i=0; i<seqs.size(); i++) {
            int len = strlen(seqs[0]);
            if (seqlen == -1)
                seqlen = len;
            else if (seqlen != len) {
                // error, sequences are not all the same length
                seqlen = -1;
                break;
            }
        }

        return seqlen;
    }

    inline char **get_seqs()
    {
        return &seqs[0];
    }


    void set_owned(bool _owned)
    {
        owned = _owned;
    }

    void extend(char **_seqs, int nseqs)
    {
        for (int i=0; i<nseqs; i++) {
            seqs.push_back(_seqs[i]);
            names.push_back("");
        }
    }

    void extend(char **_seqs, char **_names, int nseqs)
    {
        for (int i=0; i<nseqs; i++) {
            seqs.push_back(_seqs[i]);
            names.push_back(_names[i]);
        }
    }

    void append(char *seq)
    {
        seqs.push_back(seq);
        names.push_back("");
    }

    void append(string name, char *seq)
    {
        seqs.push_back(seq);
        names.push_back(name);
    }

    void clear()
    {
        if (owned) {
            const int nseqs = get_nseqs();
            for (int i=0; i<nseqs; i++)
                delete [] seqs[i];
        }
        seqs.clear();
        names.clear();
    }


    vector <char*> seqs;
    vector <string> names;

protected:
    int seqlen;
    bool owned;
};


class Sites
{
public:
    Sites(int start_coord=0, int end_coord=0) :
        start_coord(start_coord),
        end_coord(end_coord)
    {}
    ~Sites()
    {
        clear();
    }

    void append(int position, char* col, bool copy=false) {
        positions.push_back(position);
        if (copy) {
            unsigned int len = strlen(col);
            assert(len == names.size());
            char *col2 = new char [len];
            strcpy(col2, col);
            cols.push_back(col2);
        } else {
            cols.push_back(col);
        }
    }

    void clear()
    {
        for (unsigned int i=0; i<cols.size(); i++)
            delete [] cols[i];
        names.clear();
        positions.clear();
        cols.clear();
    }

    inline int length() const
    {
        return end_coord - start_coord;
    }
    
    inline int get_num_sites() const
    {
        return positions.size();
    }
    
    
    int start_coord;
    int end_coord;
    vector<string> names;
    vector<int> positions;
    vector<char*> cols;
};



class SitesMapping 
{
public:
    SitesMapping() {}
    ~SitesMapping() {}

    void init(const Sites *sites)
    {
        old_start = sites->start_coord;
        old_end = sites->end_coord;
        nsites = sites->get_num_sites();
        seqlen = sites->length();
    }

    int old_start;
    int old_end;
    int new_start;
    int new_end;
    int nsites;
    int seqlen;

    vector<int> old_sites;
    vector<int> new_sites;
    vector<int> all_sites;
};


// sequences functions
Sequences *read_fasta(FILE *infile);
Sequences *read_fasta(const char *filename);
bool write_fasta(const char *filename, Sequences *seqs);
void write_fasta(FILE *stream, Sequences *seqs);

bool check_sequences(Sequences *seqs);
bool check_sequences(int nseqs, int seqlen, char **seqs);
bool check_seq_names(Sequences *seqs);
bool check_seq_name(const char *name);
void resample_align(Sequences *aln, Sequences *aln2);

// sites functions
Sites *read_sites(FILE *infile);
Sites *read_sites(const char *filename);

Sequences *make_sequences_from_sites(Sites *sites, char default_char='A');

void find_compress_cols(const Sites *sites, int compress, 
                        SitesMapping *sites_mapping);
void compress_sites(Sites *sites, const SitesMapping *sites_mapping);


} // namespace arghmm

#endif // ARGHMM_SEQUENCES_H
