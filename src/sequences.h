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
    explicit Sequences(int seqlen=0) :
        seqlen(seqlen), owned(false)
    {}

    Sequences(char **_seqs, int nseqs, int seqlen) :
        seqlen(seqlen), owned(false)
    {
        extend(_seqs, nseqs);
    }

    // initialize from a subset of another Sequences alignment
    Sequences(const Sequences *sequences, int nseqs=-1, int _seqlen=-1, 
              int offset=0) :
        seqlen(_seqlen), owned(false)
    {
        // use same nseqs and/or seqlen by default
        if (nseqs == -1)
            nseqs = sequences->get_num_seqs();
        if (seqlen == -1)
            seqlen = sequences->length();
        
        for (int i=0; i<nseqs; i++)
            append(sequences->names[i], &sequences->seqs[i][offset]);
    }

    ~Sequences()
    {
        clear();
    }
    
    inline int get_num_seqs() const
    {
        return seqs.size();
    }

    inline int length() const
    {
        return seqlen;
    }

    inline void set_length(int _seqlen)
    {
        seqlen = _seqlen;
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
    
    bool append(string name, char *seq, int new_seqlen=-1)
    {
        // check sequence length
        if (new_seqlen > 0) {
            if (seqs.size() > 0) {
                if (new_seqlen != seqlen)
                    return false;
            } else 
                seqlen = new_seqlen;
        }
        
        seqs.push_back(seq);
        names.push_back(name);
        return true;
    }

    void clear()
    {
        if (owned) {
            const int nseqs = get_num_seqs();
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
    Sites(string chrom="", int start_coord=0, int end_coord=0) :
        chrom(chrom),
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
            char *col2 = new char [len + 1];
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

    inline int get_num_seqs() const
    {
        return names.size();
    }
    
    
    string chrom;
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

    int compress(int pos) {
        const int n = all_sites.size();
        for (int pos2 = 0; pos2<n; pos2++) {
            if (all_sites[pos2] > pos)
                return pos2;
        }
        return n - 1;
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
bool read_fasta(FILE *infile, Sequences *seqs);
bool read_fasta(const char *filename, Sequences *seqs);
bool write_fasta(const char *filename, Sequences *seqs);
void write_fasta(FILE *stream, Sequences *seqs);

bool check_sequences(Sequences *seqs);
bool check_sequences(int nseqs, int seqlen, char **seqs);
bool check_seq_names(Sequences *seqs);
bool check_seq_name(const char *name);
void resample_align(Sequences *aln, Sequences *aln2);

// sites functions
bool read_sites(FILE *infile, Sites *sites, 
                 int subregion_start=-1, int subregion_end=-1);
bool read_sites(const char *filename, Sites *sites,
                 int subregion_start=-1, int subregion_end=-1);

void make_sequences_from_sites(const Sites *sites, Sequences *sequencess, 
                               char default_char='A');

// sequence compression
void find_compress_cols(const Sites *sites, int compress, 
                        SitesMapping *sites_mapping);
void compress_sites(Sites *sites, const SitesMapping *sites_mapping);


} // namespace arghmm

#endif // ARGHMM_SEQUENCES_H
