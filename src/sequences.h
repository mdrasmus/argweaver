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
    Sequences(Sequences *sequences, int nseqs=-1, int seqlen=-1) :
        seqlen(seqlen), owned(false)
    {
        // use same nseqs and/or seqlen by default
        if (nseqs == -1)
            nseqs = sequences->get_nseqs();
        if (seqlen == -1)
            seqlen = sequences->length();
        
        for (int i=0; i<nseqs; i++)
            seqs.push_back(sequences->seqs[i]);
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


Sequences *read_fasta(FILE *infile);
Sequences *read_fasta(const char *filename);
bool write_fasta(const char *filename, Sequences *seqs);
void write_fasta(FILE *stream, Sequences *seqs);
bool check_sequences(Sequences *seqs);
bool check_sequences(int nseqs, int seqlen, char **seqs);
bool check_seq_names(Sequences *seqs);
bool check_seq_name(const char *name);
void resample_align(Sequences *aln, Sequences *aln2);


} // namespace arghmm

#endif // ARGHMM_SEQUENCES_H
