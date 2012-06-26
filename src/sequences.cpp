
#include "common.h"
#include "logging.h"
#include "parsing.h"
#include "seq.h"
#include "sequences.h"



namespace arghmm {


Sequences *read_fasta(const char *filename)
{
    // store lines until they are ready to discard
    class Discard : public vector<char*> {
    public:
        void clean() {
            for (unsigned int i=0; i<size(); i++)
                delete [] at(i);
            clear();
        }
    };


    FILE *infile = NULL;
    
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'", filename);
        return NULL;
    }
    
    char *line;
    
    Sequences *seqs = new Sequences();
    seqs->set_owned(true);
    string key;
    vector<char*> seq;
    Discard discard;
    
    while ((line = fgetline(infile))) {
        chomp(line);
        
        if (line[0] == '>') {
            if (seq.size() > 0) {  
                // add new sequence
                seqs->append(concat_strs(&seq[0], seq.size()));
                seq.clear();
                discard.clean();
            }
        
            // new key found
            key = string(&line[1]);
            delete [] line;
        } else {
            seq.push_back(trim(line));
            discard.push_back(line);
        }
    }

    fclose(infile);
    
    // add last sequence
    if (seq.size() > 0) {
        seqs->append(concat_strs(&seq[0], seq.size()));
        discard.clean();
    }

    // set sequence length
    if (seqs->set_length() < 0) {
        printError("sequences are not the same length '%s'", filename);
        return NULL;
    }
    
    return seqs;
}


bool write_fasta(const char *filename, Sequences *seqs)
{
    FILE *stream = NULL;
    
    if ((stream = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "cannot open '%s'\n", filename);
        return false;
    }

    write_fasta(stream, seqs);
    
    fclose(stream);
    return true;
}

void write_fasta(FILE *stream, Sequences *seqs)
{
    for (int i=0; i<seqs->get_nseqs(); i++) {
        fprintf(stream, ">%s\n", seqs->names[i].c_str());
        fprintf(stream, "%s\n", seqs->seqs[i]);
    }
}


bool check_sequences(Sequences *seqs)
{
    return check_sequences(
        seqs->get_nseqs(), seqs->length(), seqs->get_seqs()) &&
        check_seq_names(seqs);
}

// ensures that all characters in the alignment are sensible
// TODO: do not change alignment (keep Ns)
bool check_sequences(int nseqs, int seqlen, char **seqs)
{
    // check seqs
    // CHANGE N's to gaps
    for (int i=0; i<nseqs; i++) {
        for (int j=0; j<seqlen; j++) {
            char x = seqs[i][j];
            if (strchr("NnRrYyWwSsKkMmBbDdHhVv", x))
                // treat Ns as gaps
                x = '-';
            if (x != '-' &&
                dna2int[(int) (unsigned char) x] == -1)
            {
                // an unknown character is in the alignment
                printError("unknown char '%c' (char code %d)\n", x, x);
                return false;
            }
        }
    }
    
    return true;
}


bool check_seq_names(Sequences *seqs)
{
    for (unsigned int i=0; i<seqs->names.size(); i++) {
        if (!check_seq_name(seqs->names[i].c_str())) {
            printError("sequence name has illegal characters '%s'",
                       seqs->names[i].c_str());
            return false;
        }
    }

    return true;
}


//
// A valid gene name and species name follows these rules:
//
// 1. the first and last characters of the ID are a-z A-Z 0-9 _ - .
// 2. the middle characters can be a-z A-Z 0-9 _ - . or the space character ' '.
// 3. the ID should not be purely numerical characters 0-9
// 4. the ID should be unique within a gene tree or within a species tree
//

bool check_seq_name(const char *name)
{
    int len = strlen(name);
    
    if (len == 0) {
        printError("name is zero length");
        return false;
    }

    // check rule 1
    if (name[0] == ' ' || name[len-1] == ' ') {
        printError("name starts or ends with a space '%c'");
        return false;
    }

    // check rule 2
    for (int i=0; i<len; i++) {
        char c = name[i];
        if (!((c >= 'a' && c <= 'z') ||
              (c >= 'A' && c <= 'Z') ||
              (c >= '0' && c <= '9') ||
              c == '_' || c == '-' || c == '.' || c == ' ')) {
            printError("name contains illegal character '%c'", c);
            return false;
        }
    }

    // check rule 3
    int containsAlpha = false;
    for (int i=0; i<len; i++) {
        if (name[i] < '0' || name[i] > '9') {
            containsAlpha = true;
            break;
        }
    }
    if (!containsAlpha) {
        printError("name is purely numeric '%s'", name);
        return false;
    }

    return true;
}


void resample_align(Sequences *aln, Sequences *aln2)
{
    assert(aln->get_nseqs() == aln2->get_nseqs());
    char **seqs = aln->get_seqs();
    char **seqs2 = aln2->get_seqs();

    for (int j=0; j<aln2->length(); j++) {
        // randomly choose a column (with replacement)
        int col = irand(aln->length());
        
        // copy column
        for (int i=0; i<aln2->get_nseqs(); i++) {
            seqs2[i][j] = seqs[i][col];
        }
    }
}

} // namespace arghmm
