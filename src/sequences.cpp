
#include "common.h"
#include "logging.h"
#include "parsing.h"
#include "seq.h"
#include "sequences.h"


// TODO: add sites validation
//       - positions should be sorted and unique
//       - bases should be acceptable characters


namespace arghmm {


//=============================================================================
// input/output: FASTA

bool read_fasta(FILE *infile, Sequences *seqs)
{
    // store lines until they are ready to discard
    class Discard : public vector<char*> {
    public:
        ~Discard() { clean(); }
        void clean() {
            for (unsigned int i=0; i<size(); i++)
                delete [] at(i);
            clear();
        }
    };
    
    
    // init sequences
    seqs->clear();
    seqs->set_owned(true);

    char *line;    
    string key;
    vector<char*> seq;
    Discard discard;
    
    while ((line = fgetline(infile))) {
        chomp(line);
        
        if (line[0] == '>') {
            // parse key line

            if (seq.size() > 0) {
                // add new sequence
                char *full_seq = concat_strs(&seq[0], seq.size());
                int seqlen2 = strlen(full_seq);
                if (!seqs->append(key, full_seq, seqlen2)) {
                    printError("sequences are not the same length: %d != %d",
                               seqs->length(), seqlen2);
                    delete [] full_seq;
                    return false;
                }
                seq.clear();
                discard.clean();
            }
        
            // new key found
            key = string(&line[1]);
            delete [] line;
        } else {
            // parse sequence line

            seq.push_back(trim(line));
            discard.push_back(line);
        }
    }
    
    // add last sequence
    if (seq.size() > 0) {
        char *full_seq = concat_strs(&seq[0], seq.size());
        int seqlen2 = strlen(full_seq);
        if (!seqs->append(key, full_seq, seqlen2)) {
            printError("sequences are not the same length: %d != %d",
                       seqs->length(), seqlen2);
            delete [] full_seq;
            return false;
        }
        discard.clean();
    }
    
    return true;
}


bool read_fasta(const char *filename, Sequences *seqs)
{
    FILE *infile = NULL;    
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'", filename);
        return false;
    }
    
    bool result = read_fasta(infile, seqs);
    fclose(infile);
    
    return result;
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
    for (int i=0; i<seqs->get_num_seqs(); i++) {
        fprintf(stream, ">%s\n", seqs->names[i].c_str());
        fprintf(stream, "%s\n", seqs->seqs[i]);
    }
}



//=============================================================================
// input/output: sites file format


bool validate_site_column(char *col, int nseqs)
{
    for (int i=0; i<nseqs; i++) {
        col[i] = toupper(col[i]);
        if (col[i] != 'N' && dna2int[(int) col[i]] == -1)
            return false;
    }
    return true;
}


// Read a Sites stream
bool read_sites(FILE *infile, Sites *sites, 
                int subregion_start, int subregion_end)
{
    const char *delim = "\t";
    char *line;
    int nseqs = 0;
    
    sites->clear();
    bool error = false;
    
    // parse lines
    int lineno = 0;    
    while (!error && (line = fgetline(infile))) {
        chomp(line);
        lineno++;

        if (strncmp(line, "NAMES\t", 6) == 0) {
            // parse NAMES line
            split(&line[6], delim, sites->names);
            nseqs = sites->names.size();

        } else if (strncmp(line, "REGION\t", 7) == 0) {
            // parse RANGE line
            char chrom[51];
            if (sscanf(line, "REGION\t%50s\t%d\t%d", 
                       chrom,
                       &sites->start_coord, &sites->end_coord) != 3) {
                printError("bad REGION format");
                delete [] line;
                return false;
            }
            sites->chrom = chrom;
            sites->start_coord--;  // convert to 0-index

            // set region by subregion if specified
            if (subregion_start != -1)
                sites->start_coord = subregion_start;
            if (subregion_end != -1)
                sites->end_coord = subregion_end;
            

        } else if (strncmp(line, "RANGE\t", 6) == 0) {
            // parse RANGE line
            
            printError("deprecated RANGE line detected (use REGION instead)");
            delete [] line;
            return false;

        } else {
            // parse a site line
            
            // parse site
            int position;
            if (sscanf(line, "%d\t", &position) != 1) {
                printError("first column is not an integer (line %d)", lineno);
                delete [] line;
                return false;
            }

            // skip site if not in region            
            if (position < sites->start_coord || position >= sites->end_coord) {
                delete [] line;
                continue;
            }
            
            // skip first word
            int i=0;
            for (; line[i] != '\t'; i++) {}
            i++;
            char *col = &line[i];

            // parse bases
            unsigned int len = strlen(col);
            if (len != (unsigned int) nseqs) {
                printError("not enough bases given (line %d)", lineno);
                delete [] line;
                return false;
            }
            if (!validate_site_column(col, nseqs)) {
                printError("invalid sequence characters (line %d)", lineno);
                printError("%s\n", line);
                delete [] line;
                return false;
            }

            // convert to 0-index
            sites->append(position - 1, col, true);
        }
        
        delete [] line;
    }
    
    return true;
}


// Read a Sites alignment file
bool read_sites(const char *filename, Sites *sites, 
                int subregion_start, int subregion_end)
{
    FILE *infile;
    if ((infile = fopen(filename, "r")) == NULL) {
        printError("cannot read file '%s'", filename);
        return false;
    }
    
    bool result = read_sites(infile, sites, subregion_start, subregion_end);

    fclose(infile);
    
    return result;
}


// Converts a Sites alignment to a Sequences alignment
void make_sequences_from_sites(const Sites *sites, Sequences *sequences, 
                               char default_char)
{
    int nseqs = sites->names.size();
    int seqlen = sites->length();
    int start = sites->start_coord;
    int nsites = sites->get_num_sites();
    
    sequences->clear();
    sequences->set_owned(true);
    
    for (int i=0; i<nseqs; i++) {
        char *seq = new char [seqlen];
        
        int col = 0;
        for (int j=0; j<seqlen; j++) {
            if (col < nsites && start+j == sites->positions[col]) {
                // variant site
                seq[j] = sites->cols[col++][i];
            } else {
                seq[j] = default_char;
            }
        }

        sequences->append(sites->names[i], seq);
    }

    sequences->set_length(seqlen);
}


template<>
void apply_mask_sequences<NullValue>(Sequences *sequences, 
                                     const TrackNullValue &maskmap)
{
    const char maskchar = 'N';
    
    for (unsigned int k=0; k<maskmap.size(); k++) {
        for (int i=maskmap[k].start; i<maskmap[k].end; i++) {
            for (int j=0; j<sequences->get_num_seqs(); j++)
                sequences->seqs[j][i] = maskchar;
        }
    }
}



// Returns true if alignment column is invariant
static inline bool is_invariant_site(const char *const *seqs, 
                                     const int nseqs, const int pos)
{
    const char c = seqs[0][pos];
    for (int j=1; j<nseqs; j++) {
        if (seqs[j][pos] != c) {
            return false;
        }
    }
    return true;
}


// Converts a Sequences alignment to a Sites alignment
void make_sites_from_sequences(const Sequences *sequences, Sites *sites)
{
    int nseqs = sequences->get_num_seqs();
    int seqlen = sequences->length();
    const char * const *seqs = sequences->get_seqs();
    
    sites->clear();
    sites->start_coord = 0;
    sites->end_coord = seqlen;
    sites->names.insert(sites->names.begin(),
                        sequences->names.begin(), sequences->names.end());

    for (int i=0; i<seqlen; i++) {
        if (!is_invariant_site(seqs, nseqs, i)) {
            char *col = new char [nseqs];
            for (int j=0; j<nseqs; j++)
                col[j] = seqs[j][i];
            sites->append(i, col);
        }
    }
}



// compress the sites by a factor of 'compress'
void find_compress_cols(const Sites *sites, int compress, 
                        SitesMapping *sites_mapping)
{
    const int ncols = sites->get_num_sites();

    int blocki = 0;
    int next_block = sites->start_coord + compress;
    int half_block = compress / 2;
    
    // record old coords
    sites_mapping->init(sites);

    // special case
    if (compress == 1) {
        for (int i=sites->start_coord; i<sites->end_coord; i++) {
            sites_mapping->all_sites.push_back(i);
        }

        for (int i=0; i<ncols; i++) {
            int col = sites->positions[i];
            sites_mapping->old_sites.push_back(col);
            sites_mapping->new_sites.push_back(col - sites->start_coord);
        }

        // record new coords
        sites_mapping->new_start = 0;
        sites_mapping->new_end = sites->length();
        return;
    }
    
    // iterate through variant sites
    for (int i=0; i<ncols; i++) {
        int col = sites->positions[i];

        // find next block with variant site
        while (col >= next_block) {
            sites_mapping->all_sites.push_back(next_block - half_block);
            next_block += compress;
            blocki++;
        }

        sites_mapping->old_sites.push_back(col);
        sites_mapping->new_sites.push_back(blocki);
        sites_mapping->all_sites.push_back(col);
        next_block += compress;
        blocki++;

        // each original site should be unique
        const int n = sites_mapping->all_sites.size();
        if (n > 1)
            assert(sites_mapping->all_sites[n-1] !=
                   sites_mapping->all_sites[n-2]);
    }

    // record non-variants at end of alignment
    while (sites->end_coord >= next_block) {
        sites_mapping->all_sites.push_back(next_block - half_block);
        next_block += compress;
        blocki++;
    }

    // record new coords
    sites_mapping->new_start = 0;
    int new_end = sites->length() / compress;
    if (ncols > 0)
        sites_mapping->new_end = max(sites_mapping->new_sites[ncols-1]+1, 
                                     new_end);
    else
        sites_mapping->new_end = new_end;

    //assert(sites_mapping->all_sites.size() == sites_mapping->new_end - 
    //       sites_mapping->new_start);
}


void compress_sites(Sites *sites, const SitesMapping *sites_mapping)
{
    const int ncols = sites->cols.size();
    sites->start_coord = sites_mapping->new_start;
    sites->end_coord = sites_mapping->new_end;
    
    for (int i=0; i<ncols; i++)
        sites->positions[i] = sites_mapping->new_sites[i];
}


void uncompress_sites(Sites *sites, const SitesMapping *sites_mapping)
{
    const int ncols = sites->cols.size();
    sites->start_coord = sites_mapping->old_start;
    sites->end_coord = sites_mapping->old_end;
    
    for (int i=0; i<ncols; i++)
        sites->positions[i] = sites_mapping->old_sites[i];
}



//=============================================================================
// assert functions

bool check_sequences(Sequences *seqs)
{
    return check_sequences(
        seqs->get_num_seqs(), seqs->length(), seqs->get_seqs()) &&
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


//=============================================================================
// Misc

void resample_align(Sequences *aln, Sequences *aln2)
{
    assert(aln->get_num_seqs() == aln2->get_num_seqs());
    char **seqs = aln->get_seqs();
    char **seqs2 = aln2->get_seqs();

    for (int j=0; j<aln2->length(); j++) {
        // randomly choose a column (with replacement)
        int col = irand(aln->length());
        
        // copy column
        for (int i=0; i<aln2->get_num_seqs(); i++) {
            seqs2[i][j] = seqs[i][col];
        }
    }
}

} // namespace arghmm
