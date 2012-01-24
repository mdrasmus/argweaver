/*=============================================================================

  Matt Rasmussen
  Copyright 2007-2012

  Common sequence functions

=============================================================================*/


#ifndef SPIDIR_SEQ_H
#define SPIDIR_SEQ_H


namespace spidir {


// convert dna characters into standard numbers
extern const int dna2int[256];

// convert standard numbers to dna characters
extern const char *int2dna;

// base numbers
enum {
    DNA_A = 0,
    DNA_C = 1,
    DNA_G = 2,
    DNA_T = 3,
    DNA_PURINE,
    DNA_PRYMIDINE
};


// get the base type of a nucleotide
extern int dnatype[];



// compute background frequencies
void computeBgfreq(int nseq, char **seqs, float *bgfreq);


} // namespace spidir

#endif // SPIDIR_SEQ_H
