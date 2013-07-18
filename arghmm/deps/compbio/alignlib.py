#
# alignlib.py
# Sequence and alignment functions
#


# python libs
import sys
from collections import defaultdict

# rasmus libs
from rasmus import util

# compbio libs
from . import fasta, seqlib
from seqlib import *



#=============================================================================
# Alignment functions


def new_align(aln=None):
    """Makes a new alignment object based on the given object
        
       given      return
       -----      ------
       dict       FastaDict
       other      other
    """
    
    if aln is None:
        return fasta.FastaDict()
    elif isinstance(aln, SeqDict):
        return type(aln)()
    else:
        return fasta.FastaDict()


def mapalign(aln, keyfunc=lambda x: x, valfunc=lambda x: x):
    """Maps the keys and values of an alignment"""
    
    aln2 = new_align(aln)
    for key, val in aln.iteritems():
        aln2[keyfunc(key)] = valfunc(val)
    return aln2

            
def subalign(aln, cols):
    """Returns an alignment with a subset of the columns (cols)"""
    
    return mapalign(aln, valfunc=lambda x: "".join(util.mget(x, cols)))


def remove_empty_columns(aln):
    """
    Removes any column from an alignment 'aln' that contains only gaps
    
    A new alignment is returned
    """

    ind = []
    seqs = aln.values()
    for i in range(aln.alignlen()):
        for seq in seqs:
            if seq[i] != "-":
                ind.append(i)
                break
    
    return subalign(aln, ind)


def remove_gapped_columns(aln):
    """Removes any column form an alignment 'aln' that contains a gap
    
       A new alignment is returned
    """
    cols = zip(* aln.values())
    ind = util.find(lambda col: "-" not in col, cols)
    return subalign(aln, ind)


def require_nseqs(aln, n):
    """
    Keep only columns with atleast 'n' non gapped sequences
    """

    seqs = aln.values()
    ind = [i for i in range(aln.alignlen())
           if sum(1 for seq in seqs if seq[i] != "-") >= n]
    return subalign(aln, ind)


def get_seq_overlap(seq1, seq2):
    """
    Count number of overlapping bases between two gapped sequences
    """

    overlap = 0
    for i in range(len(seq1)):
        if seq1[i] != "-" and seq2[i] != "-":
            overlap += 1
            
    return overlap


def calc_conservation_string(aln):
    """Returns a string of stars representing the conservation of an alignment"""
    
    percids = calc_conservation(aln)
    
    # find identity positions
    identity = ""
    for pid in percids:
        if pid == 1:
            identity += "*"
        elif pid > .5:
            identity += "."
        else:
            identity += " "
    
    return identity


def calc_conservation(aln):
    """Returns a list of percent matching in each column of an alignment"""

    length = len(aln.values()[0])
    seqs = aln.values()
    percids = []
    
    # find identity positions
    identity = ""
    for i in xrange(length):
        chars = util.hist_dict(util.cget(seqs, i))
        if "-" in chars: del chars["-"]
        
        if len(chars) == 0:
            percids.append(0.0)
        else:
            pid = max(chars.values()) / float(len(aln))
            percids.append(pid)
    return percids



def print_align(aln, seqwidth = 59, spacing=2, extra=fasta.FastaDict(), 
               out=sys.stdout, order=None):
    """Pretty print an alignment"""
               
    if order == None:
        order = aln.keys()
    
    namewidth = max(map(len, order)) + spacing
    
    def mkname(name, namewidth):
        name2 = name[:namewidth]
        name2 += " " * (namewidth - len(name2))
        return name2

    identity = calc_conservation_string(aln)
    
    # print alignment
    for i in xrange(0, len(aln.values()[0]), seqwidth):
        # print sequences
        for name in order:
            print >>out, "%s %s" % (mkname(name, namewidth), 
                                    aln[name][i:i+seqwidth])
        
        # print extra
        for name in extra.keys():
            print >>out, "%s %s" % (mkname(name, namewidth), 
                                    extra[name][i:i+seqwidth])
        
        # print identity
        print >>out, (" "*namewidth) + " " + identity[i:i+seqwidth]
        print >>out


def revtranslate_align(aaseqs, dnaseqs, check=False, trim=False):
    """Reverse translates aminoacid alignment into DNA alignment
    
       Must supply original ungapped DNA.
    """
    
    align = new_align(aaseqs)
    
    for name, seq in aaseqs.iteritems():
        try:
            dna = dnaseqs[name].upper()
            dnalen = len(dna)
            aalen = sum(int(a != "-") for a in seq)
            
            if len(dna) != aalen * 3:
                if trim:
                    # make dna a multiple of three
                    dna = dna[:(len(dna) // 3) * 3]

                    if len(dna) > aalen * 3:
                        # trim dna
                        dna = dna[:aalen*3]
                    else:
                        # trim peptide to match nucleotide
                        j = 0
                        for i in xrange(len(seq)):
                            if seq[i] != '-':
                                j += 1
                                if j > len(dna) // 3:
                                    seq = seq[:i] + "-" * (len(seq) - i)
                                    break

                    aalen2 = sum(int(a != "-") for a in seq)
                    assert len(dna) == aalen2 * 3,  (
                        len(dna), aalen2 * 3)

                    util.logger("trim dna (%d) and pep (%d)" %
                                (dnalen - len(dna), aalen - aalen2))

                else:
                    # is last residue X?
                    for i in xrange(len(seq)-1, -1, -1):
                        if seq[i] == "-":
                            continue
                        if seq[i] == "X":
                            # repair
                            seq =  seq[:i] + "-" * (len(seq)-i)
                            dna = dna[:-3] #-(len(dna) % 3)]
                        break

            
            align[name] = revtranslate(seq, dna, check=check)
        except TranslateError, e:
            raise
    
    return align




#=============================================================================
# four fold degeneracy


def mark_codon_pos(seq, pos=0):
    """
    return the codon position for each base in a gapped sequence

    codon
    ATG
    012

    gaps are given codon pos -1
    Ns are counted as bases
    """
    
    codons = []

    for base in seq:
        if base != "-":
            codons.append(pos)
            pos = (pos + 1) % 3
        else:
            codons.append(-1)

    return codons


def make_codon_pos_align(aln):
    """Get the codon position of every base in an alignment"""
    
    def func(seq):
        dct = {-1: "-",
               0: "0",
               1: "1",
               2: "2"}
        return "".join(util.mget(dct, mark_codon_pos(seq)))
    return mapalign(aln, valfunc=func)


def find_aligned_codons(aln):
    """Returns the columns indices of the alignment that represent aligned 
       codons.  
    """

    ind = range(aln.alignlen())
    
    # throw out codons with non mod 3 gaps
    ind2 = []
    for i in range(0, aln.alignlen(), 3):
        bad = False
        
        for key, val in aln.iteritems():
            codon = val[i:i+3]
            if "-" in codon and codon != "---":
                bad = True
                break

        if not bad:
            ind2.extend([i, i+1, i+2])

    return ind2



def filter_aligned_codons(aln):
    """filters an alignment for only aligned codons"""

    ind = find_align_codons(aln)
    return subalign(aln, ind)


def find_four_fold(aln):
    """Returns index of all columns in alignment that are completely 
       fourfold degenerate
       
       Assumes that columns are already filtered for aligned codons
    """
    
    # create peptide alignment
    pepAln = mapalign(aln, valfunc=translate)
    
    # find peptide conservation
    pepcons = []
    pep = []
    for i in xrange(pepAln.alignlen()):
        # get a column from the peptide alignment
        col = [seq[i] for seq in pepAln.itervalues()]
        
        # compute the histogram of the column.
        # ignore gaps '-' and non-translated 'X'
        hist = util.hist_dict(col)
        if "-" in hist:
            del hist["-"]
        if "X" in hist:
            del hist["X"]
        
        # column is conserved if only one AA appears
        if len(hist) == 1:
            pepcons.append(True)
            pep.append(hist.keys()[0])
        else:
            pepcons.append(False)
            pep.append("X")
        
    
    # find four-fold sites in conserved peptides
    ind = []
    
    for i in range(0, len(aln.values()[0]), 3):
        # process only those columns that are conserved at the peptide level
        if pepcons[i//3]:
            degen = AA_DEGEN[pep[i//3]]
            for j in range(3):
                if degen[j] == 4:
                    ind.append(i+j)
    return ind


def filter_four_fold(aln):
    """returns an alignment of only four-fold degenerate sites from an 
       alignment of coding sequences
    
       This function performs the following steps:
       
       1. remove all codon columns that don't have 0 or 3 gaps
       2. keep all codon columns that code for identical AA
       3. if the codon column codes for a 4D AA, then keep its 3rd position
    """

    aln_codons = filter_align_codons(aln)
    ind = find_four_fold(aln_codons)
    return subalign(aln_codons, ind)


def calc_four_fold_dist_matrix(aln):
    names = aln.keys()

    mat = []
    # calc upper triangular
    for i in range(len(names)):
        mat.append([0.0] * (i+1))
        for j in range(i+1, len(names)):
            ind = find_four_fold(aln.get([names[i], names[j]]))

            mismatches = 0
            for k in ind:
                if aln[names[i]][k] != aln[names[j]][k]:
                    mismatches += 1
            
            if len(ind) == 0:
                mat[-1].append(1.0)
            else:            
                mat[-1].append(mismatches / float(len(ind)))

    # make symmetric
    for j in range(len(names)):
        for i in range(j):
            mat[j][i] = mat[i][j]

    return mat


def find_degen(aln):
    """Determine the degeneracy of each column in an alignment"""

    codon_ind = find_align_codons(aln)
    aln2 = subalign(aln, codon_ind)
    
    pep_aln = mapalign(aln2, valfunc=translate)
    pep = pep_aln.values()[0]
    identies = calc_conservation(pep_aln)
    
    degens = [-1] * aln.alignlen()
    
    for i in range(0, len(codon_ind), 3):
        if pep[i/3] == "X":
            continue
        degen = AA_DEGEN[pep[i/3]]
        if identies[i/3] == 1.0:
            for j in range(3):
                degens[codon_ind[i+j]] = degen[j]
                    
    return degens


def make_degen_str(aln):
    """Returns a string containing the degeneracy for each column 
       in an alignment
    """

    degens = find_degen(aln)
    degenmap = {-1: " ",
                 0: "0",
                 1: "1",
                 2: "2",
                 3: "3",
                 4: "4"}
    
    return "".join(util.mget(degenmap, degens))
    

def print_degen(aln, **args):
    """Pretty print an alignment with its degeneracy for each column"""

    extra = fasta.FastaDict()
    extra["DEGEN"] = make_degen_str(aln)
    
    print_align(aln, extra=extra, **args)


#=============================================================================
# Position Specific Scoring Matrix (PSSM)


def align2pssm(aln, pseudocounts = {}):
    pssm = []
    denom = float(len(aln)) + sum(pseudocounts.values())
    
    for i in xrange(aln.alignlen()):
        freqs = defaultdict(lambda: 0)
        for j in xrange(len(aln)):
            freqs[aln[j][i]] += 1
        
        for key in pseudocounts:
            freqs[key] += pseudocounts[key]
        
        for key in freqs:
            freqs[key] = math.log(freqs[key] / denom, 2)
        pssm.append(freqs)
    return pssm


def pssmSeq(pssm, seq):
    score = 0.0
    for i in xrange(len(seq)):
        score += pssm[i][seq[i]]
    return score



#=============================================================================
# Coordinate conversions
#
# Coordinate systems
# 
#   1. local
#       01234567
#       ATGCTGCG
# 
#   2. align
#       012222345567
#       ATG---CTG-CG
#
#   3. global
#       coordinate on chromosome on positive strand
#
# There should only be two kinds of indexing
# 1. 0-based, end exclusive (local/align coordinates)
# 2. 1-based, end inclusive (global coordinates)
#


class CoordConverter (object):
    """Converts between coordinate systems on a gapped sequence"""
    
    def __init__(self, seq):
        self.local2alignLookup = local2align(seq)
        self.align2localLookup = align2local(seq)
    
    
    def local2align(self, i, clamp=False):
        if clamp:
            return self.local2alignLookup[int(util.clamp(i, 0, 
                                           len(self.local2alignLookup)-1))]
        else:
            return self.local2alignLookup[i]


    def align2local(self, i, clamp=False):
        if clamp:
            return self.align2localLookup[int(util.clamp(i, 0, 
                                           len(self.align2localLookup)-1))]
        else:
            return self.align2localLookup[i]


    def global2local(self, gobal_coord, start, end, strand):
        """Returns local coordinate in a global region"""
        return global2local(gobal_coord, start, end, strand)
        

    def local2global(self, local_coord, start, end, strand):
        """Return global coordinate within a region from a local coordinate"""
        local2global(local_coord, start, end, strand)


    def global2align(self, global_coord, start, end, strand):
        local_coord = global2local(global_coord, start, end, strand)
    
        # throw exception for out of bounds
        if local_coord < 0 or \
           local_coord >= len(alignLookup):
            raise Exception("coordinate outside [start, end]")

        return self.local2alignLookup[local_coord]


    def align2global(self, align_coord, start, end, strand):
        local_coord = self.align2localLookup[align_coord]
        return local2global(local_coord, start, end, strand)



def local2align(seq):
    """
    Returns list of indices of non-gap characters
    
    'ATG---CTG-CG' ==> [0,1,2,6,7,8,10,11]
    
    Used to go from local -> align space
    """
    
    lookup = []
    for i in xrange(len(seq)):
        if seq[i] == "-": continue
        lookup.append(i)
    return lookup


def align2local(seq):
    """
    Returns list such that 
    
    'ATG---CTG-CG' ==> [0,1,2,2,2,3,4,5,5,6,7]
    
    Used to go from align -> local space
    """

    i = -1
    lookup = []
    for c in seq:
        if c != "-":
            i += 1
        lookup.append(i)
    return lookup



def global2local(gobal_coord, start, end, strand):
    """Returns local coordinate in a global region"""

    # swap if strands disagree
    if strand == 1:
        return gobal_coord - start
    else:
        return end - gobal_coord


def local2global(local_coord, start, end, strand):
    """Return global coordinate within a region from a local coordinate"""
    
    # swap if strands disagree
    if strand == 1:
        return local_coord + start
    else:
        return end - local_coord


def global2align(global_coord, start, end, strand, alignLookup):
    local_coord = global2local(global_coord, start, end, strand)
    
    # throw exception for out of bounds
    if local_coord < 0 or \
       local_coord >= len(alignLookup):
        raise Exception("coordinate outside [start, end]")
    
    return alignLookup[local_coord]


def align2global(align_coord, start, end, strand, localLookup):
    local_coord = localLookup[align_coord]
    return local2global(local_coord, start, end, strand)




#=============================================================================
# old code

'''
def findAlignCodons(aln):
    """find all columns of aligned codons"""
    
    codonAln = mapalign(aln, valfunc=mark_codon_pos)
    cols = map(util.hist_dict, zip(* codonAln.values()))

    ind = []
    codon = []
    gaps = defaultdict(lambda: 0)
    for i in range(len(cols)):

        if len(cols[i]) == 1:
            codon.append(i)
        elif len(cols[i]) == 2 and -1 in cols[i]:
            for key, val in aln.iteritems():
                if val[i] == "-":
                    gaps[key] += 1 
            codon.append(i)
        else:
            codon = []
        if len(codon) == 3:
            if len(gaps) == 0 or \
               util.unique([x % 3 for x in gaps.values()]) == [0]:
                ind.extend(codon)
            codon = []
            for key in gaps:
                gaps[key] = 0

    return ind


def findFourFold(aln):
    """Returns index of all columns in alignment that are completely 
       fourfold degenerate
    """
    
    aln = filterAlignCodons(aln)
    pepAln = mapalign(aln, valfunc=translate)
    pep = pepAln.values()[0]
    
    # pep conservation
    pepcons = []
    for i in xrange(pepAln.alignlen()):
        col = [seq[i] for seq in pepAln.itervalues()]
        hist = util.hist_dict(col)
        if "-" in hist:
            del hist["-"]
        if "X" in hist:
            del hist["X"]
        pepcons.append(len(hist) == 1)
        

    ind = []

    for i in range(0, len(aln.values()[0]), 3):
        if pepcons[i//3]:
            degen = AA_DEGEN[pep[i//3]]
            
            for j in range(3):
                if degen[j] == 4:
                    ind.append(i+j)
    return ind
'''
