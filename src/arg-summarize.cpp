
// C/C++ includes
#include <time.h>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <assert.h>
#include <vector>
#include <math.h>
#include <queue>
#include <map>

// arghmm includes
#include "ConfigParam.h"
#include "logging.h"
#include "parsing.h"
#include "track.h"
#include "Tree.h"
#include "tabix.h"
#include "compress.h"
#include "IntervalIterator.h"
#include <set>

//#include "allele_age.h"

using namespace argweaver;
using namespace spidir;

/* Initial version: just have it output bedgraph with a stat for each
   line of input.
 */


int summarize=0;
int getNumSample=0;
int getMean=0;
int getStdev=0;
int getQuantiles=0;
vector <double> quantiles;
char *tabix_dir=NULL;
int header=true;
char version[10]="0.2";

int print_help() {
    printf("usage: ./arg-summarize [OPTIONS] <file.bed.gz>\n"
           "   where <file.bed.gz> is bed file created by smc2bed\n"
           "OPTIONS:\n"
           "--region,-r CHR:START-END\n"
           "   retrieve statistics from this region (coordinates are 1-based)\n"
           "   To use this option, the bed file should be indexed by tabix\n"
           "   and the tabix command should be available\n"
           "--bed,-b bedfile\n"
           "   (alternative to region). Retrieve statistics from regions\n"
           "   overlapping those in this file. If bedfile contains overlapping\n"
           "   regions, those regions will be repeated in the output. Only the\n"
           "   first three columns of this file will be used.\n"
	   "--subset,-s sample-list-file\n"
	   "   subset by individual. The argument should contain haplotype names, one\n"
	   "   per line, which should be retained. All others will be removed before\n"
	   "   statistics are computed.\n"
	   "--snp,-s snpfile\n"
	   "   If used, compute statitics for specific SNPs. If used, extra information\n"
	   "   about each SNP will be printed on each row (derived/ancestral alleles and\n"
	   "   frequencies). If statistics are summarized across samples, there will be\n"
	   "   three estimates for each summary statistic- one taken across all samples,\n"
	   "   one taken across samples which agree with the consensus derived allele,\n"
           "   and one taken across only samples which agree with infiite sites (only one\n"
	   "   mutation required to show given site pattern.)\n"
	   "   The SNP file should have\n"
	   "   a header like this:\n"
	   "   #NAMES NA06985_1       NA06985_2       NA06994_1       NA06994_2\n"
	   "   and then each line should be tab-delimited with the following format:\n"
	   "   chr,start,end,AAAACAAAAA\n"
	   "   where the last column gives the alleles for each haplotype in the order\n"
	   "   listed in the header. The file should be sorted and index with tabix\n"
           "STATISTICS:\n"
           "--tree,-E\n"
           "  output newick tree strings (cannot use summary options with this)\n"
           "--tmrca,-T\n"
           "  time to the most recent common ancestor\n"
           "--branchlen,-B\n"
           "  total branch length\n"
	   "--recomb,-R\n"
	   "  recombination rate (per generation/bp)\n"
           "--tmrca-half,-H\n"
           "  Time for half of samples to reach a common ancestor\n"
           "--rth,-F\n"
           "  relative TMRCA Halftime\n"
           "--popsize,-P\n"
           "  estimated popsize\n"
	   "--allele-age,-A\n"
	   "  (requires --snp). Compute allele age for all CG snps in\n"
	   "  the region. Also returns a flag indicating whether SNP obeys\n"
	   "  infinite sites (1=obeys, 0=multiple mutations required to\n"
	   "  acheive observed site pattern)\n"
           "--numsample,-N\n"
           "  number of samples covering each region\n"
           "--mean,-M\n"
           "  return mean across all samples\n"
           "--stdev,-S\n"
           "  return standard deviation across all samples\n"
           "--quantile,-Q <q1,q2,...>\n"
           "  return the requested quantiles for each samples. q1,q2,...\n"
           "  should be comma-separated list of desired quantiles\n"
           "  (ie, 0.025,0.5,0.0975)\n"
           "OTHER OPTIONS:\n"
           "--no-header,-n\n"
           "--tabix-dir,-t <tabix_dir>\n"
           "  Specify the directory where the tabix executable is installed\n"
           "--help,-h\n"
           "   print this help.\n");
    return 0;
}

void checkResults(IntervalIterator<vector<double> > *results) {
    Interval<vector<double> > summary=results->next();
    vector<vector <double> > scores;
    while (summary.start != summary.end) {
        cout << summary.chrom << "\t" << summary.start << "\t"
             << summary.end;
        scores = summary.get_scores();
        if (scores.size() > 0) {
            vector<double> tmpScore(scores.size());
            int numscore = scores[0].size();
            assert(numscore > 0);
            for (int i=0; i < numscore; i++) {
                int have_mean = 0;
                double meanval=-1;
                for (unsigned int j=0; j < scores.size(); j++)
                    tmpScore[j] = scores[j][i];
                if (i==0 && getNumSample > 0) printf("\t%i", (int)scores.size());
                for (int j=1; j <= summarize; j++) {
                    if (getMean==j) {
                        meanval = compute_mean(tmpScore);
                        have_mean=1;
                        printf("\t%g", meanval);
                    } else if (getStdev==j) {
                        if (!have_mean)
                            meanval = compute_mean(tmpScore);
                        printf("\t%g", compute_stdev(tmpScore, meanval));
                    } else if (getQuantiles==j) {
                        vector<double> q = compute_quantiles(tmpScore, quantiles);
                        for (unsigned int k=0; k < quantiles.size(); k++) {
                            printf("\t%g", q[k]);
                        }
                    }
                }
            }
            printf("\n");
	    fflush(stdout);
        }
        summary = results->next();
    }
}

class BedLine {
public:
  BedLine(char *chr, int start, int end, int sample, char *nwk, 
	  Tree *orig_tree = NULL, Tree *pruned_tree = NULL):
    start(start), end(end), sample(sample), 
    orig_tree(orig_tree), pruned_tree(pruned_tree) {
	chrom = new char[strlen(chr)+1];
	strcpy(chrom, chr);
	if (nwk != NULL) {
	  newick = (char*)malloc((strlen(nwk)+1)*sizeof(char));
	  strcpy(newick, nwk);
	} else newick=NULL;
    };
  ~BedLine() {
    //      if (orig_tree != NULL) delete orig_tree;
    //      if (pruned_tree != NULL) delete pruned_tree;
      delete [] chrom;
      if (newick != NULL) free(newick);
  }
  char *chrom;
  int start;
  int end;
  int sample;
  Tree *orig_tree;
  Tree *pruned_tree;
  char *newick;
  vector<double> stats;
  char derAllele, otherAllele;
  int derFreq, otherFreq;
  int infSites;
};


void scoreBedLine(BedLine *line, vector<string> &statname, double allele_age=-1, int infsites=-1) {
  Tree *tree = line->pruned_tree != NULL ? line->pruned_tree : line->orig_tree;
  double bl=-1.0;
  if (line->stats.size() == statname.size()) return;
  //  fprintf(stderr, "scoreBedLine %i\t%s\t%i\t%i\n", line->sample,line->chrom, line->start, line->end); fflush(stderr);
  line->stats.resize(statname.size());
  for (unsigned int i=0; i < statname.size(); i++) {
    if (statname[i] == "tmrca")
      line->stats[i] = tree->tmrca();
    else if (statname[i]=="tmrca_half")
      line->stats[i] = tree->tmrca_half();
    else if (statname[i]=="branchlen") {
      if (bl < 0) {
	line->stats[i] = tree->total_branchlength();
	bl=line->stats[i];
      }
    }
    else if (statname[i]=="rth")
      line->stats[i] = tree->rth();
    else if (statname[i]=="popsize")
      line->stats[i] = tree->popsize();
    else if (statname[i]=="recomb") {
      if (bl < 0) bl = tree->total_branchlength();
      line->stats[i] = 1.0/(bl*(double)(line->end - line->start));
    }
    else if (statname[i]=="tree") {
      if (line->pruned_tree != NULL) {
	string tmp = line->pruned_tree->print_newick_to_string(false, true, 1, 1);
	//	line->newick = (char*)malloc((tmp.length()+1)*sizeof(char));
	//pruned tree will be fewer characters than whole tree
	sprintf(line->newick, "%s", tmp.c_str());
      }
    }
    else if (statname[i]=="allele_age")
      line->stats[i] = allele_age;
    else if (statname[i]=="inf_sites") 
      line->stats[i] = (double)infsites;
    else {
      fprintf(stderr, "Error: unknown stat %s\n", statname[i].c_str());
      exit(1);
    }
  }
}


struct CompareBedLineSample
{
  bool operator()(const BedLine *l1, const BedLine *l2) const
  {
    return l1->sample < l2->sample;
  }
};


struct CompareBedLineEnd
{
    bool operator()(const BedLine *l1, const BedLine *l2) const
    { 
      //could generalize and sort by start, but for this purpose we are only comparing entries with same start
      assert(l1->start == l2->start);
      if (l1->end == l2->end)
	return (l1->sample < l2->sample);
      return l1->end < l2->end;
    }
};

void processNextBedLine(BedLine *line, IntervalIterator<vector<double> > *results,
			vector<string> &statname, 
			char *region_chrom, int region_start, int region_end) {
  static int counter=0;
  static list<BedLine*> bedlist;

  if (line != NULL) {
    if (line->stats.size() == 0) 
      scoreBedLine(line, statname);
    if (region_chrom != NULL) {
      assert(strcmp(region_chrom, line->chrom)==0);
      if (line->end > region_end) line->end = region_end;
      if (line->start < region_start) line->start = region_start;
      assert(line->start < line->end);
    }
  }
  if (!summarize) {
    // this little bit of code ensures that output is sorted. The summarizeRegion
    // code should ensure that it comes here sorted by start coordinate, but not
    // necessarily by end coordinate (at least not when subsetting by individuals).
    // so, this stores BedLine elements just long enough until a new start 
    // coordinate is encountered, then sorts them all by end coordinate (and 
    // sample), then outputs them all.  Need to call this function one last 
    // time with line==NULL to output the final set of rows.
    if (bedlist.size() > 0 && 
	(line==NULL || bedlist.front()->start < line->start)) {
      bedlist.sort(CompareBedLineEnd());
      for (list<BedLine*>::iterator it=bedlist.begin(); it != bedlist.end(); ++it) {
	BedLine *l = *it;
	printf("%s\t%i\t%i\t%i", l->chrom, l->start, l->end, l->sample);
	for (unsigned int i=0; i < statname.size(); i++) {
	  if (statname[i]=="tree") {
	    printf("\t");
	    printf("%s", l->newick);
	  } else {
	    printf("\t%g", l->stats[i]);
	  }
	}
	printf("\n"); fflush(stdout);
	delete l;
      }
      bedlist.clear();
    }
    if (line != NULL) bedlist.push_back(line);
  } else {
    if (line != NULL) {
      results->append(line->chrom, line->start, line->end, line->stats);
      counter++;
      if (counter%100==0) {
	checkResults(results);
      }
      delete line;
    }
  }
}

class SnpStream {
public:
  SnpStream(TabixStream *snp_in) : snp_in(snp_in) {
    char tmp[1000], c;
    string str;
    assert(1==fscanf(snp_in->stream, "%s", tmp));
    assert(strcmp(tmp, "#NAMES")==0);
    done=0;
    while ('\n' != (c=fgetc(snp_in->stream)) && c!=EOF) {
      assert(c=='\t');
      assert(1==fscanf(snp_in->stream, "%s", tmp));
      str = string(tmp);
      inds.push_back(str);
    }
    if (c==EOF) done=1;
  }

  int readNext() {
    int tmpStart;
    char a;
    if (done) return 1;
    if (EOF==fscanf(snp_in->stream, "%s %i %i", chr, &tmpStart, &coord)) {
      done=1;
      return 1;
    }
    assert(tmpStart==coord-1);
    assert('\t' == fgetc(snp_in->stream));
    allele1=allele2='N';
    allele1_inds.clear();
    allele2_inds.clear();
    for (unsigned int i=0; i < inds.size(); i++) {
      a=fgetc(snp_in->stream);
      if (a=='N') continue;
      a = toupper(a);
      assert(a=='A' || a=='C' || a=='G' || a=='T');
      if (allele1=='N') {
	allele1=a;
      } 
      if (a==allele1) {
	allele1_inds.insert(inds[i]);
      } else {
	if (allele2=='N')
	  allele2=a;
	else assert(a==allele2);
	allele2_inds.insert(inds[i]);
      }
    }
    //make sure that allele1 is always minor allele
    if (allele1_inds.size() > allele2_inds.size()) {
      set<string> tmp;
      char tmpch;
      tmp = allele1_inds;
      allele1_inds = allele2_inds;
      allele2_inds = tmp;
      tmpch=allele1;
      allele1=allele2;
      allele2=tmpch;
    }
    //    printf("Read snp %i allele1=%c allele2=%c n1=%i n2=%i\n",
    //	   coord, allele1, allele2, (int)allele1_inds.size(), (int)allele2_inds.size());
    return 0;
  }


  void scoreAlleleAge(BedLine *l, vector<string> statname) {
    int num_derived, total;
    assert(l->start < coord);
    assert(l->end >= coord);
    Tree *t;
    if (l->pruned_tree != NULL) 
      t = l->pruned_tree;
    else t = l->orig_tree;
    
    set<string> prune;
    set<string> derived_in_tree;
    for (map<string,int>::iterator it=t->nodename_map.begin(); 
	 it != t->nodename_map.end(); ++it) {
      if (t->nodes[it->second]->nchildren != 0) continue;
      if (allele1_inds.find(it->first) != allele1_inds.end())
	derived_in_tree.insert(it->first);
      else if (allele2_inds.find(it->first) == allele2_inds.end())
	prune.insert(it->first);
    }
    if (derived_in_tree.size() == 0 ||
	derived_in_tree.size() + prune.size() == (unsigned int)(t->nnodes+1)/2)
      return;

    if (prune.size() > 0) {
      Tree *t2 = t->copy();
      t2->prune(prune, false);
      t = t2;
    }

    set<Node*> derived;
    for (set<string>::iterator it=derived_in_tree.begin(); it != derived_in_tree.end(); ++it) {
      map<string,int>::iterator it2 = t->nodename_map.find(*it);
      assert(it2 != t->nodename_map.end());
      derived.insert(t->nodes[it2->second]);
    }
    num_derived = (int)derived.size();
    total = (t->nnodes+1)/2;
    
    set<Node*>lca = t->lca(derived);
    set<Node*>lca2;
    int major_is_derived=0;
    if (lca.size() > 1) {
      set<Node*> derived2;
      for (map<string,int>::iterator it=t->nodename_map.begin();
	   it != t->nodename_map.end(); ++it) {
	if (t->nodes[it->second]->nchildren != 0) continue;
	if (derived.find(t->nodes[it->second]) == derived.end()) {
	  derived2.insert(t->nodes[it->second]);
	}
      }
      assert(derived2.size() > 0);
      set <Node*>lca2 = t->lca(derived2);
      if (lca2.size() < lca.size()) {
	major_is_derived=1;
	lca = lca2;
      }
    }
    assert(lca.size() >= 1);
    double age=0.0;
    for (set<Node*>::iterator it4=lca.begin(); it4 != lca.end(); ++it4) {
      Node *n = *it4;
      assert(n != t->root);
      double tempage = n->age + (n->parent->age - n->age)/2;  //midpoint of branch
      if (tempage > age) age = tempage;
    }
    scoreBedLine(l, statname, age, lca.size()==1);
    l->derAllele = (major_is_derived ? allele2 : allele1);
    l->otherAllele = (major_is_derived ? allele1 : allele2);
    l->derFreq = (major_is_derived ? total-num_derived : num_derived);
    l->otherFreq = (major_is_derived ? num_derived : total - num_derived);
    l->infSites = (lca.size() == 1);

    /*    if (summarize==0) {
      printf("%s\t%i\t%i\t%i", 
	     chr, coord-1, coord, l->sample);
      rv.print();
      printf("\n");
      }*/
    if (prune.size() > 0) {
      delete t;
    }
  }

  TabixStream *snp_in;
  vector<string> inds;
  set<string> allele1_inds;
  set<string> allele2_inds;
  char allele1, allele2;  //minor allele, major allele
  char chr[100];
  int coord;  //1-based
  int done;
};



void print_summaries(vector<double> stat) {
  double meanval=0;
  int have_mean=0;
  for (int j=1; j <= summarize; j++) {
    if (getMean==j) {
      if (stat.size() > 0) {
	meanval = compute_mean(stat);
	have_mean=1;
	printf("\t%g", meanval);
      } else printf("\tNA");
    } else if (getStdev==j) {
      if (stat.size() > 1) {
	if (!have_mean)
	  meanval = compute_mean(stat);
	printf("\t%g", compute_stdev(stat, meanval));
      } else printf("\tNA");
    } else if (getQuantiles==j) {
      if (stat.size() > 0) {
	vector<double> q = compute_quantiles(stat, quantiles);
	for (unsigned int k=0; k < quantiles.size(); k++) {
	  printf("\t%g", q[k]);
	}
      } else {
	for (unsigned int k=0; k < quantiles.size(); k++) {
	  printf("\tNA");
	}
      }
    }
  }
}


int summarizeRegionBySnp(char *snpFilename, char *filename,
			 const char *region,
			 set<string> inds, vector<string> statname,
			 vector<double> times) {
  TabixStream *snp_infile;
  TabixStream *infile;
  vector<string> token;
  map<int,BedLine*> last_entry;
  map<int,BedLine*>::iterator it;
  char chrom[1000], c;
  int start, end, sample;
  BedLine *l=NULL;

  snp_infile = new TabixStream(snpFilename, region, tabix_dir);
  if (snp_infile->stream == NULL) return 1;

  infile = new TabixStream(filename, region, tabix_dir);
  if (infile->stream == NULL) return 1;
  while (EOF != (c=fgetc(infile->stream))) {
    ungetc(c, infile->stream);
    if (c != '#') break;
    while ('\n' != (c=fgetc(infile->stream))) {
      if (c==EOF) return 0;
    }
  }
  SnpStream snpStream = SnpStream(snp_infile);
  assert(4==fscanf(infile->stream, "%s %i %i %i", 
		   chrom, &start, &end, &sample));
  assert('\t' == fgetc(infile->stream));
  char *newick = fgetline(infile->stream);
  chomp(newick);
  
  while (1) {
    list<BedLine*> bedlist;
    bedlist.clear();
    snpStream.readNext();
    if (snpStream.done) break;
    // first check already-parsed BedLines and score any that overlap SNP
    for (it=last_entry.begin(); it != last_entry.end(); it++) {
      l = it->second;
      if (l->start < snpStream.coord && l->end >= snpStream.coord) {
	snpStream.scoreAlleleAge(l, statname);
	bedlist.push_back(l);
      }
    }
    //now look through bed file until we get to one that starts after SNP
    while (start != -1 && snpStream.coord > start) {
      it = last_entry.find(sample);
      if (it == last_entry.end() || it->second->orig_tree->recomb_node == NULL) {
	Tree *orig_tree, *pruned_tree=NULL;
	if (it != last_entry.end()) {
	  l = it->second;
	  delete l->orig_tree;
	  if (inds.size() > 0) delete l->pruned_tree;
	  delete &*l;
	}
	orig_tree = new Tree(string(newick), times);
	if (inds.size() > 0) {
	  pruned_tree = orig_tree->copy();
	  orig_tree->node_map = pruned_tree->prune(inds, true);
	}
	l = new BedLine(chrom, start, end, sample, newick, orig_tree, pruned_tree);
	last_entry[sample] = l;
      } else {
	l = it->second;
	l->orig_tree->apply_spr();
	l->orig_tree->update_spr(newick, times);
	free(l->newick);
	l->newick = (char*)malloc((strlen(newick)+1)*sizeof(char));
	strcpy(l->newick, newick);
	if (inds.size() > 0) {
	  l->pruned_tree->apply_spr();
	  l->pruned_tree->update_spr_pruned(l->orig_tree);
	}
	l->start = start;
	l->end = end;
      }
      if (snpStream.coord <= end) {
	snpStream.scoreAlleleAge(l, statname);
	bedlist.push_back(l);
      }
      if (4 != fscanf(infile->stream, "%s %i %i %i", chrom, &start, &end, &sample))
	start = -1;
      else {
	assert('\t' == fgetc(infile->stream));
	delete [] newick;
	newick = fgetline(infile->stream);
	chomp(newick);
      }
    }
    if (bedlist.size() > 0) {
      if (summarize == 0) {
	bedlist.sort(CompareBedLineSample());
	for (list<BedLine*>::iterator it=bedlist.begin(); it != bedlist.end(); ++it) {
	  BedLine *l = *it;
	  printf("%s\t%i\t%i\t%i\t%c\t%c\t%i\t%i", l->chrom, snpStream.coord-1, snpStream.coord, l->sample, l->derAllele, l->otherAllele, l->derFreq, l->otherFreq); fflush(stdout);
	  for (unsigned int i=0; i < statname.size(); i++) {
	    if (statname[i]=="tree") {
	      printf("\t%s", l->newick);
	    } else if (statname[i]=="infSites") {
		printf("\t%i", (int)(l->stats[i]==1));
	    } else {
		printf("\t%g", l->stats[i]);
	    }
	  }
	  printf("\n"); fflush(stdout);
	  l->stats.clear();
	}
      } else {
	printf("%s\t%i\t%i", l->chrom, snpStream.coord-1, snpStream.coord); fflush(stdout);
	//now output three versions- one for all samples, one for same derived allele, one for infinite sites
	BedLine* first = *(bedlist.begin());
	int same=0, diff=0, infsites=0, derConstCount, derFreq, otherFreq;
	char derAllele, otherAllele;
	for (list<BedLine*>::iterator it=bedlist.begin(); it != bedlist.end(); ++it) {
	  BedLine *l = *it;
	  if (l->derAllele == first->derAllele) same++; else diff++;
	  infsites += l->infSites;
	}
	if (same >= diff) {
	  derAllele = first->derAllele;
	  otherAllele = first->otherAllele;
	  derConstCount = same;
	  derFreq = first->derFreq;
	  otherFreq = first->otherFreq;
	} else {
	  derAllele=first->otherAllele;
	  otherAllele = first->derAllele;
	  derConstCount = diff;
	  derFreq = first->otherFreq;
	  otherFreq = first->derFreq;
	}
	printf("\t%c\t%c\t%i\t%i\t%i\t%i",
	       derAllele, otherAllele, derFreq, otherFreq,
	       (int)bedlist.size(), infsites);
	for (unsigned int i=0; i < statname.size(); i++) {
	    if (statname[i] != "inf_sites") {
	  // first compute stats across all
	  vector<double> stat;
	  for (list<BedLine*>::iterator it=bedlist.begin(); it != bedlist.end(); ++it) {
	    BedLine *l = *it;
	    stat.push_back(l->stats[i]);
	  }
	  print_summaries(stat);

	  /*	  stat.clear();
	  // now stats across derived const set
	  for (list<BedLine*>::iterator it=bedlist.begin(); it != bedlist.end(); ++it) {
	    BedLine *l = *it;
	    if (l->derAllele == derAllele)
	      stat.push_back(l->stats[i]);
	  }
	  print_summaries(stat);*/

	  stat.clear();
	  //now stats for infinite sites set
	  for (list<BedLine*>::iterator it=bedlist.begin(); it != bedlist.end(); ++it) {
	    BedLine *l = *it;
	    if (l->infSites)
	      stat.push_back(l->stats[i]);
	    l->stats.clear();
	  }
	  print_summaries(stat);
	}
	}
	printf("\n"); fflush(stdout);
      }
    }
  }
  delete snp_infile;
  delete infile;
  delete [] newick;

  for (map<int,BedLine*>::iterator it=last_entry.begin(); it != last_entry.end(); ++it) {
    BedLine *l = it->second;
    if (l->pruned_tree != NULL) delete l->pruned_tree;
    if (l->orig_tree != NULL) delete l->orig_tree;
    delete &*l;
  }
  return 0;
}

int summarizeRegionNoSnp(char *filename, const char *region, 
			 set<string> inds, vector<string>statname,
			 vector<double> times) {
    TabixStream *infile;
    char c;
    char *region_chrom = NULL;
    char chrom[1000];
    vector<string> token;
    int region_start=-1, region_end=-1, start, end, sample;
    IntervalIterator<vector<double> > results;
    queue<BedLine*> bedlineQueue;
    map<int,BedLine*> bedlineMap;
    map<int,Tree*>orig_trees;
    map<int,Tree*>pruned_trees;
    map<int,Tree*>::iterator it;
    /* 
       Class BedLine contains chr,start,end, newick tree, parsed tree.
         parsed tree may be NULL if not parsing trees but otherwise will
         be populated, either by parsing the newick or an SPR operation 
         on previous tree.
       Parsed tree has recomb_node, recomb_time, coal_node, coal_time set
         (recomb_node==NULL => no recomb. Only happens in full tree at end
	 of regions analyzed by arg-sample)

       Queue bedlineQueue contains pointers to this class, will be output to
         results in order (first in, first out).
       bedlineMap<int,bedlineQueue> maps samples to pointers of the most recently 
         read instance of bedline for each sample. It points to the same objects
         as bedlineQueue (not copies).

       For each line of file:
       Read chr, start, end, sample, tree string.
       Look up bedlineMap<sample> = lastSample
       If (lastSample == NULL) {
         parse tree. Make new bedline object, add it to bedlineMap and end of bedlineQueue.
       } else if (lastSample->recomb_node != NULL) {
         apply SPR to lastSample->tree to create new parsed tree. Use this tree
         to create new bedline object, add it to bedlineMap and end of bedlineQueue.
       } else { //lastSample does not end in recomb
         assert that lastSample->end==start and lastSample->chr==chr
         update lastSample->end=end
         determine if tree string ends in recomb, update recomb_node, recomb_time, etc if so. (This is a tricky part esp if there is pruning involved).
       }
       while (first element of queue ends in recomb) {
         compute statistics for first element of queue
         add to intervaliterator
         if bedlineMap<sample>==(first element in queue), set bedlineMap<sample>=NULL
	 pop from queue
       }

      After reading all lines:
        go through queue and dump everything to intervalIterator...

     */


    infile = new TabixStream(filename, region, tabix_dir);
    if (infile->stream == NULL) return 1;

    //parse region to get region_chrom, region_start, region_end.
    // these are only needed to truncate results which fall outside
    // of the boundaries (tabix returns anything that overlaps)
    if (region != NULL) {
        split(region, "[:-]", token);
        if (token.size() != 3) {
            fprintf(stderr, "Error: bad region format; should be chr:start-end\n");
            return 1;
        }
        region_chrom = new char[token[0].size()+1];
        strcpy(region_chrom, token[0].c_str());
        region_start = atoi(token[1].c_str())-1;
        region_end = atoi(token[2].c_str());
    }
    while (EOF != (c=fgetc(infile->stream))) {
        ungetc(c, infile->stream);
        if (c!='#') break;
        while ('\n' != (c=fgetc(infile->stream))) {
           if (c==EOF) return 0;
        }
    }
    int parse_tree = (inds.size() > 0);
    if (!parse_tree) {
      for (unsigned int i=0; i < statname.size(); i++) {
         if (statname[i]!="tree") {
           parse_tree=1;
           break;
         }
      }
    }

    while (EOF != fscanf(infile->stream, "%s %i %i %i",
                         chrom, &start, &end, &sample)) {
      //      printf("got line %s %i %i %i\n", chrom, start, end, sample); fflush(stdout);
      assert('\t'==fgetc(infile->stream));
      char* newick = fgetline(infile->stream);
      chomp(newick);
      Tree *orig_tree=NULL, *pruned_tree=NULL;
      //      printf("newick: %s\n", newick);
      it = orig_trees.find(sample);
      if (it == orig_trees.end()) {  //first tree from this sample
	orig_tree = new Tree(string(newick), times);
	orig_trees[sample] = orig_tree;
	if (inds.size() > 0) {
	  pruned_tree = orig_tree->copy();
	  orig_tree->node_map = pruned_tree->prune(inds, true);
	  pruned_trees[sample] = pruned_tree;
	}
      } else {
	int parse_tree = 0;
	orig_tree = it->second;
	if (orig_tree->recomb_node == NULL) {
	  parse_tree = 1;
	  delete orig_tree;
	  orig_tree = new Tree(string(newick), times);
	  orig_trees[sample] = orig_tree;
	} else orig_tree->apply_spr();

	// set recomb_node and coal_node to next spr events indicated in newick string
        orig_tree->update_spr(newick, times);
	if (inds.size() > 0) {
	  it = pruned_trees.find(sample);
	  assert(it != pruned_trees.end());
	  pruned_tree = it->second;
	  if (parse_tree) {
	    pruned_tree = orig_tree->copy();
	    orig_tree->node_map = pruned_tree->prune(inds, true);
	    pruned_trees[sample] = pruned_tree;
	  } else pruned_tree->apply_spr();
	  pruned_tree->update_spr_pruned(orig_tree);
	}
      }

      map<int,BedLine*>::iterator it3 = bedlineMap.find(sample);
      BedLine *currline;
      if (it3 == bedlineMap.end()) {
	currline = new BedLine(chrom, start, end, sample, newick, orig_tree,
			       pruned_tree);
	bedlineMap[sample] = currline;
	bedlineQueue.push(currline);
      } else {
	currline = it3->second;
	assert(strcmp(currline->chrom, chrom)==0);
	assert(currline->end == start);
	currline->end = end;
      }

	//assume orig_tree->recomb_node == NULL is a rare occurrence that happens
	// at the boundaries of regions analyzed by arg-sample; treat these as
	// recombination events
      if (orig_tree->recomb_node == NULL ||
	  pruned_tree == NULL ||
	  pruned_tree->recomb_node != NULL) {
	scoreBedLine(currline, statname);
	bedlineMap.erase(sample);
      }

      while (bedlineQueue.size() > 0) {
	BedLine *firstline = bedlineQueue.front();
	if (firstline->stats.size() == statname.size()) {
	  processNextBedLine(firstline, &results, statname,
			     region_chrom, region_start, region_end);
	  it3 = bedlineMap.find(firstline->sample);
	  if (it3 != bedlineMap.end() && it3->second == firstline) {
	    assert(0);
	    bedlineMap.erase(firstline->sample);
	  }
	  bedlineQueue.pop();
	} else break;
      }
      delete [] newick;
    }
    infile->close();
    delete infile;

    while (bedlineQueue.size() > 0) {
      BedLine *firstline = bedlineQueue.front();
      processNextBedLine(firstline, &results, statname, 
			 region_chrom, region_start, region_end);
      delete firstline;
      bedlineQueue.pop();
    }

    if (summarize) {
        results.finish();
        checkResults(&results);
    } else {
      processNextBedLine(NULL, &results, statname, region_chrom, 
			 region_start, region_end);
    }

    it = orig_trees.begin();
    while (it != orig_trees.end()) {
      delete it->second;
      advance(it, 1);
    }
    it = pruned_trees.begin();
    while (it != pruned_trees.end()) {
      delete it->second;
      advance(it, 1);
    }
    if (region_chrom != NULL) delete[] region_chrom;
    return 0;
}

int summarizeRegion(char *snp_file, char *filename, const char *region, 
		    set<string> inds, vector<string>statname,
		    vector<double> times) {
  if (snp_file != NULL)
    return summarizeRegionBySnp(snp_file, filename, region,
				inds, statname, times);
  else return summarizeRegionNoSnp(filename, region,
				   inds, statname, times);
}


int main(int argc, char *argv[]) {
  string chr, newick;
  int opt_idx;
  vector <string>tokens;
  set <string>inds;
  char c, *region=NULL, *filename = NULL, *bedfile = NULL, *indfile = NULL,
      *timesfile=NULL;
  char *snp_file=NULL;
  int rawtrees=0, recomb=0, allele_age=0;
  vector<string> statname;
 // map<string,double> times;
  vector<double> times;
  struct option long_opts[] = {
      {"region", 1, 0, 'r'},
      {"times", 1, 0, 'm'},
      {"bed", 1, 0, 'b'},
      {"subset", 1, 0, 's'},
      {"tree", 0, 0, 'E'},
      {"tmrca", 0, 0, 'T'},
      {"recomb", 0, 0, 'R'},
      {"tmrca-half",0,0,'H'},
      {"branchlen", 0, 0, 'B'},
      {"rth", 0, 0, 'F'},
      {"popsize", 0, 0, 'P'},
      {"numsample", 0, 0, 'N'},
      {"allele-age", 0, 0, 'A'},
      {"snp-file", 1, 0, 'f'},
      {"mean", 0, 0, 'M'},
      {"stdev", 0, 0, 'S'},
      {"quantile", 1, 0, 'Q'},
      {"tabix-dir", 1, 0, 't'},
      {"no-header",0,0,'n'},
      {"help", 0, 0, 'h'},
      {0,0,0,0}};

  while (( c = (char)getopt_long(argc, argv, "r:m:b:s:TEHBRFPNAf:MSQ:t:nh",
                                 long_opts, &opt_idx)) != -1) {
      switch(c) {
      case 'r':
          region = optarg;
          break;
      case 'm':
	  timesfile = optarg;
	  break;
      case 'b':
          bedfile = optarg;
          break;
      case 's':
	  indfile = optarg;
	  break;
      case 'f':
	snp_file = optarg;
	break;
      case 'E':
          rawtrees=1;
	  statname.push_back(string("tree"));
          break;
      case 'T':
          statname.push_back(string("tmrca"));
          break;
      case 'H':
          statname.push_back(string("tmrca_half"));
          break;
      case 'B':
          statname.push_back(string("branchlen"));
          break;
      case 'F':
          statname.push_back(string("rth"));
          break;
      case 'R':
	  statname.push_back(string("recomb"));
	  recomb=1;
	  break;
      case 'P':
          statname.push_back(string("popsize"));
          break;
      case 'A':
	  allele_age=1;
          statname.push_back(string("allele_age"));
	  statname.push_back(string("inf_sites"));
          break;
      case 'N':
          getNumSample=++summarize;
          break;
      case 'M':
          getMean=++summarize;
          break;
      case 'S':
          getStdev=++summarize;;
          break;
      case 'Q': {
          getQuantiles=++summarize;
          vector<string> tokens;
          split(optarg, ',', tokens);
          for (unsigned int i=0; i < tokens.size(); i++) {
              double q=atof(tokens[i].c_str());
              //              fprintf(stderr, "getting quantile %lf\n",q);
              quantiles.push_back(q);
          }
          break;
      }
      case 't':
          tabix_dir=optarg;
          chomp(tabix_dir);
          break;
      case 'n':
          header=false;
          break;
      case 'h':
          print_help();
          return 0;
      default:
          fprintf(stderr, "Unknown option. Try ./arg-summarize --help\n");
          return 1;
      }
  }
  if (region != NULL && bedfile != NULL) {
      fprintf(stderr, "Error: --bed and --region cannot be used together.\n");
      return 1;
  }
  if (statname.size() == 0) {
     fprintf(stderr, "Error: need to specify a tree statistic\n");
     return 1;
  }

  if (summarize && statname.size()==0) {
      fprintf(stderr, "Error: need to specify a tree statistic (e.g., --tmrca, --popsize, --recomb, --allele-age, etc)\n");
      return 1;
  }
  if (summarize && rawtrees) {
     fprintf(stderr, "Error: --trees not compatible with summary statistics (--mean, --quantile, --stdev, --numsample)\n");
     return 1;
  }
  if (recomb && snp_file != NULL) {
    fprintf(stderr, "Error: cannot use --recomb and --allele-age together\n");
    return 1;
  }
  if (allele_age && snp_file == NULL) {
      fprintf(stderr, "Error: need to specify snp file with --snp to use --allele-age\n");
      return 1;
  }

  if (optind != argc-1) {
      fprintf(stderr, "Incorrect number of arguments. Try ./arg-summarize --help\n");
      return 1;
  }
  filename = argv[optind];
  if (header) {
      printf("## arg-summarize v%s\n", version);
      printf("##");
      for (int i=0; i < argc; i++) printf(" %s", argv[i]);
      printf("\n");
      printf("#chrom\tchromStart\tchromEnd");
      if (summarize==0)
          printf("\tMCMC_sample");
      if (snp_file != NULL) {
	printf("\tderAllele\tancAllele\tderFreq\tancFreq");
      }
      if (snp_file == NULL && getNumSample > 0) 
	printf("\tnumsample");
      if (summarize && snp_file) {
	printf("\tnumsample-all\tnumsample-infsites");
      }
      vector<string> stattype;
      if (snp_file == NULL) {
	stattype.push_back("");
      } else {
	stattype.push_back("-all");
	//	stattype.push_back("-derConsensus");
	stattype.push_back("-infsites");
      }

      for (unsigned int j=0; j < statname.size(); j++) {
          if (summarize==0) {
	    printf("\t%s", statname[j].c_str());
          }
	  if (statname[j] != "inf_sites") {
	  for (unsigned int k=0; k < stattype.size(); k++) {
	    for (int i=1; i <= summarize; i++) {
	      if (getMean==i) {
		printf("\t%s%s_mean", statname[j].c_str(), stattype[k].c_str());
	      } else if (getStdev==i) {
		printf("\t%s%s_stdev", statname[j].c_str(), stattype[k].c_str());
	      } else if (getQuantiles==i) {
		for (unsigned int l=0; l < quantiles.size(); l++) {
		  printf("\t%s%s_quantile_%.3f", statname[j].c_str(),
			 stattype[k].c_str(),
			 quantiles[l]);
		}
	      } 
	    }
	  }
	  }
      }
      printf("\n");
  }

  if (timesfile != NULL) {
      FILE *infile = fopen(timesfile, "r");
      double t;
      if (infile == NULL) {
	  fprintf(stderr, "Error opening %s.\n", timesfile);
	  return 1;
      }
      while (EOF != fscanf(infile, "%lf", &t)) 
	   times.push_back(t);
      std::sort(times.begin(), times.end());
      fclose(infile);
      //      fprintf(stderr, "read %i times\n", (int)times.size());
  }

  if (indfile != NULL) {
    ifstream in(indfile);
    string line;
    if (in.is_open()) {
      while ( getline(in, line) ) {
	inds.insert(line);
      }
      in.close();
      //      cerr << "subsetting to " << inds.size() << " samples from " <<
      //	indfile << "\n";
    } else {
      fprintf(stderr, "Error opening %s.\n", indfile);
      return 1;
    }
  }

  if (bedfile == NULL) {
    summarizeRegion(snp_file, filename,
		    region, inds, statname, times);
  } else {
      CompressStream bedstream(bedfile);
      char *line;
      vector<string> token;
      char *regionStr;
      if (!bedstream.stream) {
          fprintf(stderr, "error reading %s\n", bedfile);
          return 1;
      }
      while ((line = fgetline(bedstream.stream))) {
          split(line, '\t', token);
          if (token.size() < 3) {
              fprintf(stderr, "expected at least 3 files in %s\n", bedfile);
              return 1;
          }
          regionStr = new char[token[0].size()+token[1].size()+token[2].size()+3];
          int start = atoi(token[1].c_str());
          int end = atoi(token[2].c_str());
          sprintf(regionStr, "%s:%i-%i", token[0].c_str(), start+1, end);
          summarizeRegion(snp_file, filename,
                          regionStr, inds, statname, times);
          delete [] regionStr;
      }
      bedstream.close();
  }

  return 0;
}
