
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
	   //           "--allele-age,-A <snp file>\n"
	   //           "  Compute allele age for all CG snps in the region. snp_file\n"
	   //           "  needs to contain phased SNP information in the format\n"
	   //           " chr,start,end,AAAAACAAAA\n"
	   //           " a header like the following needs to identify haplotype order:\n"
	   //           " #NAMES NA06985_1       NA06985_2       NA06994_1       NA06994_2.."
           "SUMMARY OPTIONS:\n"
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
	newick = (char*)malloc((strlen(nwk)+1)*sizeof(char));
	strcpy(newick, nwk);
    };
  ~BedLine() {
    //      if (orig_tree != NULL) delete orig_tree;
    //      if (pruned_tree != NULL) delete pruned_tree;
      delete [] chrom;
      free(newick);
  }
  char *chrom;
  int start;
  int end;
  int sample;
  Tree *orig_tree;
  Tree *pruned_tree;
  char *newick;
  vector<double> stats;
};


void scoreBedLine(BedLine *line, vector<string> statname) {
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
    else if (statname[i]=="allele_age") {
      fprintf(stderr, "Error: allele age not yet implemented\n");
      exit(1);
    } else {
      fprintf(stderr, "Error: unknown stat %s\n", statname[i].c_str());
      exit(1);
    }
  }
}



void processNextBedLine(BedLine *line, IntervalIterator<vector<double> > *results,
			vector<string> statname, 
			char *region_chrom, int region_start, int region_end) {
  static int counter=0;

  if (line->stats.size() == 0) 
    scoreBedLine(line, statname);
  if (region_chrom != NULL) {
    assert(strcmp(region_chrom, line->chrom)==0);
    if (line->end > region_end) line->end = region_end;
    if (line->start < region_start) line->start = region_start;
    assert(line->start < line->end);
  }
  if (!summarize) {
    printf("%s\t%i\t%i\t%i", line->chrom, line->start, line->end, line->sample);
    for (unsigned int i=0; i < statname.size(); i++) {
      if (statname[i]=="tree") {
	printf("\t");
	printf("%s", line->newick);
      } else {
	printf("\t%g", line->stats[i]);
      }
    }
    printf("\n");
  } else {
    results->append(line->chrom, line->start, line->end, line->stats);
    counter++;
    if (counter%100==0) {
      checkResults(results);
    }
  }
}

int summarizeRegion(char *filename, const char *region, 
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
         (recomb_node==NULL => no recomb at end of region)

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
    if (infile == NULL) {
        fprintf(stderr, "Error opening %s, region=%s\n",
                filename, region==NULL ? "NULL" : region);
        return 1;
    }
    if (infile->stream == NULL) return 1;

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
	//	printf("ORIG: ");
	//	orig_tree->print_newick(stdout);
	//	printf("\n");
	if (inds.size() > 0) {
	  pruned_tree = orig_tree->copy();
	  //	  orig_tree->node_map = *pruned_tree->prune(inds, true);
	  orig_tree->node_map = pruned_tree->prune(inds, true);
	  pruned_trees[sample] = pruned_tree;
	  /*	  printf("PRUN: ");
	  pruned_tree->print_newick(stdout);
	  printf("\n");
	  fflush(stdout);
	  for (int i=0; i < orig_tree->nnodes; i++)
	    printf("node_map[%i]=%i\n", i, orig_tree->node_map.nm[i]);
	    fflush(stdout);*/
	}
      } else {
	int parse_tree = 0;
	orig_tree = it->second;
	if (orig_tree->recomb_node == NULL) {
	  //	  printf("recomb_node is NULL; reading new tree\n"); fflush(stdout);
	  assert(0);  // TODO: REMOVE THIS! IT IS OK TO BE HERE.
	  parse_tree = 1;
	  delete orig_tree;
	  orig_tree = new Tree(string(newick), times);
	  orig_trees[sample] = orig_tree;
	} else orig_tree->apply_spr();

	// set recomb_node and coal_node to next spr events indicated in newick string
        orig_tree->update_spr(newick, times);
	/*	printf("SPR1: ");
	orig_tree->print_newick(stdout);
	printf("\n");
	printf("done apply_update_spr\n");
	printf("update spr recomb_node=%i coal_node=%i\n",
	       orig_tree->recomb_node==NULL ? -1 : orig_tree->recomb_node->name,
	       orig_tree->coal_node==NULL ? -1 : orig_tree->coal_node->name);*/
	if (inds.size() > 0) {
	  it = pruned_trees.find(sample);
	  assert(it != pruned_trees.end());
	  pruned_tree = it->second;
	  if (parse_tree) {
	    assert(0); // TODO: DELETE THIS TOO. IT IS USUALLY OK TO BE HERE.
	    delete pruned_tree;
	    pruned_tree = orig_tree->copy();
	    //	    orig_tree->node_map = *pruned_tree->prune(inds, true);
	    orig_tree->node_map = pruned_tree->prune(inds, true);
	    pruned_trees[sample] = pruned_tree;
	  } else if (pruned_tree->recomb_node != NULL) {
	    pruned_tree->apply_spr();
	  }
	  /*	  printf("SPR2: ");
	  pruned_tree->print_newick(stdout); printf("\n");
	  fflush(stdout);*/

	  /*	  for (int f=0; f < orig_tree->nnodes; f++) {
            printf("map[%i]=%i\n", f, orig_tree->node_map.nm[f]);
          } 
	  fflush(stdout);*/

	  //map recomb_node and coal_node onto pruned tree
	  if (orig_tree->recomb_node == NULL) {
	    pruned_tree->recomb_node = pruned_tree->coal_node = NULL;
	  } else {
	    int num=orig_tree->node_map.nm[orig_tree->recomb_node->name];
	    //	    printf("num=%i\n", num);
	    if (num == -1) {
	      pruned_tree->recomb_node = pruned_tree->coal_node = NULL;
	    } else {
	      assert(num>=0);
 	      pruned_tree->recomb_node = pruned_tree->nodes[num];
	      pruned_tree->recomb_time = orig_tree->recomb_time;
	      num = orig_tree->node_map.nm[orig_tree->coal_node->name];
	      //	      printf("num2=%i\n", num);
	      if (num == -1) {  
		// coal node does not map; need to trace back until it does
		Node *n = orig_tree->coal_node;
		while (orig_tree->node_map.nm[n->name] == -1) {
		  //should never be root here; root should always map to pruned tree
		  assert(n->parent != NULL);
		  n = n->parent;
		}
		assert(orig_tree->coal_time-1 <= n->age);
		pruned_tree->coal_time = n->age;
		pruned_tree->coal_node = pruned_tree->nodes[orig_tree->node_map.nm[n->name]];
	      } else {
		assert(num >= 0);
		pruned_tree->coal_node = pruned_tree->nodes[num];
		pruned_tree->coal_time = orig_tree->coal_time;
	      }
	    }
	  }
	  /*	  printf("update spr2 recomb_node=%i coal_node=%i\n",
		 pruned_tree->recomb_node==NULL ? -1 : pruned_tree->recomb_node->name,
		 pruned_tree->coal_node==NULL ? -1 : pruned_tree->coal_node->name);*/
	  if (pruned_tree->recomb_node != NULL)
	    assert(pruned_tree->coal_node != NULL);
	}
      }

      map<int,BedLine*>::iterator it3 = bedlineMap.find(sample);
      BedLine *currline;
      if (it3 == bedlineMap.end()) {
	//	printf("new bedline map for sample %i (%s\t%i\t%i)\n", sample, chrom, start-offset, end-offset); fflush(stdout);
	currline = new BedLine(chrom, start, end, sample, newick, orig_tree,
			       pruned_tree);
	bedlineMap[sample] = currline;
	bedlineQueue.push(currline);
      } else {
	currline = it3->second;
	//	printf("found entry for sample %i\n", sample); fflush(stdout);
	//	printf("entry has sample %i start=%i end=%i\n", currline->sample, currline->start-offset, currline->end-offset); fflush(stdout);
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
	//	  printf("replaced bedline map for sample %i (old=(%i,%i) new=(%i,%i)\n", sample, oldstart-offset, oldend-offset, currline->start-offset, currline->end-offset); fflush(stdout);
      }

      while (bedlineQueue.size() > 0) {
	BedLine *firstline = bedlineQueue.front();
	if (firstline->stats.size() == statname.size()) {
	  processNextBedLine(firstline, &results, statname,
			     region_chrom, region_start, region_end);
	  //	  printf("processing sample %i\t(%s\t%i\t%i)\n", firstline->sample,
	  //		 firstline->chrom, firstline->start-offset, firstline->end-offset);
	  it3 = bedlineMap.find(firstline->sample);
	  if (it3 != bedlineMap.end() && it3->second == firstline) {
	    assert(0);
	    //	    printf("erasing bedlineMap entry for sample %i (%i,%i)\n", firstline->sample, firstline->start-offset, firstline->end-offset);
	    bedlineMap.erase(firstline->sample);
	  }
	  delete &*firstline;
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
      delete &*firstline;
      bedlineQueue.pop();
    }

    if (summarize) {
        results.finish();
        checkResults(&results);
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


int main(int argc, char *argv[]) {
  string chr, newick;
  int opt_idx;
  vector <string>tokens;
  set <string>inds;
  char c, *region=NULL, *filename = NULL, *bedfile = NULL, *indfile = NULL,
      *timesfile=NULL;
  char *allele_age_file=NULL;
  int numstat=0, rawtrees=0, recomb=0;
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
      {"allele-age", 1, 0, 'A'},
      {"mean", 0, 0, 'M'},
      {"stdev", 0, 0, 'S'},
      {"quantile", 1, 0, 'Q'},
      {"tabix-dir", 1, 0, 't'},
      {"no-header",0,0,'n'},
      {"help", 0, 0, 'h'},
      {0,0,0,0}};

  while (( c = (char)getopt_long(argc, argv, "r:m:b:s:TEHBRFPNA:MSQ:t:nh",
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
          numstat++;
          if (numstat > 1) {
              fprintf(stderr, "Error: can only compute one statistic at a time\n");
              exit(-1);
          }
          break;
      case 'A':
          allele_age_file = optarg;
          statname.push_back(string("allele_age"));
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
  if (summarize && (statname.size()==0 && allele_age_file==NULL)) {
      fprintf(stderr, "Error: need to specify a tree statistic\n");
      return 1;
  }
  if (summarize && rawtrees) {
     fprintf(stderr, "Error: --trees not compatible with summary statistics (--mean, --quantile, --stdev, --numsample)\n");
     return 1;
  }
  if (recomb && allele_age_file != NULL) {
    fprintf(stderr, "Error: cannot use --recomb and --allele-age together\n");
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
      /*      if (allele_age_file != NULL) {
          printf("\tder\tanc\tnsample_clear\tnsample_switch\tnsample_multmut");
	  }*/
      if (getNumSample > 0) printf("\tnumsample");
      for (unsigned int j=0; j < statname.size(); j++) {
          if (summarize==0) {
              printf("\t%s", statname[j].c_str());
          }
          for (int i=1; i <= summarize; i++) {
              if (getMean==i) {
                  printf("\t%s_mean", statname[j].c_str());
              } else if (getStdev==i) {
                  printf("\t%s_stdev", statname[j].c_str());
              } else if (getQuantiles==i) {
                  for (unsigned int k=0; k < quantiles.size(); k++) {
                      printf("\t%s_quantile_%.3f", statname[j].c_str(),
                             quantiles[k]);
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
      cerr << "subsetting to " << inds.size() << " samples from " <<
	indfile << "\n";
    } else {
      fprintf(stderr, "Error opening %s.\n", indfile);
      return 1;
    }
  }

  if (bedfile == NULL) {
      summarizeRegion(filename, //allele_age_file,
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
          summarizeRegion(filename, //allele_age_file,
                          regionStr, inds, statname, times);
          delete [] regionStr;
      }
      bedstream.close();
  }

  return 0;
}
