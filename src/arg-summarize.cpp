
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
char version[10]="0.1";

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
           "--tmrca,-T\n"
           "  time to the most recent common ancestor\n"
           "--branchlen,-B\n"
           "  total branch length\n"
           "--tmrca-half,-H\n"
           "  Time for half of samples to reach a common ancestor\n"
           "--rth,-R\n"
           "  relative TMRCA Halftime\n"
           "--popsize,-P\n"
           "  estimated popsize\n"
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

void checkResults(IntervalIterator *results) {
    Interval summary=results->next();
    while (summary.start != summary.end) {
        cout << summary.chrom << "\t" << summary.start << "\t"
             << summary.end;

        for (int i=1; i <= summarize; i++) {
            if (getNumSample==i) {
                printf("\t%i", summary.num_score());
            } else if (getMean==i) {
                printf("\t%f", summary.mean());
            } else if (getStdev==i) {
                printf("\t%f", summary.stdev());
            } else if (getQuantiles==i) {
                vector<double> q = summary.quantiles(quantiles);
                for (unsigned int j=0; j < quantiles.size(); j++) {
                    printf("\t%f", q[j]);
                }
            }
        }
        cout << "\n";
        summary = results->next();
    }
}

int summarizeRegion(char *filename, const char *region,
                    set<string> inds,   double (Tree::*stat)()) {
    TabixStream *infile;
    char *line;
    char *region_chrom = NULL;
    const char *chrom;
    vector<string> token;
    int region_start=-1, region_end=-1, start, end, sample;
    Tree *tree;
    double val;
    int counter=0;
    IntervalIterator results;

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
    while ((line = fgetline(infile->stream))) {
        token.clear();
        split(line, '\t', token);
        if (token.size() != 5) {
            fprintf(stderr, "Bad bedfile format; expected 5 columns per row"
                    "with chrom,start,end,sample,tree\n");
            return 1;
        }
        chrom = token[0].c_str();
        start = atoi(token[1].c_str());
        end = atoi(token[2].c_str());
        sample = atoi(token[3].c_str());
        tree = new Tree(token[4]);

	if (inds.size() > 0) {
            tree->prune(inds, true);
	}

        if (region_chrom != NULL) {
            assert(strcmp(region_chrom, chrom)==0);
            if (end > region_end) end = region_end;
            if (start < region_start) start = region_start;
            assert(start < end);
        }
        if (stat != NULL) {
            val = (tree->*stat)();

            if (!summarize) {
                printf("%s\t%i\t%i\t%i\t%f\n", chrom, start, end, sample, val);
            } else {
                results.append(chrom, start, end, val);
                counter++;
                if (counter%100==0) {
                    checkResults(&results);
                }
            }
            //            cout << "\t" << val;
        } else {
            printf("%s\t%i\t%i\t%i\t\t%s", chrom, start, end, sample,
                   token[4].c_str());
                //            cout << "\t" << token[4];
        }
        delete tree;
        delete [] line;
    }
    infile->close();
    delete infile;

    if (summarize) {
        results.finish();
        checkResults(&results);
    }
    if (region_chrom != NULL) delete[] region_chrom;
    return 0;
}


int main(int argc, char *argv[]) {
  string chr, newick;
  int opt_idx;
  //  Tree *tree;
  vector <string>tokens;
  set <string>inds;
  char c, *region=NULL, *filename = NULL, *bedfile = NULL, *indfile = NULL;
  int numstat=0;
  double (Tree::*stat)() = NULL;
  char statname[100];
  struct option long_opts[] = {
      {"region", 1, 0, 'r'},
      {"bed", 1, 0, 'b'},
      {"subset", 1, 0, 's'},
      {"tmrca", 0, 0, 'T'},
      {"tmrca-half",0,0,'H'},
      {"branchlen", 0, 0, 'B'},
      {"rth", 0, 0, 'R'},
      {"popsize", 0, 0, 'P'},
      {"numsample", 0, 0, 'N'},
      {"mean", 0, 0, 'M'},
      {"stdev", 0, 0, 'S'},
      {"quantile", 1, 0, 'Q'},
      {"tabix-dir", 1, 0, 't'},
      {"no-header",0,0,'n'},
      {"help", 0, 0, 'h'},
      {0,0,0,0}};

  while (( c = (char)getopt_long(argc, argv, "r:b:THBRPNMSQ:t:nh",
                                 long_opts, &opt_idx)) != -1) {
      switch(c) {
      case 'r':
          region = optarg;
          break;
      case 'b':
          bedfile = optarg;
          break;
      case 's':
	  indfile = optarg;
	  break;
      case 'T':
          if (stat == NULL) {
            stat = &Tree::tmrca;
            strcpy(statname, "tmrca");
          }
      case 'H':
          if (stat == NULL) {
            stat = &Tree::tmrca_half;
            strcpy(statname, "tmrca_half");
          }
      case 'B':
          if (stat == NULL) {
            stat = &Tree::total_branchlength;
            strcpy(statname, "branchlen");
          }
      case 'R':
          if (stat == NULL) {
            stat = &Tree::rth;
            strcpy(statname, "rth");
          }
      case 'P':
          if (stat == NULL) {
            stat = &Tree::popsize;
            strcpy(statname, "popsize");
          }
          numstat++;
          if (numstat > 1) {
              fprintf(stderr, "Error: can only compute one statistic at a time\n");
              exit(-1);
          }
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
              //              fprintf(stderr, "getting quantile %f\n",q);
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

  if (summarize && stat==NULL) {
      fprintf(stderr, "Error: need to specify a tree statistic\n");
      return 1;
  }
  //TODO: note that --prune removes NHX tags, cannot be used for recomb rate (or can it?)

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
      if (stat==NULL)
          printf("\ttree");
      else {
          if (summarize==0) {
              printf("\t%s", statname);
          }
          for (int i=1; i <= summarize; i++) {
              if (getNumSample==i) {
                  printf("\t%s_numsample", statname);
              } else if (getMean==i) {
                printf("\t%s_mean", statname);
              } else if (getStdev==i) {
                  printf("\t%s_stdev", statname);
              } else if (getQuantiles==i) {
                  for (unsigned int j=0; j < quantiles.size(); j++) {
                      printf("\t%s_quantile_%.3f", statname, quantiles[j]);
                  }
              }
          }
      }
      printf("\n");
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
      summarizeRegion(filename, region, inds, stat);
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
          summarizeRegion(filename, regionStr, inds, stat);
          delete [] regionStr;
      }
  }

  return 0;
}
