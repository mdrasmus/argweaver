
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


int summarizeRegion(char *filename, const char *region, 
                    set<string> inds, vector<string>statname,
                    vector<float> times) {
    TabixStream *infile;
    char *line, c;
    char *region_chrom = NULL;
    char chrom[1000];
    vector<string> token;
    int region_start=-1, region_end=-1, start, end, sample;
    Tree *tree=NULL;
    int counter=0;
    IntervalIterator<vector<double> > results;
    map<int,string> indmap;

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
    //    fprintf(stderr, "parse_tree=%i\n", parse_tree); fflush(stderr);
    while (EOF != fscanf(infile->stream, "%s %i %i %i",
                         chrom, &start, &end, &sample)) {
	double bl=-1.0;
	int orig_len = end-start;
        assert('\t' == fgetc(infile->stream));
        line = fgetline(infile->stream);
        chomp(line);
        if (parse_tree) {
          tree = new Tree(string(line));
	  if (times.size() > 0)
	      tree->correct_times(times);
        }
	if (inds.size() > 0) {
            tree->prune(inds, true);
	}

        if (region_chrom != NULL) {
            assert(strcmp(region_chrom, chrom)==0);
            if (end > region_end) end = region_end;
            if (start < region_start) start = region_start;
            assert(start < end);
        }
        vector<double>stats(statname.size());
        for (unsigned int i=0; i < statname.size(); i++) {
            if (statname[i] == "tmrca")
                stats[i] = tree->tmrca();
            else if (statname[i]=="tmrca_half")
                stats[i] = tree->tmrca_half();
	    else if (statname[i]=="branchlen") {
		if (bl < 0) {
		    stats[i] = tree->total_branchlength();
		    bl=stats[i];
		}
	    }
            else if (statname[i]=="rth")
                stats[i] = tree->rth();
            else if (statname[i]=="popsize")
                stats[i] = tree->popsize();
	    else if (statname[i]=="recomb") {
		if (bl < 0) bl = tree->total_branchlength();
		stats[i] = 1.0/(bl*(double)orig_len);
	    }
	    else if (statname[i]=="tree");
            else if (statname[i]=="allele_age") {
                fprintf(stderr, "Error: allele age not yet implemented\n");
                exit(1);
            } else {
		fprintf(stderr, "Error: unknown stat %s\n", statname[i].c_str());
		exit(1);
	    }
        }
        if (!summarize) {
            printf("%s\t%i\t%i\t%i", chrom, start, end, sample);
            for (unsigned int i=0; i < statname.size(); i++) {
               if (statname[i]=="tree") {
                  if (parse_tree) {
                     printf("\t");
                           //last arg means only print NHX if nothing pruned
                     tree->print_newick(stdout, false, true, 1, inds.size()==0);
                   } else {
                      printf("\t%s", line);
                   }
                } else {
                   printf("\t%g", stats[i]);
                }
            }
            printf("\n");
        } else {
            results.append(chrom, start, end, stats);
            counter++;
            if (counter%100==0) {
                checkResults(&results);
            }
        }
        if (parse_tree) delete tree;
        delete[] line;
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
  vector <string>tokens;
  set <string>inds;
  char c, *region=NULL, *filename = NULL, *bedfile = NULL, *indfile = NULL,
      *timesfile=NULL;
  char *allele_age_file=NULL;
  int numstat=0, rawtrees=0;
  vector<string> statname;
 // map<string,float> times;
  vector<float> times;

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
      {"recomb", 0, 0, 'R'},
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
      float t;
      if (infile == NULL) {
	  fprintf(stderr, "Error opening %s.\n", timesfile);
	  return 1;
      }
      while (EOF != fscanf(infile, "%f", &t)) 
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
