#include <iostream>
#include <fstream>
#include <assert.h>

#include "Tree.h"
#include "compress.h"
#include "getopt.h"
#include "parsing.h"

using namespace spidir;
using namespace argweaver;

void print_usage() {
    printf("smc2bed: This program converts a single smc file into a\n"
           "  bed file. The bed file format is chrom,start,end,sample,tree.\n"
           "  The tree nodes are labelled with NHX-style comments indicating\n"
           "  the nodes and times of the recombination event which leads to\n"
           "  the next tree.\n\n"
           "This program is intended for use combining multiple SMC files\n"
           "  from different MCMC samples; the pipeline for doing this is to\n"
           "  run smc2bed on each file, piping the results to sort-bed, then\n"
           "  bgzip. The resulting file can be indexed using tabix.\n\n");
    printf("Usage: ./smc2bed [OPTIONS] <smc-file>\n"
           "  smc-file can be gzipped\n"
           " OPTIONS:\n"
           " --region START-END\n"
           "   Process only these coordinates (1-based)\n"
           " --sample <sample>\n"
           "   Give the sample number for this file; this is important\n"
           "   when combining multiple smc files.\n"
	   " --times <times.txt>\n"
	   "   File giving the discrete times used; will help avoid"
           "   rounding error\n"
	   "   in branch lengths/ages\n");
}




int main(int argc, char *argv[]) {
    char c;
    int region[2]={-1,-1}, orig_start, orig_end, start, end;
    vector<string> names;
    string currstr;
    char *line = NULL;
    char chrom[1000];
    char *newick=NULL;
    Tree *tree=NULL;
    NodeSpr *spr;
    char *timesfile=NULL;
    int recomb_node, coal_node, sample=0, opt_idx;
    vector<double> times;
    struct option long_opts[] = {
        {"region", 1, 0, 'r'},
        {"sample", 1, 0, 's'},
	{"times", 1, 0, 't'},
        {"help", 0, 0, 'h'},
        {0,0,0,0}};
    while ((c = (char)getopt_long(argc, argv, "r:s:t:h", long_opts, &opt_idx))
           != -1) {
        switch (c) {
        case 'r':
            if (2 != (sscanf(optarg, "%d-%d", &region[0], &region[1]))) {
                fprintf(stderr, "error parsing region %s\n", optarg);
                return 1;
            }
            region[0]--;  //convert to 0-based
            break;
        case 's':
            sample = atoi(optarg);
            break;
	case 't':
	    timesfile = optarg;
	    break;
        case 'h':
            print_usage();
            return 0;
        case '?':
            fprintf(stderr, "unknown option. Try --help\n");
            return 1;
        }
    }
    if (optind != argc - 1) {
        fprintf(stderr, "Bad arguments. Try --help\n");
        return 1;
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
    CompressStream instream(argv[optind], "r");
    fprintf(stderr, "opening %s\n", argv[optind]);
    line = fgetline(instream.stream);
    chomp(line);
    if (strncmp(line, "NAMES", 5) != 0) {
        fprintf(stderr, "error: Expected first line of input to be NAMES");
        return 1;
    }
    split(&line[6], "\t", names);
    line= fgetline(instream.stream);
    chomp(line);
    if (strncmp(line, "REGION", 6) != 0) {
        fprintf(stderr, "error: Expected second line of input to be REGION");
    }
    if (sscanf(&line[7], "%1000s\t%d\t%d", chrom, &orig_start, &orig_end) != 3) {
        fprintf(stderr, "error parsing REGION string in second line\n");
        return 1;
    }
    while ((line = fgetline(instream.stream))) {
        chomp(line);
        if (strncmp(line, "TREE", 4)==0) {
            if (2 != sscanf(&line[5], "%d\t%d", &start, &end)) {
                fprintf(stderr, "error processing TREE line\n");
                return 1;
            }
            start--;  //0-based
            if (region[1] >= 0 && start >= region[1])
                break;
            if (region[0] >= 0 && end <= region[0]) {
                delete [] line;
                line = fgetline(instream.stream);
                if (line==NULL) break;
                if (strncmp(line, "SPR", 3)==0) {
                    delete [] line;
                    continue;
                }
                fprintf(stderr, "error: expected SPR after TREE line\n");
                return 1;
            }

            char *newick_end = line + strlen(line);
            newick = find(line+5, newick_end, '\t')+1;
            newick = find(newick, newick_end, '\t')+1;
	    if (tree == NULL || spr->recomb_node == NULL) {
	      if (tree != NULL) delete tree;
	      tree = new Tree(string(newick), times);
              spr = new NodeSpr(tree, newick, times);

	      //have to rename all leaf nodes and remove NHX comments
	      for (int i=0; i < tree->nnodes; i++) {
                int nodenum = atoi(tree->nodes[i]->longname.c_str());
                if (tree->nodes[i]->nchildren == 0) {
		  assert(nodenum >= 0 && (unsigned int)nodenum < names.size());
		  tree->nodes[i]->longname = names[nodenum];
                }
	      }
	    } else {
	      tree->apply_spr(spr);
	    }
            delete [] line;

            line=fgetline(instream.stream);
            if (line == NULL) {
	      spr->recomb_node = NULL;
	      spr->coal_node = NULL;
            } else if (strncmp(line, "SPR", 3)==0) {
	      int tempend;
	      if (5 != sscanf(&line[4], "%d\t%d\t%lf\t%d\t%lf",
			      &tempend, &recomb_node, &(spr->recomb_time),
			      &coal_node, &(spr->coal_time))) {
		fprintf(stderr, "error parsing SPR line\n");
		return 1;
	      }
	      if (tempend != end) {
		fprintf(stderr, "error: SPR pos does not equal TREE end\n");
		return 1;
	      }
	      if (region[1] >= 0 && tempend >= region[1]) {
		coal_node = -1;
		recomb_node = -1;
	      }
            } else {
	      fprintf(stderr, "error: expected SPR after TREE line\n");
	      return 1;
            }

            //now have to rename all leaf nodes and remove all NHX comments
	    char tmpStr[1000];
	    if (recomb_node >= 0) {
	      sprintf(tmpStr, "%i", recomb_node);
	      //	    char *tmpStr = itoa(recomb_node);
	      spr->recomb_node =
                  tree->nodes[tree->nodename_map[string(tmpStr)]];
	      sprintf(tmpStr, "%i", coal_node);
	      spr->coal_node = tree->nodes[tree->nodename_map[string(tmpStr)]];

	      if (spr->recomb_node->age-1 > spr->recomb_time)
		assert(0);
	      if (spr->recomb_node != tree->root) {
		if (spr->recomb_node->parent->age+1 < spr->recomb_time)
		  assert(0);
	      }
	      if (spr->coal_node->age-1 > spr->coal_time)
		assert(0);
	      if (spr->coal_node != tree->root) {
		if (spr->coal_node->parent->age+1 < spr->coal_time) {
		    assert(0);
		  }
	       }
	      if (times.size() > 0) spr->correct_recomb_times(times);
	    }
            if (region[0] >= 0 && start < region[0])
                start = region[0];
            if (region[1] >= 0 && end >= region[1]) {
                end = region[1];
		spr->recomb_node = spr->coal_node = NULL;
	    }
            printf("%s\t%i\t%i\t%i\t", chrom, start, end, sample);
            tree->write_newick(stdout, false, true, 1, spr);
            printf("\n");
            if (line == NULL) break;
            delete [] line;
        }
    }
    instream.close();
    return 0;
}
