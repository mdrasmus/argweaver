
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

#define VERSION_INFO "arg-summarize 0.3"

int summarize=0;
int getNumSample=0;
int getMean=0;
int getStdev=0;
int getQuantiles=0;
vector <double> quantiles;
vector<string> node_dist_leaf1;
vector<string> node_dist_leaf2;

const int EXIT_ERROR = 1;

class Config
{
public:
    Config()
    {
        make_parser();
    }
    void make_parser() {
        config.clear();

        config.add(new ConfigParamComment("Input files"));
        config.add(new ConfigParam<string>
                   ("-a", "--arg-file", "<file.bed.gz>", &argfile,
                    "Bed file containing args sampled by ARGweaver. Should"
                    " be created with smc2bed and sorted with sort-bed. If"
                    " using --region or --bedfile, also needs to be gzipped"
                    " and tabix'd"));
        config.add(new ConfigParam<string>
                   ("-r", "--region", "<chr:start-end>", &region,
                    "region to retrieve statistics from (1-based coords)"));
        config.add(new ConfigParam<string>
                   ("-b", "--bed-file", "<file.bed>", &bedfile,
                    "regions to retrieve statistics from (alternative to "
                    "region)"));
        config.add(new ConfigParam<string>
                   ("-s", "--subset", "<hap_list.txt>", &indfile,
                    "file with list of leafs to keep (rest will be pruned)"));
        config.add(new ConfigParam<string>
                   ("-f", "--snp-file", "<snp_file.bed>", &snpfile,
                    "Compute statistics for specific SNPs. Each output row will"
                    " contain information about each SNP (alleles and"
                    " frequencies), and statistics will be computed in two ways-"
                    " once across all MCMC samples, and once only for samples"
                    " where the site pattern for the SNP agrees with an infinite"
                    " sites model."
                    " The SNP file should have a header like this:\n"
                    " #NAMES NA06985_1       NA06985_2       NA06994_1 ...\n"
                    " and then each line should be tab-delimited with the format"
                    " chr,start,end,AAAACAAAAA, where the last column gives the"
                    " alleles for each haplotype in the order indicated in the"
                    " header. The file needs to be sorted and indexed with"
                    " tabix."));
        config.add(new ConfigParam<string>
                   ("-m", "--time-file", "<times.txt>", &timefile,
                    "File with list of discretized times used in trees. If"
                    " given, tree events will be set to nearest time (eliminates"
                    " some rounding error)."));

        config.add(new ConfigParamComment("Statistics to retrieve"));
        config.add(new ConfigSwitch
                   ("-E", "--tree", &rawtrees,
                    "output newick tree strings (cannot use summary options"
                    " with this)"));
        config.add(new ConfigSwitch
                   ("-T", "--tmrca", &tmrca,
                    "time to the most recent common ancestor"));
        config.add(new ConfigSwitch
                   ("-B", "--branchlen", &branchlen, "total branch length"));
        config.add(new ConfigSwitch
                   ("-R", "--recomb", &recomb,
                    "recombination rate per generation/bp (NOTE this is not"
                    " particularly well estimated by the SMC)"));
        config.add(new ConfigSwitch
                   ("-K", "--breaks", &breaks,
                    "recombinations per bp (not normalized by tree size"));
        config.add(new ConfigSwitch
                   ("-H", "--tmrca-half", &tmrca_half,
                    "time for half of samples to reach a common ancestor"));
        config.add(new ConfigSwitch
                   ("-F", "--rth", &rth,
                    "relative TMRCA Halftime (tmrca_half/tmrca)"));
        config.add(new ConfigSwitch
                   ("-P", "--popsize", &popsize,
                    "popsize estimated by coal rates in local tree"
                    " (tends to be very noisy)"));
        config.add(new ConfigSwitch
                   ("-A", "--allele-age", &allele_age,
                    "(requires --snp). Compute allele age for SNPs"));
        config.add(new ConfigParam<string>
                   ("-D", "--node-dist", "<leaf1,leaf2;leaf1,leaf3;...>",
                    &node_dist,
                    "distance between pairs of leafs. In this example return"
                    " distance between leaf1->leaf2, and leaf1->leaf3"));
        config.add(new ConfigSwitch
                   ("-Z", "--zero-len", &zero,
                    "number of branches of length zero"));
        config.add(new ConfigSwitch
                   ("-C", "--coalcounts", &coalcounts,
                    "number of coal events at each discretized time point"
                    " (requires --timefile)"));
        config.add(new ConfigSwitch
                   ("-N", "--numsample", &numsample,
                    "number of MCMC samples covering each region"));

        config.add(new ConfigParamComment("Summary options (if not given,"
                                          " statistics will be output for each"
                                          " MCMC sample)"));
        config.add(new ConfigSwitch
                   ("-M", "--mean", &mean,
                    "return mean across all MCMC samples"));
        config.add(new ConfigSwitch
                   ("-S", "--stdev", &stdev,
                    "return standard deviation across all NCMC samples"));
        config.add(new ConfigParam<string>
                   ("-Q", "--quantile", "<q1,q2,q3,...>", &quantile,
                    "return the requested quantiles for each samples"));

        config.add(new ConfigParamComment("Misceallaneous"));
        config.add(new ConfigSwitch
                   ("-n", "--no-header", &noheader, "Do not output header"));
        config.add(new ConfigParam<string>
                   ("-t", "--tabix-dir", "<tabix dir>", &tabix_dir,
                    "Specify the directory of the tabix executable"));
        config.add(new ConfigSwitch
                   ("-v", "--version", &version, "display version information"));
        config.add(new ConfigSwitch
                   ("-h", "--help", &help, "display help information"));
    }

    int parse_args(int argc, char **argv) {
        if (!config.parse(argc, (const char**) argv)) {
            if (argc < 2)
                config.printHelp();
            return EXIT_ERROR;
        }
        if (help) {
            config.printHelp();
            return EXIT_ERROR;
        }
        if (version) {
            printf(VERSION_INFO);
            return EXIT_ERROR;
        }
        return 0;
    }
    ConfigParser config;

    string argfile;
    string region;
    string bedfile;
    string indfile;
    string snpfile;
    string timefile;

    bool rawtrees;
    bool tmrca;
    bool branchlen;
    bool recomb;
    bool breaks;
    bool tmrca_half;
    bool rth;
    bool popsize;
    bool allele_age;
    string node_dist;
    bool zero;
    bool coalcounts;
    bool numsample;

    bool mean;
    bool stdev;
    string quantile;

    bool noheader;
    string tabix_dir;
    bool version;
    bool help;
};

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
                        vector<double> q =
                            compute_quantiles(tmpScore, quantiles);
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
            SprPruned *trees=NULL) :
        start(start), end(end), sample(sample),
        trees(trees) {
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
        if (newick != NULL)
            free(newick);
        //        delete trees;
    }
    char *chrom;
    int start;
    int end;
    int sample;
    SprPruned *trees;
    char *newick;
    vector<double> stats;
    char derAllele, otherAllele;
    int derFreq, otherFreq;
    int infSites;
};


void scoreBedLine(BedLine *line, vector<string> &statname, vector<double> times,
                  double allele_age=-1, int infsites=-1) {
    Tree * tree = (line->trees->pruned_tree != NULL ?
                   line->trees->pruned_tree :
                   line->trees->orig_tree);
    double bl=-1.0;
    int node_dist_idx=0;
    if (line->stats.size() == statname.size()) return;
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
        else if (statname[i]=="breaks") {
            line->stats[i] = 1.0/((double)(line->end - line->start));
        }
        else if (statname[i]=="zero_len") {
            line->stats[i] = tree->num_zero_branches();
        }
        else if (statname[i]=="tree") {
            if (line->trees->pruned_tree != NULL) {
                string tmp =
                    line->trees->pruned_tree->format_newick(false, true, 1,
                                                     &line->trees->pruned_spr);
                //pruned tree will be fewer characters than whole tree
                sprintf(line->newick, "%s", tmp.c_str());
            }
        }
        else if (statname[i]=="allele_age")
            line->stats[i] = allele_age;
        else if (statname[i]=="inf_sites")
            line->stats[i] = (double)infsites;
        else if (statname[i].substr(0, 9)=="node_dist") {
            line->stats[i] =
                tree->distBetweenLeaves(node_dist_leaf1[node_dist_idx],
                                        node_dist_leaf2[node_dist_idx]);
            node_dist_idx++;
        }
        else if (statname[i].substr(0, 10)=="coalcount.") {
            vector<double>coal_counts = tree->coalCounts(times);
            for (unsigned int j=0; j < coal_counts.size(); j++) {
                assert(i+j < statname.size() &&
                       statname[i+j].substr(0,10)=="coalcount.");
                line->stats[i+j] = coal_counts[j];
            }
            i += coal_counts.size()-1;
        }
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
        //could generalize and sort by start, but for this purpose we are
        //only comparing entries with same start
        assert(l1->start == l2->start);
        if (l1->end == l2->end)
            return (l1->sample < l2->sample);
        return l1->end < l2->end;
    }
};

void processNextBedLine(BedLine *line,
                        IntervalIterator<vector<double> > *results,
                        vector<string> &statname,
                        char *region_chrom, int region_start, int region_end,
                        vector<double> times) {
    static int counter=0;
    static list<BedLine*> bedlist;

    if (line != NULL) {
        if (line->stats.size() == 0)
            scoreBedLine(line, statname, times);
        if (region_chrom != NULL) {
            assert(strcmp(region_chrom, line->chrom)==0);
            if (line->end > region_end) line->end = region_end;
            if (line->start < region_start) line->start = region_start;
            assert(line->start < line->end);
        }
    }
    if (!summarize) {
        // this little bit of code ensures that output is sorted. The
        // summarizeRegion code should ensure that it comes here sorted by
        // start coordinate, but not necessarily by end coordinate (at least
        // not when subsetting by individuals).
        // so, this stores BedLine elements just long enough until a new start
        // coordinate is encountered, then sorts them all by end coordinate (and
        // sample), then outputs them all.  Need to call this function one last
        // time with line==NULL to output the final set of rows.
        if (bedlist.size() > 0 &&
            (line==NULL || bedlist.front()->start < line->start)) {
            bedlist.sort(CompareBedLineEnd());
            for (list<BedLine*>::iterator it=bedlist.begin();
                 it != bedlist.end(); ++it) {
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
                printf("\n");
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
        return 0;
    }


    void scoreAlleleAge(BedLine *l, vector<string> statname,
                        vector<double> times) {
        int num_derived, total;
        assert(l->start < coord);
        assert(l->end >= coord);
        Tree *t;
        if (l->trees->pruned_tree != NULL)
            t = l->trees->pruned_tree;
        else t = l->trees->orig_tree;

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
        set<Node*> derived;
        for (set<string>::iterator it=derived_in_tree.begin();
             it != derived_in_tree.end(); ++it) {
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
            set <Node*>lca2 = t->lca(derived2);
            if (lca2.size() < lca.size()) {
                major_is_derived=1;
                lca = lca2;
            }
        }
        double age=0.0;
        for (set<Node*>::iterator it4=lca.begin(); it4 != lca.end(); ++it4) {
            Node *n = *it4;
            assert(n != t->root);
            double tempage = n->age + (n->parent->age - n->age)/2;  //midpoint
            if (tempage > age) age = tempage;
        }
        if (num_derived == 0 || total-num_derived == 0) age = -1;
        scoreBedLine(l, statname, times, age, lca.size()==1);
        l->derAllele = (major_is_derived ? allele2 : allele1);
        l->otherAllele = (major_is_derived ? allele1 : allele2);
        l->derFreq = (major_is_derived ? total-num_derived : num_derived);
        l->otherFreq = (major_is_derived ? num_derived : total - num_derived);
        l->infSites = (lca.size() == 1);

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



void print_summaries(vector<double> &stat) {
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


int summarizeRegionBySnp(Config *config, const char *region,
                         set<string> inds, vector<string> statname,
                         vector<double> times) {
    TabixStream snp_infile(config->snpfile, region, config->tabix_dir);
    TabixStream infile(config->argfile, region, config->tabix_dir);
    vector<string> token;
    map<int,BedLine*> last_entry;
    map<int,BedLine*>::iterator it;
    char chrom[1000], c;
    int start, end, sample;
    BedLine *l=NULL;

    if (snp_infile.stream == NULL) return 1;
    if (infile.stream == NULL) return 1;
    while (EOF != (c=fgetc(infile.stream))) {
        ungetc(c, infile.stream);
        if (c != '#') break;
        while ('\n' != (c=fgetc(infile.stream))) {
            if (c==EOF) return 0;
        }
    }
    SnpStream snpStream = SnpStream(&snp_infile);
    if (EOF==fscanf(infile.stream, "%s %i %i %i",
                    chrom, &start, &end, &sample)) return 0;
    assert('\t' == fgetc(infile.stream));
    char *newick = fgetline(infile.stream);
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
                snpStream.scoreAlleleAge(l, statname, times);
                bedlist.push_back(l);
            }
        }
        //now look through bed file until we get to one that starts after SNP
        while (start != -1 && snpStream.coord > start) {
            it = last_entry.find(sample);
            if (it == last_entry.end() ||
                it->second->trees->orig_spr.recomb_node == NULL) {
                SprPruned *trees;
                if (it != last_entry.end()) {
                    l = it->second;
                    delete l->trees;
                    delete &*l;
                }
                trees = new SprPruned(newick, inds, times);
                l = new BedLine(chrom, start, end, sample, newick, trees);
                last_entry[sample] = l;
            } else {
                l = it->second;
                l->trees->update(newick, inds, times);
                free(l->newick);
                l->newick = (char*)malloc((strlen(newick)+1)*sizeof(char));
                strcpy(l->newick, newick);
                l->start = start;
                l->end = end;
            }
            if (snpStream.coord <= end) {
                snpStream.scoreAlleleAge(l, statname, times);
                bedlist.push_back(l);
            }
            if (4 != fscanf(infile.stream, "%s %i %i %i",
                            chrom, &start, &end, &sample))
                start = -1;
            else {
                assert('\t' == fgetc(infile.stream));
                delete [] newick;
                newick = fgetline(infile.stream);
                chomp(newick);
            }
        }
        if (bedlist.size() > 0) {
            if (summarize == 0) {
                bedlist.sort(CompareBedLineSample());
                for (list<BedLine*>::iterator it=bedlist.begin();
                     it != bedlist.end(); ++it) {
                    BedLine *l = *it;
                    printf("%s\t%i\t%i\t%i\t%c\t%c\t%i\t%i", l->chrom,
                           snpStream.coord-1, snpStream.coord, l->sample,
                           l->derAllele, l->otherAllele, l->derFreq,
                           l->otherFreq);
                    for (unsigned int i=0; i < statname.size(); i++) {
                        if (statname[i]=="tree") {
                            printf("\t%s", l->newick);
                        } else if (statname[i]=="infSites") {
                            printf("\t%i", (int)(l->stats[i]==1));
                        } else {
                            printf("\t%g", l->stats[i]);
                        }
                    }
                    printf("\n");
                    l->stats.clear();
                }
            } else {
                //now output three versions- one for all samples,
                //one for same derived allele, one for infinite sites
                BedLine* first = *(bedlist.begin());
                int same=0, diff=0, infsites=0, derConstCount,
                    derFreq, otherFreq;
                char derAllele, otherAllele;
                for (list<BedLine*>::iterator it=bedlist.begin();
                     it != bedlist.end(); ++it) {
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
                printf("%s\t%i\t%i\t%c\t%c\t%i\t%i\t%i\t%i",
                       l->chrom, snpStream.coord-1, snpStream.coord,
                       derAllele, otherAllele, derFreq, otherFreq,
                       (int)bedlist.size(), infsites);
                for (unsigned int i=0; i < statname.size(); i++) {
                    if (statname[i] != "inf_sites") {
                        // first compute stats across all
                        vector<double> stat;
                        for (list<BedLine*>::iterator it=bedlist.begin();
                             it != bedlist.end(); ++it) {
                            BedLine *l = *it;
                            stat.push_back(l->stats[i]);
                        }
                        print_summaries(stat);

                        stat.clear();
                        //now stats for infinite sites set
                        for (list<BedLine*>::iterator it=bedlist.begin();
                             it != bedlist.end(); ++it) {
                            BedLine *l = *it;
                            if (l->infSites)
                                stat.push_back(l->stats[i]);
                            l->stats.clear();
                        }
                        print_summaries(stat);
                    }
                }
                printf("\n");
            }
        }
    }
    delete [] newick;

    for (map<int,BedLine*>::iterator it=last_entry.begin();
         it != last_entry.end(); ++it) {
        BedLine *l = it->second;
        delete(l->trees);
        delete l;
    }
    return 0;
}


int summarizeRegionNoSnp(Config *config, const char *region,
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
    map<int,SprPruned*> trees;
    map<int,SprPruned*>::iterator it;
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
         parse tree. Make new bedline object, add it to bedlineMap and end
        of bedlineQueue.
      } else if (lastSample->recomb_node != NULL) {
        apply SPR to lastSample->tree to create new parsed tree. Use this tree
        to create new bedline object, add it to bedlineMap and end of
        bedlineQueue.
      } else { //lastSample does not end in recomb
        assert that lastSample->end==start and lastSample->chr==chr
        update lastSample->end=end
        determine if tree string ends in recomb, update recomb_node,
        recomb_time, etc if so.
        (This is a tricky part esp if there is pruning involved).
      }
      while (first element of queue ends in recomb) {
        compute statistics for first element of queue
        add to intervaliterator
        if bedlineMap<sample>==(first element in queue),
          set bedlineMap<sample>=NULL
        pop from queue
      }

      After reading all lines:
      go through queue and dump everything to intervalIterator...

    */


    infile = new TabixStream(config->argfile, region, config->tabix_dir);
    if (infile->stream == NULL) return 1;

    //parse region to get region_chrom, region_start, region_end.
    // these are only needed to truncate results which fall outside
    // of the boundaries (tabix returns anything that overlaps)
    if (region != NULL) {
        split(region, "[:-]", token);
        if (token.size() != 3) {
            fprintf(stderr,
                    "Error: bad region format (%s); should be chr:start-end\n",
                    region);
            return 1;
        }
        region_chrom = new char[token[0].size()+1];
        //remove commas from integer coordinates in case they are
        // copied from browser
        token[1].erase(std::remove(token[1].begin(), token[1].end(), ','),
                       token[1].end());
        token[2].erase(std::remove(token[2].begin(), token[2].end(), ','),
                       token[2].end());
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
        assert('\t'==fgetc(infile->stream));
        char* newick = fgetline(infile->stream);
        chomp(newick);
        it = trees.find(sample);
        if (it == trees.end())   //first tree from this sample
            trees[sample] = new SprPruned(newick, inds, times);
        else trees[sample]->update(newick, inds, times);

        map<int,BedLine*>::iterator it3 = bedlineMap.find(sample);
        BedLine *currline;
        if (it3 == bedlineMap.end()) {
            currline = new BedLine(chrom, start, end, sample, newick,
                                   trees[sample]);
            bedlineMap[sample] = currline;
            bedlineQueue.push(currline);
        } else {
            currline = it3->second;
            assert(strcmp(currline->chrom, chrom)==0);
            assert(currline->end == start);
            currline->end = end;
        }

        //assume orig_spr.recomb_node == NULL is a rare occurrence that happens
        // at the boundaries of regions analyzed by arg-sample; treat these as
        // recombination events
        if (trees[sample]->orig_spr.recomb_node == NULL ||
            trees[sample]->pruned_tree == NULL ||
            trees[sample]->pruned_spr.recomb_node != NULL) {
            scoreBedLine(currline, statname, times);
            bedlineMap.erase(sample);
        }

        while (bedlineQueue.size() > 0) {
            BedLine *firstline = bedlineQueue.front();
            if (firstline->stats.size() == statname.size()) {
                processNextBedLine(firstline, &results, statname,
                                   region_chrom, region_start, region_end,
                                   times);
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
                           region_chrom, region_start, region_end, times);
        //        delete firstline;
        bedlineQueue.pop();
    }

    if (summarize) {
        results.finish();
        checkResults(&results);
    } else {
        processNextBedLine(NULL, &results, statname, region_chrom,
                           region_start, region_end, times);
    }

    it = trees.begin();
    while (it != trees.end()) {
        delete it->second;
        advance(it, 1);
    }
    if (region_chrom != NULL) delete[] region_chrom;
    return 0;
}

int summarizeRegion(Config *config, const char *region,
                    set<string> inds, vector<string>statname,
                    vector<double> times) {
    if (config->snpfile.empty())
        return summarizeRegionNoSnp(config, region, inds, statname, times);
    else
        return summarizeRegionBySnp(config, region,
                                    inds, statname, times);
}


int main(int argc, char *argv[]) {
    Config c;
    int ret = c.parse_args(argc, argv);
    if (ret)
        return ret;

    if (c.argfile.empty()) {
        fprintf(stderr, "Error: must specify argfile\n");
        return 1;
    }

    vector<double> times;
    if (!c.timefile.empty()) {
        FILE *infile = fopen(c.timefile.c_str(), "r");
        double t;
        if (infile == NULL) {
            fprintf(stderr, "Error opening %s.\n", c.timefile.c_str());
            return 1;
        }
        while (EOF != fscanf(infile, "%lf", &t))
            times.push_back(t);
        std::sort(times.begin(), times.end());
        fclose(infile);
    }

    vector<string> statname;
    if (c.tmrca)
        statname.push_back(string("tmrca"));
    if (c.branchlen)
        statname.push_back(string("branchlen"));
    if (c.recomb)
        statname.push_back(string("recomb"));
    if (c.breaks)
        statname.push_back(string("breaks"));
    if (c.tmrca_half)
        statname.push_back(string("tmrca_half"));
    if (c.rth)
        statname.push_back(string("rth"));
    if (c.popsize)
        statname.push_back(string("popsize"));
    if (c.allele_age) {
        statname.push_back(string("allele_age"));
        statname.push_back(string("inf_sites"));
    }
    if (c.zero)
        statname.push_back(string("zero_len"));
    if (c.coalcounts) {
        if (c.timefile.empty()) {
            fprintf(stderr, "Error: --times required with --coalcounts\n");
            return 1;
        }
        for (unsigned int i=0; i < times.size(); i++) {
            char tmp[1000];
            sprintf(tmp, "coalcount.%i", i);
            statname.push_back(string(tmp));
        }
    }
    if (!c.node_dist.empty()) {
        vector<string> tokens1, tokens2;
        split(c.node_dist.c_str(), ';', tokens1);
        for (unsigned int i=0; i < tokens1.size(); i++) {
            split(tokens1[i].c_str(), ',', tokens2);
            if (tokens2.size() != 2) {
                fprintf(stderr, "Bad format to --node-dist argument; expect"
                        " two leaf names separated by comma");
                return 1;
            }
            statname.push_back(string("node_dist-") + tokens2[0]
                               + string(",") + tokens2[1]);
            node_dist_leaf1.push_back(tokens2[0]);
            node_dist_leaf2.push_back(tokens2[1]);
        }
    }
    if (c.rawtrees)
        statname.push_back(string("tree"));

    if (c.numsample)
        getNumSample=++summarize;
    if (c.mean)
        getMean=++summarize;
    if (c.stdev)
        getStdev=++summarize;;
    if (!c.quantile.empty()) {
        getQuantiles=++summarize;
        vector<string> tokens;
        split(c.quantile.c_str(), ',', tokens);
        for (unsigned int i=0; i < tokens.size(); i++) {
            double q=atof(tokens[i].c_str());
            //              fprintf(stderr, "getting quantile %lf\n",q);
            quantiles.push_back(q);
        }
    }

    if ((!c.region.empty()) && (!c.bedfile.empty())) {
        fprintf(stderr, "Error: --bed and --region cannot be used together.\n");
        return 1;
    }
    if (statname.size() == 0) {
        fprintf(stderr, "Error: need to specify a tree statistic\n");
        return 1;
    }

    if (summarize && statname.size()==0) {
        fprintf(stderr,
                "Error: need to specify a tree statistic (e.g., --tmrca,"
                " --popsize, --recomb, --allele-age, etc)\n");
        return 1;
    }
    if (summarize && c.rawtrees) {
        fprintf(stderr, "Error: --trees not compatible with summary statistics"
                " (--mean, --quantile, --stdev, --numsample)\n");
        return 1;
    }
    if ((c.recomb || c.breaks) && !c.snpfile.empty()) {
        fprintf(stderr, "Error: cannot use --recomb or --breaks with"
                " --allele-age\n");
        return 1;
    }
    if (c.allele_age && c.snpfile.empty()) {
        fprintf(stderr, "Error: need to specify snp file with --snp to use"
                " --allele-age\n");
        return 1;
    }

    if (!c.noheader) {
        printf("## %s\n", VERSION_INFO);
        printf("##");
        for (int i=0; i < argc; i++) printf(" %s", argv[i]);
        printf("\n");
        printf("#chrom\tchromStart\tchromEnd");
        if (summarize==0)
            printf("\tMCMC_sample");
        if (!c.snpfile.empty()) {
            printf("\tderAllele\tancAllele\tderFreq\tancFreq");
        }
        if (c.snpfile.empty() && getNumSample > 0)
            printf("\tnumsample");
        if (summarize && !c.snpfile.empty()) {
            printf("\tnumsample-all\tnumsample-infsites");
        }
        vector<string> stattype;
        if (c.snpfile.empty()) {
            stattype.push_back("");
        } else {
            stattype.push_back("-all");
            //    stattype.push_back("-derConsensus");
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
                            printf("\t%s%s_mean",
                                   statname[j].c_str(), stattype[k].c_str());
                        } else if (getStdev==i) {
                            printf("\t%s%s_stdev",
                                   statname[j].c_str(), stattype[k].c_str());
                        } else if (getQuantiles==i) {
                            for (unsigned int l=0; l < quantiles.size(); l++) {
                                printf("\t%s%s_quantile_%.3f",
                                       statname[j].c_str(),
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


    set <string>inds;
    if (!c.indfile.empty()) {
        ifstream in(c.indfile.c_str());
        string line;
        if (in.is_open()) {
            while ( getline(in, line) ) {
                inds.insert(line);
            }
            in.close();
        } else {
            fprintf(stderr, "Error opening %s.\n", c.indfile.c_str());
            return 1;
        }
    }

    if (c.bedfile.empty()) {
        summarizeRegion(&c, c.region.empty() ? NULL : c.region.c_str(),
                        inds, statname, times);
    } else {
        CompressStream bedstream(c.bedfile.c_str());
        char *line;
        vector<string> token;
        char *regionStr;
        if (!bedstream.stream) {
            fprintf(stderr, "error reading %s\n", c.bedfile.c_str());
            return 1;
        }
        while ((line = fgetline(bedstream.stream))) {
            split(line, '\t', token);
            if (token.size() < 3) {
                fprintf(stderr, "expected at least 3 files in %s\n",
                        c.bedfile.c_str());
                return 1;
            }
            regionStr = new char[token[0].size()+token[1].size()+
                                 token[2].size()+3];
            int start = atoi(token[1].c_str());
            int end = atoi(token[2].c_str());
            sprintf(regionStr, "%s:%i-%i", token[0].c_str(), start+1, end);
            summarizeRegion(&c, regionStr, inds, statname, times);
            delete [] regionStr;
        }
        bedstream.close();
    }

    return 0;
}
