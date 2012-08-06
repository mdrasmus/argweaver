
#include "ConfigParam.h"
#include "emit.h"
#include "logging.h"
#include "sample_arg.h"
#include "sequences.h"
#include "total_prob.h"


using namespace arghmm;


// version info
#define VERSION_TEXT "1.0"
#define VERSION_INFO  "\
ArgHmm " VERSION_TEXT " \n\
Matt Rasmussen\n\
Gibbs sampler for ancestral recombination graphs\n\
"

// file extensions
const char *SMC_SUFFIX = ".smc";
const char *STATS_SUFFIX = ".stats";



// parsing command-line options
class Config
{
public:

    Config() 
    {
        make_parser();
    }

    void make_parser()
    {
        config.clear();

        // input/output
	config.add(new ConfigParam<string>
		   ("-f", "--fasta", "<fasta alignment>", &fastafile, 
		    "sequence alignment in fasta format"));
	config.add(new ConfigParam<string>
		   ("-s", "--sites", "<sites alignment>", &sitesfile, 
		    "sequence alignment in sites format"));
	config.add(new ConfigParam<string>
		   ("-o", "--out", "<output prefix>", &out_prefix, 
                    "arghmm",
                    "prefix for all output filenames (default='arghmm')"));
	config.add(new ConfigParam<int>
		   ("-c", "--compress", "<compression>", &compress, 1,
                    "alignment compression factor (default=1)"));
        config.add(new ConfigParam<string>
                   ("-a", "--arg", "<SMC file>", &argfile, "",
                    "initial ARG file (*.smc) for resampling (optional)"));

        // model parameters
	config.add(new ConfigParamComment("Model parameters"));
	config.add(new ConfigParam<double>
		   ("-N", "--popsize", "<population size>", &popsize, 1e4,
                    "effective population size (default=1e4)"));
	config.add(new ConfigParam<double>
		   ("-m", "--mu", "<mutation rate>", &mu, 2.5e-8,
                    "mutations per site per generation (default=2.5e-8)"));
	config.add(new ConfigParam<double>
		   ("-r", "--rho", "<recombination rate>", &rho, 1.5e-8,
                    "recombination per site per generation (default=1.5e-8)"));
	config.add(new ConfigParam<int>
		   ("-t", "--ntimes", "<ntimes>", &ntimes, 20,
                    "number of time points (default=20)"));
	config.add(new ConfigParam<double>
		   ("", "--maxtime", "<maxtime>", &maxtime, 200e3,
                    "maximum time point in generations (default=200e3)"));

        // search
	config.add(new ConfigParamComment("Search"));
	config.add(new ConfigParam<int>
		   ("", "--climb", "<# of climb iterations>", &nclimb, 50,
                    "(default=50)"));
	config.add(new ConfigParam<int>
		   ("-n", "--iters", "<# of iterations>", &niters, 1000,
                    "(default=1000)"));
        config.add(new ConfigParam<string>
                   ("", "--region", "<start>-<end>", &region_str, "",
                    "region to resample (optional)"));
        config.add(new ConfigParam<int>
                   ("-x", "--randseed", "<random seed>", &randseed, 0,
                    "seed for random number generator (default: time)"));


        // help information
	config.add(new ConfigParamComment("Information"));
	config.add(new ConfigParam<int>
		   ("-V", "--verbose", "<verbosity level>", 
		    &verbose, LOG_LOW, 
		    "verbosity level 0=quiet, 1=low, 2=medium, 3=high"));
	//config.add(new ConfigParam<string>
	//	   ("", "--log", "<log filename>", &logfile, "", 
	//	    "log filename.  Use '-' to display on stdout."));
	config.add(new ConfigSwitch
		   ("-v", "--version", &version, "display version information"));
	config.add(new ConfigSwitch
		   ("-h", "--help", &help, 
		    "display help information"));
    }

    int parse_args(int argc, char **argv)
    {
	// parse arguments
	if (!config.parse(argc, (const char**) argv)) {
	    if (argc < 2)
		config.printHelp();
	    return 1;
	}
    
	// display help
	if (help) {
	    config.printHelp();
	    return 1;
	}
    
	// display version info
	if (version) {
	    printf(VERSION_INFO);
	    return 1;
	}
        
	return 0;
    }

    ConfigParser config;

    // input/output
    string fastafile;
    string sitesfile;
    string out_prefix;
    string argfile;

    // parameters
    double popsize;
    double mu;
    double rho;
    int ntimes;
    double maxtime;
    int compress;

    // search
    int nclimb;
    int niters;
    string region_str;
    int region[2];
    int randseed;
    
    // help/information
    int verbose;
    bool version;
    bool help;

    // logging
    FILE *stats_file;

};



bool parse_region(const char *region, int *start, int *end)
{
    return sscanf(region, "%d-%d", start, end) == 2;
}

//=============================================================================
// logging

void logProgCommands(int level, int argc, char **argv)
{
    printLog(level, "command:");
    for (int i=0; i<argc; i++) {
        printLog(level, " %s", argv[i]);
    }
    printLog(level, "\n");
}


//=============================================================================
// statistics output

void print_stats_header(FILE *stats_file)
{
    fprintf(stats_file, "prior\tlikelihood\tjoint\trecombs\tnoncompats\n");
}


void print_stats(FILE *stats_file, const ArgModel *model, 
                 const Sequences *sequences, const LocalTrees *trees)
{
    double prior = calc_arg_prior(model, trees);
    double likelihood = calc_arg_likelihood(model, sequences, trees);
    double joint = prior + likelihood;
    int nrecombs = trees->get_num_trees() - 1;

    int nseqs = sequences->get_num_seqs();
    char *seqs[nseqs];
    for (int i=0; i<nseqs; i++)
        seqs[i] = sequences->seqs[trees->seqids[i]];
    
    int noncompats = count_noncompat(trees, seqs, nseqs, sequences->length());


    fprintf(stats_file, "%f\t%f\t%f\t%d\t%d\n",
            prior, likelihood, joint, nrecombs, noncompats);

    printLog(LOG_LOW, "\n"
             "prior:      %f\n"
             "likelihood: %f\n"
             "joint:      %f\n"
             "nrecombs:   %d\n"
             "noncompats: %d\n\n",
             prior, likelihood, joint, nrecombs, noncompats);
}

//=============================================================================
// sample output

void log_local_trees(
    const ArgModel *model, const Sequences *sequences, LocalTrees *trees,
    const SitesMapping* sites_mapping, const Config *config, int iter)
{
    char iterstr[10];
    snprintf(iterstr, 10, ".%d", iter);
    string out_argfile = config->out_prefix + iterstr + SMC_SUFFIX;
    
    // write local trees uncompressed
    if (sites_mapping)
        uncompress_local_trees(trees, sites_mapping);
    write_local_trees(out_argfile.c_str(), trees, sequences, model->times);    
    if (sites_mapping)
        compress_local_trees(trees, sites_mapping);
}


//=============================================================================
// sampling methods

// build initial arg by sequential sampling
void seq_sample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                    SitesMapping* sites_mapping, Config *config)
{
    if (trees->get_num_leaves() < sequences->get_num_seqs()) {
        printLog(LOG_LOW, "Sequentially Sample Initial ARG (%d sequences)\n",
                 sequences->get_num_seqs());
        printLog(LOG_LOW, "------------------------------------------------\n");
        sample_arg_seq(model, sequences, trees);
        print_stats(config->stats_file, model, sequences, trees);
    }
}


void climb_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
               SitesMapping* sites_mapping, Config *config)
{
    printLog(LOG_LOW, "Climb Search (%d iterations)\n", config->nclimb);
    printLog(LOG_LOW, "-----------------------------\n");
    double recomb_preference = .9;
    for (int i=0; i<config->nclimb; i++) {
        printLog(LOG_LOW, "climb %d\n", i+1);
        resample_arg_climb(model, sequences, trees, recomb_preference);
        print_stats(config->stats_file, model, sequences, trees);
    }
    printLog(LOG_LOW, "\n");
}


void resample_arg_all(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                      SitesMapping* sites_mapping, Config *config)
{
    printLog(LOG_LOW, "Resample All Branches (%d iterations)\n", 
             config->niters);
    printLog(LOG_LOW, "--------------------------------------\n");
    for (int i=0; i<config->niters; i++) {
        printLog(LOG_LOW, "sample %d\n", i+1);
        resample_arg_all(model, sequences, trees);

        // logging
        print_stats(config->stats_file, model, sequences, trees);
        log_local_trees(model, sequences, trees, sites_mapping, config, i);
    }
    printLog(LOG_LOW, "\n");
}


void sample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                SitesMapping* sites_mapping, Config *config)
{
    print_stats_header(config->stats_file);

    // build initial arg by sequential sampling
    seq_sample_arg(model, sequences, trees, sites_mapping, config);
    
    if (config->region[0] != -1) {
        // region sampling
        printLog(LOG_LOW, "Resample Region (%d-%d, %d iterations)\n",
                 config->region[0], config->region[1], config->niters);
        printLog(LOG_LOW, "--------------------------------------------\n");
        resample_arg_all_region(model, sequences, trees, 
                                config->region[0], config->region[1], 
                                config->niters);

        // logging
        print_stats(config->stats_file, model, sequences, trees);
        log_local_trees(model, sequences, trees, sites_mapping, config, 0);
        
    } else{
        // climb sampling
        climb_arg(model, sequences, trees, sites_mapping, config);    
        // resample all branches
        resample_arg_all(model, sequences, trees, sites_mapping, config);
    }
    
    
    // final stats
    int nrecombs = trees->get_num_trees() - 1;
    printLog(LOG_LOW, "nrecombs %d\n", nrecombs);
}



//=============================================================================

int main(int argc, char **argv)
{
    Config c;
    int ret = c.parse_args(argc, argv);
    if (ret)
	return ret;

    // setup logging
    setLogLevel(c.verbose);
    /*
    string log_filename = c.out_prefix + LOG_SUFFIX;
    if (!(c.stats_file = fopen(stats_filename.c_str(), "w"))) {
        printError("could not open stats file '%s'", stats_filename.c_str());
        return 1;
    }    
    */
    logProgCommands(LOG_LOW, argc, argv);

    // init stats file
    string stats_filename = c.out_prefix + STATS_SUFFIX;
    if (!(c.stats_file = fopen(stats_filename.c_str(), "w"))) {
        printError("could not open stats file '%s'", stats_filename.c_str());
        return 1;
    }

    // init random number generator
    if (c.randseed == 0)
        c.randseed = time(NULL);
    srand(c.randseed);
    printLog(LOG_LOW, "random seed: %d\n\n", c.randseed);


    // setup model
    c.rho *= c.compress;
    c.mu *= c.compress;
    ArgModel model(c.ntimes, c.maxtime, c.popsize, c.rho, c.mu);


    // read sequences
    Sequences *sequences = NULL;
    Sites *sites = NULL;
    SitesMapping *sites_mapping = NULL;
    if (c.fastafile != "") {
        sequences = new Sequences();
        if (!read_fasta(c.fastafile.c_str(), sequences)) {
            delete sequences;
            sequences = NULL;
            printError("could not read fasta file");
            return 1;
        }
    }
    else if (c.sitesfile != "") {
        sites = new Sites();
        if (read_sites(c.sitesfile.c_str(), sites)) {
            if (c.compress > 1) {
                sites_mapping = new SitesMapping();
                find_compress_cols(sites, c.compress, sites_mapping);
                compress_sites(sites, sites_mapping);
            }
            
            sequences = new Sequences();
            make_sequences_from_sites(sites, sequences);
        } else {
            printError("could not read sites file");
            delete sites;
            sites = NULL;
            return 1;

        }
    } else {
        printError("must specify sequences (use --fasta or --sites)");
        return 1;
    }

    printLog(LOG_LOW, "read input sequences (nseqs=%d, length=%d)\n",
             sequences->get_num_seqs(), sequences->length());


    // get coordinates
    int start = 0;
    int end = sequences->length();
    if (sites) {
        start = sites->start_coord;
        end = sites->end_coord;
    }
    
    // setup init ARG
    LocalTrees *trees = NULL;
    if (c.argfile != "") {
        // init ARG from file
        
        trees = new LocalTrees();
        vector<string> seqnames;
        if (!read_local_trees(c.argfile.c_str(), model.times, model.ntimes,
                              trees, seqnames)) {
            printError("could not read ARG");
            return 1;
        }

        // may need to adjust start and end
        // check ARG matches sites/sequences
        if (trees->start_coord != start || trees->end_coord != end) {
            printError("trees range does not match sites: tree(start=%d, end=%d), sites(start=%d, end=%d)", trees->start_coord, trees->end_coord, start, end);
            return 1;
        }

        if (!trees->set_seqids(seqnames, sequences->names)) {
            printError("input ARG's sequence names do not match input sequences");
            return 1;
        }
        
        printLog(LOG_LOW, "read input ARG (nseqs=%d, start=%d, end=%d)\n",
                 trees->start_coord, trees->end_coord, trees->get_num_leaves());

    } else {
        // create new init ARG
        trees = new LocalTrees(start, end);
    }
    
    
    // setup coordinates for sequences
    Sequences sequences2(sequences, -1, start + sequences->length(), -start);


    // check for region sample
    if (c.region_str != "") {
        if (!parse_region(c.region_str.c_str(), &c.region[0], &c.region[1])) {
            printError("region is not specified as 'start-end'");
            return 1;
        }
    } else {
        c.region[0] = -1;
        c.region[1] = -1;
    }

    // sample ARG
    printLog(LOG_LOW, "\n");
    sample_arg(&model, &sequences2, trees, sites_mapping, &c);
    


    // clean up
    fclose(c.stats_file);
    
    if (sites)
        delete sites;
    if (sites_mapping)
        delete sites_mapping;
    delete trees;
    delete sequences;

    return 0;
}
