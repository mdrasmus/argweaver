
#include "ConfigParam.h"
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
const char *LOG_SUFFIX = ".log";



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
    
    // help/information
    int verbose;
    bool version;
    bool help;

    // logging
    FILE *logfile;

};


//=============================================================================

void print_log_header(FILE *logfile)
{
    fprintf(logfile, "prior\tlikelihood\tjoint\trecombs\n");
}


void print_log_iter(FILE *logfile, ArgModel *model, Sequences *sequences, 
                    LocalTrees *trees)
{
    double prior = calc_arg_prior(model, trees);
    double likelihood = calc_arg_likelihood(model, sequences, trees);
    double joint = prior + likelihood;
    int nrecombs = trees->get_num_trees() - 1;
    
    fprintf(logfile, "%f\t%f\t%f\t%d\n",
            prior, likelihood, joint, nrecombs);
}


void log_local_trees(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                     SitesMapping* sites_mapping, Config *config, int iter)
{
    char iterstr[10];
    snprintf(iterstr, 10, ".%d", iter);
    string out_argfile = config->out_prefix + iterstr + SMC_SUFFIX;

    //assert_uncompress_local_trees(trees, sites_mapping);

    uncompress_local_trees(trees, sites_mapping);
    write_local_trees(out_argfile.c_str(), trees, sequences, model->times);    
    compress_local_trees(trees, sites_mapping);
}


void sample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                SitesMapping* sites_mapping, Config *config)
{

    // build initial arg by sequential sampling
    print_log_header(config->logfile);
    printLog(LOG_LOW, "sequentially sample initial ARG\n");
    sample_arg_seq(model, sequences, trees);
    printLog(LOG_LOW, "\n");
    print_log_iter(config->logfile, model, sequences, trees);
    

    // climb sampling
    double recomb_preference = .9;
    for (int i=0; i<config->nclimb; i++) {
        printLog(LOG_LOW, "climb %d\n", i+1);
        resample_arg_climb(model, sequences, trees, recomb_preference);
        print_log_iter(config->logfile, model, sequences, trees);
    }
    printLog(LOG_LOW, "\n");
    

    // resample all branches
    for (int i=0; i<config->niters; i++) {
        printLog(LOG_LOW, "sample %d\n", i+1);
        resample_arg_all(model, sequences, trees);

        // logging
        print_log_iter(config->logfile, model, sequences, trees);
        log_local_trees(model, sequences, trees, sites_mapping, config, i);
    }
    printLog(LOG_LOW, "\n");
    
    
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
    string logfilename = c.out_prefix + LOG_SUFFIX;
    if (!(c.logfile = fopen(logfilename.c_str(), "w"))) {
        printError("could not open log file '%s'", logfilename.c_str());
        return 1;
    }


    // setup model
    c.rho *= c.compress;
    c.mu *= c.compress;
    ArgModel model(c.ntimes, c.maxtime, c.popsize, c.rho, c.mu);


    // read sequences
    if (c.fastafile == "" && c.sitesfile == "") {
        printError("alignment file is required (--fasta or --sites)");
        return 1;
    }
    
    Sequences *sequences = NULL;
    Sites *sites = NULL;
    SitesMapping *sites_mapping = NULL;
    if (c.fastafile != "") {
        sequences = read_fasta(c.fastafile.c_str());
    }
    else if (c.sitesfile != "") {
        sites = read_sites(c.sitesfile.c_str());
        if (sites) {
            if (c.compress > 1) {
                sites_mapping = new SitesMapping();
                find_compress_cols(sites, c.compress, sites_mapping);
                compress_sites(sites, sites_mapping);
            }
            
            sequences = make_sequences_from_sites(sites);
        }
    }

    if (!sequences) {
        printError("could not read alignment file");
        return 1;
    }
    
    // init ARG and sequences for correct range
    int start = 0;
    int end = sequences->length();
    if (sites) {
        start = sites->start_coord;
        end = sites->end_coord;
    }
    LocalTrees *trees = new LocalTrees(start, end);
    Sequences sequences2(sequences, -1, start + sequences->length(), -start);
    
    //for (unsigned int i=0; i<sequences2.names.size(); i++)
    //    printf(">> %s\n", sequences2.names[i].c_str());
    
    // sample ARG
    sample_arg(&model, &sequences2, trees, sites_mapping, &c);
    


    // clean up
    fclose(c.logfile);
    
    if (sites_mapping)
        delete sites_mapping;
    delete trees;
    delete sequences;

    return 0;
}
