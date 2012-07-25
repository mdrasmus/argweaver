
#include "logging.h"
#include "sample_arg.h"
#include "sequences.h"
#include "ConfigParam.h"

using namespace arghmm;


#define VERSION_TEXT "1.0"
#define VERSION_INFO  "\
ArgHmm " VERSION_TEXT " \n\
Matt Rasmussen\n\
Gibbs sampler for ancestral recombination graphs\n\
"


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
	config.add(new ConfigParam<double>
		   ("-c", "--compress", "<compression>", &compress, 1.0,
                    "alignment compression factor (default=1.0)"));
	config.add(new ConfigParam<string>
		   ("-o", "--out-arg", "<ARG output file>", &out_argfile, 
                    "arghmm.smc",
                    "output file for the sampled ARG (default='arghmm.sprs')"));


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
	config.add(new ConfigParam<string>
		   ("", "--log", "<log filename>", &logfile, "", 
		    "log filename.  Use '-' to display on stdout."));
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
    string out_argfile;

    // parameters
    double popsize;
    double mu;
    double rho;
    int ntimes;
    double maxtime;
    double compress;

    // search
    int nclimb;
    int niters;
    
    // help/information
    int verbose;
    bool version;
    bool help;
    string logfile;

};



void sample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                int nclimb, int niters)
{
    printLog(LOG_LOW, "sequentially sample initial ARG\n");
    sample_arg_seq(model, sequences, trees);
    printLog(LOG_LOW, "\n");
    
    // climb sampling
    double recomb_preference = .9;
    for (int i=0; i<nclimb; i++) {
        printLog(LOG_LOW, "climb %d\n", i+1);
        resample_arg_climb(model, sequences, trees, recomb_preference);
    }
    printLog(LOG_LOW, "\n");

    for (int i=0; i<niters; i++) {
        printLog(LOG_LOW, "sample %d\n", i+1);
        resample_arg_all(model, sequences, trees);
    }
    printLog(LOG_LOW, "\n");

    // final stats
    int nrecombs = trees->get_num_trees() - 1;
    printLog(LOG_LOW, "nrecombs %d\n", nrecombs);
}



int main(int argc, char **argv)
{
    Config c;
    int ret = c.parse_args(argc, argv);
    if (ret)
	return ret;

    // setup logging
    setLogLevel(c.verbose);


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
    if (c.fastafile != "") {
        sequences = read_fasta(c.fastafile.c_str());
    }
    else if (c.sitesfile != "") {
        sites = read_sites(c.sitesfile.c_str());
        if (sites)
            sequences = make_sequences_from_sites(sites);
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

    
    // sample ARG
    sample_arg(&model, &sequences2, trees, c.nclimb, c.niters);
    
    // output
    write_local_trees(c.out_argfile.c_str(), trees, sequences2, model.times);


    delete trees;
    delete sequences;

    return 0;
}
