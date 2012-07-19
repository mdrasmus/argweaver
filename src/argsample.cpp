
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
		   ("-a", "--align", "<alignment fasta>", &alignfile, 
		    "sequence alignment in fasta format"));
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
                    "arghmm.sprs",
                    "output file for the sampled ARG (default='arghmm.sprs')"));


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
    string alignfile;
    string out_argfile;

    // parameters
    double popsize;
    double mu;
    double rho;
    int ntimes;
    double maxtime;
    double compress;
    
    // help/information
    int verbose;
    bool version;
    bool help;
    string logfile;

};




int main(int argc, char **argv)
{
    Config c;
    int ret = c.parse_args(argc, argv);
    if (ret)
	return ret;


    // setup model
    c.rho *= c.compress;
    c.mu *= c.compress;
    ArgModel model(c.ntimes, c.maxtime, c.popsize, c.rho, c.mu);


    // read sequences
    if (c.alignfile == "") {
        printError("alignment file is required (--align)");
        return 1;
    }
    Sequences *sequences = read_fasta(c.alignfile.c_str());
    if (!sequences) {
        printError("could not read alignment file");
        return 1;
    }
    
    // sample arg
    LocalTrees *trees = new LocalTrees();
    sample_arg_seq(&model, sequences, trees);
    printf("length %d\n", sequences->length());
    printf("ntrees %d\n", trees->get_num_trees());


    //LocalTree *tree = trees->front().tree;
    //write_newick_tree(stdout, tree, NULL, model.times, 0, true);

    //write_local_trees(c.out_argfile.c_str(), trees, NULL, model.times);

    printf("%s\n", sequences->names[0].c_str());
    write_local_trees(c.out_argfile.c_str(), trees, *sequences, model.times);


    delete trees;
    delete sequences;
}
