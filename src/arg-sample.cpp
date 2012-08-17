

// C/C++ includes
#include <time.h>
#include <memory>

// arghmm includes
#include "compress.h"
#include "ConfigParam.h"
#include "emit.h"
#include "fs.h"
#include "logging.h"
#include "mem.h"
#include "parsing.h"
#include "sample_arg.h"
#include "sequences.h"
#include "total_prob.h"
#include <sys/stat.h>


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
const char *LOG_SUFFIX = ".log";


// parsing command-line options
class Config
{
public:

    Config()
    {
        make_parser();

        resample_region[0] = -1;
        resample_region[1] = -1;
    }

    void make_parser()
    {
        config.clear();

        // input/output
	config.add(new ConfigParam<string>
		   ("-s", "--sites", "<sites alignment>", &sitesfile, 
		    "sequence alignment in sites format"));
	config.add(new ConfigParam<string>
		   ("-f", "--fasta", "<fasta alignment>", &fastafile, 
		    "sequence alignment in FASTA format"));
	config.add(new ConfigParam<string>
		   ("-o", "--output", "<output prefix>", &out_prefix, 
                    "arg-sample",
                    "prefix for all output filenames (default='arg-sample')"));
        config.add(new ConfigParam<string>
                   ("-a", "--arg", "<SMC file>", &argfile, "",
                    "initial ARG file (*.smc) for resampling (optional)"));
        config.add(new ConfigParam<string>
                   ("", "--region", "<start>-<end>", 
                    &subregion_str, "",
                    "sample ARG for only a region of the sites (optional)"));

        // model parameters
	config.add(new ConfigParamComment("Model parameters"));
	config.add(new ConfigParam<double>
		   ("-N", "--popsize", "<population size>", &popsize, 1e4,
                    "effective population size (default=1e4)"));
	config.add(new ConfigParam<double>
		   ("-m", "--mutrate", "<mutation rate>", &mu, 2.5e-8,
                    "mutations per site per generation (default=2.5e-8)"));
	config.add(new ConfigParam<double>
		   ("-r", "--recombrate", "<recombination rate>", &rho, 1.5e-8,
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
                   ("", "--resample-region", "<start>-<end>", 
                    &resample_region_str, "",
                    "region to resample of input ARG (optional)"));
        config.add(new ConfigSwitch
		   ("", "--resume", &resume, "resume a previous run"));
        
        // misc
	config.add(new ConfigParamComment("Miscellaneous"));
 	config.add(new ConfigParam<int>
		   ("-c", "--compress-seq", "<compression factor>", 
                    &compress_seq, 1,
                    "alignment compression factor (default=1)"));
        config.add(new ConfigParam<int>
		   ("", "--sample-step", "<sample step size>", &sample_step, 
                    10, "number of iterations between steps (default=10)"));
 	config.add(new ConfigSwitch
		   ("", "--no-compress-output", &no_compress_output, 
                    "do not use compressed output"));
        config.add(new ConfigParam<int>
                   ("-x", "--randseed", "<random seed>", &randseed, 0,
                    "seed for random number generator (default: current time)"));

        // help information
	config.add(new ConfigParamComment("Information"));
	config.add(new ConfigParam<int>
		   ("-V", "--verbose", "<verbosity level>", 
		    &verbose, LOG_LOW, 
		    "verbosity level 0=quiet, 1=low, 2=medium, 3=high"));
	config.add(new ConfigSwitch
		   ("-q", "--quiet", &quiet, "suppress logging to stderr"));
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
    string subregion_str;

    // parameters
    double popsize;
    double mu;
    double rho;
    int ntimes;
    double maxtime;

    // search
    int nclimb;
    int niters;
    string resample_region_str;
    int resample_region[2];
    bool resume;
    string resume_stage;
    int resume_iter;

    // misc
    int compress_seq;
    int sample_step;
    bool no_compress_output;
    int randseed;
    
    // help/information
    bool quiet;
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

void log_intro(int level)
{
    time_t t = time(NULL);
    
    printLog(level, "argsample " VERSION_TEXT "\n");
    printLog(level, "start time: %s\n", ctime(&t));
}

void log_prog_commands(int level, int argc, char **argv)
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
    fprintf(stats_file, "stage\titer\tprior\tlikelihood\tjoint\trecombs\tnoncompats\n");
}


void print_stats(FILE *stats_file, const char *stage, int iter,
                 const ArgModel *model, 
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

    // get memory usage in MB
    double maxrss = get_max_memory_usage() / 1000.0;
    

    fprintf(stats_file, "%s\t%d\t%f\t%f\t%f\t%d\t%d\n",
            stage, iter,
            prior, likelihood, joint, nrecombs, noncompats);
    fflush(stats_file);

    printLog(LOG_LOW, "\n"
             "prior:      %f\n"
             "likelihood: %f\n"
             "joint:      %f\n"
             "nrecombs:   %d\n"
             "noncompats: %d\n"
             "max memory: %.1f MB\n\n",
             prior, likelihood, joint, nrecombs, noncompats,
             maxrss);
}

//=============================================================================
// sample output

string get_out_arg_file(const Config &config, int iter) 
{
    char iterstr[10];
    snprintf(iterstr, 10, ".%d", iter);
    return config.out_prefix + iterstr + SMC_SUFFIX;
}


bool log_local_trees(
    const ArgModel *model, const Sequences *sequences, LocalTrees *trees,
    const SitesMapping* sites_mapping, const Config *config, int iter)
{    
    //char iterstr[10];
    //snprintf(iterstr, 10, ".%d", iter);
    //string out_argfile = config->out_prefix + iterstr + SMC_SUFFIX;
    string out_argfile = get_out_arg_file(*config, iter);
    
    // write local trees uncompressed
    if (sites_mapping)
        uncompress_local_trees(trees, sites_mapping);

    // setup output stream
    FILE *out = NULL;
    if (config->no_compress_output)
        out = fopen(out_argfile.c_str(), "w");
    else
        out = open_compress((out_argfile + ".gz").c_str(), "w");
    if (!out) {
        printError("could not open '%s' for output", out_argfile.c_str());
        return false;
    }

    write_local_trees(out, trees, sequences, model->times);
    
    if (sites_mapping)
        compress_local_trees(trees, sites_mapping);

    if (config->no_compress_output)
        fclose(out);
    else
        close_compress(out);

    return true;
}


//=============================================================================

bool read_init_arg(const char *argfile, const ArgModel *model, 
                   LocalTrees *trees, vector<string> &seqnames)
{
    int len = strlen(argfile);
    bool compress = false;
    FILE *infile;
    if (len > 3 && strcmp(&argfile[len - 3], ".gz") == 0) {
        compress = true;
        infile = read_compress(argfile);
    } else {
        infile = fopen(argfile, "r");
    }
    if (!infile)
        return false;


    bool result = read_local_trees(infile, model->times, model->ntimes,
                                   trees, seqnames);

    if (compress)
        close_compress(infile);
    else
        fclose(infile);

    return result;
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
        print_stats(config->stats_file, "seq", trees->get_num_leaves(),
                    model, sequences, trees);
    }
}


void climb_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
               SitesMapping* sites_mapping, Config *config)
{
    if (config->resume)
        return;

    printLog(LOG_LOW, "Climb Search (%d iterations)\n", config->nclimb);
    printLog(LOG_LOW, "-----------------------------\n");
    double recomb_preference = .9;
    for (int i=0; i<config->nclimb; i++) {
        printLog(LOG_LOW, "climb %d\n", i+1);
        resample_arg_climb(model, sequences, trees, recomb_preference);
        print_stats(config->stats_file, "climb", i, model, sequences, trees);
    }
    printLog(LOG_LOW, "\n");
}


void resample_arg_all(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                      SitesMapping* sites_mapping, Config *config)
{
    int iter = 0;
    if (config->resume)
        iter = config->resume_iter;


    printLog(LOG_LOW, "Resample All Branches (%d iterations)\n", 
             config->niters);
    printLog(LOG_LOW, "--------------------------------------\n");
    for (int i=iter; i<config->niters; i++) {
        printLog(LOG_LOW, "sample %d\n", i+1);
        resample_arg_all(model, sequences, trees);

        // logging
        print_stats(config->stats_file, "resample", i, model, sequences, trees);

        // sample saving
        if (i % config->sample_step == 0)
            log_local_trees(model, sequences, trees, sites_mapping, config, i);
    }
    printLog(LOG_LOW, "\n");
}


// overall sampling workflow
void sample_arg(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                SitesMapping* sites_mapping, Config *config)
{
    if (!config->resume)
        print_stats_header(config->stats_file);

    // build initial arg by sequential sampling
    seq_sample_arg(model, sequences, trees, sites_mapping, config);
    
    if (config->resample_region[0] != -1) {
        // region sampling
        printLog(LOG_LOW, "Resample Region (%d-%d, %d iterations)\n",
                 config->resample_region[0], config->resample_region[1], 
                 config->niters);
        printLog(LOG_LOW, "--------------------------------------------\n");

        print_stats(config->stats_file, "resample_region", 0, 
                    model, sequences, trees);

        resample_arg_all_region(model, sequences, trees, 
                                config->resample_region[0], 
                                config->resample_region[1], 
                                config->niters);

        // logging
        print_stats(config->stats_file, "resample_region", config->niters,
                    model, sequences, trees);
        log_local_trees(model, sequences, trees, sites_mapping, config, 0);
        
    } else{
        // climb sampling
        climb_arg(model, sequences, trees, sites_mapping, config);
        // resample all branches
        resample_arg_all(model, sequences, trees, sites_mapping, config);
    }
}


//=============================================================================

bool parse_status_line(const char* line, const Config &config,
                       string &stage, int &iter, string &arg_file)
{    
    // parse stage and last iter
    vector<string> tokens;
    split(line, "\t", tokens);
    if (tokens.size() < 2) {
        printError("incomplete line in status file");
        return false;
    }
        
    string stage2 = tokens[0];
    int iter2;
    if (sscanf(tokens[1].c_str(), "%d", &iter2) != 1) {
        printError("iter column is not an integer");
        return false;
    }

    // NOTE: only resume resample stage for now
    if (stage2 != "resample")
        return true;

    // see if ARG file exists
    string out_argfile = get_out_arg_file(config, iter2);
    printf(">> %s\n", out_argfile.c_str());
    struct stat st;
    if (stat(out_argfile.c_str(), &st) == 0) {
        stage = stage2;
        iter = iter2;
        arg_file = out_argfile;
    }

    // try compress output
    out_argfile += ".gz";
    if (stat(out_argfile.c_str(), &st) == 0) {
        stage = stage2;
        iter = iter2;
        arg_file = out_argfile;
    }


    return true;
}


bool setup_resume(Config &config)
{
    if (!config.resume)
        return true;

    printLog(LOG_LOW, "Resuming previous run\n");

    // open stats file    
    string stats_filename = config.out_prefix + STATS_SUFFIX;
    printLog(LOG_LOW, "Checking previous run from stats file: %s\n", 
             stats_filename.c_str());

    FILE *stats_file;
    if (!(stats_file = fopen(stats_filename.c_str(), "r"))) {
        printError("could not open stats file '%s'",
                   stats_filename.c_str());
        return false;
    }

    // find last line of stats file that has a written ARG
    char *line = NULL;

    // skip header line
    line = fgetline(stats_file);
    if (!line) {
        printError("status file is empty");
        return false;
    }
    delete [] line;
    
    // loop through status lines
    string arg_file = "";
    while ((line = fgetline(stats_file))) {
        if (!parse_status_line(line, config, 
                               config.resume_stage, config.resume_iter, 
                               arg_file)) 
        {
            delete [] line;
            return false;
        }
        delete [] line;
    }

    if (arg_file == "") {
        printLog(LOG_LOW, "Could not find any previously writen ARG files. Try disabling resume");
        return false;
    }
    config.argfile = arg_file;
    
    printLog(LOG_LOW, "resuming at stage=%s, iter=%d, arg=%s\n", 
             config.resume_stage.c_str(), config.resume_iter,
             config.argfile.c_str());

    // clean up
    fclose(stats_file);

    return true;
}


//=============================================================================

int main(int argc, char **argv)
{
    Config c;
    int ret = c.parse_args(argc, argv);
    if (ret)
	return ret;

    // ensure output dir
    char *path = strdup(c.out_prefix.c_str());
    char *dir = dirname(path);
    if (!makedirs(dir)) {
        printError("could not make directory for output files '%s'",
                   dir);
        return 1;
    }
    free(path);
    
    // setup logging
    setLogLevel(c.verbose);
    string log_filename = c.out_prefix + LOG_SUFFIX;
    Logger *logger;
    if (c.quiet)
        logger = &g_logger;
    else
        logger = new Logger(NULL, c.verbose);

    const char *log_mode = (c.resume ? "a" : "w");
    if (!logger->openLogFile(log_filename.c_str(), log_mode)) {
        printError("could not open log file '%s'", log_filename.c_str());
        return 1;
    }
    if (!c.quiet)
        g_logger.setChain(logger);
    
    
    // log intro
    log_intro(LOG_LOW);
    log_prog_commands(LOG_LOW, argc, argv);
    Timer timer;


    // init random number generator
    if (c.randseed == 0)
        c.randseed = time(NULL);
    srand(c.randseed);
    printLog(LOG_LOW, "random seed: %d\n\n", c.randseed);


    // setup resuming
    if (!setup_resume(c)) {
        printError("resume failed.");
        return 1;
    }


    // setup model
    c.rho *= c.compress_seq;
    c.mu *= c.compress_seq;
    ArgModel model(c.ntimes, c.maxtime, c.popsize, c.rho, c.mu);


    // read sequences
    Sequences sequences;
    Sites *sites = NULL;
    auto_ptr<Sites> sites_ptr;
    SitesMapping *sites_mapping = NULL;
    auto_ptr<SitesMapping> sites_mapping_ptr;
    
    if (c.fastafile != "") {
        // read FASTA file
        
        if (!read_fasta(c.fastafile.c_str(), &sequences)) {
            printError("could not read fasta file");
            return 1;
        }

        printLog(LOG_LOW, 
                 "read input sequences (nseqs=%d, length=%d)\n",
                 sequences.get_num_seqs(), sequences.length());
    }
    else if (c.sitesfile != "") {
        // read sites file
        
        // parse subregion if given
        int subregion[2] = {-1, -1};
        if (c.subregion_str != "") {
            if (!parse_region(c.subregion_str.c_str(), 
                              &subregion[0], &subregion[1])) {
                printError("subregion is not specified as 'start-end'");
                return 1;
            }
        }

        // read sites
        sites = new Sites();
        sites_ptr = auto_ptr<Sites>(sites);
        if (!read_sites(c.sitesfile.c_str(), sites, 
                        subregion[0], subregion[1])) {
            printError("could not read sites file");
            return 1;
        }
        printLog(LOG_LOW, 
                 "read input sites (chrom=%s, start=%d, end=%d, length=%d, nseqs=%d, nsites=%d)\n",
                 sites->chrom.c_str(), sites->start_coord, sites->end_coord,
                 sites->length(), sites->get_num_seqs(),
                 sites->get_num_sites());

        // sanity check for sites
        if (sites->get_num_sites() == 0) {
            printLog(LOG_LOW, "no sites given.  terminating.\n");
            return 1;
        }

        if (c.compress_seq > 1) {
            // sequence compress requested
            sites_mapping = new SitesMapping();
            sites_mapping_ptr = auto_ptr<SitesMapping>(sites_mapping);
            find_compress_cols(sites, c.compress_seq, sites_mapping);
            compress_sites(sites, sites_mapping);
        }
            
        make_sequences_from_sites(sites, &sequences);
            
    } else {
        // no input sequence specified
        printError("must specify sequences (use --fasta or --sites)");
        return 1;
    }



    // get coordinates
    int start = 0;
    int end = sequences.length();
    if (sites) {
        start = sites->start_coord;
        end = sites->end_coord;
    }
    
    // setup init ARG
    LocalTrees *trees = NULL;
    auto_ptr<LocalTrees> trees_ptr;
    if (c.argfile != "") {
        // init ARG from file
        
        trees = new LocalTrees();
        trees_ptr = auto_ptr<LocalTrees>(trees);
        vector<string> seqnames;
        if (!read_init_arg(c.argfile.c_str(), &model, trees, seqnames)) {
            printError("could not read ARG");
            return 1;
        }

        if (!trees->set_seqids(seqnames, sequences.names)) {
            printError("input ARG's sequence names do not match input sequences");
            return 1;
        }
        
        printLog(LOG_LOW, "read input ARG (chrom=%s, start=%d, end=%d, nseqs=%d)\n",
                 trees->chrom.c_str(), trees->start_coord, trees->end_coord, 
                 trees->get_num_leaves());

        // compress input tree if compression is requested
        if (sites_mapping)
            compress_local_trees(trees, sites_mapping, true);

        // TODO: may need to adjust start and end
        // check ARG matches sites/sequences
        if (trees->start_coord != start || trees->end_coord != end) {
            printError("trees range does not match sites: tree(start=%d, end=%d), sites(start=%d, end=%d) [compressed coordinates]", 
                       trees->start_coord, trees->end_coord, start, end);
            return 1;
        }
    } else {
        // create new init ARG
        trees = new LocalTrees(start, end);
        trees_ptr = auto_ptr<LocalTrees>(trees);
    }
    
    // set chromosome name
    if (sites)
        trees->chrom = sites->chrom;
    
    
    // setup coordinates for sequences
    Sequences sequences2(&sequences, -1, start + sequences.length(), -start);


    // check for region sample
    if (c.resample_region_str != "") {
        if (!parse_region(c.resample_region_str.c_str(), 
                          &c.resample_region[0], &c.resample_region[1])) {
            printError("region is not specified as 'start-end'");
            return 1;
        }

        if (sites_mapping) {
            c.resample_region[0] = sites_mapping->compress(c.resample_region[0]);
            c.resample_region[1] = sites_mapping->compress(c.resample_region[1]);
        }
        
    }
    
    
    // init stats file
    string stats_filename = c.out_prefix + STATS_SUFFIX;
    const char *stats_mode = (c.resume ? "a" : "w");
    if (!(c.stats_file = fopen(stats_filename.c_str(), stats_mode))) {
        printError("could not open stats file '%s'", stats_filename.c_str());
        return 1;
    }

    // sample ARG
    printLog(LOG_LOW, "\n");
    sample_arg(&model, &sequences2, trees, sites_mapping, &c);
    
    // final log message
    double maxrss = get_max_memory_usage() / 1000.0;
    printTimerLog(timer, LOG_LOW, "sampling time: ");
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);
    printLog(LOG_LOW, "FINISH\n");

    // clean up
    fclose(c.stats_file);
    
    
    return 0;
}
