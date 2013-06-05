

// C/C++ includes
#include <time.h>
#include <memory>
#include <sys/stat.h>

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
#include "track.h"



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


// debug options level
const int DEBUG_OPT = 1;


const int EXIT_ERROR = 1;



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
		   ("-s", "--sites", "<sites alignment>", &sites_file, 
		    "sequence alignment in sites format"));
	config.add(new ConfigParam<string>
		   ("-f", "--fasta", "<fasta alignment>", &fasta_file, 
		    "sequence alignment in FASTA format"));
	config.add(new ConfigParam<string>
		   ("-o", "--output", "<output prefix>", &out_prefix, 
                    "arg-sample",
                    "prefix for all output filenames (default='arg-sample')"));
        config.add(new ConfigParam<string>
                   ("-a", "--arg", "<SMC file>", &arg_file, "",
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
        config.add(new ConfigParam<double>
                   ("", "--time-step", "<time>", &time_step, 0,
                    "linear time step in generations (optional)"));
        config.add(new ConfigParam<string>
                   ("", "--times-file", "<times filename>", &times_file, "",
                    "file containing time points (optional)"));
	config.add(new ConfigParam<string>
		   ("-M", "--mutmap", "<mutation rate map file>", &mutmap, "",
                    "mutation map file (optional)"));
	config.add(new ConfigParam<string>
		   ("-R", "--recombmap", "<recombination rate map file>", 
                    &recombmap, "",
                    "recombination map file (optional)"));
	config.add(new ConfigParam<string>
		   ("", "--maskmap", "<sites mask>", 
                    &maskmap, "",
                    "mask map file (optional)"));


        // search
	config.add(new ConfigParamComment("Search"));
	config.add(new ConfigParam<int>
		   ("-n", "--iters", "<# of iterations>", &niters, 1000,
                    "(default=1000)"));
        config.add(new ConfigParam<string>
                   ("", "--resample-region", "<start>-<end>", 
                    &resample_region_str, "",
                    "region to resample of input ARG (optional)"));
        config.add(new ConfigSwitch
		   ("", "--resume", &resume, "resume a previous run"));
        config.add(new ConfigSwitch
		   ("", "--overwrite", &overwrite, 
                    "force an overwrite of a previous run"));
        config.add(new ConfigSwitch
		   ("", "--gibbs", &gibbs, 
                    "use Gibbs sampling"));
        
        // misc
	config.add(new ConfigParamComment("Miscellaneous"));
 	config.add(new ConfigParam<int>
		   ("-c", "--compress-seq", "<compression factor>", 
                    &compress_seq, 1,
                    "alignment compression factor (default=1)"));
	config.add(new ConfigParam<int>
		   ("", "--climb", "<# of climb iterations>", &nclimb, 0,
                    "(default=0)"));
        config.add(new ConfigParam<int>
		   ("", "--sample-step", "<sample step size>", &sample_step, 
                    10, "number of iterations between steps (default=10)"));
 	config.add(new ConfigSwitch
		   ("", "--no-compress-output", &no_compress_output, 
                    "do not use compressed output"));
        config.add(new ConfigParam<int>
                   ("-x", "--randseed", "<random seed>", &randseed, 0,
                    "seed for random number generator (default=current time)"));

        // advance options
        config.add(new ConfigParamComment("Advanced Options", DEBUG_OPT));
        config.add(new ConfigParam<double>
                   ("", "--prob-path-switch", "<probability>", 
                    &prob_path_switch, .1,
                    "removal path switch (default=.1)", DEBUG_OPT));
        config.add(new ConfigSwitch
                   ("", "--infsites", &infsites, 
                    "assume infinite sites model (at most one mutation per site)", 
                    DEBUG_OPT));
        config.add(new ConfigParam<int>
                   ("", "--resample-window", "<window size>", 
                    &resample_window, 100000,
                    "sliding window for resampling (default=100000)", DEBUG_OPT));
        config.add(new ConfigParam<int>
                   ("", "--resample-window-iters", "<iterations>", 
                    &resample_window_iters, 10,
                    "number of iterations per sliding window for resampling (default=10)", DEBUG_OPT));


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
        config.add(new ConfigSwitch
                   ("", "--help-advanced", &help_debug, 
                    "display help information about advanced options"));
    }

    int parse_args(int argc, char **argv)
    {
	// parse arguments
	if (!config.parse(argc, (const char**) argv)) {
	    if (argc < 2)
		config.printHelp();
	    return EXIT_ERROR;
	}
    
	// display help
	if (help) {
	    config.printHelp();
	    return EXIT_ERROR;
	}

        // display debug help
        if (help_debug) {
            config.printHelp(stderr, DEBUG_OPT);
            return EXIT_ERROR;
        }
    
	// display version info
	if (version) {
	    printf(VERSION_INFO);
	    return EXIT_ERROR;
	}
        
	return 0;
    }

    ConfigParser config;

    // input/output
    string fasta_file;
    string sites_file;
    string out_prefix;
    string arg_file;
    string subregion_str;

    // model parameters
    double popsize;
    double mu;
    double rho;
    int ntimes;
    double maxtime;
    double time_step;
    string times_file;
    string mutmap;
    string recombmap;
    string maskmap;
    ArgModel model;

    // search
    int nclimb;
    int niters;
    string resample_region_str;
    int resample_region[2];
    bool resume;
    bool overwrite;
    string resume_stage;
    int resume_iter;
    int resample_window;
    int resample_window_iters;
    bool gibbs;

    // misc
    int compress_seq;
    int sample_step;
    bool no_compress_output;
    int randseed;
    double prob_path_switch;
    bool infsites;
    
    // help/information
    bool quiet;
    int verbose;
    bool version;
    bool help;
    bool help_debug;

    // logging
    FILE *stats_file;
};



bool parse_region(const char *region, int *start, int *end)
{
    return sscanf(region, "%d-%d", start, end) == 2;
}

//=============================================================================
// logging

// log the program version and start time
void log_intro(int level)
{
    time_t t = time(NULL);
    
    printLog(level, "arg-sample " VERSION_TEXT "\n");
    printLog(level, "start time: %s", ctime(&t));  // newline is in ctime
}

// log the command used
void log_prog_commands(int level, int argc, char **argv)
{
    printLog(level, "command:");
    for (int i=0; i<argc; i++) {
        printLog(level, " %s", argv[i]);
    }
    printLog(level, "\n");
}


// log the model used
void log_model(const ArgModel &model)
{
    printLog(LOG_LOW, "\n");
    printLog(LOG_LOW, "model: \n");
    printLog(LOG_LOW, "  mu = %e\n", model.mu);
    printLog(LOG_LOW, "  rho = %e\n", model.rho);
    printLog(LOG_LOW, "  ntimes = %d\n", model.ntimes);
    printLog(LOG_LOW, "  times = [");
    for (int i=0; i<model.ntimes-1; i++)
        printLog(LOG_LOW, "%f,", model.times[i]);
    printLog(LOG_LOW, "%f]\n", model.times[model.ntimes-1]);
    printLog(LOG_LOW, "  popsizes = [");
    for (int i=0; i<model.ntimes-1; i++)
        printLog(LOG_LOW, "%f,", model.popsizes[i]);
    printLog(LOG_LOW, "%f]\n", model.popsizes[model.ntimes-1]);

    if (isLogLevel(LOG_HIGH)) {
        printLog(LOG_HIGH, "mutmap = [\n");
        for (unsigned int i=0; i<model.mutmap.size(); i++) {
            printLog(LOG_HIGH, "%d\t%d\t%e\n", 
                     model.mutmap[i].start, model.mutmap[i].end,
                     model.mutmap[i].value);
        }
        printLog(LOG_HIGH, "]\n");

        printLog(LOG_HIGH, "recombmap = [\n");
        for (unsigned int i=0; i<model.recombmap.size(); i++) {
            printLog(LOG_HIGH, "%d\t%d\t%e\n", 
                     model.recombmap[i].start, model.recombmap[i].end,
                     model.recombmap[i].value);
        }
        printLog(LOG_HIGH, "]\n");
    }

    printLog(LOG_LOW, "\n");
}


//=============================================================================
// alignment compression

template<class T>
void compress_track(Track<T> &track, SitesMapping *sites_mapping,
                    double compress_seq, bool is_rate)
{
    Track<T> track2;
    
    if (sites_mapping) {
        // get block lengths
        vector<int> blocks;
        for (unsigned int i=0; i<track.size(); i++)
            blocks.push_back(track[i].length());
        
        // compress block lengths
        vector<int> blocks2;
        sites_mapping->compress_blocks(blocks, blocks2);

        // build compress track
        int start = sites_mapping->new_start;
        for (unsigned int i=0; i<track.size(); i++) {
            int end = start + blocks2[i];
            if (end > start)
                track2.append(track[i].chrom, start, end, track[i].value);
            start = end;
        }

        // replace track
        track.clear();
        track.insert(track.begin(), track2.begin(), track2.end());
    }


    // compress rate
    if (is_rate) {
        for (unsigned int i=0; i<track.size(); i++)
            track[i].value *= compress_seq;
    }
}


template<class T>
void compress_mask(Track<T> &track, SitesMapping *sites_mapping)
{
    if (sites_mapping) {
        for (unsigned int i=0; i<track.size(); i++) {
            track[i].start = sites_mapping->compress(track[i].start);
            track[i].end = sites_mapping->compress(track[i].end);
        }
    }
}


void compress_model(ArgModel *model, SitesMapping *sites_mapping,
                    double compress_seq)
{
    model->rho *= compress_seq;
    model->mu *= compress_seq;

    compress_track(model->mutmap, sites_mapping, compress_seq, true);
    compress_track(model->recombmap, sites_mapping, compress_seq, true);
}




//=============================================================================
// statistics output

void print_stats_header(FILE *stats_file)
{
    fprintf(stats_file, "stage\titer\tprior\tlikelihood\tjoint\trecombs\tnoncompats\targlen\n");
}


void print_stats(FILE *stats_file, const char *stage, int iter,
                 ArgModel *model, 
                 const Sequences *sequences, LocalTrees *trees,
                 const SitesMapping* sites_mapping, const Config *config)
{
    // calculate number of recombinations
    int nrecombs = trees->get_num_trees() - 1;

    // calculate number of non-compatiable sites
    int nseqs = trees->get_num_leaves();
    char *seqs[nseqs];
    for (int i=0; i<nseqs; i++)
        seqs[i] = sequences->seqs[trees->seqids[i]];
    int noncompats = count_noncompat(trees, seqs, nseqs, sequences->length());

    // get memory usage in MB
    double maxrss = get_max_memory_usage() / 1000.0;


    // calculate likelihood, prior, and joint probabilities
    // uncompressed local trees 
    if (sites_mapping)
        uncompress_local_trees(trees, sites_mapping);

    double prior = calc_arg_prior(&config->model, trees);
    double likelihood = calc_arg_likelihood(&config->model, sequences, trees,
                                            sites_mapping);
    double joint = prior + likelihood;
    double arglen = get_arglen(trees, config->model.times);

    // recompress local trees
    if (sites_mapping)
        compress_local_trees(trees, sites_mapping);
    
    // output stats
    fprintf(stats_file, "%s\t%d\t%f\t%f\t%f\t%d\t%d\t%f\n",
            stage, iter,
            prior, likelihood, joint, nrecombs, noncompats, arglen);
    fflush(stats_file);

    printLog(LOG_LOW, "\n"
             "prior:      %f\n"
             "likelihood: %f\n"
             "joint:      %f\n"
             "nrecombs:   %d\n"
             "noncompats: %d\n"
             "arglen:     %f\n"
             "max memory: %.1f MB\n\n",
             prior, likelihood, joint, nrecombs, noncompats, arglen, maxrss);

}

//=============================================================================
// sample output


// Returns the iteration-specific ARG filename
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
    string out_arg_file = get_out_arg_file(*config, iter);
    if (!config->no_compress_output)
        out_arg_file += ".gz";

    // write local trees uncompressed
    if (sites_mapping)
        uncompress_local_trees(trees, sites_mapping);

    // setup output stream
    CompressStream stream(out_arg_file.c_str(), "w");
    if (!stream.stream) {
        printError("cannot write '%s'", out_arg_file.c_str());
        return false;
    }
        

    write_local_trees(stream.stream, trees, sequences, model->times);
    
    if (sites_mapping)
        compress_local_trees(trees, sites_mapping);
    
    return true;
}


//=============================================================================

bool read_init_arg(const char *arg_file, const ArgModel *model, 
                   LocalTrees *trees, vector<string> &seqnames)
{
    CompressStream stream(arg_file, "r");
    if (!stream.stream) {
        printError("cannot read '%s'", arg_file);
        return false;
    }

    return read_local_trees(stream.stream, model->times, model->ntimes,
                            trees, seqnames);
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
        sample_arg_seq(model, sequences, trees, true);
        print_stats(config->stats_file, "seq", trees->get_num_leaves(),
                    model, sequences, trees, sites_mapping, config);
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
        print_stats(config->stats_file, "climb", i, model, sequences, trees,
                    sites_mapping, config);
    }
    printLog(LOG_LOW, "\n");
}


void resample_arg_all(ArgModel *model, Sequences *sequences, LocalTrees *trees,
                      SitesMapping* sites_mapping, Config *config)
{
    // setup search options
    double frac_leaf = .5;
    //int window = 100000;
    //int niters = 10;
    int window = config->resample_window;
    int niters = config->resample_window_iters;
    window /= config->compress_seq;
    int step = window / 2;

    // set iteration counter
    int iter = 1;
    if (config->resume)
        iter = config->resume_iter + 1;
    else {
        // save first ARG (iter=0)
        print_stats(config->stats_file, "resample", 0, model, sequences, trees,
                    sites_mapping, config);
        log_local_trees(model, sequences, trees, sites_mapping, config, 0);
    }


    printLog(LOG_LOW, "Resample All Branches (%d iterations)\n", 
             config->niters);
    printLog(LOG_LOW, "--------------------------------------\n");
    for (int i=iter; i<=config->niters; i++) {
        printLog(LOG_LOW, "sample %d\n", i);
        Timer timer;
        if (config->gibbs)
            resample_arg(model, sequences, trees);
        else
            resample_arg_mcmc_all(model, sequences, trees, frac_leaf,
                                  window, step, niters);
        printTimerLog(timer, LOG_LOW, "sample time:");

        
        // logging
        print_stats(config->stats_file, "resample", i, model, sequences, trees,
                    sites_mapping, config);

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
                    model, sequences, trees, sites_mapping, config);

        resample_arg_region(model, sequences, trees, 
                            config->resample_region[0], 
                            config->resample_region[1], 
                            config->niters);

        // logging
        print_stats(config->stats_file, "resample_region", config->niters,
                    model, sequences, trees, sites_mapping, config);
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
    string out_arg_file = get_out_arg_file(config, iter2);
    struct stat st;
    if (stat(out_arg_file.c_str(), &st) == 0) {
        stage = stage2;
        iter = iter2;
        arg_file = out_arg_file;
    }

    // try compress output
    out_arg_file += ".gz";
    if (stat(out_arg_file.c_str(), &st) == 0) {
        stage = stage2;
        iter = iter2;
        arg_file = out_arg_file;
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
        if (!parse_status_line(
            line, config, config.resume_stage, config.resume_iter, arg_file)) {
            delete [] line;
            return false;
        }
        delete [] line;
    }

    if (arg_file == "") {
        printLog(LOG_LOW, "Could not find any previously written ARG files. Try disabling resume\n");
        return false;
    }
    config.arg_file = arg_file;
    
    printLog(LOG_LOW, "resuming at stage=%s, iter=%d, arg=%s\n", 
             config.resume_stage.c_str(), config.resume_iter,
             config.arg_file.c_str());

    // clean up
    fclose(stats_file);

    return true;
}


bool check_overwrite(Config &config)
{
    // check for stats file    
    string stats_filename = config.out_prefix + STATS_SUFFIX;
    bool exists = !access(stats_filename.c_str(), F_OK);
    if (config.resume || config.overwrite || !exists) {
        return true;
    } else {
        printError("Stats file already exists: %s", stats_filename.c_str());
        printError("To force overwrite use --overwrite option");
        return false;
    }
}


bool ensure_output_dir(const char *outdir)
{
    char *path = strdup(outdir);
    char *dir = dirname(path);
    bool result = true;
    if (!makedirs(dir)) {
        printError("could not make directory for output files '%s'", dir);
        result = false;
    }
    free(path);
    return result;
}


//=============================================================================


int main(int argc, char **argv)
{
    // parse command line arguments
    Config c;
    int ret = c.parse_args(argc, argv);
    if (ret)
	return ret;

    // ensure output dir
    if (!ensure_output_dir(c.out_prefix.c_str()))
        return EXIT_ERROR;
    
    // check overwriting
    if (!check_overwrite(c))
        return EXIT_ERROR;
    
    // setup logging
    setLogLevel(c.verbose);
    string log_filename = c.out_prefix + LOG_SUFFIX;
    Logger *logger;
    if (c.quiet) {
        // log only to file
        logger = &g_logger;
    } else {
        // log to both stdout and file
        logger = new Logger(NULL, c.verbose);
        g_logger.setChain(logger);
    }
    const char *log_mode = (c.resume ? "a" : "w");
    if (!logger->openLogFile(log_filename.c_str(), log_mode)) {
        printError("could not open log file '%s'", log_filename.c_str());
        return EXIT_ERROR;
    }    
    
    // log intro
    if (c.resume)
        printLog(LOG_LOW, "RESUME\n");
    log_intro(LOG_LOW);
    log_prog_commands(LOG_LOW, argc, argv);
    Timer timer;


    // init random number generator
    if (c.randseed == 0)
        c.randseed = time(NULL);
    srand(c.randseed);
    printLog(LOG_LOW, "random seed: %d\n", c.randseed);


    // try to resume a previous run
    if (!setup_resume(c)) {
        printError("resume failed.");
        if (c.overwrite) {
            c.resume = false;
            printLog(LOG_LOW, "Resume failed.  Sampling will start from scratch since overwrite is enabled.\n");
        } else {
            return EXIT_ERROR;
        }
    }


    // read sequences
    Sites sites;
    Sequences sequences;
    SitesMapping *sites_mapping = NULL;
    auto_ptr<SitesMapping> sites_mapping_ptr;
    Region seq_region;
    Region seq_region_compress;

    
    if (c.fasta_file != "") {
        // read FASTA file
        
        if (!read_fasta(c.fasta_file.c_str(), &sequences)) {
            printError("could not read fasta file");
            return EXIT_ERROR;
        }
        seq_region.set("chr", 0, sequences.length());

        printLog(LOG_LOW, "read input sequences (nseqs=%d, length=%d)\n",
                 sequences.get_num_seqs(), sequences.length());
        
        // if compress requested, make sites object
        if (c.compress_seq > 1)
            make_sites_from_sequences(&sequences, &sites);
        
    } else if (c.sites_file != "") {
        // read sites file
        
        // parse subregion if given
        int subregion[2] = {-1, -1};
        if (c.subregion_str != "") {
            if (!parse_region(c.subregion_str.c_str(), 
                              &subregion[0], &subregion[1])) {
                printError("subregion is not specified as 'start-end'");
                return EXIT_ERROR;
            }
            subregion[0] -= 1; // convert to 0-index
        }

        // read sites
        CompressStream stream(c.sites_file.c_str());
        if (!stream.stream || 
            !read_sites(stream.stream, &sites, subregion[0], subregion[1])) {
            printError("could not read sites file");
            return EXIT_ERROR;
        }
        stream.close();

        printLog(LOG_LOW, "read input sites (chrom=%s, start=%d, end=%d, length=%d, nseqs=%d, nsites=%d)\n",
                 sites.chrom.c_str(), sites.start_coord, sites.end_coord,
                 sites.length(), sites.get_num_seqs(),
                 sites.get_num_sites());

        // sanity check for sites
        if (sites.get_num_sites() == 0) {
            printLog(LOG_LOW, "no sites given.  terminating.\n");
            return EXIT_ERROR;
        }
        seq_region.set(sites.chrom, sites.start_coord, sites.end_coord);
            
    } else {
        // no input sequence specified
        printError("must specify sequences (use --fasta or --sites)");
        return EXIT_ERROR;
    }

    // compress sequences
    if (sites.get_num_sites() > 0) {
        sites_mapping = new SitesMapping();
        sites_mapping_ptr = auto_ptr<SitesMapping>(sites_mapping);

        find_compress_cols(&sites, c.compress_seq, sites_mapping);
        compress_sites(&sites, sites_mapping);
        make_sequences_from_sites(&sites, &sequences);
    }
    seq_region_compress.set(seq_region.chrom, 0, sequences.length());



    // mask sequences
    TrackNullValue maskmap;
    if (c.maskmap != "") {
        // read mask
        CompressStream stream(c.maskmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &maskmap, seq_region)) {
            printError("cannot read mask map '%s'", 
                       c.maskmap.c_str());
            return EXIT_ERROR;
        }

        // apply mask
        if (sites_mapping)
            compress_mask(maskmap, sites_mapping);
        apply_mask_sequences(&sequences, maskmap);

        // report number of masked sites
        bool *masked = new bool [sequences.length()];
        find_masked_sites(sequences.get_seqs(), sequences.get_num_seqs(), 
                          sequences.length(), masked);
        int nmasked = 0;
        for (int i=0; i<sequences.length(); i++)
            nmasked += int(masked[i]);
        delete [] masked;
        printLog(LOG_LOW, "masked %d (%.1f%%) sites\n", nmasked, 
                 100.0 * nmasked / double(sequences.length()));
    }


    // setup model parameters
    if (c.times_file != "") {
        printError("not implemented yet");
        return EXIT_ERROR;
    } else if (c.time_step)
        c.model.set_linear_times(c.time_step, c.ntimes);
    else
        c.model.set_log_times(c.maxtime, c.ntimes);
    c.model.rho = c.rho;
    c.model.mu = c.mu;
    const double infsites_penalty = 1e-100; // TODO: make configurable
    if (c.infsites)
        c.model.infsites_penalty = infsites_penalty;
    c.model.set_popsizes(c.popsize, c.model.ntimes);

    // read model parameter maps if given
    if (c.mutmap != "") {
        CompressStream stream(c.mutmap.c_str(), "r");
        if (!stream.stream || 
            !read_track_filter(stream.stream, &c.model.mutmap, seq_region)) {
            printError("cannot read mutation rate map '%s'", c.mutmap.c_str());
            return EXIT_ERROR;
        }
    }
    if (c.recombmap != "") {
        CompressStream stream(c.recombmap.c_str(), "r");
        if (!stream.stream ||
            !read_track_filter(stream.stream, &c.model.recombmap, seq_region)) {
            printError("cannot read recombination rate map '%s'", 
                       c.recombmap.c_str());
            return EXIT_ERROR;
        }
    }
    
    // make compressed model
    ArgModel model(c.model);
    if (!model.setup_maps(seq_region.chrom, seq_region.start, seq_region.end))
        return EXIT_ERROR;
    compress_model(&model, sites_mapping, c.compress_seq);
    
    // log original model
    log_model(model);

    
    // setup init ARG
    LocalTrees *trees = NULL;
    auto_ptr<LocalTrees> trees_ptr;
    if (c.arg_file != "") {
        // init ARG from file
        
        trees = new LocalTrees();
        trees_ptr = auto_ptr<LocalTrees>(trees);
        vector<string> seqnames;
        if (!read_init_arg(c.arg_file.c_str(), &c.model, trees, seqnames)) {
            printError("could not read ARG");
            return EXIT_ERROR;
        }

        if (!trees->set_seqids(seqnames, sequences.names)) {
            printError("input ARG's sequence names do not match input sequences");
            return EXIT_ERROR;
        }
        
        printLog(LOG_LOW, "read input ARG (chrom=%s, start=%d, end=%d, nseqs=%d)\n",
                 trees->chrom.c_str(), trees->start_coord, trees->end_coord, 
                 trees->get_num_leaves());

        // check ARG matches sites/sequences
        if (trees->start_coord != seq_region.start || 
            trees->end_coord != seq_region.end) {
            printError("trees range does not match sites: tree(start=%d, end=%d), sites(start=%d, end=%d) [compressed coordinates]", 
                       trees->start_coord, trees->end_coord, 
                       seq_region.start, seq_region.end);
            return EXIT_ERROR;
        }

        // compress input ARG if compression is requested
        if (sites_mapping)
            compress_local_trees(trees, sites_mapping, true);
        
    } else {
        // create new init ARG
        trees = new LocalTrees(seq_region_compress.start, 
                               seq_region_compress.end);
        trees->chrom = seq_region.chrom;
        trees_ptr = auto_ptr<LocalTrees>(trees);
    }
    

    // check for region sample
    if (c.resample_region_str != "") {
        if (!parse_region(c.resample_region_str.c_str(), 
                          &c.resample_region[0], &c.resample_region[1]))
        {
            printError("--resample-region is not specified as 'start-end'");
            return EXIT_ERROR;
        }
        c.resample_region[0] -= 1; // convert to 0-index

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
        return EXIT_ERROR;
    }

    // get memory usage in MB
    double maxrss = get_max_memory_usage() / 1000.0;
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);

    // sample ARG
    printLog(LOG_LOW, "\n");
    sample_arg(&model, &sequences, trees, sites_mapping, &c);
    
    // final log message
    maxrss = get_max_memory_usage() / 1000.0;
    printTimerLog(timer, LOG_LOW, "sampling time: ");
    printLog(LOG_LOW, "max memory usage: %.1f MB\n", maxrss);
    printLog(LOG_LOW, "FINISH\n");

    // clean up
    fclose(c.stats_file);
    
    
    return 0;
}
