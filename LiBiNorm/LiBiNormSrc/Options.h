// ***************************************************************************
// Options.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 3 May 2018
// ---------------------------------------------------------------------------
// Top level code options
// ***************************************************************************

#ifndef OPTIONS_H
#define OPTIONS_H


#define LIBINORM_VERSION "2.0"
//	Bam/gff file reading
#define DEFAULT_FEATURE_TYPE_EXON "exon" 
#define DEFAULT_GTF_ID_ATTRIBUTE "gene_id"
#define DEFAULT_GFF_ID_ATTRIBUTE "gene"
#define DEFAULT_COUNT_MODE intersect_strict
#define DEFAULT_COUNT_MODE_HTSEQ_COMPATIBLE intersect_union

#define DEF_THREADS 3 //-p

//	The transcript types that are ignored when parsing gtf files.
//	These are included in htseq_compatible mode
#define IGNORED_GTF_TRANSCRIPT_TYPES  "retained_intron" 

//	Put reads into cache file when number of reads exceed READ_CACHE_SIZE
#ifdef _DEBUG
#define READ_CACHE_SIZE 50000
#else
#define READ_CACHE_SIZE 2000000
#endif

//	Report progress every REP_LEN entries
#define REP_LEN 100000

//	When matching forward and reverse reads the position information can also used to pair the reads
//	anabling this tends to worsen the mismatch between htseq-count and LiBinomrm, particularly when data is
//	cached
#define MATCH_USING_BOTH_POSITIONS

//	Read selection for paremeter estimation
#define DEF_MAX_READS_FOR_PARAM_ESTIMATION 100000000  // -d
#define MAX_READS_GENE 100

//	Nelder Mead
#define USE_NELDER_MEAD_FOR_INITIAL_VALUES
#define NELDER_MEAD_ITERATIONS 2000
#define EDGE_PENALTY_MULTIPLIER 1000

//	Prior: Nelder Mead and MCMC slope for parameters.  By setting a multiplier to zero the code is disabled
//	These provide common priors for all models, as implemented in priorFunc in ModelData.cpp  
// #define D_PRIOR_TARGET 0.6
// #define D_PRIOR_MULTIPLIER 0
// #define H_PRIOR_TARGET 100
// #define H_PRIOR_MULTIPLIER 0
//	These provides an additional priors for parameter H model E, as implemented in priorFuncE in ModelData.cpp  
#define E_H_PRIOR_TARGET 20
#define E_H_PRIOR_MULTIPLIER 2
// #define E_D_PRIOR_TARGET 0.5
// #define E_D_PRIOR_MULTIPLIER 2

//	MCMC operation
#define NUMBER_OF_MCMC_RUNS 10  // -r
#define MCMC_ITERATIONS 2000  // -s     // MCMC iterations if no Nelder Mead initialisation
#define NELDER_MCMC_ITERATIONS 200  // -s   MCMC iterations if we pre-initialise using Nelder Mead
#define MCMC_JUMP_SIZE 0.01
#define DEFAULT_MODEL ModelBD // -n
#define MCMC_SIGMA 1	

//	Outputting results
#define DEFAULT_NORMALISATION_GENE_LENGTH 1000
#define MAX_GENE_LENGTH_FOR_NORM_PLOT 20000 

//	Allows a maximum gene length to be specified
#define DEF_LENGTH_OF_GENE_FOR_PARAM_ESTIMATION 20000 // -e

//	Parameters for printing bias plots
#define N_BIAS_BINS 100
#define N_BIAS_GENE_SEGMENTS 500
#define READ_THRESHOLD 0	//Genes with fewer reads are ignored
#define MAX_VALUE 0.04		//Output values clipped to this


//	The insert size for matching pairs should be A and -A.  In some datasets they are A and A.  By using the absolute value of the
//	insert size we can ensure that the pairs are still matched up
#define USE_ABS_INSERT_TO_MATCH_READS

//	Prints a file containing the predicted read distrubtion for different mRNA lengths
#define PRINT_DISTRIBUTION

//  Output gtf file as bed data
#define OUTPUT_BED_DATA

//	Output selected genes from gtf file as FASTA file
#define OUTPUT_FASTA_FILE

#ifdef OUTPUT_FASTA_FILE
#define HISAT2
#endif

//	******************************************************************************
//	The following are options that would normally be disabled in the release version of the code


//	Use the original option wehere 'best' selected the model with the lowest log liklihood.
// #define SELECT_BY_LL

//	The seed for the random number generator used to select reads.  If undefined then a random seed is generated
// #define SELECT_READS_SEED 2017

//	Use this to add the mode which creates fastq files based on bam files with artificial problems
// #define MAKE_FASTQ_MODE

//	Use this to enable various addditional tools for exploring LiBiNorm data.  This currently
//	provides the additional land, land2 and gene run modes
// #define LIBITOOLS

//	Use this to add the mode where duplicates in bam files can be removed
// #define DEDUP_MODE

//  Output detailed results of interpreting the Feature file
// #define OUTPUT_FEATURE_DATA

//	Output info on how each read is identified
// #define OUTPUT_READ_MAPPING_INFO

//	Adds option of just use selected genes (cont mode)
// #define USE_GENES_FROM_GENELIST

//	Adds a menu option that allows the seed to be specified
// #define SELECT_READ_SEED

//	Use this mode to run a model with specific parameters which are loaded in using the -i
//	command.  This affects how the values are setin ModelParameters.cpp
// #define INITIAL_VALUES

//	Data from longer and shorter genes are given additional weightings.  The current values differ from the original
//	Use this to reset to the original
// #define ORIGINAL_WEIGHTING

//	For consistent selection of reads and reproducing the MATLAB algorithms
// #define REPRODUCE_MATLAB	

//	Adds -x option where the program will pause at the end rather than simply exiting.  Useful for debugging
// #define PAUSE_AT_END_OPTION	

// Adds -w option which results in the output of additional debug messages
// #define OUTPUT_DEBUG_MESSAGES

//	Some of the code in ModelData.cpp has also been writtent using vectors which is slower but the code
//	more closely matches the MATLAB code
// #define VECTOR_MATHS

//	The original MATLAB code had an error in setting the initial values for mcmc runs which this 
//	simulates (ModelData.cpp)
// #define SIMULATE_MATLAB_BUG

//	Output landscape files in the original format
// #define LANDSCAPE_FORMAT_1

//	Enable this to produce a heatmap data file that can be used by the mathematica script
// #define MATHMATICA_FILE

//  Print full MCMC run data
// #define PRINT_MCMC_RUN_DATA

//	Create a landscape file using the reads that were actually used for parameter discovery
// #define MAKE_LANDSCAPE_OF_USED_READS

//	Just use some genes for parameter discovery.  The -g paremeter in model mode is used to group the genes into blocks of N such
//	that only alternate blocks are used for parameter discovery.  The names of genes that are not used are output as XXX in the count
//	file.  -N also divides the genes into blocks of N but swaps which blocks are used ano not used for parameter discoverey
// #define USE_GROUPS_OF_GENES_FOR_DISCOVERY

//	Convert the input data into flattened data with the same number of reads in each gene, but distributed evenly through the gene
//	but with spaces at either end.  Genes with fewer than N genes are excluded
// #define CREATE_FLAT_DATA

//	Use absolute rather than log h value when doing Nelder Mead and MCMC
// #define ABS_H_PARAM

//	Causes the code to halt waiting for user input before exiting after a failure
#ifdef _DEBUG
#define PAUSE_ON_EXIT_FAILURE
#endif

//	Use this option to use the parameters associated with the most likly parameter set
//	https://sciencehouse.wordpress.com/2010/06/23/mcmc-and-fitting-models-to-data/
//	rather than the median values
// #define USE_PARAMS_FROM_LOWEST_LL


#endif