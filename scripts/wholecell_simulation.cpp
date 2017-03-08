// Source code for whole-cell simulation of E. coli protein synthesis.
// Author: Arvind R. Subramaniam
//
// Source code was adapted from following references:
//
// Subramaniam AR, Zid BM, O'Shea EK.
// Cell 159 (5): 1200—1211 (2014)
//
// Shah P, Ding Y, Niemczyk M, Kudla G, and Plotkin JB.
// Cell 153 (7): 1589-1601 (2013)

///////////////////////////////////////////////////////////////////////////////
// DECLARING HEADER FILES

using namespace std;
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <sys/time.h>
#include <sys/stat.h>
#include <cstdlib>
#include <string.h>
#include <boost/regex.hpp>
#include <unistd.h>
#include <fstream>

///////////////////////////////////////////////////////////////////////////////
// FIXED PARAMETERS

#define AVOGADRO_NUMBER 6.023e23
// Max gene length required for initialization of mRNA data structures
#define MAX_GENE_LEN_ALLW 5000 

///////////////////////////////////////////////////////////////////////////////
// DEFAULT GLOBAL VARIABLES

// Seed for random number generator.
int seed = 111;
// Includes peptidyl-transfer and translocation time,
// but not codon-dependent tRNA selection time
// which is added separately below,
// command-line option provided.
double codonMinElongationTime = 0.045;
// Number of genes,
// required for initiatization of data structure.
int totalNumberOfGenes = 1067;
// Total ribosomes.  Command line option provided.
int totalNumberOfRibosomes = 7.3e4;
// Fraction of ribosomes bound to mRNA
double boundRibosomeFraction = 0.85;
// Total tRNAs
int totalNumberOfTrna = 6.8e5;        
// Total mRNAs, only initialization here
// Actual value calculated from mRNA copy number file.
int totalNumberOfMrna = 0;
// Total space within a cell
double cellVolume = 1.0e-18;
// codons covered by a single ribosome
int ribosomeFootprintSize = 10;
// Rate at which ribosomes prematurely terminate during normal elongation cycle.
// command-line option provided.
double backgroundPrematureTerminationRate = 0;
// Rate (s^-1) at which ribosomes selectively terminate when 
// the elongation rate falls below threshold tRNA accommodation rate (s^-1)
// command-line option provided.
double thresholdTrnaAccommodationRate = 0;
double selectivePrematureTerminationRate = 0;
// Rate of premature termination for ribosomes that have
// been hit by another ribosome from 5' side (s^-1)
// command-line option provided.
double prematureTerminationRateUpon5PrimeHit = 0;
// Rate of premature termination for ribosomes that have
// been hit by another ribosome from 3' side (s^-1)
// command-line option provided.
double prematureTerminationRateUpon3PrimeHit = 0;
// kcat/km for ternary complex ribosome association,
// commandline option provided.
// Command line option provided.
double ribosomeTrnaKcatByKm = 2.0e7;
// starvationfactor less than 1 implies that
// the cell is starving for an amino acid.
// Command line option provided.
double starvationFactor = 1;
// These factors multiply
// the aminoacylation rate constants for the tRNAs
// corresponding to the liming amino acid.
// Command line option provided.
double relativeTrnaAminoacylationRates[5]
= {1,1,1,1,1};
// Run options
// Controls what outputs are printed to files.
// Look at end of main function for file details.
int printOptions[5] = {1,1,1,1,1};
// Prefix for output file names.
// Common-line option provided.
char outfilePrefix[350] = "";
// Variable to temporarily store
// the names of different output file names.
char outfile[350] = "x";    
// This file will contain the initiation rate of mRNAs,
// mRNA copy number in the cell,a numerical codon sequence
// for all mRNAs that are simulated.
// Set from command-line.
char mrnaFile[350] = "";           
// Input Files specifying parameters for different tRNA species.
// tRNA concentrations.
char trnaConcentrationFile[150] = "../annotations/simulations/wholecell.trna.concentrations.tsv";
// Cognate tRNA for each codon.
char cognateTrnaFile[150] = "../annotations/simulations/wholecell.cognate.pairs.tsv";
// Codon-tRNA interaction strength.
char pairingStrengthFile[150] = "../annotations/simulations/wholecell.cognate.weights.tsv";
// Aminoacylation rates of all tRNAs.
char aminoacylationRateFile[150] = "../annotations/simulations/wholecell.trna.aminoacylation.rate.tsv";
// Total time for simulation.
// Command line option provided.
double totalTime = 15;
// Threshold time for starting analysis.
// Command line option provided.
double thresholdTime = 10;

////////////////////////////////////////////////////////////////////////////////
// CREATE CLASSES for different molecule species - ribosome, mrna, gene, trna,
// codon
class Ribosome {
  // Note that new ribosome objects are created only when
  // they initiate on an mRNA. This saves space, but might
  // be modified in the future to keep track of both free
  // and mRNA-bound ribosomes.
  public:
    // Bound to which mRNA
    int mRNA;
    // Position on mRNA
    int pos;
    // tRNA at the A-site
    int AtRNA;
    // tRNA at the P-site
    int PtRNA;
    // Variable that keeps track of whether the ribosome has empty A site
    bool AsiteEmpty;
    // Variables that keeps track of whether the ribosome was
    // hit from 5' and/or 3' sides
    bool hitFrom5Prime;
    bool hitFrom3Prime;
    // The variables below are not biological (meaning, a real ribosome
    // are not expected to have these state variables), but they are kept
    // track during the simulation
    // Time of translation initiation
    double timeOfLastInitiation;
    // Time of arrival at current codon, used for calculating ribosome density.
    double timeOfLastElongation;
    // Position in the list A site empty ribosomes
    int aempty_list_pos;
    // Position in the list of translocation-competent ribosomes
    int tl_list_pos;
    // Positions in the list of 5' or 3' hit ribosomes
    int hit5_list_pos;
    int hit3_list_pos;
};

class transcript {
  public:
    // Gene id
    int gene;
    // Number of initiation events
    int numberOfInitiations;
    // Number of translation events
    int numberOfFullProteins;
    // Number of translation events
    int numberOfAbortedProteins;
    // Time of last initiation event
    double timeOfLastMrnaInitiation;
    // Average time to initiation
    double averageInitiationTime;
    // Average translation time
    double averageElongationTime;
    // Average translation time
    double averageAbortionTime;
};

class gene {
  public:
    int seq[MAX_GENE_LEN_ALLW];// Codon sequence of the gene
    int length;// Length of the gene
    int exp;// Gene expression level
    double ini_prob;// Initiation probability of the mRNA
};

class trna {
  public:
    int id;// tRNA id, ranges between 1 and 40.
    char anticodon[4];// unmodified anticodon sequence
    char aminoacid[4];// three letter amino acid code
    double concn;// tRNA concentration
    double trnaRibosomeRateConstant;// kcat/Km for the trna-ternary complex to
    // bind to the ribosome at the A-site.
    double aminoacylationRateConstant;// this can be decreased to account for charging.
};

class codon {
  public:
    int id;// codon id, ranges from 0 to 60.
    char triplet[4];// codon seuquence.
    int trna[2];// cognate trna ids, range from 0 to 40.
    // 0 implies no second cognate trna.
    double weight[2];// codon-trna pairing weight.
    // 0 implies no second cognate trna.
};

////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS FOR READING INPUT FILES AND COMMAND-LINE PARAMETERS

// Reading the processed sequence file
int Read_sequence_File(char *filename, gene *Genes) {
   FILE *fh;
  int c1=0,c2=0;
  char curr_char;
  char buffer[5000];

  fh=fopen(filename, "r");

  fscanf(fh,"%lf",&Genes[c1].ini_prob);

  do {
    fscanf(fh,"%d",&Genes[c1].exp);

    c2 = 0;
    do {
      fscanf(fh,"%d",&Genes[c1].seq[c2]);
      c2++;
      curr_char = fgetc(fh);
    }while(curr_char != '\n');
    Genes[c1].length = c2;
    c1++;
  }while(fscanf(fh,"%lf",&Genes[c1].ini_prob) ==1);
}

// Read in tRNA concentration file
int Read_tRNA_File(char *trnaConcnFile, char *trnaAminoacylationRateFile, trna *trna_list) {
   FILE *fh1, *fh2;
  int loop1=0;
  int loop2=0;
  char buffer[256];
  int id;
  double rate;

  fh1=fopen(trnaConcnFile, "r");

  fgets( buffer, sizeof(buffer), fh1);// Read the header line first.

  while(fscanf(fh1,"%d",&trna_list[loop1].id)==1) {
    fscanf(fh1,"%s%lf%s",
        &trna_list[loop1].anticodon,
        &trna_list[loop1].concn,
        &trna_list[loop1].aminoacid);
    loop1++;
  }

  buffer[0] = 0;

  fh2=fopen(trnaAminoacylationRateFile, "r");

  fgets( buffer, sizeof(buffer), fh2);// Read the header line first.

  while(fscanf(fh2,"%d",&id)==1) {
    char anticodon[4];// I do not store the anticodon.
    fscanf(fh2,"%s%lf", &anticodon, &trna_list[loop2].aminoacylationRateConstant);
    loop2++;
  }

}

// Read in codon-trna pairs and codon-trna weights.
int Read_Codon_Files(char *cognateTrnaFile, char *pairingStrengthFile, codon *Codon) {
   FILE *pair_file;
  FILE *weight_file;
  int loop1=0;
  char buffer[256];

  pair_file=fopen(cognateTrnaFile, "r");

  if (!pair_file)                 // Check if file exists
  {   printf("\nPair file doesn't exist.\n");
    exit(1);
  }

  weight_file=fopen(pairingStrengthFile, "r");

  if (!weight_file)                 // Check if file exists
  {   printf("\nWeight file doesn't exist.\n");
    exit(1);
  }

  fgets( buffer, sizeof(buffer), pair_file);// Read the header line first.
  fgets( buffer, sizeof(buffer), weight_file);// Read the header line first.

  while(fscanf(pair_file,"%d",&Codon[loop1].id)==1) {
    char codon1[4], codon2[4];// I do not store the codon.
    int temp_id = 0;
    fscanf(pair_file,"%s%d%d", &Codon[loop1].triplet, &Codon[loop1].trna[0], &Codon[loop1].trna[1]);
    fscanf(weight_file,"%d%s%lf%lf", &temp_id, &codon2, &Codon[loop1].weight[0], &Codon[loop1].weight[1]);
    loop1++;
  }
}


// Read in commandline arguments
void Read_Commandline_Args(int argc, char *argv[], trna *trna_list) {
  int argIndex;
  FILE *fh;
  int c1;
  int trnaloop;
  for (argIndex=1;argIndex<argc;argIndex++) {
    if ( !strcmp(argv[argIndex], "--output-prefix") ) {
      strcat( outfilePrefix, argv[argIndex+1] );
    }
    if ( !strcmp(argv[argIndex], "--input-genes") ) {
      strcat( mrnaFile, argv[argIndex+1] );
      strcat( outfilePrefix, "mrnafile_" );
      string s(argv[argIndex+1]);
      boost::smatch m;
      boost::regex e("/([^/]*).csv");

      if (boost::regex_search(s,m,e)) {
        string temp(m[1]);
        strcat( outfilePrefix, temp.c_str() );
      }
      
    }
  }
  for (argIndex=1;argIndex<argc;argIndex++) {
    if ( !strcmp(argv[argIndex], "--total-ribosomes") ) {
      totalNumberOfRibosomes = atoi( argv[argIndex+1] );
      strcat( outfilePrefix, "_ribosomes_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }
    if ( !strcmp(argv[argIndex], "--total-genes") ) {
      totalNumberOfGenes = atoi( argv[argIndex+1] );
      strcat( outfilePrefix, "_numberofgenes_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }    
    if ( !strcmp(argv[argIndex], "--ribotrnakcatkm") ) {
      ribosomeTrnaKcatByKm = atof( argv[argIndex+1] );
      strcat( outfilePrefix, "_ribotrnakcatkm_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }
    if ( !strcmp(argv[argIndex], "--min-elong-time") ) {
      codonMinElongationTime = atof( argv[argIndex+1] );
      strcat( outfilePrefix, "_minelongtime_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }
    if ( !strcmp(argv[argIndex], "--aminoacid") ) {
      strcat( outfilePrefix, "_" );
      strcat( outfilePrefix, argv[argIndex+1] );
      if ( !strcmp(argv[argIndex+2], "starvationfactor" ) ) {
        starvationFactor = atof(argv[argIndex+3]);
        strcat( outfilePrefix, "_starvationfactor_" );
        strcat( outfilePrefix, argv[argIndex+3] );
        for (trnaloop=0;trnaloop<61;trnaloop++) {
          if ( !strcmp( trna_list[trnaloop].aminoacid, argv[argIndex+1] ) ) {
            trna_list[trnaloop].aminoacylationRateConstant *= starvationFactor;
          }
        }
      }

      if ( !strcmp(argv[argIndex+2], "totalcognatetrnafactor" ) ) {
        starvationFactor = atof(argv[argIndex+3]);
        strcat( outfilePrefix, "_totalcognatetrnafactor_" );
        strcat( outfilePrefix, argv[argIndex+3] );
        for (trnaloop=0;trnaloop<61;trnaloop++) {
          if (!strcmp( trna_list[trnaloop].aminoacid, argv[argIndex+1] )) {
            trna_list[trnaloop].concn *= starvationFactor;
          }
        }
      }

      if ( !strcmp(argv[argIndex+2], "relativetrnaaminoacylationrates" ) ) {
        strcat( outfilePrefix, "_relativetrnaaminoacylationrates" );
        c1 = 0;
        while (argv[c1+argIndex+3][0] != '-' ) {
          if (c1+argIndex+3 == argc-1) {
            break;
          }
          relativeTrnaAminoacylationRates[c1] = atof(argv[c1+argIndex+3]);
          strcat( outfilePrefix, "_" );
          strcat( outfilePrefix, argv[c1+argIndex+3] );
          c1 += 1;
        }
        c1 = 0;
        for (trnaloop=0;trnaloop<61;trnaloop++) {
          if ( !strcmp( trna_list[trnaloop].aminoacid, argv[argIndex+1] ) ) {
            trna_list[trnaloop].aminoacylationRateConstant *= relativeTrnaAminoacylationRates[c1];
            c1 += 1;
          }
        }
      }
    }

    if ( !strcmp(argv[argIndex], "--background-preterm-rate") ) {
      backgroundPrematureTerminationRate = atof(argv[argIndex+1]);
      strcat( outfilePrefix, "_bkgdpreterm_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }

    if ( !strcmp(argv[argIndex], "--selective-preterm-rate") ) {
      selectivePrematureTerminationRate = atof(argv[argIndex+1]);
      strcat( outfilePrefix, "_selpreterm_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }

    if ( !strcmp(argv[argIndex], "--threshold-accommodation-rate") ) {
      thresholdTrnaAccommodationRate = atof(argv[argIndex+1]);
      strcat( outfilePrefix, "_thresholdtrnarate_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }

    if ( !strcmp(argv[argIndex], "--5prime-preterm-rate") ) {
      prematureTerminationRateUpon5PrimeHit = atof(argv[argIndex+1]);
      strcat( outfilePrefix, "_5primepreterm_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }

    if ( !strcmp(argv[argIndex], "--3prime-preterm-rate") ) {
      prematureTerminationRateUpon3PrimeHit = atof(argv[argIndex+1]);
      strcat( outfilePrefix, "_3primepreterm_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }

    if ( !strcmp(argv[argIndex], "--threshold-time") ) {
      thresholdTime = atof(argv[argIndex+1]);
      strcat( outfilePrefix, "_thresholdtime_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }

    if ( !strcmp(argv[argIndex], "--total-time") ) {
      totalTime = atof(argv[argIndex+1]);
      strcat( outfilePrefix, "_totaltime_" );
      strcat( outfilePrefix, argv[argIndex+1] );
    }
  }

  // This makes a new directory with an output prefix based on all
  // command-line arguments
  mkdir(outfilePrefix, 0700);
}

// This function updates all the state variables when
// a ribosome leaves an mRNA due to either normal or
// premature termination.
int TerminateRibosome(
    //Long list of arguments starts here:
    int r_id,
    int m_id,
    int c_id,
    Ribosome* ribosome_list,
    trna* trna_list,
    transcript* mRNA,
    gene* Genes,
    int* ribosomeBoundTrna, int* deacylatedFreeTrna,
    double &totalAminoacylationRate,
    int** R_grid,
    int** mrna_initiable,
    int &Rf,
    int pos_in_aempty_list,
    int pos_in_tl_list,
    int pos_in_hit5_list,
    int pos_in_hit3_list,
    int* n_Rb_ae, int &n_Rb_ao, int &n_Rb_hit5, int &n_Rb_hit3,
    int** Rb_ae, int* Rb_ao, int* Rb_hit5, int* Rb_hit3,
    int &total_n_mrna_free,
    int* n_mrna_free,
    double* initiationScaledMrnaCopyNumber,
    double &totalInitiationScaledMrnaCopyNumber,
    double &currentTime,
    double* summedElongationTimes,
    double** avgRibosomeDensityGrid,
    int** numberOfElongationsGrid,
    int* numberOfElongations,
    int* numberOfAbortions,
    int &next_avail_ribo
    ) {
    int fivePrimeRibosomePosition,
        fivePrimeRibosomeId,
        threePrimeRibosomePosition,
        threePrimeRibosomeId;
    int g_id;

  // Release P-site tRNA
  // We do not track the initiator tRNA at the start codon
  if (ribosome_list[r_id].pos > 1) {
    ribosomeBoundTrna[ribosome_list[r_id].PtRNA]--;
    deacylatedFreeTrna[ribosome_list[r_id].PtRNA]++;
    // Update aminoacylation rate.
    totalAminoacylationRate += trna_list[ribosome_list[r_id].PtRNA].aminoacylationRateConstant;
  }
  // Release A-site tRNA if it was occupied
  if (ribosome_list[r_id].AsiteEmpty == false) {
    ribosomeBoundTrna[ribosome_list[r_id].AtRNA]--;
    deacylatedFreeTrna[ribosome_list[r_id].AtRNA]++;
    // Update aminoacylation rate.
    totalAminoacylationRate += trna_list[ribosome_list[r_id].AtRNA].aminoacylationRateConstant;
  }

  // Update the corresponding mRNA position to empty.
  R_grid[m_id][ribosome_list[r_id].pos] = totalNumberOfRibosomes;

  // Free a ribosome upon termination
  Rf++;

  // Update number and list of ribosomes with empty A site
  if (ribosome_list[r_id].AsiteEmpty == true) {
    pos_in_aempty_list = ribosome_list[r_id].aempty_list_pos;
    n_Rb_ae[c_id]--;
    Rb_ae[c_id][pos_in_aempty_list] = Rb_ae[c_id][n_Rb_ae[c_id]];
    ribosome_list[Rb_ae[c_id][n_Rb_ae[c_id]]].aempty_list_pos = pos_in_aempty_list;
  }

  // Update number and list of 5' hit ribosomes
  if (ribosome_list[r_id].hitFrom5Prime == true) {
    pos_in_hit5_list = ribosome_list[r_id].hit5_list_pos;
    n_Rb_hit5--;
    Rb_hit5[pos_in_hit5_list] = Rb_hit5[n_Rb_hit5];
    ribosome_list[Rb_hit5[n_Rb_hit5]].hit5_list_pos = pos_in_hit5_list;
  }

  // Update number and list of 3' hit ribosomes
  if (ribosome_list[r_id].hitFrom3Prime == true) {
    pos_in_hit3_list = ribosome_list[r_id].hit3_list_pos;
    n_Rb_hit3--;
    Rb_hit3[pos_in_hit3_list] = Rb_hit3[n_Rb_hit3];
    ribosome_list[Rb_hit3[n_Rb_hit3]].hit3_list_pos = pos_in_hit3_list;
  }

  // Update number and list of translocation competent ribosomes
  // In our current model, this situation will never be satisfied.
  if (ribosome_list[r_id].AsiteEmpty == false &&
      ribosome_list[r_id].hitFrom3Prime == false) {
    pos_in_tl_list = ribosome_list[r_id].tl_list_pos;
    n_Rb_ao--;
    Rb_ao[pos_in_tl_list] = Rb_ao[n_Rb_ao];
    ribosome_list[Rb_ao[n_Rb_ao]].tl_list_pos = pos_in_tl_list;
  }

  // If there is an immediate 5' ribosome, it is not hit from 3' anymore,
  // and it is now translocation competent if its A-site is occupied.
  fivePrimeRibosomePosition = ribosome_list[r_id].pos - ribosomeFootprintSize;
  fivePrimeRibosomeId = totalNumberOfRibosomes;
  if (fivePrimeRibosomePosition > 0) {
    fivePrimeRibosomeId = R_grid[m_id][fivePrimeRibosomePosition];
  }
  if (fivePrimeRibosomeId < totalNumberOfRibosomes &&
      ribosome_list[fivePrimeRibosomeId].hitFrom3Prime == true) {
    pos_in_hit3_list = ribosome_list[fivePrimeRibosomeId].hit3_list_pos;
    ribosome_list[fivePrimeRibosomeId].hitFrom3Prime = false;
    n_Rb_hit3--;
    Rb_hit3[pos_in_hit3_list] = Rb_hit3[n_Rb_hit3];
    ribosome_list[Rb_hit3[n_Rb_hit3]].hit3_list_pos = pos_in_hit3_list;

    ribosome_list[fivePrimeRibosomeId].tl_list_pos = n_Rb_ao;
    Rb_ao[n_Rb_ao] = fivePrimeRibosomeId;
    n_Rb_ao++;
  }

  // If there is an immediate 3' ribosome, it is not hit from 5' anymore.
  threePrimeRibosomePosition = ribosome_list[r_id].pos + ribosomeFootprintSize;
  threePrimeRibosomeId = totalNumberOfRibosomes;
  if (threePrimeRibosomePosition <= Genes[mRNA[m_id].gene].length - 1) {
    threePrimeRibosomeId = R_grid[m_id][threePrimeRibosomePosition];
  }
  if (threePrimeRibosomeId < totalNumberOfRibosomes &&
      ribosome_list[threePrimeRibosomeId].hitFrom5Prime == true) {
    pos_in_hit5_list = ribosome_list[threePrimeRibosomeId].hit5_list_pos;
    ribosome_list[threePrimeRibosomeId].hitFrom5Prime = false;
    n_Rb_hit5--;
    Rb_hit5[pos_in_hit5_list] = Rb_hit5[n_Rb_hit5];
    ribosome_list[Rb_hit5[n_Rb_hit5]].hit5_list_pos = pos_in_hit5_list;
  }

  // The mRNA is now initiation-competent if a 'ribosomefootprintsize'
  // region got empty 
  if (ribosome_list[r_id].pos <= ribosomeFootprintSize) {
    // find gene id
    g_id = mRNA[m_id].gene;
    // add the mrna to the list of free mrnas for that gene
    mrna_initiable[g_id][n_mrna_free[g_id]] = m_id;
    // increase the number of free mrnas for that gene
    n_mrna_free[g_id]++;
    // increase total number of free mrnas
    total_n_mrna_free++;
    // increase the initiation scaled mrna copy number for that gene
    initiationScaledMrnaCopyNumber[g_id] += Genes[g_id].ini_prob;
    // increase the total initiation scaled mrna copy number
    totalInitiationScaledMrnaCopyNumber += Genes[g_id].ini_prob;
  }

  // Update variables that are tracked during the simulation
  if (currentTime>thresholdTime) {
    // For estimation of avg elongation times of codons
    summedElongationTimes[c_id] += currentTime-ribosome_list[r_id].timeOfLastElongation;
    // Keep track of the ribosome density at the corresponding codon position.
    avgRibosomeDensityGrid[m_id][ribosome_list[r_id].pos] += currentTime-ribosome_list[r_id].timeOfLastElongation;
    // Keep track of the number of times that this position on this mRNA was elongated / aborted.
    numberOfElongationsGrid[m_id][ribosome_list[r_id].pos] ++;
    if (ribosome_list[r_id].pos >= Genes[mRNA[m_id].gene].length - 1) {
      mRNA[m_id].numberOfFullProteins++;
      mRNA[m_id].averageElongationTime += currentTime - ribosome_list[r_id].timeOfLastInitiation;
      numberOfElongations[c_id]++;
    } else {
      mRNA[m_id].numberOfAbortedProteins++;
      mRNA[m_id].averageAbortionTime += currentTime - ribosome_list[r_id].timeOfLastInitiation;
      numberOfAbortions[c_id]++;
    }
  }

  // The data of the terminated ribosome is replaced by that of 
  // the most recent initiated ribosome.
  next_avail_ribo--;
  if (r_id!=next_avail_ribo) {
    ribosome_list[r_id].mRNA = ribosome_list[next_avail_ribo].mRNA;
    ribosome_list[r_id].pos = ribosome_list[next_avail_ribo].pos;
    ribosome_list[r_id].AtRNA = ribosome_list[next_avail_ribo].AtRNA;
    ribosome_list[r_id].PtRNA = ribosome_list[next_avail_ribo].PtRNA;
    ribosome_list[r_id].hitFrom5Prime = ribosome_list[next_avail_ribo].hitFrom5Prime;
    ribosome_list[r_id].hitFrom3Prime = ribosome_list[next_avail_ribo].hitFrom3Prime;
    ribosome_list[r_id].AsiteEmpty = ribosome_list[next_avail_ribo].AsiteEmpty;
    ribosome_list[r_id].timeOfLastInitiation = ribosome_list[next_avail_ribo].timeOfLastInitiation;
    ribosome_list[r_id].timeOfLastElongation = ribosome_list[next_avail_ribo].timeOfLastElongation;
    ribosome_list[r_id].aempty_list_pos = ribosome_list[next_avail_ribo].aempty_list_pos;
    ribosome_list[r_id].tl_list_pos = ribosome_list[next_avail_ribo].tl_list_pos;
    ribosome_list[r_id].hit5_list_pos = ribosome_list[next_avail_ribo].hit5_list_pos;
    ribosome_list[r_id].hit3_list_pos = ribosome_list[next_avail_ribo].hit3_list_pos;
    R_grid[ribosome_list[r_id].mRNA][ribosome_list[r_id].pos] = r_id;

    // If the last initiated ribosome had empty A-site,
    // update the ribosome id in the A site empty list accordingly.
    if (ribosome_list[r_id].AsiteEmpty == true) {
      c_id = Genes[mRNA[ribosome_list[r_id].mRNA].gene].seq[ribosome_list[r_id].pos];
      Rb_ae[c_id][ribosome_list[r_id].aempty_list_pos] = r_id;
    } 
    // If the last initiated ribosome was hit from 5' end,
    // update the ribosome id in the 5' hit list accordingly.
    if (ribosome_list[r_id].hitFrom5Prime == true) {
      Rb_hit5[ribosome_list[r_id].hit5_list_pos] = r_id;
    }
    // If the last initiated ribosome was hit from 3' end,
    // update the ribosome id in the 3' hit list accordingly.
    if (ribosome_list[r_id].hitFrom3Prime == true) {
      Rb_hit3[ribosome_list[r_id].hit3_list_pos] = r_id;
    } 
    // If the last initiated ribosome was translocation competent,
    // update the ribosome id in the tl list accordingly.
    if (ribosome_list[r_id].AsiteEmpty == false &&
        ribosome_list[r_id].hitFrom3Prime == false) {
      Rb_ao[ribosome_list[r_id].tl_list_pos] = r_id;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
// MAIN FUNCTION

int main(int argc, char *argv[]) {
////////////////////////////////////////////////////////////////////////////////
// READ INPUT FILES AND PARAMETERS 

  int c1, c2, c3, c4;
  struct timeval time1;
  struct timeval time2;
  FILE *f1, *f2, *f3, *f4, *f5, *f6, *f6a, *f6b, *f7, *f8;
  Ribosome* ribosome_list = new Ribosome[totalNumberOfRibosomes];
  trna* trna_list = new trna[61];
  codon* Codon = new codon[61];

  gettimeofday(&time1, NULL);

  srand(seed);
  // Read in the trna concentrations and aminoacylation rates;
  Read_tRNA_File(trnaConcentrationFile, aminoacylationRateFile, trna_list);

  // Read in the codon trna pairings and their corresponding weights;
  Read_Codon_Files(cognateTrnaFile, pairingStrengthFile, Codon);

  Read_Commandline_Args( argc, argv, trna_list );

  // Number of genes can be modifed from command line
  gene* Genes = new gene[totalNumberOfGenes];
  
  // Read in the numeric seq
  Read_sequence_File(mrnaFile, Genes);

  printf("\nRead input files succcessfully.\n"
      "Total run time = %0.0fs.\n",totalTime);

///////////////////////////////////////////////////////////////////////////////
  // Calculate the total number of mRNAs based on the second col in the fasta file.
  totalNumberOfMrna = 0;
  int geneloop;
  for (geneloop=0;geneloop<totalNumberOfGenes;geneloop++) {
    totalNumberOfMrna += Genes[geneloop].exp;
  }

  // initialize an array of mRNAs
  transcript* mRNA = new transcript[totalNumberOfMrna];
  int mrnaloop=0;
  int exprloop;
  for (geneloop=0;geneloop<totalNumberOfGenes;geneloop++) {
    for (exprloop=0;exprloop<Genes[geneloop].exp;exprloop++) {
      mRNA[mrnaloop].gene = geneloop;
      mRNA[mrnaloop].numberOfInitiations = 0;
      mRNA[mrnaloop].timeOfLastMrnaInitiation = 0.0;
      mRNA[mrnaloop].numberOfFullProteins = 0;
      mRNA[mrnaloop].numberOfAbortedProteins = 0;
      mRNA[mrnaloop].averageInitiationTime = 0.0;
      mRNA[mrnaloop].averageElongationTime = 0.0;
      mRNA[mrnaloop].averageAbortionTime = 0.0;
      mrnaloop++;
    }
  }

  // Variables to track during the simuation
  int deacylatedFreeTrna[61];// Number of uncharged free trna of each type.
  int acylatedFreeTrna[61];// Number of chaged free tRNAs of each type.
  int ribosomeBoundTrna[61];// Number of ribosome bound tRNAs of each type.
  int totalTrna[61];// Total number of tRNAs of type.
  int Rf = totalNumberOfRibosomes;// Total number of free ribosomes
  // (at start all ribosomes are free)
  int Rb[61];// Number of ribosomes bound to each codon
  int n_mrna_free[totalNumberOfGenes];// Number of initiable mRNAs of each gene
  // (first 'ribosomeFootprintSize' codons unbound by ribosomes)
  int total_n_mrna_free;// Total number of initiable mRNAs
  int x;
  // Number of bound ribosomes at each codon that empty A site
  int n_Rb_ae[61];
  // Number of bound ribosomes that can translocate
  int n_Rb_ao = 0;
  // Number of bound ribosomes that have been 5' or 3' hit by a following ribosome
  int n_Rb_hit5 = 0;
  int n_Rb_hit3 = 0;
  // The number of times 5' hit and 3' hit ribosomes are formed
  int numberOfHit5[totalNumberOfMrna]; 
  int numberOfHit3[totalNumberOfMrna]; 
  int r_id;// Current ribosome id.
  // id of the ribosome that is just 5' or 3' to the current ribosome.
  int threePrimeRibosomeId;
  int fivePrimeRibosomeId;
  // location of the ribosome that is just 5' or 3' to the current ribosome.
  int fivePrimeRibosomePosition;
  int threePrimeRibosomePosition;
  // variables for storing position of current ribosome in various lists.
  int pos_in_aempty_list, pos_in_tl_list, pos_in_hit5_list, pos_in_hit3_list;
  int m_id;// Current mRNA id.
  int c_id;// Current codon id.
  int cognate_trna_id = 100;// Current cognate tRNA id.
  int c2_id;
  int g_id;
  int pos;
  double coin, prob_g;
  int obs_max_len = 0;// Observed max gene length
  int obs_max_exp = 0;// Observed max gene expression
  int next_avail_ribo = 0;
  // keep track of whether the ribosome is terminating from end of mRNA
  bool term_flag = false;
  // check if the ribosome is prematurely terminating.
  bool abort_flag = false;
  int numberOfElongations[61];// Number of times a codon type is elongated
  int numberOfAbortions[61];// Number of times a codon type is elongated
  double summedElongationTimes[61];
  double initiationScaledMrnaCopyNumber[totalNumberOfGenes];
  double totalInitiationScaledMrnaCopyNumber = 0.0;
  double totalInitiationRate,
         totalTranslocationRate,
         totalTrnaAccommodationRate,
         totalBackgroundPrematureTerminationRate,
         totalSelectivePrematureTerminationRate,
         totalPrematureTerminationRateUpon5PrimeHit,
         totalPrematureTerminationRateUpon3PrimeHit,
         totalAminoacylationRate;
  double trnaAccommodationRate[61];
  double avgAcylatedTrnaConcn[61],
         avgDeacylatedTrnaConcn[61],
         avgRibosomeBoundTrnaConcn[61];// Average number of free tRNAs of each type
  // (averaged by time)
  double avg_Rf = 0.0;// Average number of free ribosomes
  double avg_n_Rb_ae = 0.0;// Average number of A-empty ribosomes
  double avg_n_Rb_ao = 0.0;// Average number of Translocation competent ribosomes
  double avg_n_Rb_hit3 = 0.0;// Average number of 3' hit ribosomes
  double avg_n_Rb_hit5 = 0.0;// Average number of 5' hit ribosomes
  // (averaged by time)
  int tracking_total_n_Rb_e = 0,
      tracking_Rf = 0,
      tracking_n_Rb_hit5 = 0,
      tracking_n_Rb_hit3 = 0;
  double trackingISMCN = 0.0;
  double currentTime = 0.0;
  double previousRate, currentRate;
  double totalRate = 0.0;
  double randomRate;
  double randomAbortionRate = 0.0;
  int t_print = floor(thresholdTime);
  char tmp_t[10];
  int numberOfIterations = 0;
  double singleTrnaAccommodationRate[61];
  int cognatetrnaloop;
  int temp_trna_total, temp_trna;
  // Bound ribosomes that can translocate
  int* Rb_ao = new int[totalNumberOfRibosomes];
  // Bound ribosomes to each codon that have A site empty 
  int** Rb_ae = new int*[62];
  // Bound ribosomes that have been hit from 5' side.
  int* Rb_hit5 = new int[totalNumberOfRibosomes];
  // Bound ribosomes that have been hit from 3' side.
  int* Rb_hit3 = new int[totalNumberOfRibosomes];
  for (c1=0;c1<62;c1++) {
    Rb_ae[c1] = new int[totalNumberOfRibosomes];
  }

  // Scale initiation rate by mRNA copy number to find
  // total initiation rate easily in each cycle.
  for (c1=0;c1<totalNumberOfGenes;c1++) {
    if (Genes[c1].exp > obs_max_exp) {
      obs_max_exp = Genes[c1].exp;
    }
    if (Genes[c1].length > obs_max_len) {
      obs_max_len = Genes[c1].length;
    }
    n_mrna_free[c1] = Genes[c1].exp;
    total_n_mrna_free += n_mrna_free[c1];
    initiationScaledMrnaCopyNumber[c1] = n_mrna_free[c1]*Genes[c1].ini_prob;
    totalInitiationScaledMrnaCopyNumber += initiationScaledMrnaCopyNumber[c1];
    numberOfHit5[c1] = 0;
    numberOfHit3[c1] = 0;
  }

  // List to figure out which mRNAs can be initiated
  // based on no bound ribosomes from pos=0 to pos>=ribosomeFootprintSize
  int **mrna_initiable = new int*[totalNumberOfGenes];
  for (c1=0;c1<totalNumberOfGenes;c1++) {
    mrna_initiable[c1] = new int[obs_max_exp];
  }

  // Initialize R_grid and avgRibosomeDensityGrid.
  // R_grid now contains the id of ribosome at each mRNA position.
  // If there is no ribosome then that position get the value totalNumberOfRibosomes instead of 0 as 0 is a ribosome id
  // The state of the system with respect to mRNAs and bound ribosomes
  // avgRibosomeDensityGrid keeps track of the average ribosome density at each position of each mRNA after threshold time has elapsed.
  int** R_grid = new int*[totalNumberOfMrna];
  int** numberOfElongationsGrid = new int*[totalNumberOfMrna];
  double** avgRibosomeDensityGrid = new double*[totalNumberOfMrna];
  for (c1=0;c1<totalNumberOfMrna;c1++) {
    R_grid[c1] = new int[obs_max_len];
    numberOfElongationsGrid[c1] = new int[obs_max_len];
    avgRibosomeDensityGrid[c1] = new double[obs_max_len];
    for (c2=0;c2<Genes[mRNA[c1].gene].length;c2++) {
      R_grid[c1][c2] = totalNumberOfRibosomes;
      avgRibosomeDensityGrid[c1][c2] = 0.0;
      numberOfElongationsGrid[c1][c2] = 0;
    }
  }

  // Initialize free mRNAs
  c3=0;
  for (c2=0;c2<totalNumberOfGenes;c2++) {
    for (c1=0;c1<Genes[c2].exp;c1++) {
      mrna_initiable[c2][c1]=c3;// At the begininning all mRNAs can be initiated (store their ids here)
      c3++;
    }
  }

  // Initialize elongation time vectors and number of A site empty ribosomes to each codon.
  int codonloop;
  for (codonloop=0;codonloop<61;codonloop++) {
    summedElongationTimes[codonloop] = 0.0;
    numberOfElongations[codonloop] = 0;
    numberOfAbortions[codonloop] = 0;
    n_Rb_ae[codonloop] = 0;
  }

  // Scale tRNA concentration to match the total number of tRNAs at that growth rate.
  // Scale tRNA aminoacylation rate and ribosome associate rate to change units.
  int trnaloop;
  double temp_trnatotal = 0;
  for (trnaloop=0;trnaloop<61;trnaloop++) {
    temp_trnatotal += trna_list[trnaloop].concn;
  }
  for (trnaloop=0;trnaloop<61;trnaloop++) {
    trna_list[trnaloop].concn *= (totalNumberOfTrna / temp_trnatotal);
    trna_list[trnaloop].aminoacylationRateConstant = trna_list[trnaloop].aminoacylationRateConstant / (cellVolume * 1000 * AVOGADRO_NUMBER);// converting from (moles per litre)-1s-1 to (molecules per cell)-1s-1 units.
    trna_list[trnaloop].trnaRibosomeRateConstant = ribosomeTrnaKcatByKm / (cellVolume * 1000 * AVOGADRO_NUMBER);// converting from (moles per litre)-1s-1 to (molecules per cell)-1s-1 units.
    // Initialize the pools of deacylated, acylated and ribosome bound tRNAs.
    // All trnas are initially free and deacylated.
    deacylatedFreeTrna[trnaloop] = trna_list[trnaloop].concn;
    acylatedFreeTrna[trnaloop] = 0;
    ribosomeBoundTrna[trnaloop] = 0;
    totalTrna[trnaloop] = deacylatedFreeTrna[trnaloop] + acylatedFreeTrna[trnaloop] + ribosomeBoundTrna[trnaloop];
    avgAcylatedTrnaConcn[trnaloop] = 0.0;
    avgDeacylatedTrnaConcn[trnaloop] = 0.0;
    avgRibosomeBoundTrnaConcn[trnaloop] = 0.0;
  }

///////////////////////////////////////////////////////////////////////////////
// SIMULATION LOOP
// Key Steps:
// * Calculate rates of all reaction channels using state matrix
// * Pick a reaction channel
// * Update states of ribosomes and all other state variables 
// (tRNA concentrations, mRNA states, lists of ribosomes in different states)

  // Total Aminoacylation Rate calculated  initially and updated whenever deacylatedFreeTrna changes.
  totalAminoacylationRate = 0.0;
  for (trnaloop=0;trnaloop<61;trnaloop++) {
    totalAminoacylationRate += ( trna_list[trnaloop].aminoacylationRateConstant * deacylatedFreeTrna[trnaloop] );
  }

  // Till current time is less than max time
  while(currentTime<totalTime)
  {
    // RATE CALCULATION
    // Calculate the total rate of Gillespie simulation at this step as the 
    // sum of the rates of initiation, elongation, aminoacylation and abortion.
    totalRate = 0.0;
    totalRate += totalAminoacylationRate;

    // Total Initiation Rate
    totalInitiationRate = totalInitiationScaledMrnaCopyNumber*Rf/( (1-boundRibosomeFraction) * totalNumberOfRibosomes );
    totalRate += totalInitiationRate;

    // Total trna Accommodation Rate
    totalTrnaAccommodationRate = 0.0;
    for (codonloop=0;codonloop<61;codonloop++) {
      singleTrnaAccommodationRate[codonloop] = 0.0;
      for (cognatetrnaloop=0; cognatetrnaloop<2; cognatetrnaloop++) {
        singleTrnaAccommodationRate[codonloop] += ( trna_list[ Codon[codonloop].trna[cognatetrnaloop] ].trnaRibosomeRateConstant * acylatedFreeTrna[ Codon[codonloop].trna[cognatetrnaloop] ] * Codon[codonloop].weight[cognatetrnaloop] );
      }
  trnaAccommodationRate[codonloop] = singleTrnaAccommodationRate[codonloop] * n_Rb_ae[codonloop];
      totalTrnaAccommodationRate += trnaAccommodationRate[codonloop];
    }
    totalRate += totalTrnaAccommodationRate;

    // Total Translocation Rate
    totalTranslocationRate = n_Rb_ao / codonMinElongationTime;
    totalRate += totalTranslocationRate;

    // Total Background Premature Termination Rate (Asite empty ribosomes only)
    totalBackgroundPrematureTerminationRate = 0.0;
    for (codonloop=0;codonloop<61;codonloop++) {
      totalBackgroundPrematureTerminationRate += backgroundPrematureTerminationRate * n_Rb_ae[codonloop];
    }
    totalRate += totalBackgroundPrematureTerminationRate;

    // Total Selective Premature Termination Rate (Asite empty ribosomes only)
    totalSelectivePrematureTerminationRate = 0.0;
    for (codonloop=0;codonloop<61;codonloop++) {
      if (singleTrnaAccommodationRate[codonloop] < thresholdTrnaAccommodationRate) {
        totalSelectivePrematureTerminationRate += selectivePrematureTerminationRate * n_Rb_ae[codonloop];
      } else {
        continue;
      }
    }
    totalRate += totalSelectivePrematureTerminationRate;

    // Premature Termination Rate For Ribosomes that have been hit by another
    // Ribosome from the 5' side and 3' side
    totalPrematureTerminationRateUpon5PrimeHit = prematureTerminationRateUpon5PrimeHit * n_Rb_hit5;
    totalRate += totalPrematureTerminationRateUpon5PrimeHit;

    totalPrematureTerminationRateUpon3PrimeHit = prematureTerminationRateUpon3PrimeHit * n_Rb_hit3;
    totalRate += totalPrematureTerminationRateUpon3PrimeHit;

    // Increment time. Following Shah et al, we have not implemented the strict
    //  Gillespie prescription of random time drawn from an exponential distribution. 
    //  Since all calculated quantities are steady-state values,
    //  this difference does not matter.
    currentTime +=  ( 1 / totalRate );

///////////////////////////////////////////////////////////////////////////////

    // Tracking time-averaged free ribosomes and tRNAs for final output.
    if (currentTime>thresholdTime) {
      if (printOptions[3]==1) {
        for (trnaloop=0;trnaloop<61;trnaloop++) {
          avgAcylatedTrnaConcn[trnaloop] += (double)acylatedFreeTrna[trnaloop] / totalRate;
          avgDeacylatedTrnaConcn[trnaloop] += (double)deacylatedFreeTrna[trnaloop] / totalRate;
          avgRibosomeBoundTrnaConcn[trnaloop] += (double)ribosomeBoundTrna[trnaloop] / totalRate;
        }
        avg_Rf += (double)Rf / totalRate;
        avg_n_Rb_hit3 += (double)n_Rb_hit3 / totalRate;
        avg_n_Rb_hit5 += (double)n_Rb_hit5 / totalRate;
        for (codonloop=0;codonloop<61;codonloop++) {
          avg_n_Rb_ae += (double)n_Rb_ae[codonloop] / totalRate;
        }
        avg_n_Rb_ao += (double)n_Rb_ao / totalRate;
      }
    }

    // Draw a random reaction channel by choosing between aminoacylation, initiation, trna accommodation, translocation and various termination reactions and then picking the sub channel according to the respective rate.

    // Find a random rate that is used for comparison below
    randomRate = totalRate * rand()/((double)RAND_MAX);
    previousRate = 0.0;
    currentRate = 0.0;

    // Aminoacylation
    previousRate = currentRate;
    currentRate += totalAminoacylationRate;
    // Pick tRNA that is aminoacylated with a probability proportional to its aminoacylation rate.
    if (randomRate > previousRate && randomRate < currentRate) {
      for (trnaloop=0;trnaloop<61;trnaloop++) {
        previousRate += ( trna_list[trnaloop].aminoacylationRateConstant * deacylatedFreeTrna[trnaloop] );
        if (randomRate < previousRate) {
          deacylatedFreeTrna[trnaloop]--;// Decrease deacylated trna by one.
          // Update aminoacylation rate.
          totalAminoacylationRate -= trna_list[trnaloop].aminoacylationRateConstant;
          acylatedFreeTrna[trnaloop]++;// Increase acylated trna by one.
          break;
        }
      }
    }

    // Translation initiation
    previousRate = currentRate;
    currentRate += totalInitiationRate;

    if (randomRate > previousRate && randomRate < currentRate) {
      /// Choose various reaction channels: ///
      // - Gene based on initiation probability
      // - Random mRNA once the gene is picked
      // - First codon on mRNA since this is initiation
      // - Free Ribosome sequentially from the pool of free ribosomes

      // Pick the gene for mRNA based on their initiation rates.
      for (geneloop=0;geneloop<totalNumberOfGenes;geneloop++) {
        previousRate += initiationScaledMrnaCopyNumber[geneloop] *Rf / ( (1-boundRibosomeFraction) * totalNumberOfRibosomes );
        if (randomRate < previousRate) {
          // Once a gene is selected pick a random mRNA
          c2 = rand()%n_mrna_free[geneloop];
          m_id = mrna_initiable[geneloop][c2];
          // Update free mRNA list for subsequent rounds of initiation
          n_mrna_free[geneloop]--;
          total_n_mrna_free--;
          initiationScaledMrnaCopyNumber[geneloop] -= Genes[geneloop].ini_prob;
          totalInitiationScaledMrnaCopyNumber -= Genes[geneloop].ini_prob;
          if (c2!=n_mrna_free[geneloop]) {
            mrna_initiable[geneloop][c2] = mrna_initiable[geneloop][n_mrna_free[geneloop]];
          }
          break;
        }
      }

      // The start codon is at the ribosome P-site upon initiation
      // Hence asite codon is 1 and not 0 below.
      c_id = Genes[mRNA[m_id].gene].seq[1];

      // Ribosomes are picked sequentially starting at r_id = 0
      r_id = next_avail_ribo;
      next_avail_ribo++;

      // Update the state of the initiating ribosome.
      ribosome_list[r_id].mRNA = m_id;
      ribosome_list[r_id].pos = 1;
      ribosome_list[r_id].hitFrom5Prime = false;
      ribosome_list[r_id].hitFrom3Prime = false;
      ribosome_list[r_id].timeOfLastInitiation = currentTime;
      ribosome_list[r_id].timeOfLastElongation = currentTime;

      // Update number of free ribosomes
      Rf--;

      // Update number and list of ribosomes with empty A site
      ribosome_list[r_id].AsiteEmpty = true;
      ribosome_list[r_id].aempty_list_pos = n_Rb_ae[c_id];
      Rb_ae[c_id][n_Rb_ae[c_id]] = r_id;
      n_Rb_ae[c_id]++;

      // Update the ribosome grid upon initiation
      R_grid[m_id][ribosome_list[r_id].pos] = r_id;

      // Update variables that are tracked during the course of the simulation.
      if (currentTime>thresholdTime) {
        mRNA[m_id].numberOfInitiations++;
        // The average initation time is finally calculated at the end of the
        // simulation by dividing the number below by the number of initiations.
        mRNA[m_id].averageInitiationTime += currentTime-mRNA[m_id].timeOfLastMrnaInitiation;
      }
      // Update mRNA and mRNA grid
      mRNA[m_id].timeOfLastMrnaInitiation = currentTime;
    }

    // tRNA Accommodation
    previousRate = currentRate;
    currentRate += totalTrnaAccommodationRate;

    if (randomRate > previousRate && randomRate < currentRate) {

      /// Choose various reaction channels: ///
      // - Codon, Ribosome, mRNA, tRNA

      // Pick a random codon based on trna accommodation rate
      for (codonloop=0;codonloop<61;codonloop++) {
        previousRate += trnaAccommodationRate[codonloop];
        if (randomRate < previousRate) {
          c_id = codonloop;
          break;
        }
      }

      // Randomly pick an A-site empty ribosome bound to codon c_id
      pos_in_aempty_list = rand()%n_Rb_ae[c_id];
      r_id = Rb_ae[c_id][pos_in_aempty_list];
      // The mRNA is whichever mRNA the picked ribosome was bound to.
      m_id = ribosome_list[r_id].mRNA;

      // Pick a random A-site cognate tRNA based on respective elongation probabilities.
      coin = rand()/((double)RAND_MAX);
      prob_g = 0.0;

      for (cognatetrnaloop=0;cognatetrnaloop<2;cognatetrnaloop++) {
        prob_g += ( trna_list[ Codon[c_id].trna[cognatetrnaloop] ].trnaRibosomeRateConstant
            * acylatedFreeTrna[ Codon[c_id].trna[cognatetrnaloop] ]
		    * Codon[c_id].weight[cognatetrnaloop] ) / singleTrnaAccommodationRate[c_id];
        if (coin<prob_g) {
          cognate_trna_id = Codon[c_id].trna[cognatetrnaloop];
          break;
        }
      }

      // Ribosome mRNA ID, position, codon ID, P site tRNA are unchanged
      // Set the A-site tRNA to be the one picked above. 
      // matters only for A-site empty ribosomes
      ribosome_list[r_id].AtRNA = cognate_trna_id;
      ribosome_list[r_id].AsiteEmpty = false;

      // Update list of Asite empty ribosomes
      n_Rb_ae[c_id]--;
      if (pos_in_aempty_list != n_Rb_ae[c_id]) {
        Rb_ae[c_id][pos_in_aempty_list] = Rb_ae[c_id][n_Rb_ae[c_id]];
        ribosome_list[Rb_ae[c_id][n_Rb_ae[c_id]]].aempty_list_pos = pos_in_aempty_list;
      }


      // Update list of 5' hit ribosomes if necessary
      // Upon A-site accommodation, ribosome is assumed to go out of the
      // hit from 5' state.
      if (ribosome_list[r_id].hitFrom5Prime == true) {
        ribosome_list[r_id].hitFrom5Prime = false;
        pos_in_hit5_list = ribosome_list[r_id].hit5_list_pos;
        n_Rb_hit5--;
        Rb_hit5[pos_in_hit5_list] = Rb_hit5[n_Rb_hit5];
        ribosome_list[Rb_hit5[n_Rb_hit5]].hit5_list_pos = pos_in_hit5_list;
      }

      // If there is an immediate 3' ribosome, the accommodated ribosome goes
      // into a 3' hit state, or else the accommodated ribosome goes
      // into a translocation competent state.
      threePrimeRibosomePosition = ribosome_list[r_id].pos + ribosomeFootprintSize;
      threePrimeRibosomeId = totalNumberOfRibosomes;
      if (threePrimeRibosomePosition <= Genes[mRNA[m_id].gene].length - 1) {
        threePrimeRibosomeId = R_grid[m_id][threePrimeRibosomePosition];
      }
      if (threePrimeRibosomeId < totalNumberOfRibosomes) {
        ribosome_list[r_id].hitFrom3Prime = true;
        ribosome_list[r_id].hit3_list_pos = n_Rb_hit3;
        Rb_hit3[n_Rb_hit3] = r_id;
        n_Rb_hit3++;
        numberOfHit3[m_id]++;
      } else {
        Rb_ao[n_Rb_ao] = r_id;
        ribosome_list[r_id].tl_list_pos = n_Rb_ao;
        n_Rb_ao++;
      }

      // If there is an immediate 3' ribosome and if the A-site of the 
      // downstream ribosome 
      // is also empty, the downstream ribosome goes into a 5' hit state.
      if (threePrimeRibosomeId < totalNumberOfRibosomes &&
          ribosome_list[threePrimeRibosomeId].AsiteEmpty == true) {
        ribosome_list[threePrimeRibosomeId].hitFrom5Prime = true;
        ribosome_list[threePrimeRibosomeId].hit5_list_pos = n_Rb_hit5;
        Rb_hit5[n_Rb_hit5] = threePrimeRibosomeId;
        n_Rb_hit5++;
        numberOfHit5[m_id]++;
      }
  
      // Decrease the number of free acylated tRNA and increase the number of 
      // ribosome bound tRNA by 1.
      acylatedFreeTrna[ribosome_list[r_id].AtRNA]--;
      ribosomeBoundTrna[ribosome_list[r_id].AtRNA]++;
    }

    // Ribosome translocation
    previousRate = currentRate;
    currentRate += totalTranslocationRate;

    if (randomRate > previousRate && randomRate < currentRate) {
      // Randomly pick a ribosome that is translocation competent
      pos_in_tl_list = rand()%n_Rb_ao;
      r_id = Rb_ao[pos_in_tl_list];
      m_id = ribosome_list[r_id].mRNA;
      c_id = Genes[mRNA[m_id].gene].seq[ribosome_list[r_id].pos];
      c2_id = Genes[mRNA[m_id].gene].seq[ribosome_list[r_id].pos+1];

      // Release P-site tRNA
      // We do not track the initiator tRNA at the start codon
      if (ribosome_list[r_id].pos > 1) {
        ribosomeBoundTrna[ribosome_list[r_id].PtRNA]--;
        deacylatedFreeTrna[ribosome_list[r_id].PtRNA]++;
        // Update aminoacylation rate.
        totalAminoacylationRate += trna_list[ribosome_list[r_id].PtRNA].aminoacylationRateConstant;
      }

      // Move the A-site tRNA to the P-site tRNA
      ribosome_list[r_id].PtRNA = ribosome_list[r_id].AtRNA;
      
      // Update the corresponding mRNA positions
      R_grid[m_id][ribosome_list[r_id].pos] = totalNumberOfRibosomes;
      R_grid[m_id][ribosome_list[r_id].pos + 1] = r_id;
  
      // Update number and list of ribosomes with empty A site
      if (ribosome_list[r_id].AsiteEmpty == false) {
        ribosome_list[r_id].AsiteEmpty = true;
        ribosome_list[r_id].aempty_list_pos = n_Rb_ae[c2_id];
        Rb_ae[c2_id][n_Rb_ae[c2_id]] = r_id;
        n_Rb_ae[c2_id]++;
      }
  
      // Update number and list of translocation competent ribosomes
      n_Rb_ao--;
      Rb_ao[pos_in_tl_list] = Rb_ao[n_Rb_ao];
      ribosome_list[Rb_ao[n_Rb_ao]].tl_list_pos = pos_in_tl_list;
  
      // If there is an immediate 5' ribosome that is hit from 3',
      // it is not hit from 3' anymore and is translocation competent.
      fivePrimeRibosomePosition = ribosome_list[r_id].pos - ribosomeFootprintSize;
      fivePrimeRibosomeId = totalNumberOfRibosomes;
      if (fivePrimeRibosomePosition > 0) {
        fivePrimeRibosomeId = R_grid[m_id][fivePrimeRibosomePosition];
      }
      if (fivePrimeRibosomeId < totalNumberOfRibosomes &&
          ribosome_list[fivePrimeRibosomeId].hitFrom3Prime == true) {
        pos_in_hit3_list = ribosome_list[fivePrimeRibosomeId].hit3_list_pos;
        ribosome_list[fivePrimeRibosomeId].hitFrom3Prime = false;
        n_Rb_hit3--;
        Rb_hit3[pos_in_hit3_list] = Rb_hit3[n_Rb_hit3];
        ribosome_list[Rb_hit3[n_Rb_hit3]].hit3_list_pos = pos_in_hit3_list;
  
        ribosome_list[fivePrimeRibosomeId].tl_list_pos = n_Rb_ao;
        Rb_ao[n_Rb_ao] = fivePrimeRibosomeId;
        n_Rb_ao++;
      }
  
      // The mRNA is now initiation-competent if a 'ribosomefootprintsize'
      // region got empty 
      if (ribosome_list[r_id].pos == ribosomeFootprintSize) {
        // find gene id
        g_id = mRNA[m_id].gene;
        // add the mrna to the list of free mrnas for that gene
        mrna_initiable[g_id][n_mrna_free[g_id]] = m_id;
        // increase the number of free mrnas for that gene
        n_mrna_free[g_id]++;
        // increase total number of free mrnas
        total_n_mrna_free++;
        // increase the initiation scaled mrna copy number for that gene
        initiationScaledMrnaCopyNumber[g_id] += Genes[g_id].ini_prob;
        // increase the total initiation scaled mrna copy number
        totalInitiationScaledMrnaCopyNumber += Genes[g_id].ini_prob;
      }
  
      // Update variables that are tracked during the simulation
      if (currentTime>thresholdTime) {
        // For estimation of avg elongation times of codons
        summedElongationTimes[c_id] += currentTime-ribosome_list[r_id].timeOfLastElongation;
        // Keep track of the ribosome density at the corresponding codon position.
        avgRibosomeDensityGrid[m_id][ribosome_list[r_id].pos] += currentTime-ribosome_list[r_id].timeOfLastElongation;
        // Keep track of the number of times that this position on this mRNA was elongated / aborted.
        numberOfElongationsGrid[m_id][ribosome_list[r_id].pos] ++;
        numberOfElongations[c_id]++;
      }

      // Update ribosome parameters
      ribosome_list[r_id].pos++;
      ribosome_list[r_id].timeOfLastElongation = currentTime;

      // If ribosome reached the end of the message, terminate translation.
      if (ribosome_list[r_id].pos == Genes[mRNA[m_id].gene].length -1) {
        TerminateRibosome(
          r_id,
          m_id,
          c2_id,
          ribosome_list,
          trna_list,
          mRNA,
          Genes,
          ribosomeBoundTrna, deacylatedFreeTrna,
          totalAminoacylationRate,
          R_grid,
          mrna_initiable,
          Rf,
          pos_in_aempty_list,
          pos_in_tl_list,
          pos_in_hit5_list,
          pos_in_hit3_list,
          n_Rb_ae, n_Rb_ao, n_Rb_hit5, n_Rb_hit3,
          Rb_ae, Rb_ao, Rb_hit5, Rb_hit3,
          total_n_mrna_free,
          n_mrna_free,
          initiationScaledMrnaCopyNumber,
          totalInitiationScaledMrnaCopyNumber,
          currentTime,
          summedElongationTimes,
          avgRibosomeDensityGrid,
          numberOfElongationsGrid,
          numberOfElongations,
          numberOfAbortions,
          next_avail_ribo
      );
      }
    }

    // Premature termination upon collision from the 5' side
    previousRate = currentRate;
    currentRate += totalPrematureTerminationRateUpon5PrimeHit;
    if (randomRate > previousRate && randomRate < currentRate) {
      // Randomly pick a ribosome that has been hit from the 5' side
      pos_in_hit5_list = rand()%n_Rb_hit5;
      r_id = Rb_hit5[pos_in_hit5_list];
      m_id = ribosome_list[r_id].mRNA;
      c_id = Genes[mRNA[m_id].gene].seq[ribosome_list[r_id].pos];
      TerminateRibosome(
        r_id,
        m_id,
        c_id,
        ribosome_list,
        trna_list,
        mRNA,
        Genes,
        ribosomeBoundTrna, deacylatedFreeTrna,
        totalAminoacylationRate,
        R_grid,
        mrna_initiable,
        Rf,
        pos_in_aempty_list,
        pos_in_tl_list,
        pos_in_hit5_list,
        pos_in_hit3_list,
        n_Rb_ae, n_Rb_ao, n_Rb_hit5, n_Rb_hit3,
        Rb_ae, Rb_ao, Rb_hit5, Rb_hit3,
        total_n_mrna_free,
        n_mrna_free,
        initiationScaledMrnaCopyNumber,
        totalInitiationScaledMrnaCopyNumber,
        currentTime,
        summedElongationTimes,
        avgRibosomeDensityGrid,
        numberOfElongationsGrid,
        numberOfElongations,
        numberOfAbortions,
        next_avail_ribo
    );
    }

    // Premature termination upon collision from the 3' side
    previousRate = currentRate;
    currentRate += totalPrematureTerminationRateUpon3PrimeHit;
    if (randomRate > previousRate && randomRate < currentRate) {
      // Randomly pick a ribosome that has been hit from the 3' side
      pos_in_hit3_list = rand()%n_Rb_hit3;
      r_id = Rb_hit3[pos_in_hit3_list];
      m_id = ribosome_list[r_id].mRNA;
      c_id = Genes[mRNA[m_id].gene].seq[ribosome_list[r_id].pos];
      TerminateRibosome(
        r_id,
        m_id,
        c_id,
        ribosome_list,
        trna_list,
        mRNA,
        Genes,
        ribosomeBoundTrna, deacylatedFreeTrna,
        totalAminoacylationRate,
        R_grid,
        mrna_initiable,
        Rf,
        pos_in_aempty_list,
        pos_in_tl_list,
        pos_in_hit5_list,
        pos_in_hit3_list,
        n_Rb_ae, n_Rb_ao, n_Rb_hit5, n_Rb_hit3,
        Rb_ae, Rb_ao, Rb_hit5, Rb_hit3,
        total_n_mrna_free,
        n_mrna_free,
        initiationScaledMrnaCopyNumber,
        totalInitiationScaledMrnaCopyNumber,
        currentTime,
        summedElongationTimes,
        avgRibosomeDensityGrid,
        numberOfElongationsGrid,
        numberOfElongations,
        numberOfAbortions,
        next_avail_ribo
    );
    }

    // Premature background termination of A-site empty ribosomes
    previousRate = currentRate;
    currentRate += totalBackgroundPrematureTerminationRate;
    if (randomRate > previousRate && randomRate < currentRate) {
      // Pick a random codon based on background premature  terminatino rate
      for (codonloop=0;codonloop<61;codonloop++) {
        previousRate += backgroundPrematureTerminationRate * n_Rb_ae[codonloop];
        if (randomRate < previousRate) {
          c_id = codonloop;
          break;
        }
      }
      // Randomly pick an A-site empty ribosome bound to codon c_id
      pos_in_aempty_list = rand()%n_Rb_ae[c_id];
      r_id = Rb_ae[c_id][pos_in_aempty_list];
      m_id = ribosome_list[r_id].mRNA;
      TerminateRibosome(
        r_id,
        m_id,
        c_id,
        ribosome_list,
        trna_list,
        mRNA,
        Genes,
        ribosomeBoundTrna, deacylatedFreeTrna,
        totalAminoacylationRate,
        R_grid,
        mrna_initiable,
        Rf,
        pos_in_aempty_list,
        pos_in_tl_list,
        pos_in_hit5_list,
        pos_in_hit3_list,
        n_Rb_ae, n_Rb_ao, n_Rb_hit5, n_Rb_hit3,
        Rb_ae, Rb_ao, Rb_hit5, Rb_hit3,
        total_n_mrna_free,
        n_mrna_free,
        initiationScaledMrnaCopyNumber,
        totalInitiationScaledMrnaCopyNumber,
        currentTime,
        summedElongationTimes,
        avgRibosomeDensityGrid,
        numberOfElongationsGrid,
        numberOfElongations,
        numberOfAbortions,
        next_avail_ribo
      );
    }

    // Premature selective termination of A-site empty ribosomes
    previousRate = currentRate;
    currentRate += totalSelectivePrematureTerminationRate;
    if (randomRate > previousRate && randomRate < currentRate) {
      // Pick a codon based on selective premature  termination rate
      for (codonloop=0;codonloop<61;codonloop++) {
        if (singleTrnaAccommodationRate[codonloop] < thresholdTrnaAccommodationRate) {
          previousRate += selectivePrematureTerminationRate * 
                            n_Rb_ae[codonloop];
          if (randomRate < previousRate) {
            c_id = codonloop;
            break;
          }
        } else {
          continue;
        }
      }
      // Randomly pick an A-site empty ribosome bound to codon c_id
      pos_in_aempty_list = rand()%n_Rb_ae[c_id];
      r_id = Rb_ae[c_id][pos_in_aempty_list];
      m_id = ribosome_list[r_id].mRNA;
      TerminateRibosome(
        r_id,
        m_id,
        c_id,
        ribosome_list,
        trna_list,
        mRNA,
        Genes,
        ribosomeBoundTrna, deacylatedFreeTrna,
        totalAminoacylationRate,
        R_grid,
        mrna_initiable,
        Rf,
        pos_in_aempty_list,
        pos_in_tl_list,
        pos_in_hit5_list,
        pos_in_hit3_list,
        n_Rb_ae, n_Rb_ao, n_Rb_hit5, n_Rb_hit3,
        Rb_ae, Rb_ao, Rb_hit5, Rb_hit3,
        total_n_mrna_free,
        n_mrna_free,
        initiationScaledMrnaCopyNumber,
        totalInitiationScaledMrnaCopyNumber,
        currentTime,
        summedElongationTimes,
        avgRibosomeDensityGrid,
        numberOfElongationsGrid,
        numberOfElongations,
        numberOfAbortions,
        next_avail_ribo
      );
    }

    // Print state regularly during simulation
    int total_n_Rb_ae = 0;
    for (int codon=0;codon<62;codon++) {
      total_n_Rb_ae += n_Rb_ae[codon];
    }

    if (
        //tracking_total_n_Rb_e + tracking_n_Rb_hit3 != total_n_Rb_e + n_Rb_hit3 ||
        numberOfIterations%10000000 == 0 ||
        //tracking_n_Rb_hit5 != n_Rb_hit5 ||
        //tracking_n_Rb_hit3 != n_Rb_hit3
        //m_id != -1
        // mRNA[m_id].gene != -1 // ||
        //initiationScaledMrnaCopyNumber[2] != trackingISMCN
        false
       ) {
      usleep(200000);
      for (int pos=0; pos<Genes[mRNA[m_id].gene].length; pos++) {
        if (R_grid[m_id][pos] != totalNumberOfRibosomes) {
          if (ribosome_list[R_grid[m_id][pos]].AsiteEmpty == true) {
            printf("%d ", pos);
          } else {
            printf("%d* ", pos);
          }
        }
      }
      printf("\n");
      printf("iterations = %d\n"
          "time = %0.2fs\n"
          "totalRate = %0.2fs-1\n"
          "totalInitiationRate = %0.2fs-1\n"
          "totalTrnaAccommodationRate = %0.2fs-1\n"
          "totalTranslocationRate = %0.2fs-1\n"
          "totalAminoacylationRate = %0.2fs-1\n"
          "totalPrematureTerminationRateUpon5PrimeHit = %0.2fs-1\n"
          "totalPrematureTerminationRateUpon3PrimeHit = %0.2fs-1\n"
          "numberOfAsiteEmptyRibosomes= %d\n"
          "numberOfTranslocationCompetentRibosomes= %d\n"
          "numberOf5primeHitRibosomes = %d\n"
          "numberOf3primeHitRibosomes = %d\n"
          "numberOfFreeRibosomes = %d\n"
          "totalNumberOfRibosomes = %d\n"
          "c_id = %d\nm_id = %d\ng_id = %d\nr_id = %d\npos = %d\n"
          "ribosome_list.hitFrom5Prime = %s\n"
          "ribosome_list.hitFrom3Prime = %s\n"
          "ribosome_list.AsiteEmpty = %s\n"
          "TrackingISMCN = %0.5f, ISMCN[2]=%0.5f\n"
          "///////////////////////////////////////////\n",
          numberOfIterations, currentTime, totalRate,
          totalInitiationRate, totalTrnaAccommodationRate, totalTranslocationRate, totalAminoacylationRate,
          totalPrematureTerminationRateUpon5PrimeHit,
          totalPrematureTerminationRateUpon3PrimeHit,
          total_n_Rb_ae, n_Rb_ao, n_Rb_hit5, n_Rb_hit3, Rf, total_n_Rb_ae + n_Rb_ao + n_Rb_hit3 + Rf,
          c_id, m_id, mRNA[m_id].gene, r_id, ribosome_list[r_id].pos,
          ribosome_list[r_id].hitFrom5Prime ? "true" : "false",
          ribosome_list[r_id].hitFrom3Prime ? "true" : "false",
          ribosome_list[r_id].AsiteEmpty ? "true" : "false",
          trackingISMCN, initiationScaledMrnaCopyNumber[2]
          );
      tracking_Rf = Rf;
      tracking_n_Rb_hit5 = n_Rb_hit5;
      tracking_n_Rb_hit3 = n_Rb_hit3;
      trackingISMCN = initiationScaledMrnaCopyNumber[2];
    }
    numberOfIterations++;

  }

  delete [] R_grid;// Free up the memory occupied by this huge matrix.

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Process output for printing

  int n_trans[totalNumberOfGenes];
  int n_aborts[totalNumberOfGenes];
  int n_tot=0;
  int n_abort_total=0;
  int n_tot_all=0;
  double g_etimes[totalNumberOfGenes];
  double g_atimes[totalNumberOfGenes];
  double g_ini[totalNumberOfGenes];
  double dCST=0.0;

  if(printOptions[1]==1 || printOptions[2]==1 || printOptions[3]==1) {
    c3 = 0;
    for (c1=0;c1<totalNumberOfGenes;c1++) {
      n_trans[c1] = 0;
      g_ini[c1]=0.0;
      g_etimes[c1]=0.0;

      for (c2=0;c2<Genes[c1].exp;c2++) {
        n_trans[c1] += mRNA[c3].numberOfFullProteins;
        n_aborts[c1] += mRNA[c3].numberOfAbortedProteins;
        g_ini[c1] += mRNA[c3].averageInitiationTime;
        g_etimes[c1] += mRNA[c3].averageElongationTime;
        g_atimes[c1] += mRNA[c3].averageAbortionTime;
        n_tot += mRNA[c3].numberOfFullProteins;
        n_abort_total += mRNA[c3].numberOfAbortedProteins;
        c3++;
      }
    }
  }

  // Print the output
  // Elongation times of all codons
  if(printOptions[0]==1) {
    strcpy(outfile,outfilePrefix);
    f2 = fopen(strcat(outfile,"/etimes.out"),"w");
    fprintf(f2,"Codon\tNum_of_events\tAvg_elong_time(sec)\tNum_of_abortions\n");
    for (c1=0;c1<61;c1++) {
      dCST = summedElongationTimes[c1]/(double)numberOfElongations[c1];
      fprintf(f2,"%d\t%s\t%d\t%0.4f\t%d\n",c1, Codon[c1].triplet,numberOfElongations[c1],dCST, numberOfAbortions[c1] );
    }
    fclose(f2);
  }

  // Average total elongation times of all genes
  if(printOptions[1]==1) {
    strcpy(outfile,outfilePrefix);
    f3 = fopen(strcat(outfile,"/gene_totetimes.out"),"w");

    fprintf(f3,"Gene\tNum_of_full_events\tAvg_total_elong_time(sec)\tNum_of_aborts\tAvg_total_abort_time(sec)\n");
    for (c1=0;c1<totalNumberOfGenes;c1++) {
      dCST = g_etimes[c1]/(double)n_trans[c1];
      fprintf(f3,"%d\t%d\t%g\t%d\t%g\n",c1,n_trans[c1],dCST, n_aborts[c1], g_atimes[c1]/(double)n_aborts[c1]);
    }
    fclose(f3);
  }

  // Average time between initiation of all genes
  if(printOptions[2]==1) {
    strcpy(outfile,outfilePrefix);
    f4 = fopen(strcat(outfile,"/gene_initimes.out"),"w");

    fprintf(f4,"Gene\tNum_of_events\tAvg_initiation_time(sec)\n");
    for (c1=0;c1<totalNumberOfGenes;c1++) {
      dCST = g_ini[c1]/(double)n_trans[c1];
      fprintf(f4,"%d\t%d\t%g\n",c1,n_trans[c1],dCST);
    }
    fclose(f4);
  }

  // Average number of free ribosomes and tRNAs at equilibrium
  if(printOptions[3]==1) {
    strcpy(outfile,outfilePrefix);
    f5 = fopen(strcat(outfile,"/avg_ribo_tRNA.out"),"w");
    fprintf( f5, "Time-averaged number of ribosomes and trnas in different states\n" );
    avg_Rf = avg_Rf/(currentTime-thresholdTime);
    fprintf(f5,"free_ribosomes\t%0f\n",avg_Rf);
    avg_n_Rb_ao = avg_n_Rb_ao/(currentTime-thresholdTime);
    fprintf(f5,"asiteoccupied_ribosomes\t%0f\n",avg_n_Rb_ao);
    avg_n_Rb_ae = avg_n_Rb_ae/(currentTime-thresholdTime);
    fprintf(f5,"asiteempty_ribosomes\t%0f\n",avg_n_Rb_ae);
    avg_n_Rb_hit5 = avg_n_Rb_hit5/(currentTime-thresholdTime);
    fprintf(f5,"hit5_ribosomes\t%0f\n",avg_n_Rb_hit5);
    avg_n_Rb_hit3 = avg_n_Rb_hit3/(currentTime-thresholdTime);
    fprintf(f5,"hit3_ribosomes\t%0f\n",avg_n_Rb_hit3);
    for (trnaloop=0;trnaloop<61;trnaloop++) {
      if (totalTrna[trnaloop] > 0) {
        avgAcylatedTrnaConcn[trnaloop] = avgAcylatedTrnaConcn[trnaloop]/(currentTime-thresholdTime);
        avgDeacylatedTrnaConcn[trnaloop] = avgDeacylatedTrnaConcn[trnaloop]/(currentTime-thresholdTime);
        avgRibosomeBoundTrnaConcn[trnaloop] = avgRibosomeBoundTrnaConcn[trnaloop]/(currentTime-thresholdTime);
        fprintf(f5,"acylated_trna_%d_%s\t%0f\n", trnaloop, trna_list[trnaloop].anticodon, avgAcylatedTrnaConcn[trnaloop] );
        fprintf(f5,"deacylated_trna_%d_%s\t%0f\n", trnaloop, trna_list[trnaloop].anticodon, avgDeacylatedTrnaConcn[trnaloop] );
        fprintf(f5,"ribosomebound_trna_%d_%s\t%0f\n", trnaloop, trna_list[trnaloop].anticodon, avgRibosomeBoundTrnaConcn[trnaloop] );
      }
    }
    fclose(f5);
  }

  // Time-averaged ribosome density on each gene that was tracked after threshold time was crossed.
  if(printOptions[4]==1) {
    strcpy(outfile,outfilePrefix);
    f6 = fopen(strcat(outfile,"/gene_ribo_density.out"),"w");
    strcpy(outfile,outfilePrefix);
    f6a = fopen(strcat(outfile,"/gene_number_elongs.out"),"w");
    strcpy(outfile,outfilePrefix);
    f6b = fopen(strcat(outfile,"/gene_number_collisions.out"),"w");
    fprintf(f6b,"Gene\tNum_of_5'_hits\tNumber_of_3'_hits\n");

    double** geneSteadyStateRibosomeDensity;// Sum up the ribosome density on all mRNAs for each gene.
    int** geneTotalNumberOfElongations;// Sum up the total number of elongations on all mRNAs for each gene.
    int** geneTotalNumberOfHit5; // Sump up the number of 5' hits on all mRNAs for each gene.
    int** geneTotalNumberOfHit3; // Sump up the number of 3' hits on all mRNAs for each gene.
    geneSteadyStateRibosomeDensity = (double**) malloc(sizeof(double *) * totalNumberOfGenes);
    geneTotalNumberOfElongations = (int**) malloc(sizeof(int *) * totalNumberOfGenes);
    geneTotalNumberOfHit5 = (int**) malloc(sizeof(int *) * totalNumberOfGenes);
    geneTotalNumberOfHit3 = (int**) malloc(sizeof(int *) * totalNumberOfGenes);

    m_id = 0;
    int start_mid;
    for (geneloop=0;geneloop<totalNumberOfGenes;geneloop++) {
      geneSteadyStateRibosomeDensity[geneloop] = (double*) malloc(sizeof(double *) * obs_max_len);
      geneTotalNumberOfElongations[geneloop] = (int*) malloc(sizeof(int *) * obs_max_len);
      geneTotalNumberOfHit5[geneloop] = 0;
      geneTotalNumberOfHit3[geneloop] = 0;
      start_mid = m_id;
      for (c1=0; c1<Genes[geneloop].length; c1++) {
        geneSteadyStateRibosomeDensity[geneloop][c1] = 0.0;
        geneTotalNumberOfElongations[geneloop][c1] = 0;
        for (c2=0; c2<Genes[geneloop].exp; c2++) {
          geneSteadyStateRibosomeDensity[geneloop][c1] += avgRibosomeDensityGrid[start_mid+c2][c1];
          geneTotalNumberOfElongations[geneloop][c1] += numberOfElongationsGrid[start_mid+c2][c1];
          if (c1 == 0) {
            m_id++;
            geneTotalNumberOfHit5[geneloop] += numberOfHit5[start_mid+c2];
            geneTotalNumberOfHit3[geneloop] += numberOfHit3[start_mid+c2];
          }
        }
      }
    }

    for (geneloop=0; geneloop<totalNumberOfGenes; geneloop++) {
      for (c2=0;c2<Genes[geneloop].length;c2++) {
        fprintf(f6,"%0.2f ", geneSteadyStateRibosomeDensity[geneloop][c2]);
        fprintf(f6a,"%d ", geneTotalNumberOfElongations[geneloop][c2]);
      }
      fprintf(f6,"\n");
      fprintf(f6a,"\n");
      fprintf(f6b,"%d\t%d\t%d\n",geneloop,geneTotalNumberOfHit5[geneloop],geneTotalNumberOfHit3[geneloop]);
    }
    fclose(f6);
    fclose(f6a);
    fclose(f6b);
    free( geneSteadyStateRibosomeDensity );
    free( geneTotalNumberOfElongations );
    free( geneTotalNumberOfHit5 );
    free( geneTotalNumberOfHit3 );
  }

  strcpy( outfile, outfilePrefix);
  FILE *fh;
  fh = fopen( strcat( outfile, "/gene_data.in" ), "w" );
  fprintf(fh, "gene_id\tlength\texpression\tinit_rate\tfirst_10_codons\n" );
  for (c1=0; c1<totalNumberOfGenes; c1++) {
    fprintf(fh, "%d\t%d\t%d\t%0.4f", c1, Genes[c1].length, Genes[c1].exp, Genes[c1].ini_prob );
    if (Genes[c1].length >= 10) {
      for (c2=0; c2<10; c2++) {
        fprintf( fh, "\t%d", Genes[c1].seq[c2] );
      }
      fprintf( fh, "\n" );
    }
  }
  fclose(fh);

  strcpy( outfile, outfilePrefix);
  fh = fopen( strcat( outfile, "/trna_data.in" ), "w" );
  fprintf(fh, "trna_id\tanticodon\tconcn\tribosome_association_rate\taminoacylation_rate\n" );
  for (c1=0; c1<61; c1++) {
    if (trna_list[c1].concn > 0) {
      fprintf(fh, "%d\t%s\t%0.0f\t%0.3g\t%0.3g\n", c1, trna_list[c1].anticodon,
          trna_list[c1].concn, trna_list[c1].trnaRibosomeRateConstant, trna_list[c1].aminoacylationRateConstant );
    }
  }
  fclose(fh);

  strcpy( outfile, outfilePrefix);
  fh = fopen( strcat( outfile, "/codon_data.in" ), "w" );
  fprintf(fh, "codon_id\tcodon\ttrna1\ttrna2\tweight1\tweight2\n" );
  for (c1=0; c1<61; c1++) {
    fprintf(fh, "%d\t%s\t%0d\t%d\t%0.2f\t%0.2f\n", c1, Codon[c1].triplet,
        Codon[c1].trna[0], Codon[c1].trna[1], Codon[c1].weight[0], Codon[c1].weight[1] );
  }
  fclose(fh);

  strcpy( outfile, outfilePrefix);
  fh = fopen( strcat( outfile, "/simulation_parameters.c" ), "w" );
  fprintf(fh, "char command[600] = \"%s ", argv[0] );
  for (int argIndex=1;argIndex<argc;argIndex++) {
    fprintf(fh, "%s ", argv[argIndex] );
  }
  fprintf(fh, "\";\n");

  gettimeofday(&time2, NULL);
  fprintf(fh, "//Program run time: %d seconds\n", time2.tv_sec - time1.tv_sec );
  time_t result = time(NULL);
  fprintf(fh, "//Program end time: %s\n", ctime(&result) );
  fprintf(fh, "//INPUT VARIABLES:\n", ctime(&result) );
  fprintf(fh, 
    "// Seed for random number generator.\n"
    "int seed = %d;\n", seed);
  fprintf(fh, 
    "// Includes peptidyl-transfer and translocation time,\n"
    "// but not codon-dependent tRNA selection time\n"
    "// which is added separately below,\n"
    "// command-line option provided.\n"
    "double codonMinElongationTime = %0.4g;\n", codonMinElongationTime);
  fprintf(fh, 
    "// Number of genes,\n"
    "// required for initiatization of data structure.\n"
    "int totalNumberOfGenes = %d;\n", totalNumberOfGenes);
  fprintf(fh, 
    "// Total ribosomes.  Command line option provided.\n"
    "int totalNumberOfRibosomes = %d;\n", totalNumberOfRibosomes);
  fprintf(fh, 
    "// Fraction of ribosomes bound to mRNA\n"
    "double boundRibosomeFraction = %0.3g;\n", boundRibosomeFraction);
  fprintf(fh, 
    "// Total tRNAs\n"
    "int totalNumberOfTrna = %d;\n", totalNumberOfTrna);       
  fprintf(fh, 
    "// Total mRNAs\n"
    "// Actual value calculated from mRNA copy number file.\n"
    "int totalNumberOfMrna = %d;\n", totalNumberOfMrna);
  fprintf(fh, 
    "// Total space within a cell, m^-3\n"
    "double cellVolume = %0.3g;\n", cellVolume);
  fprintf(fh, 
    "// codons covered by a single ribosome\n"
    "int ribosomeFootprintSize = %d;\n", ribosomeFootprintSize);
  fprintf(fh, 
    "// Rate at which ribosomes prematurely terminate during normal elongation cycle.\n"
    "// Command line option provided.\n"
    "double backgroundPrematureTerminationRate = %0.2g;\n", backgroundPrematureTerminationRate);
  fprintf(fh, 
    "// Rate of premature termination for ribosomes that have\n"
    "// been hit by another ribosome from 5' side (s^-1)\n"
    "// Command line option provided.\n"
    "double prematureTerminationRateUpon5PrimeHit = %0.2g;\n", prematureTerminationRateUpon5PrimeHit);
  fprintf(fh, 
    "// Rate of premature termination for ribosomes that have\n"
    "// been hit by another ribosome from 3' side (s^-1)\n"
    "// Command line option provided.\n"
    "double prematureTerminationRateUpon3PrimeHit = %0.2g;\n", prematureTerminationRateUpon3PrimeHit);
  fprintf(fh, 
    "// kcat/km for ternary complex ribosome association,\n"
    "// commandline option provided.\n"
    "// Command line option provided.\n"
    "double ribosomeTrnaKcatByKm = %0.2g;\n", ribosomeTrnaKcatByKm);
  fprintf(fh, 
    "// starvationfactor less than 1 implies that\n"
    "// the cell is starving for an amino acid.\n"
    "// Command line option provided.\n"
    "double starvationFactor = %0.3g;\n", starvationFactor);
  fprintf(fh, 
    "// These factors multiply\n"
    "// the aminoacylation rate constants for the tRNAs\n"
    "// corresponding to the liming amino acid.\n"
    "// Command line option provided.\n"
    "double relativeTrnaAminoacylationRates[5] = {%0.2g,%0.2g,%0.2g,%0.2g,%0.2g};\n",
    relativeTrnaAminoacylationRates[0],
    relativeTrnaAminoacylationRates[1],
    relativeTrnaAminoacylationRates[2],
    relativeTrnaAminoacylationRates[3],
    relativeTrnaAminoacylationRates[4]);
  fprintf(fh, 
    "// Run options\n"
    "// Controls what outputs are printed to files.\n"
    "// Look at end of main function for file details.\n"
    "int printOptions[5] = {%d,%d,%d,%d,%d};\n",
    printOptions[0],
    printOptions[1],
    printOptions[2],
    printOptions[3],
    printOptions[4]);
  fprintf(fh, 
    "// Prefix for output file names.\n"
    "// Common-line option provided.\n"
    "char outfilePrefix[350] = \"%s\";\n", outfilePrefix);
  fprintf(fh, 
    "// Variable to temporarily store\n"
    "// the names of different output file names.\n"
    "char outfile[350] = \"%s\";\n", outfile);    
  fprintf(fh, 
    "// This file will contain the initiation rate of mRNAs,\n"
    "// mRNA copy number in the cell,a numerical codon sequence\n"
    "// for all mRNAs that are simulated.\n"
    "// Set from command-line.\n"
    "char mrnaFile[350] = \"%s\";\n", mrnaFile);
  fprintf(fh, 
    "// Input Files specifying parameters for different tRNA species.\n"
    "// tRNA concentrations.\n"
    "char trnaConcentrationFile[150] = \"%s\";\n", trnaConcentrationFile);
  fprintf(fh, 
    "// Cognate tRNA for each codon.\n"
    "char cognateTrnaFile[150] = \"%s\";\n", cognateTrnaFile);
  fprintf(fh, 
    "// Codon-tRNA interaction strength.\n"
    "char pairingStrengthFile[150] = \"%s\";\n", pairingStrengthFile);
  fprintf(fh, 
    "// Aminoacylation rates of all tRNAs.\n"
    "char aminoacylationRateFile[150] = \"%s\";\n", aminoacylationRateFile);
  fprintf(fh, 
    "// Total time simulated in seconds\n"
    "double totalSimulatedTime = %0.7g;\n", currentTime);
  fprintf(fh, 
    "// Threshold time for starting analysis in seconds.\n"
    "double thresholdTime = %0.2g;\n", thresholdTime);
  fclose(fh);

  // Free the pointer arrays created above
  delete [] Genes;
  delete [] ribosome_list;
  delete [] mRNA;
  delete [] Codon;
  delete [] trna_list;
  delete [] Rb_ao;
  delete [] Rb_ae;
  delete [] Rb_hit3;
  delete [] Rb_hit5;
  delete [] mrna_initiable;
  delete [] avgRibosomeDensityGrid;
  delete [] numberOfElongationsGrid;

  printf("Finished Simulation. Total simulated time: %0.8gs.\n",currentTime);
}
