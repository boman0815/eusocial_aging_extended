#ifndef parameters_h   // checks if parameters are already included
#define parameters_h

//colony dynamics
const int nColony = 1000;									//number of colonies
const int tEnd = 20;										//simulation time	
const int dataInterval = 10;								//save date every dataInterval time steps 

//life history parameters
const int ageMax = 20;                                      // maximum lifespan 
const double foragingRisk	= 0.0;							//replaces extrinsic mortality
const double survivalOffset	= 0.5;//choose values from -1:1 for initial survival of the first queen created 0:1 hyperbolic tanget function (ind.cpp l.174) 
const double productivity	= 1.0;
const double queenFecundity = 10.0;							// higher setting of fecundity makes eggs cheaper ( more eggs per resource unit) 
const double workerFecundity = 0.0;	//was 5 in the test runs// higher setting of fecundity makes eggs cheaper ( more eggs per resource unit) 

// 4 scenarios
const bool gyneReplacesQueen = true;      //true working  fine , when using worker replacemnt the code generates in error in colony.cpp l 143 (assert(offspring.size()))
const bool workerReplacesQueen = true;				//false working fine 
//
//
//


//genetics
//this should represent the mutation accumulation setting , mutations just affect the given age class not other effects
const int nGene = 40;//was set to 100 and working, now back to 40 (20queen,20 worker genes?!//number of genes; has to be even;ech gene has an effect on all age classes and the effect is summed up (ind. cpp l. 170?); more genes basically lead to a larger target for mutations 
const int nCaste = 2;									//queens and workers
const int nTrait = 1;									//survival and allocation
const double mutationRate = 1.0e-4;						//rate of mutations per gene per generation
const double sigma[nCaste][nTrait] = { { 0.08 },{ 0.08 } };//was set to 0.01 in the test runs //mutational effect size over age classes, increase= more age classes are affected by a mutation   ; *sigma[2] points to alpham, *sigma[0], *sigma[1] works
const double alpha[nTrait] = {0.0};						//correlation between subsequent age classes within trait  
const double beta[nTrait][(nCaste * (nCaste - 1)) / 2] = { { 0.0 } };   //genetic correlation within trait, between castes ; orig value 0.0

const double eta[nTrait] = {-0.2};						//mutation bias (relative to mutational effect size)  the final value is -.02*the standart derivation of the mutational effect size
//sigma, alpha and beta are important to set up the correlation marix in ind.cpp l. 50ff; 
														
	//old stuff													
//const double ggamma[nCaste][(nTrait * (nTrait - 1)) / 2] =  //antagonistic pleiotropy; correlation between trait, within caste
//    {{-0.1},{-0.1}};
//const double delta = 0.1;								//antagonistic pleiotropy delay


#endif