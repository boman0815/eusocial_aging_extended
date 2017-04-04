/*=============================================================================================================
												   simulation.cpp
===============================================================================================================

 Model implementation; entry point main()

 C++-code accompanying:

		authors
		title

 Written by:
		Boris Kramer and G. Sander van Doorn
		Groningen Institute for Evolutionary Life Sciences (Gelifes)
		University of Groningen
		the Netherlands

 Program version
		29/03/2016	: first version

 ============================================================================================================*/

#include "colony.h"
#include "random.h"
#include "utils.h"
#include "parameters.h"
#include <numeric>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iterator>
#include <fstream>
#include <sstream>

 /*=============================================================================================================
													global variables
 =============================================================================================================*/

unsigned int simulationId;
Colony * population[nColony];//? what does this do? 
std::ofstream ofs;

//std::ofstream ofSurvival; // not tested yet 25012017
//std::ofstream ofDeaths;
std::ofstream ofParameters;
//std::vector<int> deadQueenCounter(ageMax, 0);// vector of 0 length ageMax
//std::vector<int> deadWorkerCounter(ageMax, 0);

/*=============================================================================================================
												   model implementation
=============================================================================================================*/

///Checks if a file exists
/// @param filename name of the file
/// @returns if the file exists
bool is_regular_file(const std::string& filename)
{
	std::fstream f;
	f.open(filename.c_str(), std::ios::in);
	return f.is_open();
}

// construct  outfile names  (sander life history, survival patterns and death counts ) no  needed here 
std::string constructFilenameParameters(const long long simulationId) // why that complicated approach to get the name for the file to write? 
{
	std::stringstream par; // im still confused bt the usage of the names s f g ofs for the objects of class os
	par << "parameters_" << simulationId << ".csv";
	return par.str();   // add headers?? 
}

//not tested jet 
////construct other filenames
//std::string constructFilenameSurvival(const long long simulationId) // why that complicated approach to get the name for the file to write? 
//{
//	std::stringstream surv;
//	surv << "survival_" << simulationId << ".csv";
//	return surv.str();
//}
//
////construct other filenames
//std::string constructFilenameDeaths(const long long simulationId) // why that complicated approach to get the name for the file to write? 
//{
//	std::stringstream death;
//	death << "deaths_" << simulationId << ".csv";
//	return death.str();
//}
// end richel stuff 25.1.2017




void init()
// initialisation
{
	//create simulation id and use it as a seed for the random number generator
	simulationId = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
	rnd::set_seed(simulationId);

	// create mutational variance covariance matrix 
	Gene::initGeneticArchitecture();

	// create initial population
	Female *q = new Female(); //only one single female to start with? from which the whole population is drawn.
	for (int i = 0; i < nColony; ++i)
		population[i] = new Colony(new Queen(q, Male(q)));// initialize with mated queen?!
	delete q;

	// prepare data arrays
	for (int k = 0; k <= ageMax; ++k)
		Queen::ageOfDeath[k] = Female::ageOfDeath[k] = 0u; //[k] refers to the entry in the vector/ why do they need to be initialized with unsigned 0 here? they are statics of the classes

	// open file sander for : t		occupied	w/ queen	avg #workers	avg queen age	avg queen lifespan	avg worker lifespan		queen intrinsic survival	worker survival	
	std::ostringstream str;
	str << "data_" << simulationId << ".csv";
	ofs.open(str.str());
	ofs << "t" << ',' << ',' // creates the first row of the data_xx.csv file (header) 
		<< "occupied" << ','
		<< "w/ queen" << ','
		<< "avg #workers" << ','
		<< "avg queen age" << ','
		<< "avg queen lifespan" << ','
		<< "avg worker lifespan" << ','
		<< ',' << "queen intrinsic survival";
	for (int a = 0; a < ageMax; ++a)
		ofs << ',' << "age class " << a + 2;
		ofs << ',' << "worker intrinsic survival";
	for (int a = 0; a < ageMax; ++a)
		ofs << ',' << "age class " << a+2; // puts the commas after queen surv in the header 
	ofs << "\n";

	// open parameters_xxx.csv  
	const std::string fileName = constructFilenameParameters(simulationId);
	ofParameters.open(fileName);														// new richel version
	if (!ofParameters.is_open())														// If it is open we can do our writing to the file.
		std::cout << "Something went wrong with opening the parameter file!" << "\n"; 	// old--> throw std::runtime_error("Cannot open file");
		ofParameters << "Parameter" << "," << "Value" << "," << "Value" << "\n"			//header
		<< "nColony" << "," << nColony << "\n" 											
		<< "tEnd" << "," << tEnd << "\n"
		<< "dataInterval" << "," << dataInterval << "\n"
		<< "ageMax" << "," << ageMax << "\n"
		<< "foragingRisk" << "," << foragingRisk << "\n"
		<< "survivalOffset" << "," << survivalOffset << "\n"
		<< "productivity" << "," << productivity << "\n"
		<< "queenFecundity" << "," << queenFecundity << "\n"
		<< "workerFecundity" << "," << workerFecundity << "\n"
		<< "gyneReplacesQueen" << "," << gyneReplacesQueen << "\n"
		<< "workerReplacesQueen" << "," << workerReplacesQueen << "\n"
		<< "nGene" << "," << nGene << "\n"
		<< "nCaste" << "," << nCaste << "\n"
		<< "nTrait" << "," << nTrait << "\n"
		<< "mutationRate" << "," << mutationRate << "\n"
		<< "sigma[nCaste][nTrait]" << "," << *sigma[0] << "," << *sigma[1] << "\n"    
		<< "alpha[nTrait]" << "," << *alpha << "\n"           
		<< "beta[nTrait][(nCaste * (nCaste - 1)) / 2]" << "," << *beta[0] << "\n"
		<< " eta[nTrait]" << "," << *eta << "\n";     
	ofParameters.close();																// After you are done with the file always close it.

//
//
//// untested 25012016
////open file for suvival rates
//	const std::string fileName2 = constructFilenameSurvival(simulationId);
//	ofSurvival = std::ofstream(fileName2);
//	if (!ofSurvival.is_open()) throw std::runtime_error("Cannot open file"); // construct header for csv file 
//	ofSurvival << "generation ," << "  factor , " << " ," << "  caste , ";
//	for (int u = 1; u <= ageMax; ++u) {
//		ofSurvival << "age_class_" << u << " ,";
//	}
//	ofSurvival << '\n';
//	ofSurvival.fill(','); // commas between entries? how does c++ know where to write the data and even more interesting can this also be used in a loop to grab the data again.. 
//
//	//open file for death counts 
//	const std::string fileName3 = constructFilenameDeaths(simulationId);
//	ofDeaths = std::ofstream(fileName3);
//	if (!ofDeaths.is_open()) throw std::runtime_error("Cannot open file");
//	ofDeaths << "generation ," << "  factor , " << " ," << "  caste , "; //creates the header for the csv file 
//	for (int v = 1; v <= ageMax; ++v) {
//		ofDeaths << "age_class_" << v << " ,";
//	}
//	ofDeaths << '\n';
//	ofDeaths.fill(','); // commas between entries? how does c++ know where to write the data and even more interesting can this also be used in a loop to grab the data again.. 
//	//
	// other richel stuff from previous model  25.1.2017
	//



	//// write .txt file with parameters ... working
	//std::ofstream parameterFile ("parameterFile.txt");
	////parameterFile<< "parameters_" << simulationId << ".txt";   // not working
	//if (parameterFile.is_open())   // Always check to see if the file is open and for errors.
	//{
	//	parameterFile << "nColony" << "\t" << nColony << "\n" 	// If it is open we can do our writing to the file.
	//		<< "tEnd" << "\t" << tEnd << "\n"
	//		<< "dataInterval" << "\t" << dataInterval << "\n"
	//		<< "ageMax" << "\t" << ageMax << "\n"
	//		<< "foragingRisk" << "\t" << foragingRisk << "\n"
	//		<< "survivalOffset" << "\t" << survivalOffset << "\n"
	//		<< "productivity" << "\t" << productivity << "\n"
	//		<< "queenFecundity" << "\t" << queenFecundity << "\n"
	//		<< "workerFecundity" << "\t" << workerFecundity << "\n"
	//		<< "gyneReplacesQueen" << "\t" << gyneReplacesQueen << "\n"
	//		<< "workerReplacesQueen" << "\t" << workerReplacesQueen << "\n"
	//		<< "nGene" << "\t" << nGene << "\n"
	//		<< "nCaste" << "\t" << nCaste << "\n"
	//		<< "nTrait" << "\t" << nTrait << "\n"
	//		<< "mutationRate" << "\t" << mutationRate << "\n"
	//		<< "sigma[nCaste][nTrait]" << "\t" << *sigma[0] << "\t" << *sigma[1] << "\n"   // how do i obtain those values  not working 
	//		<< "alpha[nTrait]" << "\t" << *alpha << "\n"           // how do i obtain those values  not working 
	//		<< "beta[nTrait][(nCaste * (nCaste - 1)) / 2]" << "\t" << *beta[0] << "\n"  // how do i obtain those values  not working 
	//		<< " eta[nTrait]" << "\t" << *eta << "\n";     // how do i obtain those values  not working 
	//}
	//else
	//{
	//	std::cout << "Something went wrong with opening the parameter file!" << "\n"; 	// If the file isn't open something went wrong. Point that out.
	//}
	//parameterFile.close();  // After you are done with the file always close it.
}

void dataAnalysis(int t) // sander can you please explain me the details here? 
{
	std::vector<double> sumxq(ageMax, 0.0), sumxw(ageMax, 0.0);
	size_t sumw = 0u, sumq = 0u, sumn = 0u, suma = 0u;//  0u=  0 unsigned; size_t ? 

	for (int i = 0; i < nColony; ++i)
		population[i]->getData(sumn, sumq, sumw, suma, sumxq, sumxw);//what are those variables? sumn=number of colonies ; sumq number of queens, sumw number of workers , is summs everz vaue up for the whole simulation, than in the next lines the mans are calculated 

	//std::cout << "Time = "<< t << '\t' << "sumn = " <<sumn << '\t' << "sumq * 1.0 / sumn = "<< sumq * 1.0 / sumn << '\t' << "sumw * 1.0 / sumn = "<< sumw * 1.0 / sumn << '\n'; // sanders output deactivated for cluster runs sumn should be the number of occupied colonies ; sumq is the number of queens

	ofs << t << ',' << ',' // output into data_simulationxxx.csv file 
		<< sumn * 1.0 / nColony << ','//occupied colonies
		<< sumq * 1.0 / sumn << ','   // colonies with queen
		<< sumw * 1.0 / sumn << ','   // average number of workers per colony
		<< (sumq ? suma * 1.0 / sumq : nan("")) << ',';

	sumn = suma = 0u; //clear the variables after writing to file or better before using them again
	for (int k = 0; k <= ageMax; ++k) {
		double n = Queen::ageOfDeath[k];//why is the scope resolution operator used here ? k is the age classes thus n is the number of deaths in a given age class? 
		sumn += n;    //sums all the entries of the ages of death vector 
		suma += k * n;//
		Queen::ageOfDeath[k] = 0u; //after reading the entry its replaces with a 0 for the next timestep	
	}
	ofs << suma * 1.0 / sumn << ','; // calculates average queen lifespan of the queens that have died in this timestep

	sumn = suma = 0u;
	for (int k = 0; k <= ageMax; ++k) {
		double n = Female::ageOfDeath[k];
		sumn += n;
		suma += k * n;
		Female::ageOfDeath[k] = 0u;
	}
	ofs << suma * 1.0 / sumn << ',';

	for (int a = 0; a < ageMax; ++a)
		ofs << ',' << (sumq ? sumxq[a] / sumq : nan(""));// what?
	ofs << ',';
	for (int a = 0; a < ageMax; ++a)
		ofs << ',' << (sumw ? sumxw[a] / sumw : nan(""));
	ofs << '\n';

	//std::cout << "sumn = " << sumn << "\n" // boris output to check what in the sumn etc 
	//	<< "sumq = " << sumq << "\n"
	//	<< "suma = " << suma << "\n";
	ofs.flush();
}

 //create a sturcture for the measurements that are goint to be written to a file //The size of Measurment will be the number of colonies
struct Measurement
{
	///Time of measurement, in timesteps
	int t;  // integer for timestep of the simulation

	///Queen ages. A negative number means NA, no queen
	std::vector<double> queenAges; //vector to store the queen ages 

	///Colony sizes
	std::vector<int> sizes; //vector to store the colony sizes

	///Fitmessess
	std::vector<double> fitnesses; // vector to store the fitnes of the colonies 

	//ages at death   to store in another file tried but failed 
	//std::vector<int> qDeathCounts;
	//std::vector<int> wDeathCounts;
};

Measurement measure(const int t, Colony ** pop , const int nColonies)// member function? 
{
	Measurement m;
	m.t = t;
	//Measure all queen ages
	for (int i = 0; i != nColonies; ++i)
	{
		const Colony * const this_colony = pop[i];
		assert(this_colony);
		const Queen * const queen = this_colony->getQueen();//why this ?
		if (!queen)// no queen put -1 for the age
		{
			m.queenAges.push_back(-1.0);
		}
		else
		{
			// No pointer:
			// Queen q;
			// q.getAge()
			//
			// With a  pointer:
			// Queen * q;
			// q->getAge()
			// (*q).getAge()
			//
			m.queenAges.push_back(queen->getAge());// get the age and put it at the end of the vector
		}
	}
	//Measure colony sizes
	for (int i = 0; i != nColonies; ++i)
	{
		const Colony * const this_colony = pop[i];
		assert(this_colony);
		m.sizes.push_back(this_colony->size());
	}
	//Measure fitnesses
	for (int i = 0; i != nColonies; ++i)
	{
		const Colony * const this_colony = pop[i];
		assert(this_colony);
		m.fitnesses.push_back(this_colony->getFitness());
	}
	return m;// is an object of struct measure that contains t, queenages, sizes and fitnesses
}

//Measurement deaths(const int t, Queen, Female, const int ageMax)// tried but failed to put the right input 
//{
//	Measurement de;
//	de.t = t;
//
//	for (int hi = 0; hi < ageMax; ++hi) {
//		de.qDeathCounts.push_back(Queen::ageOfDeath[hi]);
//		de.wDeathCounts.push_back(Female::ageOfDeath[hi]);
//	}
//		
//	return de;// 
//}

//richel for results file 
std::ostream& operator<<(std::ostream& os, const Measurement& m) noexcept // writing to the outstream not sure why this has to be here? and when is it called? 
{
	assert(m.queenAges.size() == m.sizes.size());//checs if the size of the ages vector is similar to the sizes vector
	assert(m.queenAges.size() == m.fitnesses.size());

	const int sz = m.queenAges.size();
	for (int i = 0; i != sz; ++i)
	{
		os << m.t << ',' << m.queenAges[i] << ',' << m.sizes[i] << ',' << m.fitnesses[i] << '\n';// this basicall stores all the information i want to have in the results file at the end
	}
	return os; // this is all the data
}
//richel for results file 
void saveResults(const std::vector<Measurement>& results, const std::string& filename)
{
	std::ofstream f(filename);
	f << "t,queen_age,colony_size,fitness\n";

	for (const Measurement& m: results) f << m;
}

//write the daeths per age class into another file deaths_xxx.csv
int writeFile(int t) {  //so far this only writes one file for one line of results  the both age at death vectors, to improove just open the and use a function to write to it when needed. 
	std::ofstream aFile;
	aFile.open("deaths_time_"+std::to_string(t) +"_"+ std::to_string(simulationId) + ".csv");
	aFile << "t";
	for (int ha = 0; ha < ageMax; ++ha)
	{
		aFile << ',' << "queenAge_" << ha + 1;
	}
	for (int ha = 0; ha < ageMax; ++ha)
	{
		aFile << ',' << "workerAge_" << ha + 1;
	}
	aFile << "\n";
	aFile << t;
	for (int ha = 0; ha < ageMax; ++ha)
	{
		aFile << "," << Queen::ageOfDeath[ha];
	}
	for (int ha = 0; ha < ageMax; ++ha)
	{
		aFile << "," << Female::ageOfDeath[ha];
	}
	aFile << "\n";
	aFile.close();
	return 0;
}


/*=============================================================================================================
												   main()
=============================================================================================================*/
int main()    // the designated start of the program
{
	//preliminaries
	init();
	std::vector<Measurement> results;

	for (int t = 0; t < tEnd; ) {
		// testing output stuff 
		//std::cout << "beta stuff" << "\t" << beta[nTrait][(nCaste * (nCaste - 1)) / 2] << '\t' << nTrait << '\t' << *beta << '\t' << *beta[0] << "\t" << *beta[1] << "\t" << beta[nCaste * (nCaste - 1) / 2] << '\n';

		// production of offspring and survival
		size_t sum = 0u;
		std::vector<size_t> weights(nColony, 0u); // weights reprecesnt the reproduction, dont really know why or how this is different from offspring, should be nearly the same //creates vector of length ncolony with unsigned 0s
		for (int i = 0; i < nColony; ++i) {
			if (population[i]->size()) { // if there are individuals in the colony? 
				population[i]->reproduce();
				population[i]->survival();
				sum += weights[i] = population[i]->nrOffspring(); // count the number of offspring per colony and sum it all up 
			}
		}

		// establishment of new colonies throughout the simulation only when the colony is emppty (no workers, no queen and no offspring (weights) 
		for (int i = 0; i < nColony && sum; ) { // go through all colonies

			//20170328 get fitness and size of a colony working
			//how do i acess the queen and its member functions to ge the age of an individual queen in the colony
			//std::cout << "t = " << t << " " << "i = " << i << " Fitness = "
			//	<< population[i]->getFitness() << " colony size = " << population[i]->size()
			//	<< " queen age = ";
			//if (population[i]->getQueen()) //Does a Queen exist, i.e. is it not nullptr
			//{
			//	std::cout << population[i]->getQueen()->getAge(); /// getQueen is a function of the class colony to get the queen, then we use a function to get the age of the queen which is a memeber of the class queen
			//}
			//else
			//{
			//	std::cout << "NA";
			//}
			//std::cout <<" weights = " << weights[i]<< "\n";
			////////

			if (population[i]->size() == 0u && weights[i] == 0u) {  // i = colony we are looking at j colony the gyne is taken from // 02022017 added : && weights[i] == 0u to make sure no individual are taken from the empty offspring vector 
				size_t j = std::discrete_distribution<size_t>(weights.begin(), weights.end())(rnd::rng); // to calculate potential fitness ..expected based on fitness
				//std::cout << i << " = i " << j << " = j  " << weights[j] << " = weights are equal to: ";
				Female * tmp = population[j]->getGyne(); // get new queen
				Queen * q = new Queen(tmp, Colony::getDrone());
				delete tmp;// delete temp female
				if (q->doesSurvive()) {
					delete population[i]; //old colony killed, but shouldnt it survive for the workers? 
					population[i] = new Colony(q);// new colony established
					++i;
				}
				else delete q;
				--weights[j];
				--sum;
			}
			else ++i;
		}

		bool terminate = true;
		for (int i = 0; i < nColony; ++i) {
			population[i]->addWorkers();
			if (population[i]->size())
				terminate = false;
		}
		Colony::clearDrones();

		if (terminate) {
			std::cout << "population extinct\n";
			break;
		}
		//20170328 working, outputs all the entries of ages at death array
		//for (int hi = 0; hi < ageMax; ++hi) {
		//	std::cout << "deaths age queen "<<hi<<" = " << Queen::ageOfDeath[hi] << " , "<< " worker ="<< Female::ageOfDeath[hi];
		//	if (hi == ageMax-1 ) std::cout << "\n";
		//}
		//std::cout << "deaths age 1 = " << Queen::ageOfDeath[1] << "\n";
		if (t == tEnd - 1 || t== tEnd-2 )//output the age ofdeath vector to a file 
		{
			writeFile(t);// stores the age at death distribution but just one line so far 
		}

		++t;//update the simulation time after everything is done in the timestep

		if (t % dataInterval == 0)  // run analysis at each dataInterval
		{
			dataAnalysis(t);
			// here the measurements (queen age , col size, fitness) are made and and stored in results, the file will be written after the last timestep of the simulation
			const Measurement result = measure(t, population, nColony); //here the results are masured by creating a structure of measurement named result?!
			results.push_back(result);// adds the measures for the next timesteps 
		}  
		//const Measurement result2 = deaths (t, Queen& , &Female);// tried but failed 
		
	}

	//saving the  queenAge, colony size and fitness of the colony to a file after the simulation ends
	const std::string results_filename = "results_" + std::to_string(simulationId) + ".csv";// create filename for the size age results 
	saveResults(results, results_filename); // here the results file is saved (end of simulation)
	assert(is_regular_file(results_filename));//check if there is a regular file 

	for (int i = 0; i < nColony; ++i)
		delete population[i];//delete all colonies after the end of the simulation
	return 0;
}

