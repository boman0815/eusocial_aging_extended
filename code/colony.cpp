/*=============================================================================================================
                                                   colony.cpp
===============================================================================================================

 Implementation of the class Colony
 
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
 
 =============================================================================================================*/

#include <cassert>
#include <cmath>
#include "colony.h"
#include "random.h"

std::vector<Male> Colony::drones;

//default constructor for Colony
Colony::Colony(Queen * const q) : queen(q), fitness(0.0)
{

}

// destructor
Colony::~Colony() 
{
	if(queen != nullptr) delete queen;
	while (!workers.empty()) {
		delete workers.back();
		workers.pop_back();
	}
}

// reproduce function 
// if no queen calculate fitness and workers lay eggs? 
// includes a life history, where the queen and the workers make allocation decisions of they forage of lay eggs , f they forage then they are exposed to different levels of exrrinsiv mortality given by the mortality functions; queens trade of resources vs reproduction
//else queen lays eggs allocation limited to 1.0
//fitness: W= (eggs*resources) /(eggs + resources) 
void Colony::reproduce()  //:: scope resolution operator to define methods outside class definition
{
	double workerMortality;// creates mortality  and reproduction as thery are linkes to a trade off as resources are linkes to the risc of dieing from extrinsic sources 
	if (isQueenless()) { // workers only produce drones if there is no queen //  generally drones and queens are only poduced if needed in the simulation (queen dies; new queen needs to mate) 
		const size_t n = workers.size();
		assert(n); // d just checked if n exist,and the program terminates if not true
		double optimalWorkerAllocation = 1.0 / (1.0 + sqrt(workerFecundity / productivity));// ?investment into reproduction   // trade off between 
		workerMortality = 1.0 - pow(1.0 - foragingRisk, 1.0 - optimalWorkerAllocation);// calculates worker mortality// pow(base ,exponent)
	
		const double eggs = workerFecundity * optimalWorkerAllocation;// why this term especially since  workerFecundity is already in optimal worker allocation
		const double resources = productivity * (1.0 - optimalWorkerAllocation);//whats not invested into eggs
		fitness = n * eggs * resources / (eggs + resources);// calculates fitness of the colony 
		for (size_t i = rnd::poisson(fitness); i > 0u; --i) // ???
			drones.push_back(Male(workers[rnd::integer(n)]));// workers produce drones ?!
	}
	else {
		double optimalQueenAllocation = size() / (1.0 + sqrt(queenFecundity / productivity)); // if no workers 20 percents  energy goes into egglaying // size = worker number? // where is the decision tht the queen forages or not ??
		if (optimalQueenAllocation > 1.0)
			optimalQueenAllocation = 1.0; //  only eggs produced  n extrinsic mortality
		queen->extrinsicMortality(1.0 - pow(1.0 - foragingRisk, 1.0 - optimalQueenAllocation)); // =(1-fR)^(1-oQA)   , why not just call the function? why is the formula in here ? 
		workerMortality = foragingRisk; // why is this "calculated" here? 
		
		const double eggs = queenFecundity * optimalQueenAllocation;//
		const double resources = productivity * (1.0 - optimalQueenAllocation + workers.size());
		fitness = eggs * resources / (eggs + resources); // if no resources eggs dont survive // --> eggs survival// just number of surviving eggs
		for (size_t i = rnd::poisson(fitness); i > 0u; --i) { // get the actual number oof  eggs to calculate the fitness
			Female *tmp = new Female(queen);//tmp is a pointer that point to an object of type female; if an object is created like this arrow operator is used to access methods see line below  otherwise female tmp to create object ans tmp.extrinsicMortality to call for the method
			tmp->extrinsicMortality(0.0);
			offspring.push_back(tmp); // but the offspring is a private memebr of the class colony, or is this another offspring? eggs 
			drones.push_back(Male(queen)); //queens also produce drones?!
		}
	}
	for (Female * worker : workers)
		worker->extrinsicMortality(workerMortality);// calculates worker survival probability feeds int does survive function
	//
	
	//
}


// survival function // interesting mix of infividual / colony level functions 
void Colony::survival() 
{
	// worker survival
	for (size_t i = 0u; i < workers.size();) { // increase age of all workers
		workers[i]->incrementAge();
		if (workers[i]->doesSurvive()) ++i; // check if the individuals survives
		else {
			delete workers[i];
			workers[i] = workers.back();// whats the differecnes btween the back ans pop_back
			workers.pop_back();
		}
	}

	// queen survival 
	if (!isQueenless()) {// if queen is alive 
		queen->incrementAge();
		if (!queen->doesSurvive()) {// if the queen doesnt survive
			delete queen;
			queen = nullptr;// why is this needed? 
		}
	}
	//I
	// queen replacement     
	while (isQueenless() && gyneReplacesQueen && offspring.size() && drones.size()) {// if  there is no queen and gyne replacement is true and there is offspring and drones in the colony then :
		queen = new Queen(offspring.back(), getDrone());// new queen queen is created from the last entry in the existing offspring vector of the colony
		if (!queen->doesSurvive()) { //then she has to survive
			delete queen;  // otherwise gets deleted 
			queen = nullptr;
		}
		delete offspring.back(); // when the queen was replaced by the offspring then  the last entry gets removed from the vector (the individuals that became the new queen 
		offspring.pop_back();
	}
	//II
	if (isQueenless() && workerReplacesQueen && workers.size() && drones.size()) {
		size_t i = rnd::integer(workers.size());
		queen = new Queen(workers[i], getDrone());
		delete workers[i];
		workers[i] = workers.back();
		workers.pop_back();
	}
	//III
	// added 02022017 to remove bug that appeared when existing workers replace the queen but offspring.size was 0 ; now if there is no worker the offspring can develop into a new queen mate and continue the colony
	while (isQueenless() && workerReplacesQueen && offspring.size() && drones.size()) {// what?
		queen = new Queen(offspring.back(), getDrone());// create new queen
		if (!queen->doesSurvive()) { // why is this here again? 
			delete queen;
			queen = nullptr;
		}
		delete offspring.back();
		offspring.pop_back();
	}
	/// end new code 02022017
}

void Colony::addWorkers()
{
	while (!offspring.empty()) { // why can the vector class methods be used here? this meant if there are workers? 
		if (offspring.back()->doesSurvive()) // what is thei combination? 
			workers.push_back(offspring.back());//?
		else
			delete offspring.back(); // offspring is a vector of the individuals of a colony? 
		offspring.pop_back();
	}
}

// select new queen for the colony? 
Female* Colony::getGyne()// ??? 
{
	//std::cout << "offspring.size() =" << "\t" << offspring.size() << "\n";  // number of undeveloped offspring (wheights) ..? 
	assert(offspring.size());// test if the vector exists and has entries? // the simulation breaks here if gyne repqueen= false & workerRepQueen= true
	Female *tmp = offspring.back();
	offspring.pop_back();
	return tmp;
}

// select male for the colony? 
const Male& Colony::getDrone() 
{
	assert(drones.size());
	return drones[rnd::integer(drones.size())];// randomly chooses one drone ?
}

// get data, s how do those work in detail? in the end we get average values in the written data
void Colony::getData(size_t &sumn, size_t& sumq, size_t& sumw, size_t& suma, std::vector<double> &sumxq, std::vector<double> &sumxw) const // how does this function work? 
{
	if (size()) {			 //if size is not empty
		++sumn;				 // number of occupied colonies
		if (!isQueenless()) {// if there is a queen? 
			++sumq;			 //number of colonies with a queen
			suma += queen->getAge(); // ? queen.getAge how is that differtn?
			//std::cout << "queenage = "<< queen->getAge()<< "\n";///20170328 testing on how to get the age working when data analysis is called 
			queen->getData(sumxq); //age specific survival for each individuals
		}
		sumw += workers.size();
		for (Female * worker : workers)
			worker->getData(sumxw);
	}
}