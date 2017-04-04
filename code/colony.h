/*=============================================================================================================
                                                   colony.h
===============================================================================================================

 Definition of the class Colony
 
 C++-code accompanying:
 
        authors boris and sander 
        title
 
 Written by:
        Boris Kramer and G. Sander van Doorn
        Groningen Institute for Evolutionary Life Sciences (Gelifes)
        University of Groningen
        the Netherlands
 
 Program version
        29/03/2016	: first version
 
 =============================================================================================================*/

#ifndef colony_h
#define colony_h

#include "individual.h"
#include <vector>

 //defines the class COLONY

class Colony { 
public:
	Colony(Queen * const);// constructor why  a pointer? why const not in fron t of queen? 
	~Colony(); //destructor
	Queen const * getQueen() const { return queen; }  //from old code// const * declares a pointer , to declare what kind of data the pointer points to; why queen and not the class Queen? function and pointer?   
	void reproduce();//function prototype, new
	void survival();//function prototype, new
	void addWorkers();//function prototype, new
	size_t size() const { return workers.size() + (queen != nullptr); }// how does this work? returns colony size, but why queen nullpointer?  
	size_t nrOffspring() const { return offspring.size(); } //std::size_t is the unsigned integer type of the result of the sizeof operator as well as the sizeof... operator and the alignof operator 
	bool isQueenless() const { return queen == nullptr; }
	Female *getGyne();
	void getData(size_t&, size_t&, size_t&, size_t&, std::vector<double>&, std::vector<double>&) const; // why 4 times the same input to the function? 
	static const Male& getDrone(); // how does this work? 
	static void clearDrones() { drones.clear(); }
	// 20170328 added boris
	double getFitness() const { return fitness; } // working: returns the fitness of a colony with population[i]->getFitness()

	///Get the number of dead workers per age class after 'growColony'. Size of the vector equals the number of age classes
	//const std::vector<int>& getDeadWorkers() const { return deadWorkers; }
	///Get the number of dead queens per age class after 'growColony'. Size of the vector equals the number of age classes
	//const std::vector<int>& getDeadQueens() const { return deadQueens; }

	//29032017 generate output queen ages colony size and fitness/prodi into a new file at different times of the simulation 
	// double getFitness () const {return fitness};//
	//output : simulation time, queen id, col size, fitness at the beginning and at the end of the simulation
	// age at death distribution collected during the lastt 100 or so steps of the simulation to fit mortality models 
	//as we had before accumulated ages per death ot just a table: with ages at death for queens and workers.

	//work in progress 31012016
	///Get the number of dead workers per age class after 'growColony'. Size of the vector equals the number of age classes
	//const std::vector<int>& getDeadWorkers() const { return deadWorkers; }

	///Get the number of dead queens per age class after 'growColony'. Size of the vector equals the number of age classes
	//const std::vector<int>& getDeadQueens() const { return deadQueens; }

private:
	double fitness;
	std::vector<Female*> workers;// what do those vectors do exactly
	Queen* queen; //???
	std::vector<Female*> offspring;// what do those vectors do exactly
	static std::vector<Male> drones;// what do those vectors do exactly static means just one vector for all members of the class colony 

	//20170328
	// including death counts for individuals
	//The number of dead workers per age class (after 'growColony')
	//std::vector<int> deadWorkers;
	///The number of dead queens per age class (after 'growColony')
	//std::vector<int> deadQueens;
	
};

#endif