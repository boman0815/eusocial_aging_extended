/*=============================================================================================================
                                                   individual.h
===============================================================================================================

 Definition of the class Indidivual
 
 C++-code accompanying:
 
        authors  test test
        title
 
 Written by:
        Boris Kramer and G. Sander van Doorn
        Groningen Institute for Evolutionary Life Sciences (Gelifes)
        University of Groningen
        the Netherlands
 
 Program version
        29/03/2016	: first version
 
 =============================================================================================================*/

#ifndef individual_h
#define individual_h

#include<array>
#include<vector>
#include "parameters.h"
#include "corand.h"

class Gene {
public:
	Gene() : dx(K, 0.0) {} // constructor for the start of the simulation 
	Gene(const Gene&) = default; // alternative constructor for subsequent queen generations? 
	void mutate() { dx += geneticArchitecture->operator()(); }  // what does this code do? 
	double operator()(const int &cs, const int &tr, const int &ag) { return dx[ag + (tr + cs * nTrait) * ageMax]; } // what is the perator? ()() why double brackets
	void getData(std::vector<double>&, std::vector<double>&) const; // is this also a function prototype? why the const after
	static void initGeneticArchitecture();// static member function or method, it just exists one for an objects of that class ; to call a static methos the the scope resolution operator is used :: 
private:
	Vector dx;
	static rnd::MultiNormal const * geneticArchitecture;
	const static int K = nCaste * nTrait * ageMax; // static class member  means that there is onle one of those for objects of the class
};

class Male;
class Individual {
	friend class Male; //male is friend of individuals. 
public:
	enum Traits { survival };// is 0 
	void incrementAge() { ++age; }   // if called just increses the age of the worker or queen by 1 
	void extrinsicMortality(double deathRate) { extrinsicSurvival = 1.0 - deathRate; }// extrinsic mortality function for all individuals? 
	bool doesSurvive() const { return age < ageMax && rnd::bernoulli(extrinsicSurvival * lifeHistory(survival, age)); } // why the const in front of the function? also why extrinsic survival here and not extrinsic mortality ?
    friend std::ostream& operator<< (std::ostream&, const Individual&); // what does thi do? 
	double getLifeHistoryTrait(const Traits &tr) const { return lifeHistory(static_cast<Matrix::size_type>(tr), age); }
	const Matrix& getLifeHistory() const { return lifeHistory; }
	int getAge() const { return age; }
	void getData(std::vector<double>&) const;
protected:
    Individual() : lifeHistory(nTrait, ageMax), age(0) {}
	Individual(const Individual&) = default;
	std::array<Gene, nGene> genome;
	int age;
	double extrinsicSurvival;
    Matrix lifeHistory;
	virtual void develop() = 0; 
	void mutate(); 
};

class Queen;

class Female : public Individual { // what? 
	friend class Queen;
public:
	Female() = default;
	Female(const Female&) = delete;
	Female(Queen const * const mother);
	~Female() { ++ageOfDeath[age]; } // the destructor adds the death to the age of death vector? 
	static std::array<size_t, ageMax + 1u> ageOfDeath;
private:
	void develop();
};

class Male {
	friend class Female;
public:
	Male(Individual const * const);
	Male(const Male&) = default;
private:
	Male() = delete;
	void mutate();
	std::array<Gene, nGene / 2> genome;//haploid genome
};

class Queen: public Individual {
public:
	friend class Female;   
	Queen() = delete;// is this some kind of constructor? 
	Queen(Female const * const, const Male&);
	~Queen() { ++ageOfDeath[age]; } //how does this work? adds one at the given age class
	static std::array<size_t, ageMax + 1u> ageOfDeath; // a static is in this case one array for all individuals in the class queen 
private:
	Male sperm;
	void develop();
};

#endif