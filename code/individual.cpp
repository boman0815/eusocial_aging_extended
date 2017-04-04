/*=============================================================================================================
                                                   individual.cpp
===============================================================================================================

 Implementation of the class Indidivual
 
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

#include "individual.h"
#include "random.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

std::array<size_t, ageMax + 1u> Female::ageOfDeath;
std::array<size_t, ageMax + 1u> Queen::ageOfDeath;

/*=============================================================================================================
                                Implementation of the class Individual::Gene
=============================================================================================================*/


// sander : how doe this matrices work? ? !!!! 
rnd::MultiNormal const * Gene::geneticArchitecture;

extern unsigned int simulationId;
void Gene::initGeneticArchitecture()
{
	//create partial correlation matrix
	Matrix P(K,K, 0.0);

	//partial correlations within trait// set it up matrix p//mutational variance covariance matrix(multinormal distribution), parms alpha and beta
    for (int cs1 = 0; cs1 < nCaste; ++cs1) {
        const int offset = cs1 * nCaste - ((cs1 + 1) * (cs1 + 2)) / 2;
		for (int tr = 0; tr < nTrait; ++tr)
			for (int i = 0; i < ageMax; ++i) {
				const int k1 = i + (tr + cs1 * nTrait) * ageMax;
				P(k1, k1) = 1.0;
                if (i - 1 >= 0) P(k1, k1 - 1) = -alpha[tr];
                if (i + 1 < ageMax) P(k1, k1 + 1) = -alpha[tr];
                for (int cs2 = cs1 + 1; cs2 < nCaste; ++cs2) {
                    const int k2 = i + (tr + cs2 * nTrait) * ageMax;
                    P(k2, k1) = P(k1, k2) = -beta[tr][offset + cs2] * P(k1, k1);
                    if (i - 1 >= 0) P(k2 - 1, k1) = P(k1, k2 - 1) = beta[tr][offset + cs2] * P(k1, k1 - 1);
                    if (i + 1 < ageMax) P(k2 + 1, k1) = P(k1, k2 + 1) = beta[tr][offset + cs2] * P(k1, k1 + 1);
                }
			}
    }
	/*
	//partial correlations between traits
    for (int tr1 = 0; tr1 < nTrait; ++tr1) {
        const int offset = tr1 * nTrait - ((tr1 + 1) * (tr1 + 2)) / 2;
        for (int cs = 0; cs < nCaste; ++cs)
            for (int i = 0; i < ageMax; ++i) {
                const int k1 = i + (tr1 + cs * nTrait) * ageMax;
                for (int tr2 = tr1 + 1; tr2 < nTrait; ++tr2) {
                    double wgt = 1.0;
                    for (int j = i; j < ageMax; ++j) {
                        const int k2 = j + (tr2 + cs * nTrait) * ageMax;
						P(k1, k2) = P(k2, k1) = - wgt * ggamma[cs][offset + tr2];
                        wgt *= delta;
					}
				}
            }
    }*/

    //write partial correlation matrix to file
	std::ostringstream oss;
	oss << "architecture_" << simulationId << ".csv";  //  writes the filename ofthe file to be written
	std::ofstream ofs(oss.str().c_str());
	verify(ofs.is_open());
	ofs.fill(',');
	ofs << "partial correlation matrix\n\n" << P << '\n'; // writes header in the file and then the partial correlation matrix

	//calculate the mutational variance/covariance matrix
	P = P.inverse();
    Vector bias(K);
    for(int i = 0; i < K; ++i) {
        const double si = sigma[i / (ageMax * nTrait)][(i / ageMax) % nTrait];
		const double etai = eta[(i / ageMax) % nTrait];
        for(int j = i + 1; j < K; ++j) {
            const double sj = sigma[j / (ageMax * nTrait)][(j / ageMax) % nTrait];
            P(i , j) = P(j , i) *= si * sj / sqrt(P(i, i) * P(j, j));
        }
        P(i, i) = sqr(si);
        bias[i] = etai * si;
    }
	ofs << "\nmutational variance covariance matrix\n\n" << P << '\n'; //writes the name and the matrix into the file 

	//initialize multivariate Gaussian distribution of mutational effects
	geneticArchitecture = new rnd::MultiNormal(P, bias);
    
    ofs << "\nstandard deviations of the multivariate normal and eigenvectors of the mutational variance covariance matrix\n\n"
        << *geneticArchitecture << '\n';

    ofs << "\nsampled mutations\n\n";
    for (int i = 0; i < K; ++i) ofs << geneticArchitecture->operator()() << '\n'; // writes the mutations , what does geneticArchitecture->operator()() do? how does it work?
    
    ofs.close(); // close and save the file to the disk 
}

void Gene::getData(std::vector<double> &sumx, std::vector<double> &sumxx) const
{
	for (int i = 0; i < K; ++i){
		double tmp = dx[i];
		sumx[i] += tmp;
		sumxx[i] += tmp * tmp;
	}
}

/*=============================================================================================================
                                Implementation of the classes Individual, Sperm, Worker and Queen
=============================================================================================================*/

Male::Male(Individual const * const mother)
//constructor for production of drone by worker or queen
{
	for (int i = 0; i < nGene; i += 2)
		genome[i / 2] = rnd::uniform() < 0.5 ? mother->genome[i] : mother->genome[i + 1];
	mutate();
}

Female::Female(Queen const * const mother)
//constructor for sexual reproduction -> makes a new female (worker or gyne)
{
	for (int i = 0; i < nGene; i += 2) {
		genome[i] = rnd::uniform() < 0.5 ? mother->genome[i] : mother->genome[i + 1];
		genome[i + 1] = mother->sperm.genome[i / 2];
			}
	mutate();
	develop();
}

Queen::Queen(Female const * const female, const Male &drone) : Individual(*female), sperm(drone)
//constructor -> promotes female to queen
{
	develop();
}

void Individual::mutate()
{
	if (rnd::uniform() < nGene * ageMax * mutationRate) 
		genome[rnd::integer(nGene)].mutate();
}

void Male::mutate()
{
	if (rnd::uniform() < nGene * ageMax * mutationRate / 2) // because haploid? 
		genome[rnd::integer(nGene / 2)].mutate();
}

void Female::develop()// how does this work? 
{
	// accumulate genetic effects and compute age-specific probabilities of survival
	for (int a = 0; a < ageMax; ++a) {
		lifeHistory(survival, a) = survivalOffset; //survival offset(parms) -1:1 = 0:1 
		for (int i = 0; i < nGene; ++i)
			lifeHistory(survival, a) += genome[i](1, survival, a);
		lifeHistory(survival, a) = 0.5 * (1.0 + tanh(2.0 * lifeHistory(survival, a))); //determines age specific survival hyperbolic tangest function with a slopte of 1 at the inflection point 
	}
}

void Queen::develop()// how does this work?
{
	// accumulate genetic effects and compute age-specific probabilities of survival
	for (int a = 0; a < ageMax; ++a) {
		lifeHistory(survival, a) = survivalOffset;
		for (int i = 0; i < nGene; ++i)
			lifeHistory(survival, a) += genome[i](0, survival, a); // only difference is the 0 here to 
		lifeHistory(survival, a) = 0.5 * (1.0 + tanh(2.0 * lifeHistory(survival, a)));
	}
}

std::ostream& operator<< (std::ostream &os, const Individual &obj) // how does this work? 
{
    os << obj.lifeHistory << '\n'; 
    return os;
	
}
 
void Individual::getData(std::vector<double>& sumx) const    
{
	for (int a = 0; a < ageMax; ++a)   
		sumx[a] += lifeHistory(survival, a);  //age specific survival per individual loops through age classes
}