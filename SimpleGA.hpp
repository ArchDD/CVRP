#pragma once
using namespace std;
#include <vector>
#include "Chromosome.hpp"

class SimpleGA
{
private:
	vector<Chromosome*>* population;
	vector<Chromosome*> offsprings;
	vector<Node*>* nodes;
	int dimension;
	int capacity;

	int samples;
	int generations;
	double mutationProbability;

	Chromosome* bestSolution;
public:
	SimpleGA(vector<Chromosome*>* p, vector<Node*>* n, int d, int c);
	void generatePopulation();
	void evaluatePopulation();
	void reproduceOffspring();
	void replacePopulation();
	void evaluateSolution();
	void stepGA();
	void writeResult();
	void run();

	void pmx(Chromosome* p1, Chromosome* p2);
	void vrpCrossover(Chromosome* p1, Chromosome* p2);
	void swapGenes(int i, Chromosome* p1, Chromosome* p2, Chromosome* chromosome);
	Chromosome* swapMutation(Chromosome* ch);
};
