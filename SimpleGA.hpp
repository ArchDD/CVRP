#pragma once
using namespace std;
#include <vector>
#include "Chromosome.hpp"

class SimpleGA
{
private:
	Chromosome* bestSolution;
	vector<Chromosome*>* population;
	vector<Node*>* nodes;
	int dimension;
	int capacity;

	int samples;
	int generations;
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
};
