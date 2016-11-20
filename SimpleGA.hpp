#pragma once
using namespace std;
#include <vector>
#include "Chromosome.hpp"

class SimpleGA
{
private:
	vector<Chromosome*>* population;
	vector<Node*>* nodes;
	int dimension;
	int capacity;
public:
	SimpleGA(vector<Chromosome*>* p, vector<Node*>* n, int d, int c);
	void generatePopulation();
	void evaluatePopulation();
	void selectParents();
	void replacePopulation();
	void stepGA();
	void writeResult();
	void run();
};
