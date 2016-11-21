#pragma once
using namespace std;
#include <vector>

class Node
{
public:
	int x, y, demand;
};

class Vehicle
{
public:
	int load;
	vector<int> route;
};

class Chromosome
{
private:
	vector<int> customers;
public:
	int dimension;
	int capacity;
	vector<Node*>* nodes;
	vector<Vehicle*> genes;
	int fitness;
	double probability;
	double cost;

	Chromosome(vector<Node*>* n, int d, int c);
	Chromosome(Chromosome* chromosome);
	void initialise();
	void createSubroute();
	void appendSubroute();
	void evaluateFitness();
	void evaluateProbability();
	void evaluateLoad(Vehicle* vehicle);

	void pmx(Chromosome* p1, Chromosome* p2);
	void swapGenes(int i, Chromosome* p1, Chromosome* p2, Chromosome* chromosome);
};
