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
	int dimension, capacity;
	vector<Node*>* nodes;
	vector<int> customers;
public:
	vector<Vehicle*> genes;
	int fitness;
	float probability;

	Chromosome(vector<Node*>* n, int d, int c);
	void initialise();
	void createSubroute();
	void appendSubroute();
	void evaluateFitness();
	void evaluateProbability();
};
