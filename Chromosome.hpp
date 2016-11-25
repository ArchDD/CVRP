#pragma once
using namespace std;
#include <vector>

class Node
{
public:
	int x, y, demand;
	//vector<Node*> closest;
	float *distances;
	Node();
};

class Vehicle
{
public:
	int load;
	vector<int> route;
	Vehicle() {};
	Vehicle(Vehicle* v) { load = v->load; route = v->route; };
};

class Chromosome
{
private:
	vector<int> customers;
protected:
	int dimension;
	int capacity;
	vector<Node*>* nodes;
public:
	vector<Vehicle*> genes;
	float fitness;
	float probability;
	float cost;

	Chromosome(vector<Node*>* n, int d, int c, bool init);
	Chromosome(Chromosome* chromosome);
	void initialise();
	void createSubroute();
	void appendSubroute();
	void evaluateFitness();
	double evaluatePreciseCost();
	void evaluateProbability();
	void evaluateLoad(Vehicle* vehicle);
	int* getArrayRepresentation(int* chromosome, int size);
	void free();
};
