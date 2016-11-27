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
protected:
	vector<Node*>* nodes;
public:
	int load;
	vector<int> route;
	Vehicle(vector<Node*>* n)
	{
		nodes = n;
		load = 0;
		route.push_back(0);
	}

	Vehicle(Vehicle* v)
	{
		nodes = v->nodes;
		load = v->load;
		route = v->route;
	}

	void push(int n)
	{
		route.push_back(n);
		load += (*nodes)[n]->demand;
	}

	void pop()
	{
		load-= (*nodes)[route.back()]->demand;
		route.pop_back();
	}

	void insert(int i, int n)
	{
		route.insert(route.begin()+i, n); load += (*nodes)[n]->demand;
	}

	void erase(int i)
	{
		load-= (*nodes)[route[i]]->demand;
		route.erase(route.begin()+i);
	}

	void replace(int i, int n)
	{
		load -= (*nodes)[route[i]]->demand;
		route[i] = n;
		load += (*nodes)[n]->demand;
	}
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
	int permutation[249];
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
	void setArrayRepresentation();
	void free();
};
