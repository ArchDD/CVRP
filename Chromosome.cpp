#include "Chromosome.hpp"
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <math.h>

Node::Node()
{ 
	distances = (float*)malloc(sizeof(float) * 250);
}

Chromosome::Chromosome(vector<Node*>* n, int d, int c)
{
	nodes = n;
	dimension = d;
	capacity = c;

	initialise();
}

Chromosome::Chromosome(Chromosome* chromosome)
{
	nodes = chromosome->nodes;
	dimension = chromosome->dimension;
	capacity = chromosome->capacity;
	fitness = chromosome->fitness;
	probability = chromosome->probability;
	cost = chromosome->cost;

	for (int i = 0; i < chromosome->genes.size(); i++)
	{
		Vehicle* v = new Vehicle();
		v->route = chromosome->genes[i]->route;
		v->load = chromosome->genes[i]->load;
		genes.push_back(v);
	}
}

void Chromosome::initialise()
{
	// Create copy of customers for gene pool
	for (int i = 0; i < nodes->size(); i++)
		customers.push_back(i);

	while (customers.size() > 1)
	{
		createSubroute();
		appendSubroute();
	}

	evaluateFitness();
}

void Chromosome::clearRoute()
{
	for (int i = 0; i < genes.size(); i++)
	{
		delete(genes[i]);
	}
	genes.clear();
}

void Chromosome::createSubroute()
{
	// Random integer in range 1 to customer size
	int i = (rand() % (customers.size()-1)) + 1;
	int n = customers[i];
	Vehicle* v = new Vehicle();
	v->route.push_back(0);
	v->route.push_back(n);
	v->route.push_back(0);
	genes.push_back(v);
	customers.erase(customers.begin()+i);
}


void Chromosome::appendSubroute()
{
	// Return if only depot left
	if (customers.size() == 1) return;

	int i = genes.size()-1;
	int j = (rand() % (customers.size()-1)) + 1;
	int n = customers[j];

	// Update vehicle load then try to fit in last subroute
	evaluateLoad(genes[i]);
	if (genes[i]->load + (*nodes)[n]->demand > capacity)
	{
		createSubroute();
	}
	else
	{
		// Remove last 0, then add j and 0 again (O(1) operators)
		genes[i]->route.pop_back();
		genes[i]->route.push_back(n);
		genes[i]->route.push_back(0);
		customers.erase(customers.begin()+j);
		appendSubroute();
	}
}

void Chromosome::evaluateFitness()
{
	float c = 0.0f;
	// Evaluate through each vehicle
	for (int i = 0; i < genes.size(); i++)
	{
		Vehicle* v = genes[i];
		// Retrieve node value from route
		for (int j = 0; j < v->route.size()-1; j++)
		{
			c += (*nodes)[v->route[j]] -> distances[v->route[j+1]];
		}
	}

	cost = c;
	fitness = 100000.0f/cost;
}

void Chromosome::evaluateLoad(Vehicle* vehicle)
{
	int load = 0;
	for (int i = 1; i < vehicle->route.size()-1; i++)
	{
		int n = vehicle->route[i];
		load += (*nodes)[n]->demand;
	}
	vehicle->load = load;
}

double Chromosome::evaluatePreciseCost()
{
	double c = 0.0;
	// Evaluate through each vehicle
	for (int i = 0; i < genes.size(); i++)
	{
		Vehicle* v = genes[i];
		// Retrieve node value from route
		for (int j = 0; j < v->route.size()-1; j++)
		{
			int n1 = v->route[j], n2 = v->route[j+1];
			int x1 = (*nodes)[n1]->x;
			int y1 = (*nodes)[n1]->y;
			int x2 = (*nodes)[n2]->x;
			int y2 = (*nodes)[n2]->y;

			int x = (x1 - x2);
			int y = (y1 - y2);
			double distanceSquared = (x*x + y*y)*1.0;
			c += sqrt(distanceSquared);
		}
	}
	return c;
}

bool Chromosome::containsGene(int n)
{
	for (int i = 0; i < genes.size(); i++)
	{
		for (int j = 1; j < genes[i]->route.size()-1; j++)
		{
			if (n == genes[i]->route[j])
				return true;
		}
	}
	return false;
}

void Chromosome::free()
{
	for (int i = 0; i < genes.size();i++)
	{
		delete(genes[i]);
	}

	delete this;
}