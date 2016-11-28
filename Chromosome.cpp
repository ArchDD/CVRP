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
	closest = -1;
}

Chromosome::Chromosome(vector<Node*>* n, int d, int c, bool init)
{
	nodes = n;
	dimension = d;
	capacity = c;
	
	if (init)
		initialise();
}

Chromosome::Chromosome(Chromosome* chromosome)
{
	nodes = chromosome->nodes;
	dimension = chromosome->dimension;
	capacity = chromosome->capacity;
	cost = chromosome->cost;

	for (int i = 0; i < chromosome->genes.size(); i++)
	{
		Vehicle* v = new Vehicle(chromosome->genes[i]);
		genes.push_back(v);
	}
}

void Chromosome::initialise()
{
	// Create copy of customers for gene pool
	for (int i = 0; i < nodes->size(); i++)
		customers.push_back(i);
	createSubroute();
	while (customers.size() > 1)
	{
		appendSubroute();
	}
	if (genes.back()->route.back() != 0)
			genes.back()->push(0);
	customers.clear();
	evaluateFitness();
}

void Chromosome::createSubroute()
{
	// Random integer in range 1 to customer size
	int i = (rand() % (customers.size()-1)) + 1;
	int n = customers[i];
	Vehicle* v = new Vehicle(nodes);
	v->push(n);
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
	if (genes[i]->load + (*nodes)[n]->demand > capacity)
	{
		genes[i]->push(0);
		createSubroute();
	}
	else
	{
		genes[i]->push(n);
		customers.erase(customers.begin()+j);
		appendSubroute();
	}
}

void Chromosome::evaluateFitness()
{
	float c = 0.0f;
	// Evaluate through each vehicle
	for (int i = 0; i < genes.size(); i++)
		c += evaluateCost(genes[i]);
	cost = c;
}

float Chromosome::evaluateCost(Vehicle* v)
{
	float c = 0.0f;
	if (v->route.size() <= 2) return cost;
	for (int i = 0; i < v->route.size()-1; i++)
	{
		c += (*nodes)[v->route[i]] -> distances[v->route[i+1]];
	}
	return c;
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

void Chromosome::setArrayRepresentation()
{
	int c = 0;
	for (int i = 0; i < genes.size(); i++)
	{
		Vehicle *v = genes[i];
		for (int j = 1; j < v->route.size()-1; j++)
		{
			permutation[c] = v->route[j];
			c+=1;
		}
	}
}

void Chromosome::free()
{
	for (int i = 0; i < genes.size();i++)
	{
		delete(genes[i]);
	}

	delete this;
}