#include "Chromosome.hpp"
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <math.h>

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

void Chromosome::createSubroute()
{
	// Random integer in range 1 to customer size
	srand(time(0));
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
	srand(time(0));
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
	double c = 0.0;
	// Evaluate through each vehicle
	for (int i = 0; i < genes.size(); i++)
	{
		Vehicle* v = genes[i];
		// Retrieve node value from route
		for (int j = 0; j < v->route.size()-1; j++)
		{
			int n1 = v->route[j];
			int x1 = (*nodes)[n1]->x;
			int y1 = (*nodes)[n1]->y;

			int n2 = v->route[j+1];
			int x2 = (*nodes)[n2]->x;
			int y2 = (*nodes)[n2]->y;

			int x = (x1 - x2);
			int y = (y1 - y2);

			double distanceSquared = (x*x + y*y)*1.0;
			//printf("distance squared %f\n", distanceSquared);

			c += sqrt(distanceSquared);
		}
	}

	cost = c;
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
	//printf("load: %d\n", load);
}

// CROSSOVERS
// PMX crossover is a non-destructive operator
void Chromosome::pmx(Chromosome* p1, Chromosome* p2)
{
	// Uniformly select crossover points
	srand(time(0));
	int begin = (rand() % dimension-1) + 1;
	srand(time(0));
	int end = (rand() % dimension-1) + 1;
	if (begin > end)
	{
		int tmp = begin;
		begin = end;
		end = begin;
	}

	Chromosome* b1 = new Chromosome(p1);
	Chromosome* b2 = new Chromosome(p2);

	for (int i = begin; i < end; i++)
	{
		swapGenes(i, p1, p2, b1);
		swapGenes(i, p1, p2, b2);
	}
}

void Chromosome::swapGenes(int pos, Chromosome* p1, Chromosome* p2, Chromosome* chromosome)
{
	int g1, g2;
	int i1, i2, i3, i4, j1, j2, j3, j4;
	int c = 0;

	// Get genes at position
	for (int i = 0; i < p1->genes.size(); i++)
	{
		// Ignore depot
		for (int j = 1; j < p1->genes[i]->route.size()-1; j++)
		{
			if (c < pos)
				c++;
			else if (c==pos)
			{
				g1 = p1->genes[i]->route[j];
				i1 = i;
				j1 = j;
				j = p1->genes[i]->route.size()-1;
				i = p1->genes.size();
			}
		}
	}

	c = 0;
	for (int i = 0; i < p2->genes.size(); i++)
	{
		// Ignore depot
		for (int j = 1; j < p2->genes[i]->route.size()-1; j++)
		{
			if (c < pos)
				c++;
			else if (c==pos)
			{
				g2 = p2->genes[i]->route[j];
				i2 = i;
				j2 = j;
				j = p2->genes[i]->route.size()-1;
				i = p2->genes.size();
			}
		}
	}

	// Swap chromosome's g1 and g2
	for (int i = 0; i < chromosome->genes.size(); i++)
	{
		// Ignore depot
		for (int j = 1; j < chromosome->genes[i]->route.size()-1; j++)
		{
			if (chromosome->genes[i]->route[j] == g1)
			{
				i3 = i;
				j3 = j;
			}
			if (chromosome->genes[i]->route[j] == g2)
			{
				i4 = i;
				j4 = j;
			}
		}
	}

	int tmp = chromosome->genes[i3]->route[j3];
	chromosome->genes[i3]->route[j3] = chromosome->genes[i4]->route[j4];
	chromosome->genes[i4]->route[j4] = tmp;
}

// MUTATIONS
//
