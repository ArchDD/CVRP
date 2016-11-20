#include "Chromosome.hpp"
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

Chromosome::Chromosome(vector<Node*>* n, int d, int c)
{
	nodes = n;
	dimension = d;
	capacity = c;

	initialise();
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
}

void Chromosome::createSubroute()
{
	// Random integer in range 1 to customer size
	int i = (rand() % customers.size()-1) + 1;
	Vehicle* v = new Vehicle();
	v->route.push_back(0);
	v->route.push_back(i);
	v->route.push_back(0);
	genes.push_back(v);
	customers.erase(customers.begin()+i);
}

void Chromosome::appendSubroute()
{
	// Return if only depot left
	if (customers.size() == 1) return;

	int i = genes.size()-1;
	int j = (rand() % customers.size()-1) + 1;
	// Try to fit in last subroute
	if (genes[i]->load + (*nodes)[j]->demand > capacity)
	{
		createSubroute();
	}
	else
	{
		// Remove last 0, then add j and 0 again (O(1) operators)
		genes[i]->route.pop_back();
		genes[i]->route.push_back(j);
		genes[i]->route.push_back(0);
		customers.erase(customers.begin()+j);
		appendSubroute();
	}
}