#include "SimpleGA.hpp"
#include <cstdio>
#include <iostream>

SimpleGA::SimpleGA(vector<Chromosome*>* p, vector<Node*>* n, int d, int c)
{
	population = p;
	nodes = n;
	dimension = d;
	capacity = c;
}

void SimpleGA::generatePopulation()
{
	// Chromosome constructors will create random paths
	for (int i = 0; i < 1; i++)
	{
		Chromosome* chromosome = new Chromosome(nodes, dimension, capacity);
		population->push_back(chromosome);
	}
}

void SimpleGA::evaluatePopulation()
{
	// Chromosomes have evaluateFitness method
}

void SimpleGA::selectParents()
{
	// roulette
}

void SimpleGA::replacePopulation()
{

}

void SimpleGA::stepGA()
{
	// Evaluate population
	// Select parents
	// Reproduce using genetic operators
	// Replace population with offspring
}

void SimpleGA::writeResult()
{

}

void SimpleGA::run()
{
	printf("Starting simple genetic algorithm\n");
	// Generate initial population
	generatePopulation();
	// Loop
	stepGA();
	// End
	writeResult();
	printf("Ending simple genetic algorithm\n");
}