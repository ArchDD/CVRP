#include "SimpleGA.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

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
		bestSolution = chromosome;
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
	printf("Writing result\n");

	ofstream file;
	file.setf(ios::fixed,ios::floatfield);
	file.precision(3);

	file.open("best-solution.txt", ofstream::out | ofstream::trunc);
	file << "login dd13282 61545\n";
	file << "name Dillon Keith Diep\n";
	file << "algorithm Genetic Algorithm\n";
	file << "cost " << bestSolution->cost <<endl;

	printf("Best Cost: %f\n", bestSolution->cost);
	string n;
	stringstream convert;

	for (int i = 0; i < bestSolution->genes.size(); i++)
	{
		//printf("Vehicle: %d\n", i+1);
		//printf("Route Length: %d\n", bestSolution->genes[i]->route.size());
		string s = "";
		for (int j = 0; j < bestSolution->genes[i]->route.size(); j++)
		{
			convert << (bestSolution->genes[i]->route[j]+1);
			n = convert.str();
			convert.str("");
			convert.clear();
			s = s + n + "->";
		}
		s = s.substr(0, s.size()-2);
		file << s + "\n";
	}

	file.close();
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