#include "SimpleGA.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <algorithm>

SimpleGA::SimpleGA(vector<Chromosome*>* p, vector<Node*>* n, int d, int c)
{
	population = p;
	nodes = n;
	dimension = d;
	capacity = c;

	generations = 1;
	samples = 5;
}

void SimpleGA::generatePopulation()
{
	// Chromosome constructors will create random paths
	for (int i = 0; i < samples; i++)
	{
		Chromosome* chromosome = new Chromosome(nodes, dimension, capacity);
		population->push_back(chromosome);
	}
}

void SimpleGA::evaluatePopulation()
{
	// Chromosomes have evaluateFitness method
	for (int i = 0; i < samples; i++)
	{
		(*population)[i]->evaluateFitness();
	}

	double totalFitness = 0.0;
	for (int i = 0; i < samples; i++)
	{
		// Aggregate total fitness
		totalFitness += (*population)[i]->fitness;
	}

	// Determine probability
	for (int i = 0; i < samples; i++)
	{
		Chromosome* c = (*population)[i];
		c->probability = c->fitness / totalFitness;
	}
}

void SimpleGA::reproduceOffspring()
{
	// Select parents
	// Roulette wheel selection
	while(offsprings.size() < dimension)
	{
		Chromosome* c1;
		Chromosome* c2;
		double p = 0.0;
		double p1 = (double)rand() / RAND_MAX;
		double p2 = (double)rand() / RAND_MAX;

		for (int i = 0; i < samples; i++)
		{
			Chromosome* c = (*population)[i];
			if (p1 > p && p1 < p+c->probability)
				c1 = c;
			if (p2 > p && p2 < p+c->probability)
				c2 = c;

			p += c->probability;
		}

		if (c1 == NULL || c2 == NULL) printf("Null parents selected\n");
		// Reselect if same chromosome
		if (c1 == c2)
			reproduceOffspring();
		else
			pmx(c1, c2);
	}

	// Keep copy of best solution
	Chromosome* best = new Chromosome(bestSolution);
	// Start freeing extra chromosomes
	while(offsprings.size() != dimension-1)
	{
		Chromosome* last = offsprings[offsprings.size()-1];
		offsprings.pop_back();
		free(last);
	}
	offsprings.push_back(best);
}


void SimpleGA::replacePopulation()
{
	for (int i = samples-1; i >= 0; i--)
	{
		free((*population)[i]);
		population->pop_back();
	}
	
	for (int i = samples-1; i >= 0; i--)
	{
		population->push_back(offsprings[i]);
		offsprings.pop_back();
	}
}

void SimpleGA::evaluateSolution()
{
	double c = numeric_limits<double>::max();
	for (int i = 0; i < population->size(); i++)
	{
		Chromosome* chromosome = (*population)[i];
		if (chromosome->cost < c)
		{
			c = chromosome->cost;
			bestSolution = chromosome;
		}
	}
}

void SimpleGA::stepGA()
{
	for (int i = 0; i < generations; i++)
	{
		srand(time(0));
		reproduceOffspring();
		evaluatePopulation();
		replacePopulation();
		evaluateSolution();
	}
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
	evaluatePopulation();
	// Loop
	stepGA();
	// End
	writeResult();
	printf("Ending simple genetic algorithm\n");
}



// CROSSOVERS
// PMX crossover is a non-destructive operator
void SimpleGA::pmx(Chromosome* p1, Chromosome* p2)
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
	offsprings.push_back(b1);
	offsprings.push_back(b2);
}

void SimpleGA::swapGenes(int pos, Chromosome* p1, Chromosome* p2, Chromosome* chromosome)
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