#include "SimpleGA.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <algorithm>
#include <ctime>

SimpleGA::SimpleGA(vector<Chromosome*>* p, vector<Node*>* n, int d, int c)
{
	population = p;
	nodes = n;
	dimension = d;
	capacity = c;

	generations = 5;
	samples = 250;
	mutationProbability = 0.1;
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
		for (int j = 0; j < (*population)[i]->genes.size(); j++)
			(*population)[i]->evaluateLoad( (*population)[i]->genes[j] );
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
	while(offsprings.size() < samples)
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
		if (c1 != c2)
			pmx(c1, c2);
	}

	// Free extra chromosomes
	while(offsprings.size() > samples)
	{
		Chromosome* last = offsprings[offsprings.size()-1];
		offsprings.pop_back();
		last->free();
	}
	// Keep copy of best solution
	if (bestSolution != NULL)
	{
		Chromosome* best = new Chromosome(bestSolution);
		bestSolution = best;
		Chromosome* last = offsprings.back();
		offsprings.pop_back();
		last->free();
		offsprings.push_back(best);
	}

}


void SimpleGA::replacePopulation()
{
	for (int i = samples-1; i >=0; i--)
	{
		population->back()->free();
		population->pop_back();
	}

	for (int i = samples-1; i >=0; i--)
	{
		population->push_back(offsprings.back());
		offsprings.pop_back();

		// Mutate if not best solution
		double p = (double)rand() / RAND_MAX;
		if (p < mutationProbability && population->back() != bestSolution)
			population->back() = swapMutation(population->back());
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
	/*clock_t t = clock();
	for (int i = 0; i < generations; i++)
	{
		//clock_t d = clock();
		reproduceOffspring();
		//d = clock() - d; printf("reproduceOffspring: %f seconds\n", ((double)d)/CLOCKS_PER_SEC); d = clock();

		evaluatePopulation();
		//d = clock() - d; printf("evaluatePopulation: %f seconds\n", ((double)d)/CLOCKS_PER_SEC); d = clock();

		replacePopulation();
		//d = clock() - d; printf("replacePopulation: %f seconds\n", ((double)d)/CLOCKS_PER_SEC); d = clock();

		evaluateSolution();
		//d = clock() - d; printf("evaluateSolution: %f seconds\n", ((double)d)/CLOCKS_PER_SEC); d = clock();
	}
	t = clock() - t;
	printf("Time taken: %f seconds\n", ((double)t)/CLOCKS_PER_SEC);*/

	clock_t t1 = clock(), t2 = clock(), t3 = clock();
	float f = 0.0f;
	int i = 0, batch = 10;
	float timeLimit = 1.0f * 10.0f, batchTime = 0.0f;

	while (f < timeLimit)
	{
		reproduceOffspring();
		replacePopulation();
		evaluatePopulation();
		evaluateSolution();

		if (i % batch == 0)
		{
			t2 = clock() - t1;
			f = ((float)t2 / CLOCKS_PER_SEC);
			t3 = clock() - t3;
			batchTime = ((float)t3 / CLOCKS_PER_SEC);
		}
		i++;
	}
	t1 = clock() - t1;
	printf("Time taken: %f seconds\n", ((double)t1)/CLOCKS_PER_SEC);
	printf("Iterations completed: %d\n", i);
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
	int begin = (rand() % dimension-1) + 1;
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

	// Push only evaluated load is acceptable
	bool acceptable = true;
	for (int i = 0; i < b1->genes.size(); i++)
		if (b1->genes[i]->load > capacity)
			acceptable = false;
	if (acceptable)
		offsprings.push_back(b1);

	acceptable = true;
	for (int i = 0; i < b2->genes.size(); i++)
		if (b2->genes[i]->load > capacity)
			acceptable = false;
	if (acceptable)
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
			else if (c>=pos)
			{
				g1 = p1->genes[i]->route[j];
				i1 = i;
				j1 = j;
				j = p1->genes[i]->route.size()-1;
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
			else if (c>=pos)
			{
				g2 = p2->genes[i]->route[j];
				i2 = i;
				j2 = j;
				j = p2->genes[i]->route.size()-1;
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

	//Evaluate new load
	chromosome->evaluateLoad(chromosome->genes[i3]);
	chromosome->evaluateLoad(chromosome->genes[i4]);
}

// MUTATIONS
Chromosome* SimpleGA::swapMutation(Chromosome* ch)
{
	// Select two random customers
	int i1 = rand() % ch->genes.size();
	int j1 = (rand() % (ch->genes[i1]->route.size()-2)) + 1;
	int i2 = rand() % ch->genes.size();
	int j2 = (rand() % (ch->genes[i2]->route.size()-2)) + 1;

	// Trying to swap the customers will invert or interchange path
	Chromosome* mutation = new Chromosome(ch);
	int tmp = mutation->genes[i1]->route[j1];
	mutation->genes[i1]->route[j1] = mutation->genes[i2]->route[j2];
	mutation->genes[i2]->route[j2] = tmp;

	mutation->evaluateLoad(mutation->genes[i1]);
	mutation->evaluateLoad(mutation->genes[i2]);
	if (mutation->genes[i1]->load <= capacity && mutation->genes[i2]->load <= capacity)
	{
		ch->free();
		return mutation;
	}
	mutation->free();
	return ch;
}

