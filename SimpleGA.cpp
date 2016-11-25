#include "SimpleGA.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <limits>
#include <algorithm>
#include <ctime>
#include <map>
#include <iterator>

SimpleGA::SimpleGA(vector<Chromosome*>* p, vector<Node*>* n, int d, int c)
{
	population = p;
	nodes = n;
	dimension = d;
	capacity = c;
	samples = 100;
	mutationProbability = 0.2f;
}

void SimpleGA::generatePopulation()
{
	// Chromosome constructors will create random paths
	for (int i = 0; i < samples; i++)
	{
		Chromosome* chromosome = new Chromosome(nodes, dimension, capacity, true);
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

	float totalFitness = 0.0f;
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
	while(offsprings.size() < samples)
	{
		//vector<Chromosome*> parents = rouletteSelection();
		vector<Chromosome*> parents = tournamentSelection(10);
		Chromosome* c1 = parents[0];
		Chromosome* c2 = parents[1];

		if (c1 == NULL || c2 == NULL) printf("Null parents selected\n");
		// Reselect if same chromosome
		if (c1 != c2)
		{
			/*float p = (float)rand() / RAND_MAX;
			if (p < 0.7f)
				pmx(c1, c2);
			else if (p < 0.9f)
				scx(c1, c2);
			else
				vrpCrossover(c1, c2);*/
			ox(c1, c2);
			//pmx(c1, c2);
		}
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
		float p = (float)rand() / RAND_MAX;
		if (p < mutationProbability && population->back() != bestSolution)
		{
			p = (float)rand() / RAND_MAX;
			/*if (p < 0.5f)
				population->back() = swapMutation(population->back());
			else if (p < 0.75f)
				population->back() = inversionMutation(population->back());
			else
				population->back() = insertionMutation(population->back());*/
			//population->back() = inversionMutation(population->back());
			p = (float)rand() / RAND_MAX;
			if (p < 0.5f)
				population->back() = swapMutation(population->back());
			else
				population->back() = inversionMutation(population->back());
		}
	}
}

void SimpleGA::evaluateSolution()
{
	float c = numeric_limits<float>::max();
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

// Resets a portion of chromosomes from random position once a while
void SimpleGA::filtration()
{
	map <int, int> filtered;
	while(filtered.size() < samples/10)
	{
		int s = rand() % samples;
		if (filtered.count(s) == 0 && (*population)[s] != bestSolution)
		{
			filtered[s] = s;
			(*population)[s]->free();
			(*population)[s] = new Chromosome(nodes, dimension, capacity, true);
		}
	}
}

void SimpleGA::stepGA()
{
	clock_t t1 = clock(), t2 = clock(), t3 = clock();
	float f = 0.0f;
	int i = 0, batch = 10;
	float timeLimit = 10.0f * 60.0f, batchTime = 0.0f;

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
		if (i % 200 == 0)
			filtration();
	}
	t1 = clock() - t1;
	printf("Time taken: %f seconds\n", ((float)t1)/CLOCKS_PER_SEC);
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
	double cost = bestSolution->evaluatePreciseCost();
	file << "cost " << cost <<endl;

	printf("Best Cost: %f\n", cost);

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

// SELECTIONS
vector<Chromosome*> SimpleGA::rouletteSelection()
{
	Chromosome* c1;
	Chromosome* c2;
	float p = 0.0f;
	float p1 = (float)rand() / RAND_MAX;
	float p2 = (float)rand() / RAND_MAX;

	for (int i = 0; i < samples; i++)
	{
		Chromosome* c = (*population)[i];
		if (p1 > p && p1 < p+c->probability)
			c1 = c;
		if (p2 > p && p2 < p+c->probability)
			c2 = c;

		p += c->probability;
	}
	vector<Chromosome*> parents;
	parents.push_back(c1);
	parents.push_back(c2);
	return parents;
}

vector<Chromosome*> SimpleGA::tournamentSelection(int n)
{
	vector<Chromosome*> selection;
	
	while (selection.size() < n)
	{
		int v = rand() % population->size();
		Chromosome* c = (*population)[v];
		bool contains = false;
		for (int i = 0; i < selection.size(); i++)
			if (selection[i] == c)
				contains = true;
		if (!contains)
			selection.push_back(c);
	}
	Chromosome* c1;
	Chromosome* c2;
	int a, b;
	float cost1 = numeric_limits<float>::max();
	float cost2 = cost1;

	for (int i = 0; i < n; i++)
	{
		if (selection[i]->cost < cost2)
		{
			b = i;
			cost2 = selection[i]->cost;
		}
		if (selection[i]->cost < cost1)
		{
			int tmp = a, tmpCost = cost1;
			a = i;
			cost1 = selection[i]->cost;
			b = tmp;
			cost2 = tmpCost;
		}
	}
	vector<Chromosome*> parents;
	parents.push_back(selection[a]);
	parents.push_back(selection[b]);
	return parents;
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
	else
		b1->free();

	acceptable = true;
	for (int i = 0; i < b2->genes.size(); i++)
		if (b2->genes[i]->load > capacity)
			acceptable = false;
	if (acceptable)
		offsprings.push_back(b2);
	else
		b2->free();
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

void SimpleGA::vrpCrossover(Chromosome* p1, Chromosome* p2)
{
	Chromosome* chromosome = new Chromosome(nodes, dimension, capacity, false);

	// Create selection of customers
	map<int, int> customers;
	for (int i = 0; i < dimension; i++)
		customers[i] = i;
	// Remove depot
	customers.erase(0);

	while (customers.size() > 0)
	{
		// Select random customer
		map<int, int>::iterator it = customers.begin();
		advance(it, rand() % customers.size() );
		int c = it->first;

		// Pick subroutes including customer from parents
		Vehicle* v1;
		Vehicle* v2;
		Vehicle* v3;
		for (int i = 0; i < p1->genes.size(); i++)
		{
			for (int j = 0; j < p1->genes[i]->route.size(); j++)
			{
				if (p1->genes[i]->route[j] == c)
				{
					v1 = p1->genes[i];
				}

			}
		}
		for (int i = 0; i < p2->genes.size(); i++)
		{
			for (int j = 0; j < p2->genes[i]->route.size(); j++)
			{
				if (p2->genes[i]->route[j] == c)
				{
					v2 = p2->genes[i];
				}

			}
		}

		bool subset1 = true, subset2 = true;
		// Determine whether subroutes are subsets of remaining customers
		for (int i = 0; i < v1->route.size(); i++)
		{
		
			if (customers.count( v1->route[i] ) == 0)
				subset1 = false;
		}
		for (int i = 0; i < v2->route.size(); i++)
		{
			if (customers.count( v2->route[i] ) == 0)
				subset2 = false;
		}

		int s1 = 1, s2 = 1;
		if (subset1 && subset2)
		{
			float f1 = (float)s1;
			float f2 = (float)s2;
			float p1 = f2 / (f1 + f2);
			//float p2 = f1 / (f1 + f2);
			float p = (float)rand() / RAND_MAX;
			//printf("p %f p1 %f\n", p , p1);
			if (p < p1)
			{
				v3 = v1;
				s1+=1;
			}
			else
			{
				v3 = v2;
				s2+=1;
			}
			chromosome->genes.push_back ( new Vehicle(v3) );
			// Remove subroute from customers
			for (int i = 1; i < v3->route.size()-1; i++)
			{
				printf("v3 %x size %d i %d\n", v3, v3->route.size(), i);
				int n = v3->route[i];
				customers.erase(n);
			}
		}
		else if (subset1)
		{
			v3 = v1;
			s1+=1;
			chromosome->genes.push_back ( new Vehicle(v3) );
			// Remove subroute from customers
			for (int i = 1; i < v3->route.size()-1; i++)
			{
				customers.erase(v3->route[i]);
			}
		}
		else if (subset2)
		{
			v3 = v2;
			s2+=1;
			chromosome->genes.push_back ( new Vehicle(v3) );
			// Remove subroute from customers
			for (int i = 1; i < v3->route.size()-1; i++)
			{
				customers.erase(v3->route[i]);
			}
		}
		else if (subset1 == false || subset2 == false)
		{
			Vehicle* v = new Vehicle();
			v->route.push_back(0);
			v->route.push_back(0);
			bool overCapacity = false;
			while(customers.size() > 0 && overCapacity == false)
			{
				// Select random customer left
				it = customers.begin();
				advance(it, rand() % customers.size());
				int j = it->first;

				chromosome->evaluateLoad(v);
				if (v->load + (*nodes)[j]->demand <= capacity)
				{
					v->route.pop_back();
					v->route.push_back(j);
					v->route.push_back(0);
					customers.erase(j);
				}
				else
				{
					overCapacity = true;
				}
			}
			chromosome->genes.push_back(v);
		}
	}
	chromosome->evaluateFitness();
	for (int i = 0; i< chromosome->genes.size(); i++)
	{
		chromosome->evaluateLoad(chromosome->genes[i]);
		if (chromosome->genes[i]->load > capacity)
			for (int j = 0; j < chromosome->genes[i]->route.size(); j++)
				printf("%d->",chromosome->genes[i]->route[j]);
	}
	offsprings.push_back(chromosome);
}

// Sequential Constructive Crossover Operator (SCX)
void SimpleGA::scx(Chromosome* p1, Chromosome* p2)
{
	Chromosome* chromosome = new Chromosome(nodes, dimension, capacity, false);

	// Choose first node randomly
	int node = (rand() % dimension-1) + 1;

	// Create route map to quickly determine if nodes are already chosen
	map<int, int> route;
	route[0] = 0;
	route[node] = node;

	// Create selection of customers
	map<int, int> customers;
	for (int i = 0; i < dimension; i++)
		customers[i] = i;
	// Remove depot and first random node
	customers.erase(0);
	customers.erase(node);
	int a, b;

	Vehicle* v = new Vehicle();
	v->route.push_back(0);
	v->route.push_back(node);
	v->route.push_back(0);
	chromosome->genes.push_back(v);

	while(!customers.empty())
	{
		bool found = false, selected = false;
		for (int i = 0; i < p1->genes.size(); i++)
		{
			for (int j = 1; j < p1->genes[i]->route.size()-1; j++)
			{
				if (!found)
				{
					if (node == p1->genes[i]->route[j])
						found = true;
				}
				else if (!selected)
				{
					if (route.count(p1->genes[i]->route[j]) == 0)
					{
						a = p1->genes[i]->route[j];
						selected = true;
					}
				}
			}
		}
		if (selected == false)
		{
			for (int i = 1; i < dimension; i++)
			{
				if (route.count(i) == 0)
					a = i;
			}
		}

		found = false;
		selected = false;
		for (int i = 0; i < p2->genes.size(); i++)
		{
			for (int j = 1; j < p2->genes[i]->route.size()-1; j++)
			{
				if (!found)
				{
					if (node == p2->genes[i]->route[j])
						found = true;
				}
				else if (!selected)
				{
					if (route.count(p2->genes[i]->route[j]) == 0)
					{
						b = p2->genes[i]->route[j];
						selected = true;
					}
				}
			}
		}
		if (selected == false)
		{
			for (int i = 1; i < dimension; i++)
			{
				if (route.count(i) == 0)
					b = i;
			}
		}

		int c;
		if ((*nodes)[node]->distances[a] < (*nodes)[node]->distances[b])
			c = a;
		else
			c = b;

		chromosome->evaluateLoad(v);
		if (v->load + (*nodes)[c]->demand <= capacity)
		{
			v->route.pop_back();
			v->route.push_back(c);
			v->route.push_back(0);
		}
		else
		{
			v = new Vehicle();
			v->route.push_back(0);
			v->route.push_back(c);
			v->route.push_back(0);
			chromosome->genes.push_back(v);
		}
		route[c] = c;
		customers.erase(c);
	}
	offsprings.push_back(chromosome);
}

// A modified order crossover
void SimpleGA::ox(Chromosome* p1, Chromosome* p2)
{
	p1->evaluateFitness();
	Chromosome* chromosome = new Chromosome(nodes, dimension, capacity, false);
	// Copy contents into array of chromosomes without depot
	int size = dimension-1;
	int* chromosome1 = (int*)malloc(sizeof(float) * size);
	int* chromosome2 = (int*)malloc(sizeof(float) * size);
	int* chromosome3 = (int*)malloc(sizeof(float) * size);
	int c1 = 0, c2 = 0;

	for (int i = 0; i < p1->genes.size(); i++)
	{
		Vehicle *v = p1->genes[i];
		for (int j = 1; j < v->route.size()-1; j++)
		{
			chromosome1[c1] = v->route[j];
			c1+=1;
		}
	}
	for (int i = 0; i < p2->genes.size(); i++)
	{
		Vehicle *v = p2->genes[i];
		for (int j = 1; j < v->route.size()-1; j++)
		{
			chromosome2[c2] = v->route[j];
			c2+=1;
		}
	}
	// Select a subtour from p1
	int start = rand() % size;
	int end = rand() % size;
	if (start > end)
	{
		int tmp = start;
		start = end;
		end = tmp;
	}
	// Create a map for quick duplicate check
	map<int, int> customers;
	// Copy subtour to offspring
	for (int i = start; i < end; i++)
	{
		int n = chromosome1[i];
		chromosome3[i] = n;
		customers[n] = n;
	}
	// Fill in rest of nodes with other parent
	int k = 0, l = 0;
	while (k < start)
	{
		int val = l;
		int n = chromosome2[l];
		if (customers.count(n) != 0)
		{
			l+=1;
			n = chromosome2[l];
		}
		else
		{
			chromosome3[k] = n;
			k+=1;
			l+=1;
		}
	}
	k = end;
	while (k < size)
	{
		int n = chromosome2[l];
		if (customers.count(n) != 0)
		{
			l+=1;
			n = chromosome2[l];
		}
		else
		{
			chromosome3[k] = n;
			k+=1;
			l+=1;
		}
	}

	// Repair
	Vehicle *v = new Vehicle();
	v->route.push_back(0);
	v->route.push_back(0);
	chromosome->genes.push_back(v);
	for (int i = 0; i < size; i++)
	{
		chromosome->evaluateLoad(v);
		int n = chromosome3[i];
		if (v->load + (*nodes)[n]->demand <= capacity)
		{
			v->route.pop_back();
			v->route.push_back(n);
			v->route.push_back(0);
		}
		else
		{
			v = new Vehicle();
			v->route.push_back(0);
			v->route.push_back(n);
			v->route.push_back(0);
			chromosome->genes.push_back(v);
		}
	}
	offsprings.push_back(chromosome);
}
// Position-Based Crossover
/*Chromosome* SimpleGA::positionCrossover(Chromosome* p1, Chromosome* p2)
{
	// Select set of position from one parent at random
	// Create child by copying nodes on at positions
	// Remove set of nodes from other parent
	// Fill in rest of nodes
}*/

// Cycle (CX) Crossover


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

Chromosome* SimpleGA::inversionMutation(Chromosome* ch)
{
	Chromosome* mutation = new Chromosome(ch);
	// Select a random subroute
	int i = rand() % mutation->genes.size();
	Vehicle *v = mutation->genes[i];
	int j1 = (rand() % (v->route.size()-2))+1;
	int j2 = (rand() % (v->route.size()-2))+1;
	if (j2 < j1)
	{
		int tmp = j1;
		j1 = j2;
		j2 = tmp;
	}
	// Invert subtour of subroute
	for (int j = 0; j < ((j1+j2)/2)-j1; j++)
	{
		int tmp = v->route[j1+j];
		v->route[j1+j] = v->route[j2-j];
		v->route[j2-j] = tmp;
	}
	mutation->evaluateLoad(v);
	if (v->load <= capacity)
	{
		ch->free();
		return mutation;
	}
	mutation->free();
	return ch;
}

Chromosome* SimpleGA::insertionMutation(Chromosome* ch)
{
	Chromosome* mutation = new Chromosome(ch);
	// Selct random customer and position
	int cus = (rand() % dimension-1) + 1;
	// Remove customer
	for (int i = 0; i < mutation->genes.size(); i++)
	{
		Vehicle *v = mutation->genes[i];
		for (int j = 1; j < v->route.size()-1; j++)
		{
			if (v->route[j]== cus)
				v->route.erase(v->route.begin() + j);
		}
	}	
	int pos = (rand() % dimension-1) + 1;
	int c = 0;
	Vehicle *v;
	for (int i = 0; i < mutation->genes.size(); i++)
	{
		v = mutation->genes[i];
		for (int j = 1; j < v->route.size()-1; j++)
		{
			if (c == pos)
			{
				v->route.insert(v->route.begin()+j, cus);
				mutation->evaluateLoad(v);
				if (v->load <= capacity)
				{
					ch->free();
					return mutation;
				}
				else
				{
					mutation->free();
					return ch;
				}
			}
			c+=1;
		}
	}
	mutation->free();
	return ch;
}
