#include "SimpleGA.hpp"
#include <cstdio>
#include <iostream>
#include <limits>
#include <algorithm>
#include <ctime>
#include <map>
#include <iterator>

SimpleGA::SimpleGA(vector<Node*>* n, int d, int c, int p)
{
	nodes = n;
	dimension = d;
	capacity = c;
	samples = p;
	mutationProbability = 0.2f;
}

void SimpleGA::generatePopulation()
{
	// Chromosome constructors will create random paths
	for (int i = 0; i < samples; i++)
	{
		Chromosome* chromosome = new Chromosome(nodes, dimension, capacity, true);
		population.push_back(chromosome);
	}
}

void SimpleGA::evaluatePopulation()
{
	// Evaluate population fitness and track best 
	float c = numeric_limits<float>::max();
	for (int i = 0; i < samples; i++)
	{
		Chromosome* chromosome = population[i];
		population[i]->evaluateFitness();
		if (chromosome->cost < c)
		{
			c = chromosome->cost;
			bestSolution = chromosome;
		}
	}

	/*// Uncomment for roulette wheel selection
	float totalFitness = 0.0f;
	for (int i = 0; i < samples; i++)
	{
		// Aggregate total fitness
		totalFitness += population[i]->fitness;
	}

	// Determine probability
	for (int i = 0; i < samples; i++)
	{
		Chromosome* c = population[i];
		c->probability = c->fitness / totalFitness;
	}*/
}

void SimpleGA::reproduceOffspring()
{
	// Select parents
	while(offsprings.size() < samples)
	{
		vector<Chromosome*> parents = tournamentSelection(10);
		Chromosome* c1 = parents[0];
		Chromosome* c2 = parents[1];
		c1->setArrayRepresentation();
		c2->setArrayRepresentation();
		// Reselect if same chromosome
		if (c1 != c2)
		{
			float p = (float)rand() / RAND_MAX;
			if (p < 0.6f)
				scx(c1, c2);
			else if (p < 0.7f)
				pmx(c1, c2);
			else if (p < 0.8f)
				pb(c1, c2);
			else if (p < 0.9f)
				ox(c1, c2);
			else if (p < 0.95f)
				cx(c1, c2);
			else
				vrpCrossover(c1, c2);
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
		population.back()->free();
		population.pop_back();
	}

	for (int i = samples-1; i >=0; i--)
	{
		population.push_back(offsprings.back());
		offsprings.pop_back();

		// Mutate if not best solution
		float p = (float)rand() / RAND_MAX;
		if (p < mutationProbability && population.back() != bestSolution)
		{
			p = (float)rand() / RAND_MAX;
			if (p < 0.5f)
				population.back() = insertionMutation(population.back());
			else if (p < 0.75f)
				population.back() = inversionMutation(population.back());
			else
				population.back() = swapMutation(population.back());
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
		if (filtered.count(s) == 0 && population[s] != bestSolution)
		{
			filtered[s] = s;
			population[s]->free();
			population[s] = new Chromosome(nodes, dimension, capacity, true);
		}
	}
}

void SimpleGA::start()
{
	// Generate initial population
	generatePopulation();
	evaluatePopulation();
}

void SimpleGA::step()
{
	reproduceOffspring();
	replacePopulation();
	evaluatePopulation();
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
		Chromosome* c = population[i];
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
		int v = rand() % population.size();
		Chromosome* c = population[v];
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
	int size = dimension-1;
	int *chromosome1 = p1->permutation;
	int *chromosome2 = p2->permutation;

	// Uniformly select crossover points
	int begin = (rand() % (size-1)) + 1;
	int end = (rand() % (size-1)) + 1;
	if (begin > end)
	{
		int tmp = begin;
		begin = end;
		end = begin;
	}

	Chromosome* b1 = new Chromosome(nodes, dimension, capacity, false);
	Chromosome* b2 = new Chromosome(nodes, dimension, capacity, false);

	for (int i = begin; i < end; i++)
	{
		int g1 = chromosome1[i], g2 = chromosome2[i];
		swapGenes(g1, g2, chromosome1);
		swapGenes(g1, g2, chromosome2);
	}

	repair(b1, p1, p2, chromosome1, size);
	repair(b2, p1, p2, chromosome2, size);

	offsprings.push_back(b1);
	offsprings.push_back(b2);
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
			float p = (float)rand() / RAND_MAX;
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
			Vehicle* v = new Vehicle(nodes);
			bool overCapacity = false;
			while(customers.size() > 0 && overCapacity == false)
			{
				// Select random customer left
				it = customers.begin();
				advance(it, rand() % customers.size());
				int j = it->first;
				if (v->load + (*nodes)[j]->demand <= capacity)
				{
					v->push(j);
					customers.erase(j);
				}
				else
				{
					overCapacity = true;
				}
			}
			v->push(0);
			chromosome->genes.push_back(v);
		}
	}
	offsprings.push_back(chromosome);
}

// Sequential Constructive Crossover Operator (SCX)
void SimpleGA::scx(Chromosome* p1, Chromosome* p2)
{
	Chromosome* chromosome = new Chromosome(nodes, dimension, capacity, false);
	int size = dimension-1;
	// Choose first node of parent 1
	int *chromosome1 = p1->permutation;
	int *chromosome2 = p2->permutation;
	int chromosome3[size];
	for (int i = 0; i < size; i++)
		chromosome3[i] = -1;
	int node = chromosome1[0];
	chromosome3[0] = node;
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
	int a, b, k = 1;

	while(!customers.empty())
	{
		bool found = false, selected = false;
		for (int i = 0; i < size; i++)
		{
				if (!found)
				{
					if (node == chromosome1[i])
						found = true;
				}
				else if (!selected)
				{
					if (route.count(chromosome1[i]) == 0)
					{
						a = chromosome1[i];
						selected = true;
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

		found = false, selected = false;
		for (int i = 0; i < size; i++)
		{
				if (!found)
				{
					if (node == chromosome2[i])
						found = true;
				}
				else if (!selected)
				{
					if (route.count(chromosome2[i]) == 0)
					{
						b = chromosome2[i];
						selected = true;
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

		node = c;
		chromosome3[k] = c;
		route[c] = c;
		customers.erase(c);
		k+=1;
	}
	repair(chromosome, p1, p2, chromosome3, size);
	offsprings.push_back(chromosome);
}

// A modified order crossover
void SimpleGA::ox(Chromosome* p1, Chromosome* p2)
{
	Chromosome* chromosome = new Chromosome(nodes, dimension, capacity, false);
	// Copy contents into array of chromosomes without depot
	int size = dimension-1;
	int *chromosome1 = p1->permutation;
	int *chromosome2 = p2->permutation;
	int chromosome3[size];
	// Must initialise to non-node values
	for (int i = 0; i < size; i++)
		chromosome3[i] = -1;

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
	fillNodes(chromosome3, chromosome2, size, &customers);
	repair(chromosome, p1, p2, chromosome3, size);
	offsprings.push_back(chromosome);
}

// Position-Based Crossover
void SimpleGA::pb(Chromosome* p1, Chromosome* p2)
{
	Chromosome* chromosome = new Chromosome(nodes, dimension, capacity, false);
	int size = dimension-1;
	int *chromosome1 = p1->permutation;
	int *chromosome2 = p2->permutation;
	int chromosome3[size];
	for (int i = 0; i < size; i++)
		chromosome3[i] = -1;
	// Select set of position from one parent at random and set those to offspring
	map<int, int> customers;
	int num = rand() % size;
	for (int i = 0; i < num; i++)
	{
		int pos = rand() % size;
		int n = chromosome1[pos];
		chromosome3[pos] = n;
		customers[n] = n;
	}
	// Fill in rest of nodes
	fillNodes(chromosome3, chromosome2, size, &customers);
	repair(chromosome, p1, p2, chromosome3, size);
	offsprings.push_back(chromosome);
}

// Cycle (CX) Crossover
void SimpleGA::cx(Chromosome* p1, Chromosome* p2)
{
	Chromosome* offspring1 = new Chromosome(nodes, dimension, capacity, false);
	Chromosome* offspring2 = new Chromosome(nodes, dimension, capacity, false);
	int size = dimension-1;
	int *chromosome1 = p1->permutation;
	int *chromosome2 = p2->permutation;
	int chromosome3[size];
	int chromosome4[size];
	for (int i = 0; i < size; i++)
	{
		chromosome3[i] = -1;
		chromosome4[i] = -1;
	}
	// Create maps of value to index for fast lookup
	map<int, int> m1, m2;
	for (int i = 0; i < size; i++)
	{
		int a = chromosome1[i], b = chromosome2[i];
		m1[a] = i;
		m2[b] = i;
	}

	map<int, int> customers;
	for (int i = 1; i < dimension; i++)
		customers[i] = i;
	// Find cycles
	vector< vector<int> > cycles1;
	vector< vector<int> > cycles2;
	vector< vector<int> > positions;
	while(!customers.empty())
	{
		vector<int> cycle1;
		vector<int> cycle2;
		vector<int> position;
		int i = 0, pos = -1;
		while(pos == -1)
		{
			if (customers.count(chromosome1[i]) != 0)
				pos = i;
			i++;
		}
		int start = chromosome1[pos];
		int next1 = start;
		int next2 = chromosome2[pos];
		cycle1.push_back(next1);
		cycle2.push_back(next2);
		position.push_back(pos);
		customers.erase(next1);
		while(start != next2 && !customers.empty())
		{
			pos = m1[next2];
			next1 = chromosome1[pos];
			next2 = chromosome2[pos];
			cycle1.push_back(next1);
			cycle2.push_back(next2);
			position.push_back(pos);
			customers.erase(next1);
		}
		cycles1.push_back(cycle1);
		cycles2.push_back(cycle2);
		positions.push_back(position);
	}
	for (int i = 0; i < cycles1.size(); i++)
	{
		for (int j = 0; j < cycles1[i].size(); j++)
		{
			if (i % 2 == 0)
			{
				chromosome3[positions[i][j]] = cycles1[i][j];
				chromosome4[positions[i][j]] = cycles2[i][j];
			}
			else
			{
				chromosome3[positions[i][j]] = cycles2[i][j];
				chromosome4[positions[i][j]] = cycles1[i][j];
			}
		}
	}
	// Fill in rest of nodes
	//fillNodes(chromosome3, chromosome2, size, &customers);
	repair(offspring1, p1, p2, chromosome3, size);
	//fillNodes(chromosome4, chromosome2, size, &customers);
	repair(offspring2, p1, p2, chromosome4, size);
	offsprings.push_back(offspring1);
	offsprings.push_back(offspring2);
}

void SimpleGA::swapGenes(int g1, int g2, int chromosome[])
{
	int x, y, size = dimension-1;
	// Locate positions of genes
	for (int i = 0; i < size; i++)
	{
		if (chromosome[i] == g1)
			x = i;
		if (chromosome[i] == g2)
			y = i;
	}

	int tmp = chromosome[x];
	chromosome[x] = chromosome[y];
	chromosome[y] = tmp;
}

void SimpleGA::fillNodes(int chromosome[], int parent[], int size, map<int, int>* existing)
{
	// Fill in rest of nodes with other parent
	int i = 0, j = 0;
	while (i < size)
	{

		int a = chromosome[i];
		int b = parent[j];
		// Skip existing chromosomes
		while (existing->count(a) != 0)
		{
			i+=1;
			a = chromosome[i];
		}
		while (existing->count(b) != 0)
		{
			j+=1;
			b = parent[j];
		}
		
		chromosome[i] = b;
		i+=1;
		j+=1;
	}
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
	mutation->genes[i1]->replace(j1, mutation->genes[i2]->route[j2]);
	mutation->genes[i2]->replace(j2, tmp);
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
		v->replace(j1+j, v->route[j2-j]);
		v->replace(j2-j, tmp);
	}
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
	int cus = (rand() % (dimension-1)) + 1;
	// Remove customer
	for (int i = 0; i < mutation->genes.size(); i++)
	{
		Vehicle *v = mutation->genes[i];
		for (int j = 1; j < v->route.size()-1; j++)
		{
			if (v->route[j] == cus)
			{
				v->erase(j);
				// Remove empty routes
				if (v->route.size() == 2)
				{
					delete(v);
					mutation->genes.erase(mutation->genes.begin() + i);
				}
			}
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
				v->insert(j, cus);
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

void SimpleGA::split(Chromosome* chromosome, Vehicle* v1, int i)
{
	// Split by two or three
	int p = (float)rand()/RAND_MAX;
	if (p < 0.5f)
	{
		Vehicle* v2 = new Vehicle(nodes);
		chromosome->genes.insert(chromosome->genes.begin()+i, v2);
		int s = v1->route.size()-2;
		v1->pop();
		for (int j = s; j > s/2; j--)
		{
			int n = v1->route.back();
			v1->pop();
			v2->push(n);
		}
		v1->push(0);
		v2->push(0);
	}
	else
	{
		Vehicle* v2 = new Vehicle(nodes);
		Vehicle* v3 = new Vehicle(nodes);
		chromosome->genes.insert(chromosome->genes.begin()+i, v2);
		chromosome->genes.insert(chromosome->genes.begin()+i+1, v3);
		int s = v1->route.size()-2;
		v1->pop();
		for (int j = s; j > (2*s)/3; j--)
		{
			int n = v1->route.back();
			v1->pop();
			v2->push(n);
		}
		for (int j = (2*s)/3; j > s/3; j--)
		{
			int n = v1->route.back();
			v1->pop();
			v3->push(n);
		}
		v1->push(0);
		v2->push(0);
		v3->push(0);
	}
}

// Repair strategies
void SimpleGA::repair(Chromosome* chromosome, Chromosome* p1, Chromosome* p2, int ch[], int size)
{
	float p = (float)rand() / RAND_MAX;
	if (p < 0.6f)
		greedyRepair(chromosome, ch, size);
	else if (p < 0.8f)
		randomRepair(chromosome, ch, size);
	else
		inheritanceRepair(chromosome, p1, p2, ch);
	//greedyRepair(chromosome, ch, size);
}

void SimpleGA::greedyRepair(Chromosome* chromosome, int ch[], int size)
{
	// Repair
	Vehicle *v = new Vehicle(nodes);
	chromosome->genes.push_back(v);
	for (int i = 0; i < size; i++)
	{
		int n = ch[i];
		if (v->load + (*nodes)[n]->demand <= capacity)
		{
			v->push(n);
		}
		else
		{
			v->push(0);
			v = new Vehicle(nodes);
			v->push(n);
			chromosome->genes.push_back(v);
		}
	}
	if (v->route.back() != 0)
		v->push(0);
}

void SimpleGA::randomRepair(Chromosome* chromosome, int ch[], int size)
{
	// Repair
	Vehicle *v = new Vehicle(nodes);
	chromosome->genes.push_back(v);
	int cap = capacity/2 + ((rand() % (capacity/2))+1);
	for (int i = 0; i < size; i++)
	{
		int n = ch[i];
		if (v->load + (*nodes)[n]->demand <= cap)
		{
			v->push(n);
		}
		else
		{
			cap = capacity/2 + ((rand() % (capacity/2))+1);
			v->push(0);
			v = new Vehicle(nodes);
			v->push(n);
			chromosome->genes.push_back(v);
		}
	}
	if (v->route.back() != 0)
		v->push(0);
}

void SimpleGA::inheritanceRepair(Chromosome* chromosome, Chromosome* p1, Chromosome* p2, int ch[])
{
	int size = dimension-1;
	map<int, int> customers;
	for (int i = 1; i < dimension; i++)
			customers[i] = i;

	// Look for existing routes from parents
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < p1->genes.size(); j++)
		{
			Vehicle* v = p1->genes[j];
			int e = v->route.size()-1;
			bool route = true;
			for (int k = 1; k < e; k++)
			{	
				if (ch[i+k-1] != v->route[k] || customers.count(ch[i+k-1]) == 0)
					route = false;
			}
			if (route)
			{
				chromosome->genes.push_back( new Vehicle(v) );
				for (int k = 1; k < e; k++)
					customers.erase(v->route[k]);
			}
		}
		for (int j = 0; j < p2->genes.size(); j++)
		{
			Vehicle* v = p2->genes[j];
			int e = v->route.size()-1;
			bool route = true;
			for (int k = 1; k < e; k++)
			{
				if (ch[i+k-1] != v->route[k] || customers.count(ch[i+k-1]) == 0)
					route = false;
			}
			if (route)
			{
				chromosome->genes.push_back( new Vehicle(v) );
				for (int k = 1; k < e; k++)
					customers.erase(v->route[k]);
			}
		}
	}
	// Create subroutes of remaining
	bool create = true;
	Vehicle* v = NULL;
	for (int i = 0; i < size; i++)
	{
		if (customers.count(ch[i]) != 0)
		{
			if (create)
			{
				create = false;
				v = new Vehicle(nodes);
				chromosome->genes.push_back(v);
				v->push(ch[i]);
			}
			else
			{
				v->push(ch[i]);
			}
			customers.erase(ch[i]);
		}
		else
		{
			if (create == false)
				v->push(0);
			create = true;
		}
	}
	if (v != NULL && v->route.back() != 0)
		v->push(0);

	// Split routes that are over capacity;
	int i = 0;
	while(i < chromosome->genes.size())
	{
		Vehicle* v1 = chromosome->genes[i];
		while (v1->load > capacity)
		{
			split(chromosome, v1, i);
		}
		i++;
	}
}

void SimpleGA::free()
{
	for (int i = 0; i < population.size(); i++)
	{
		population[i]->free();
	}
}