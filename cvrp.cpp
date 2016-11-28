#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include <fstream>
#include <string>
#include <sstream>

#include "SimpleGA.hpp"
#include "Chromosome.hpp"

using namespace std;

Chromosome* bestSolution;
vector<Node*> nodes;
int dimension;
int capacity;

unsigned long long getTotalSystemMemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}

int readInputFile(string fileName)
{
	ifstream f(fileName.c_str());
	if (f.is_open())
	{
		string s;
		int x, y, d;
		// Read format until nodes
		f >> s >> s >> dimension >> s >> s >> capacity >> s;

		cout << "Reading: " << s << endl;

		for (int i = 0; i < dimension; i++)
		{
			f >> x >> x >> y;
			Node* n = new Node();
			n->x = x;
			n->y = y;
			nodes.push_back(n);
		}
		f >> s;
		cout << "Reading: " << s << endl;

		for (int i = 0; i < dimension; i++)
		{
			f >> d >> d;
			nodes[i]->demand = d;
		}

		f.close();
	}
	else
	{
		cout << "File name " << fileName << " not found." << endl;
	}
}

void preprocess()
{
	printf("Begin preprocessing stage\n");
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			int x1 = nodes[i]->x;
			int y1 = nodes[i]->y;
			int x2 = nodes[j]->x;
			int y2 = nodes[j]->y;

			int x = (x1 - x2);
			int y = (y1 - y2);
			float distanceSquared = (x*x + y*y)*1.0;
			float distance = sqrt(distanceSquared);
			nodes[i]->distances[j] = distance;
			if (nodes[i]->closest == -1 && i != j)
				nodes[i]->closest = j;
			else if (nodes[i]->distances[j] < nodes[i]->distances[nodes[i]->closest] && i != j)
				nodes[i]->closest = j;
		}
	}

	printf("Ending preprocessing stage\n\n");
}

void postprocess()
{
	printf("\nBegin post-processing stage\n");
	printf("Cost before:\t\t\t\t%.3f\n", bestSolution->evaluatePreciseCost());

	double tic, toc, timeLimit = 30.0f;
	struct timeval timstr;
	gettimeofday(&timstr,NULL);
	tic = timstr.tv_sec+(timstr.tv_usec/1000000.0);

	#pragma omp parallel for schedule(auto)
	for (int i = 0; i < bestSolution->genes.size(); i++)
	{
		Vehicle v = Vehicle(bestSolution->genes[i]);
		float c1 = bestSolution->evaluateCost(bestSolution->genes[i]);
		while( next_permutation(v.route.begin()+1, v.route.end()-1) && toc-tic < timeLimit)
		{
			float c2 = bestSolution->evaluateCost(&v);
			if (c2 < c1)
			{
				bestSolution->genes[i]->route = v.route;
				c1 = bestSolution->evaluateCost(bestSolution->genes[i]);
			}
			gettimeofday(&timstr,NULL);
			toc = timstr.tv_sec+(timstr.tv_usec/1000000.0);
		}
	}
	
	if (toc-tic >= timeLimit)
		printf("Time limit exceeded.\n");

	printf("Cost after:\t\t\t\t%.3f\n", bestSolution->evaluatePreciseCost());
	printf("Ending post-processing stage\n\n");
}

void writeResult()
{
	printf("Writing best result to file.\n\n");

	ofstream file;
	file.setf(ios::fixed,ios::floatfield);
	file.precision(3);

	file.open("best-solution.txt", ofstream::out | ofstream::trunc);
	file << "login dd13282 61545\n";
	file << "name Dillon Keith Diep\n";
	file << "algorithm Genetic Algorithm\n";
	double cost = bestSolution->evaluatePreciseCost();
	file << "cost " << cost <<endl;

	string n;
	stringstream convert;

	for (int i = 0; i < bestSolution->genes.size(); i++)
	{
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

void free()
{
	for (int i = 0; i < nodes.size(); i++)
	{
		free(nodes[i]->distances);
		free(nodes[i]);
	}
	delete(bestSolution);
}

int main()
{
	double tic, toc;
	struct timeval timstr;
	gettimeofday(&timstr,NULL);
	tic = timstr.tv_sec+(timstr.tv_usec/1000000.0);

	printf("\n[BEGIN EXECUTION]\n");
	printf("Total memory available: %lld\n", getTotalSystemMemory());
	srand(time(0));
	readInputFile("fruitybun250.vrp");

	preprocess();

	bestSolution = new Chromosome(&nodes, dimension, capacity, true);

	int NUM_THREADS = omp_get_max_threads();
	omp_set_num_threads(NUM_THREADS);
	printf("Threads available: %d\n", NUM_THREADS);

	int n = NUM_THREADS, p = 256;
	SimpleGA* simpleGA[n];

	printf("Starting %d genetic algorithm instances.\n", n);
	for (int i = 0; i < n; i++)
	{
		simpleGA[i] = new SimpleGA(&nodes, dimension, capacity, p);
		simpleGA[i]->start();
	}

	printf("Parallel processing...\n");
	double timeLimit = 119.0 * 60.0, t = 0.0, timed = 60.0;
	int iterations = 0, steps = 100;

	while(t < timeLimit)
	{
		#pragma omp parallel for schedule(auto)
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < steps; j++)
				simpleGA[i]->step();
			if (iterations % (steps*5) == 0)
				simpleGA[i]->filtration();
		}

		iterations += steps;

		gettimeofday(&timstr,NULL);
		toc = timstr.tv_sec+(timstr.tv_usec/1000000.0);
		t = toc-tic;
		if (t > timed)
		{
			printf("%.0f minute has passed.\n", t/60.0);
			timed+=60.0;
		}
	}

	for (int i = 0; i < n; i++)
	{
		if (simpleGA[i]->bestSolution->cost < bestSolution->cost)
		{
			bestSolution->free();
			bestSolution = new Chromosome(simpleGA[i]->bestSolution);
		}
		simpleGA[i]->free();
	}
	printf("End genetic algorithms.\n");

	postprocess();
	
	gettimeofday(&timstr,NULL);
	toc = timstr.tv_sec+(timstr.tv_usec/1000000.0);
	writeResult();
	printf("Elapsed time:\t\t\t\t%.2f (s)\n", toc-tic);
	printf("Number of parallel GAs ran:\t\t%d\n", n);
	printf("Population size:\t\t\t%d\n", p);
	printf("Iterations completed:\t\t\t%d\n", iterations);
	printf("Best solution cost:\t\t\t%.3f\n", bestSolution->evaluatePreciseCost());
	
	free();
	printf("[END EXECUTION]\n\n");
	return 0;
}
