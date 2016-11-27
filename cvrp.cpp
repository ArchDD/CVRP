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
		}
	}

	printf("Ending preprocessing stage\n\n");
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


	int n = 8;
	SimpleGA* simpleGA[n];

	printf("Starting %d genetic algorithms.\n", n);
	for (int i = 0; i < n; i++)
	{
		simpleGA[i] = new SimpleGA(&nodes, dimension, capacity);
		simpleGA[i]->start();
	}

	printf("Parallel processing...\n");
	double timeLimit = 1.0 * 10.0, t = 0.0;
	int iterations = 0, steps = 20;
	while(t < timeLimit)
	{
		#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			simpleGA[i]->step(steps);
		}

		iterations += 200;
		gettimeofday(&timstr,NULL);
		toc = timstr.tv_sec+(timstr.tv_usec/1000000.0);
		t = toc-tic;
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
	
	gettimeofday(&timstr,NULL);
	toc = timstr.tv_sec+(timstr.tv_usec/1000000.0);
	writeResult();
	printf("Elapsed time:\t\t\t\t%.2f (s)\n", toc-tic);
	printf("Number of parallel GAs ran:\t\t%d\n", n);
	printf("Iterations completed:\t\t\t%d\n", iterations);
	printf("Best solution cost:\t\t\t%.3f\n", bestSolution->evaluatePreciseCost());
	
	free();
	printf("[END EXECUTION]\n\n");
	return 0;
}
