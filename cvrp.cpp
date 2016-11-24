#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include "SimpleGA.hpp"
#include "Chromosome.hpp"
#include <unistd.h>

using namespace std;

vector<Node*> nodes;
int dimension;
int capacity;
vector<Chromosome*> population;

unsigned long long getTotalSystemMemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}

int readInputFile(string fileName)
{
	printf("Total memory available: %lld\n", getTotalSystemMemory());
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

	printf("Ending preprocessing stage\n");
}

void free()
{
	for (int i = 0; i < nodes.size(); i++)
	{
		free(nodes[i]->distances);
		free(nodes[i]);
	}
}

int main()
{
	srand(time(0));
	readInputFile("fruitybun250.vrp");

	preprocess();

	SimpleGA* simpleGA = new SimpleGA(&population, &nodes, dimension, capacity);
	simpleGA->run();

	free(simpleGA);
	free();

	return 0;
}
