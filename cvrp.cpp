#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "SimpleGA.hpp"
#include "Chromosome.hpp"
using namespace std;

vector<Node*> nodes;
unsigned int dimension;
unsigned int capacity;
unsigned int iterations = 5;
unsigned int trucks = 1;
unsigned int cost = 0;
vector<Chromosome*> population;

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

void free()
{
	for (int i = 0; i < nodes.size(); i++)
	{
		free(nodes[i]);
	}
}

int main()
{
	readInputFile("fruitybun250.vrp");

	SimpleGA* simpleGA = new SimpleGA(&population, &nodes, dimension, capacity);
	simpleGA->run();

	free();

	return 0;
}
