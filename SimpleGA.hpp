#pragma once

#include <cstdio>
#include <iostream>
using namespace std;

class SimpleGA
{
private:
public:
	SimpleGA();

	void generatePopulation();

	void evaluatePopulation();

	void selectParents();

	void replacePopulation();

	void stepGA();

	void writeResult();

	void run();
};
