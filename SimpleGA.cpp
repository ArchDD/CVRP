#include "SimpleGA.hpp"

SimpleGA::SimpleGA()
{

};

void SimpleGA::generatePopulation()
{

};

void SimpleGA::evaluatePopulation()
{

}

void SimpleGA::selectParents()
{

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