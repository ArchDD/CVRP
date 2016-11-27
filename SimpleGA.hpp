#pragma once
using namespace std;
#include <vector>
#include <map>
#include "Chromosome.hpp"

class SimpleGA
{
private:
	vector<Chromosome*> population;
	vector<Chromosome*> offsprings;
	vector<Node*>* nodes;
	int dimension;
	int capacity;

	int samples;
	int generations;
	float mutationProbability;
public:
	Chromosome* bestSolution;
	SimpleGA(vector<Node*>* n, int d, int c, int p);
	void generatePopulation();
	void evaluatePopulation();
	void reproduceOffspring();
	void replacePopulation();
	void filtration();
	void start();
	void step();


	vector<Chromosome*> rouletteSelection();
	vector<Chromosome*> tournamentSelection(int n);
	void pmx(Chromosome* p1, Chromosome* p2);
	void vrpCrossover(Chromosome* p1, Chromosome* p2);
	void scx(Chromosome* p1, Chromosome* p2);
	void ox(Chromosome* p1, Chromosome* p2);
	void pb(Chromosome* p1, Chromosome* p2);
	void ob(Chromosome* p1, Chromosome* p2);
	void cx(Chromosome* p1, Chromosome* p2);
	void swapGenes(int g1, int g2, int chromosome[]);
	void fillNodes(int chromosome[], int parent[], int size, map<int, int>* existing);

	Chromosome* swapMutation(Chromosome* ch);
	Chromosome* inversionMutation(Chromosome* ch);
	Chromosome* insertionMutation(Chromosome* ch);
	Chromosome* shuffleMutation(Chromosome* ch);
	Chromosome* splitMutation(Chromosome* ch);
	void split(Chromosome* chromosome, Vehicle* v1, int i);

	void repair(Chromosome* chromosome, Chromosome* p1, Chromosome* p2, int ch[], int size);
	void greedyRepair(Chromosome* chromosome, int ch[], int size);
	void inheritanceRepair(Chromosome* chromosome, Chromosome* p1, Chromosome* p2, int ch[]);
	void randomRepair(Chromosome* chromosome, int ch[], int size);

	void free();
};
