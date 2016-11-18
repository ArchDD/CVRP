#include "cvrp.hpp"
using namespace std;

unsigned int iterations = 5;

int main()
{
	SimpleGA* simpleGA = new SimpleGA();
	simpleGA->run();
	return 0;
}
