#include <cstdio>
#include <iostream>

using namespace std;

class SimpleGA
{
private:
public:
	SimpleGA()
	{

	};

	void run()
	{
		printf("Starting simple genetic algorithm\n");
		printf("Ending simple genetic algorithm\n");
	};
};

int main()
{
	SimpleGA* simpleGA = new SimpleGA();
	simpleGA->run();
	return 0;
}
