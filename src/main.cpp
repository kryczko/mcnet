#include <iostream>
#include "config.h"
#include "system.h"

int main(int argc, char ** argv) {
	Config config;
	config.printConfig();

	System system(config);
	// std::cout << "Min energy is: " << system.minEnergy() << " Max energy is: " << system.maxEnergy() << std::endl;
	system.run();
	return 0;
}