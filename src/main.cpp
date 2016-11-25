#include <iostream>
#include "config.h"
#include "system.h"

int main(int argc, char ** argv) {
	Config config;
	config.printConfig();

	System system(config);
	system.run();
	return 0;
}