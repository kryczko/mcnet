#ifndef _RNG_H_
#define _RNG_H_ 

#include <random>
#include <chrono>


class Rng {
	private:
		std::mt19937 gen;
	public:
		double random() {
			return (double) this->gen() / this->gen.max();
		}

		int randint(int min, int max) {
			return (int) (max - min) * this->random() + min;
		}

		Rng(int seed=std::chrono::system_clock::now().time_since_epoch().count()) {
			std::mt19937 gen(seed);
			this->gen = gen;
		}
};

#endif