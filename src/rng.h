#ifndef _RNG_H_
#define _RNG_H_ 

#include <random>


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

		Rng(int seed=1) {
			std::mt19937 gen(seed);
			this->gen = gen;
		}
};

#endif