#ifndef _STATS_H_
#define _STATS_H_

class Stats {
	int n_tries = 0;
	int n_accepts = 0;
	int r_adjust = 1000;
public:
	int rAdjust() {
		return this->r_adjust;
	}

	void accept() {
		this->n_accepts ++;
	}

	void increment() {
		this->n_tries ++;
	}

	double acceptRatio() {
		return (double) (this->n_accepts) / (this->r_adjust);
	}

	int nTries() {
		return this->n_tries;
	}

	void reset() {
		this->n_accepts = 0;
	}

	void globalReset() {
		this->n_tries = 0;
		this->n_accepts = 0;
	}

	int nAccepts() {
		return this->n_accepts;
	}
};

#endif