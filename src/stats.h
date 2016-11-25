#ifndef _STATS_H_
#define _STATS_H_

class Stats {
	int n_tries = 0;
	int n_accepts = 0;
public:
	void accept() {
		this->n_accepts ++;
	}

	void increment() {
		this->n_tries ++;
	}

	double rejectRatio() {
		return (double) (this->n_tries - this->n_accepts) / (this->n_tries);
	}

	int nTries() {
		return this->n_tries;
	}
};

#endif