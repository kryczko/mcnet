#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <sstream>

#include "rng.h"
#include "stats.h"
#include "config.h"

class System {
    
    // Max radius an atom can move - this will be adjusted such that the number of rejects is 50%
    double rmax = 0.5;
    double beta = 0.0005;
    double current_energy;

    int n_atoms = 0;
    int r_adjust = 10;
    int max_trials;
    int seed;
    
    std::vector<std::string> atom_types;
    std::vector<std::vector<double>> atom_positions;
    
    Rng rng;
    Stats stats;
    Config config;


public:

    void allocateRAMDisk() {
        system("sudo mkdir -p memory/");
        system("sudo mount -t tmpfs -o size=2048M tmpfs memory/");
    }

    void deallocateRAMDisk() {
        system("sudo umount memory/");
        system("sudo rm -r memory/");
    }

    void split(const std::string &s, char delim, std::vector<std::string> &elems) {
        std::stringstream ss;
        ss.str(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
    }

    double aenetEnergy(std::vector<std::vector<double>> coords) {
        std::ofstream trial("memory/trial.xsf");
        trial << "CRYSTAL\nPRIMVEC\n";
        trial << "  " << this->config.xDim() << "  0.000 0.000\n  0.000  " 
              << this->config.yDim() << "  0.000\n  0.000  0.000  " << this->config.xDim() << std::endl;
        trial << "PRIMCOORD\n";
        trial << this->n_atoms << " 1\n";
        for (int i = 0; i < this->n_atoms; i ++) {
            trial << this->atom_types[i] << " ";
            for (auto val : this->atom_positions[i]) {
                trial << val << "  ";
            }
            trial << std::endl;
        }
        trial.close();

        system("cd aenet && /opt/aenet-1.0.0/bin/predict.x-1.0.0-gfortran_mpi predict.in ../memory/trial.xsf > ../memory/aenet_out.txt 2>&1 && cd ..");
        std::ifstream energy_file("memory/aenet_out.txt");
        std::string stuff;
        double energy;
        for (std::string line; getline(energy_file, line);) {
            if (line.find("Total energy") != std::string::npos) {
                // we found the total energy
                std::vector<std::string> items;
                split(line, ' ', items);
                energy = stof(items[items.size() - 2]);
            }
        }
        return energy;
    }

    double pbcWrapDist(double input) {
		int i  = input;
		if (std::abs(input) >= 0.5) {
			if (input > 0.0) {i += 1;}
			if (input < 0.0) {i -= 1;}
		}
		return i;
	}

    double r(std::vector<double> a, std::vector<double> b) {
    	double dx = a[0] - b[0];
    	double dy = a[1] - b[1];
    	double dz = a[2] - b[2];
    	dx -= this->config.xDim() * this->pbcWrapDist(dx / this->config.xDim());
    	dy -= this->config.yDim() * this->pbcWrapDist(dy / this->config.yDim());
    	dz -= this->config.zDim() * this->pbcWrapDist(dz / this->config.zDim());
    	return sqrt( dx*dx + dy*dy + dz*dz );
    }

    double lennardJonesEnergy(std::vector<std::vector<double>> coords) {
    	// v(r) = 4 epsilon (( sigma / r) ^12 - (sigma / r) ^ 6)
    	double toten = 0.0;
    	double epsilon = 1e4; // boltzmann constant times 300 K
    	double sigma = 1e-10;
    	double J_to_eV = 1.6e-19; 
    	for (int i = 0; i < coords.size(); i ++) {
    		for (int j = 0; j < coords.size(); j ++) {
    			if (i != j) {
    				double dr = this->r(coords[i], coords[j]);
    				toten += pow(1 / dr, 12) - pow(1 / dr, 6);
    			}

    		}
    	}
    	return toten * 2.0 * epsilon;

    }

    void pbc_wrap(double& val, int index) {
        double lat;
        if (index == 0) {
            // x-axis
            lat = this->config.xDim();
        } else if (index == 1) {
            // y-axis
            lat = this->config.yDim();
        } else if (index == 2) {
            // z-axis
            lat = this->config.zDim();
        } else {
            std::cout << "ERROR: The index being passed into the pbc_wrap function is incorrect. Exiting...\n";
            exit(-1);
        }

        if (val > lat) {
            val -= lat;
        } else if (val < 0.0) {
            val += lat;
        }
    }

    void run() {
    	while (this->stats.nTries() < this->max_trials) {
    		this->trialMove();
    		if (this->stats.nTries() % this->r_adjust == 0) {
    			if (this->stats.rejectRatio() > 0.5) {
    				this->rmax *= 1.05;
    			} else {
    				this->rmax *= 0.95;
    			}
    			std::cout << "Iteration: " << this->stats.nTries() << " Reject ratio: " << this->stats.rejectRatio() << " Rmax: " << this->rmax << "  Energy: " << this->current_energy << std::endl;

    		}
    	}
    }

    void trialMove() {
        int atom_index = this->rng.randint(0, this->n_atoms);
        std::vector<std::vector<double>> trial_atom_positions = this->atom_positions;
        for (int i = 0; i < 3; i ++) {
            trial_atom_positions[atom_index][i] += (2.0 * this->rng.random() - 1) * this->rmax;
            pbc_wrap(trial_atom_positions[atom_index][i], i);
        }

        double next_energy = this->lennardJonesEnergy(trial_atom_positions);
        double delta_e = next_energy - this->current_energy;
        double delta_e_beta = delta_e * this->beta;
        if (delta_e_beta < 75.0) {
        	if (delta_e_beta < 0.0) {
	        	this->current_energy = next_energy;
	        	this->atom_positions = trial_atom_positions;
	        	this->stats.accept();
        	} else if (exp( - delta_e_beta ) > this->rng.random() ) {
	        	this->current_energy = next_energy;
	        	this->atom_positions = trial_atom_positions;
	        	this->stats.accept();
        	}
        }
        this->stats.increment();
    }

    void readXyzFile(std::string filename) {
        std::ifstream xyzfile(filename.c_str());
        if (!xyzfile) {
            std::cout << "ERROR: We could not open your XYZ file. Please check the input file 'input.yaml' and verify your XYZ file. Exiting...\n";
            exit(-1);
        }
        xyzfile >> n_atoms;
        std::string aname;
        double x, y, z;
        while (!xyzfile.eof()) {
            xyzfile >> aname >> x >> y >> z;
            this->atom_types.push_back(aname);
            std::vector<double> coords{x, y, z};
            this->atom_positions.push_back(coords);
        }
        std::cout << "INFO: Read in XYZ file.\n";
        this->current_energy = this->lennardJonesEnergy(this->atom_positions);
    }

    System(Config config) {
        this->config = config;
        this->seed = config.getSeed();
        this->rng = Rng(seed);
        this->max_trials = config.getMaxTrials();
        this->allocateRAMDisk();
        this->readXyzFile(config.xyzFile());
    }

    ~System() {
        this->deallocateRAMDisk();
    }

};

#endif