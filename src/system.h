#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "rng.h"
#include "stats.h"

class System {
    
    // Max radius an atom can move - this will be adjusted such that the number of rejects is 50%
    double rmax = 0.5;

    int n_atoms = 0;
    int seed;
    
    std::vector<std::string> atom_types;
    std::vector<std::vector<double>> atom_positions;
    
    Rng rng;
    Stats stats;


public:

    void allocateRAMDisk() {
        system("sudo mkdir -p memory/");
        system("sudo mount -t tmpfs -o size=2048M tmpfs memory/");
    }

    void deallocateRAMDisk() {
        system("sudo umount memory/");
        system("sudo rm -r memory/");
    }

    double energy(std::vector<std::vector<double>> coords) {
        ofstream trial("memory/trial.xsf");
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

        system("/opt/aenet/")
    }

    void pbc_wrap(double& val, int index) {
        double lat = -1.0;
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

    void trialMove() {
        int atom_index = this->rng.randint(0, this->n_atoms);
        std::vector<std::vector<double>> trial_atom_positions = this->atom_positions;
        for (int i = 0; i < 3; i ++) {
            trial_atom_positions[atom_index][i] += (2.0 * this->rng.random() - 1) * this->rmax;
            pbc_wrap(trial_atom_positions[atom_index][i], i);
        }


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
    }

    System(Config config) {
        this->config = config;
        this->seed = config.getSeed();
        this->rng = Rng(seed);
        this->allocateRAMDisk();
    }

    ~System() {
        this->deallocateRAMDisk();
    }

};

#endif