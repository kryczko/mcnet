#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <sstream>
#include <stdio.h>
#include <unordered_map>

#include "rng.h"
#include "stats.h"
#include "config.h"

using namespace std;

typedef std::unordered_map<std::string, int> StrToIntDict;
typedef std::unordered_map<string, vector<double>> StrToVecDict;

class System {
    
    // Max radius an atom can move - this will be adjusted such that the number of rejects is 50%
    double rmax = 2.0;
    double beta = 1.0;
    double current_energy;

    int n_atoms = 0;
    int max_trials;
    int seed;
    
    std::vector<std::string> atom_types;
    std::vector<std::vector<double>> atom_positions;
    vector<vector<double>> lj_forces;

    StrToIntDict atom_counts;
    StrToVecDict rdfs;

    Rng rng;
    Stats stats;
    Config config;

public:

    void writeXSF() {
        this->ljForces();
        char outputfile[100];
        sprintf(outputfile, "output/xsf/structure%05d.xsf", this->stats.nTries());
        system("mkdir -p output/xsf");
        ofstream output(outputfile);
        output << "# total energy = " << this->current_energy << " eV\n\nCRYSTAL\nPRIMVEC\n";
        output << to_string(this->config.xDim()) + " 0.000 0.000\n0.000 " + to_string(this->config.xDim()) + " 0.000\n0.000 0.000 " + to_string(this->config.zDim()) + "\n";
        output  << "PRIMCOORD\n";
        output << this->n_atoms << "  1\n";
        for (int i = 0; i < this->atom_types.size(); i ++) {
            output << this->atom_types[i] << " " << this->atom_positions[i][0] << " " << this->atom_positions[i][1] << " " << this->atom_positions[i][2] << "  " << this->lj_forces[i][0] << " " << this->lj_forces[i][1] << " " << this->lj_forces[i][2] << endl;
        }
        output.close();
    }

    void ljForces() {
        double epsilon = 1.0;
        this->lj_forces.resize(this->n_atoms);
        for (int i = 0; i < this->n_atoms; i ++) {
            double fx(0), fy(0), fz(0);
            for (int j = 0; j < this->n_atoms; j ++) {
                if (i != j) {
                    double dr = this->r(this->atom_positions[i], this->atom_positions[j]);
                    double dx = this->atom_positions[i][0] - this->atom_positions[j][0];
                    double dy = this->atom_positions[i][1] - this->atom_positions[j][1];
                    double dz = this->atom_positions[i][2] - this->atom_positions[j][2];
                    // std::cout << i << " " << j << " " << dr << " " << pow(1 / dr, 12) - pow(1 / dr, 6) << std::endl;
                    double val = pow(1 / dr, 6) - pow(1 / dr, 12) / ( dr * dr );
                    fx += val * dx;
                    fy += val * dy;
                    fz += val * dz;
                }
            }
            vector<double> dummy(3);
            dummy[0] = fx * 12 * epsilon;
            dummy[1] = fy * 12 * epsilon;
            dummy[2] = fz * 12 * epsilon;
            this->lj_forces[i] = dummy;
        }
    }


    void rdf(string type1, string type2) {
        for (int i = 0; i < this->n_atoms; i ++) {
            if (this->atom_types[i] == type1) {
                for (int j = 0; j < this->atom_types.size(); j ++) {
                    if (this->atom_types[j] == type2 && i != j) {
                        double dist = this->r(this->atom_positions[i], this->atom_positions[j]);
                        if (dist <= (this->config.rdfMaxDist() + this->config.rdfIncrement())) {
                            int bin = dist / this->config.rdfIncrement();
                            // cout << this->stats.nTries() << " ";
                            // cout << type1 << " " << type2 << " ";
                            // cout << "Before: " << this->rdfs[type1 + type2][bin];  
                            this->rdfs[type1 + type2][bin] ++;
                            // cout << " After: " << this->rdfs[type1 + type2][bin] << endl;
                        }
                    }
                }
            }
        }
    }

    void updateRdfs() {
        for (auto& elem : this->rdfs) {
            string name = elem.first;
            this->rdf( string(1, name[0]), string(1, name[1]));
        }
        ofstream output("output/rdfs.dat");
        output << "#  ";
        int n_bins = this->config.rdfMaxDist() / this->config.rdfIncrement();
        for (auto& elem : this->rdfs) {
            output << elem.first << "  ";
        }
        output << "\n";
        for (int i = 1; i < n_bins + 1; i ++) {
            double r_3 = pow((i+1)*this->config.rdfIncrement(), 3) - pow((i)*this->config.rdfIncrement(), 3);
            double shell = 4.0 * 3.14159 * r_3 / 3.0;
            double Vcell = this->config.xDim() * this->config.yDim() * this->config.zDim();
            output << (i)*this->config.rdfIncrement() + this->config.rdfIncrement()*0.5 << "  ";
            for (auto& elem : this->rdfs) {
                int atom_count1 = this->atom_counts[string(1, elem.first[0])];
                int atom_count2 = this->atom_counts[string(1, elem.first[1])];
                output <<  (Vcell * elem.second[i]) / (this->stats.nTries() * shell * atom_count1 * atom_count2)  << "  ";
            }
            output << "\n";
        } 
        output.close();
    }
    
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
        double epsilon = 1.0; 
        double sigma = 1e-10;
        double J_to_eV = 1.6e-19; 
        for (int i = 0; i < coords.size(); i ++) {
            for (int j = 0; j < coords.size(); j ++) {
                if (i != j) {
                    double dr = this->r(coords[i], coords[j]);
                    // std::cout << i << " " << j << " " << dr << " " << pow(1 / dr, 12) - pow(1 / dr, 6) << std::endl;
                    toten += pow(1 / dr, 12) - pow(1 / dr, 6);
                }
            }
        }
        // std::cout << toten * 2.0 * epsilon << "\n";
        // std::string stuff;
        // std::cin >> stuff;
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
        double ratio = abs(val / lat);
        if (val > lat) {
            val -= lat * ratio;
        } else if (val < 0.0) {
            val += lat * ratio;
        }
    }

    void run() {
        ofstream output("output/trajectory.xyz");
        while (this->stats.nTries() < this->max_trials) {
            this->trialMove();
            if (this->stats.nTries() % this->stats.rAdjust() == 0) {
                if (this->stats.rejectRatio() > 0.5) {
                    this->rmax *= 1.05;
                } else {
                    this->rmax *= 0.95;
                }
                std::cout << "Iteration: " << this->stats.nTries() << " Reject ratio: " << this->stats.rejectRatio() << " Rmax: " << this->rmax << "  Energy: " << this->current_energy << std::endl;
                this->stats.reset();
            }
            output << this->n_atoms << "\n\n";
            for (int i = 0; i < this->n_atoms; i ++) {
                output << this->atom_types[i] << "  " << this->atom_positions[i][0] << "  " << this->atom_positions[i][1] << "  " << this->atom_positions[i][2] << endl;
            }
        this->updateRdfs();
        }
        output.close();
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
        this->writeXSF();
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
            aname = "";
            xyzfile >> aname >> x >> y >> z;
            if (aname != "") {
                this->atom_types.push_back(aname);
                std::vector<double> coords{x, y, z};
                this->atom_positions.push_back(coords);
            }
        }
        
        for (int i = 0; i < this->atom_types.size(); i ++) {
            string atom_name = this->atom_types[i];
            StrToIntDict::const_iterator got = this->atom_counts.find(atom_name);
            if (got == this->atom_counts.end()) {
                this->atom_counts.insert({atom_name, 1});
            } else {
                this->atom_counts.at(atom_name) ++;
            }
        }

        std::cout << "INFO: Read in XYZ file.\n";
        this->current_energy = this->lennardJonesEnergy(this->atom_positions);
    }

    void initDistros() {
        // init the rdfs
        system("mkdir -p output");
        for (auto& elem : this->atom_counts) {
            for (auto& elem2 : this->atom_counts) {
                vector<double> distro((int) this->config.rdfMaxDist() / this->config.rdfIncrement(), 0.0);
                StrToVecDict::const_iterator got = this->rdfs.find(elem2.first + elem.first);
                if (got == this->rdfs.end()) {
                    this->rdfs.insert({elem.first + elem2.first, distro});
                }
            }
        }
    }

    System(Config config) {
        this->config = config;
        this->seed = config.getSeed();
        this->rng = Rng(seed);
        this->max_trials = config.getMaxTrials();
        this->allocateRAMDisk();
        this->readXyzFile(config.xyzFile());
        this->initDistros();
    }

    ~System() {
        this->deallocateRAMDisk();
    }

};

#endif