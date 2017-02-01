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
#include <iomanip>
#include <unordered_map>
#include <set>
#include <thread>
#include <future>

#include "rng.h"
#include "stats.h"
#include "config.h"

using namespace std;

typedef std::unordered_map<std::string, int> StrToIntDict;
typedef std::unordered_map<string, vector<double>> StrToVecDict;

class System {
    
    // Max radius an atom can move - this will be adjusted such that the number of rejects is 50%
    double rmax = 0.01;
    double beta = 35.1652183807007;
    // double beta = 11.60452205665357;
    double current_energy;
    double min_energy = 1e30;
    double max_energy = -1e30;
    double avg_distance;

    int n_atoms = 0;
    int n_outputs = 0;
    int max_trials;
    int seed;
    
    vector<double> msd;
    vector<double> energies;
    vector<double> energy_distro;
    std::vector<std::string> atom_types;
    std::vector<std::vector<double>> atom_positions;
    vector<vector<double>> original_positions;
    vector<vector<double>> lj_forces;
    vector<vector<double>> aenet_forces;

    StrToIntDict atom_counts;
    StrToVecDict rdfs;


    Rng rng;
    Stats stats;
    Config config;

    ofstream summary_output;
    ofstream msd_output;

public:

    double min(vector<double> list) {
        double min_val = 1e30;
        for (auto elem : list) {
            if (elem < min_val) {
                min_val = elem;
            }
        }
        return min_val;
    }

    double max(vector<double> list) {
        double max_val = -1e30;
        for (auto elem : list) {
            if (elem > max_val) {
                max_val = elem;
            }
        }
        return max_val;
    }

    double sum(vector<double> list) {
        double sum_val = 0;
        for (auto elem : list) {
            sum_val += elem;
        }
        return sum_val;
    }

    void updateMsd() {
        double dist_sum = 0;
        for (int i = 0; i < this->n_atoms; i ++) {
            double dx = this->original_positions[i][0] - this->atom_positions[i][0];
            double dy = this->original_positions[i][1] - this->atom_positions[i][1];
            double dz = this->original_positions[i][2] - this->atom_positions[i][2];
            dist_sum += dx * dx + dy * dy + dz * dz;
        }
        this->msd.push_back(dist_sum / this->n_atoms);
        this->msd_output << this->msd[this->msd.size() - 1] << "\n";
        this->msd_output.flush();
    }

    void updateEnergyDistro() {
        int bins = 100;
        this->energy_distro.clear();
        this->energy_distro.resize(bins, 0.0);
        // double min_val = this->min(this->energies);
        // double diff = this->max(this->energies) - min_val;
        double min_val = -1.0;
        double diff = 2.0;
        if (!diff) {
            return;
        }
        for (int i = 0; i < this->energies.size(); i ++) {
            int bin = bins * (this->energies[i] + abs(min_val)) / diff;
            if (bin >= 0 && bin < bins) {
            // cout << this->energies[i] << " " << min_val << " " << diff << " " << bin << endl;
            this->energy_distro[bin] ++;
        }
        }
        ofstream output("output/energy_distro.dat");
        double incr = diff / bins;
        for (int i = 0; i < bins; i ++) {
            output << min_val + i*incr << "\t" << this->energy_distro[i]  / this->sum(this->energy_distro) << endl;
            output << min_val + (i+1)*incr << "\t" << this->energy_distro[i] / this->sum(this->energy_distro) << endl;
        }
        output.close();
    }

    void writeXSF() {
        this->ljForces();
        char outputfile[100];
        sprintf(outputfile, "output/xsf/structure%010d.xsf", this->stats.nTries());
        system("mkdir -p output/xsf");
        ofstream output(outputfile);
        output << "# total energy = " << this->current_energy << " eV\n\nCRYSTAL\nPRIMVEC\n";
        output << to_string(this->config.xDim()) + " 0.000 0.000\n0.000 " + to_string(this->config.xDim()) + " 0.000\n0.000 0.000 " + to_string(this->config.zDim()) + "\n";
        output  << "PRIMCOORD\n";
        output << this->n_atoms << "  1\n";
        for (int i = 0; i < this->atom_types.size(); i ++) {
            output << fixed << setprecision(5) << this->atom_types[i] << " " << this->atom_positions[i][0] << " " << this->atom_positions[i][1] << " " << this->atom_positions[i][2] << "  " << this->lj_forces[i][0] << " " << this->lj_forces[i][1] << " " << this->lj_forces[i][2] << endl;
        }
        output.close();
    }

    void writeSummary() {
        if (this->stats.nTries() == 0) {
            this->summary_output << "# Iteration\tenergy\tavg distance\n\n";
        }
        this->summary_output << this->stats.nTries() << "\t" << this->current_energy << "\t" << this->avg_distance << endl;
    }

    void ljForces() {
        double epsilon = 1.0;
        this->lj_forces.resize(this->n_atoms);
        for (int i = 0; i < this->n_atoms; i ++) {
            double fx(0), fy(0), fz(0);
            for (int j = 0; j < this->n_atoms; j ++) {
                if (i != j) {
                    double dr = this->r(this->atom_positions[i], this->atom_positions[j]);
                    if (dr < 2.5) {
                        double dx = this->atom_positions[i][0] - this->atom_positions[j][0];
                        double dy = this->atom_positions[i][1] - this->atom_positions[j][1];
                        double dz = this->atom_positions[i][2] - this->atom_positions[j][2];
                        dx -= this->config.xDim() * this->pbcWrapDist(dx / this->config.xDim());
                        dy -= this->config.xDim() * this->pbcWrapDist(dy / this->config.xDim());
                        dz -= this->config.xDim() * this->pbcWrapDist(dz / this->config.xDim());
                        // std::cout << i << " " << j << " " << dr << " " << pow(1 / dr, 12) - pow(1 / dr, 6) << std::endl;
                        double val = pow(0.69 / dr, 6) - pow(0.69 / dr, 12) / ( dr * dr );
                        fx += val * dx;
                        fy += val * dy;
                        fz += val * dz;
                    }
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
            this->rdf( string(name.size(), name[0]), string(name.size(), name[1]));
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
                string type1 = string(1, elem.first[0]);
                string type2 = string(1, elem.first[1]);
                int atom_count1 = this->atom_counts[type1];
                int atom_count2 = this->atom_counts[type2];
                if (type1 == type2) {
                    atom_count1 --;
                }
                output <<  (Vcell * elem.second[i]) / (this->n_outputs * shell * atom_count1 * atom_count2)  << "  ";
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

    double cutOffEnergy(double r) {
        double max = 1e10;
        double min = 1e2;
        return - (max - min) * r + max;
    }

    void getForces() {
        this->aenet_forces.clear();
        ifstream aenet_input("memory/aenet_out.txt");
        string stuff, name;
        bool collect = false;
        double x, y, z, fx, fy, fz;
        int n_atoms = 0;
        while(!aenet_input.eof()) {
            if (collect == false) {
                aenet_input >> stuff;
            }
            if (stuff == "--------------------------------------------------------------------------------------") {
                collect = true;
            }
            if (collect) {
                aenet_input >> name >> x >> y >> z >> fx >> fy >> fz;
                vector<double> f = {fx, fy, fz};
                this->aenet_forces.push_back(f);
                n_atoms ++;
            }
            if (n_atoms == this->n_atoms) {
                collect = false;
            }
        }
    }

    void split(const std::string &s, char delim, std::vector<std::string> &elems) {
        std::stringstream ss;
        ss.str(s);
        std::string item;
        while (std::getline(ss, item, delim)) {
            elems.push_back(item);
        }
    }

    vector<int> nearestNeighbors(int id, vector<vector<double>> coords) {
        vector<int> nns;
        double cut_off = 2.5;
        for (int i = 0; i < coords.size(); i ++) {
            if (i != id) {
                if (this->r(coords[id], coords[i]) < cut_off) {
                    nns.push_back(i);
                }
            }
        }
        return nns;
    }

    double aenetEnergy(std::vector<std::vector<double>> coords, int thread=0) {
        std::ofstream trial(("memory/trial" + to_string(thread) + ".xsf").c_str());
        trial << "CRYSTAL\nPRIMVEC\n";
        trial << "  " << this->config.xDim() << "  0.000 0.000\n  0.000  " 
              << this->config.yDim() << "  0.000\n  0.000  0.000  " << this->config.xDim() << std::endl;
        trial << "PRIMCOORD\n";
        trial << coords.size() << " 1\n";
        for (int i = 0; i < coords.size(); i ++) {
            trial << this->atom_types[i] << "  ";
            for (auto val : coords[i]) {
                trial << val << "  ";
            }
            trial << std::endl;
        }
        trial.close();
        system(("cd aenet && /opt/aenet-1.0.0/bin/predict.x-1.0.0-gfortran_mpi predict.in ../memory/trial" + to_string(thread) + ".xsf > ../memory/aenet_out" + to_string(thread) + ".txt 2>&1 && cd ..").c_str());

        std::ifstream energy_file(("memory/aenet_out" + to_string(thread) + ".txt").c_str());
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
        // system("rm -f memory/*");
        return energy;
    }

    double pbcWrapDist(double input) {
        int i  = input;
        if (std::abs(input - i) >= 0.5) {
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
        double epsilon = 2.0; 
        double shift = (pow(1 / 2.5, 12) - pow(1 / 2.5, 6));
        for (int i = 0; i < coords.size(); i ++) {
            // vector<int> nns = this->nearestNeighbors(i, coords);
            for (int j = 0; j < coords.size(); j ++) {
                if (i != j) {
                    double dr = this->r(coords[i], coords[j]);
                    // std::cout << i << " " << j << " " << dr << " " << pow(1 / dr, 12) - pow(1 / dr, 6) << std::endl;
                    if (dr < 2.5) { 
                        toten += (pow(0.69 / dr, 12) - pow(0.69 / dr, 6) - shift);
                    }
                }
            }
        }
        // std::cout << toten * 2.0 * epsilon << "\n";
        // std::string stuff;
        // std::cin >> stuff;
        return toten * 2.0 * epsilon;

    }

    double localLennardJonesEnergy(std::vector<vector<double>> coords, int index) {
        double toten = 0.0;
        double epsilon = 2.0;
        double shift = (pow(0.69 / 2.5, 12) - pow(0.69 / 2.5, 6));
        for (int i = 0; i < coords.size(); i ++) {
            if (i != index) {
                double dr = this->r(coords[index], coords[i]);
                if (dr < 2.5) {
                    toten += (pow(0.69 / dr, 12) - pow(0.69 / dr, 6) - shift);
                } 
            }
        }
        return 4.0 * epsilon * toten;
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

    double maxEnergy() {
        double max_energy = -1e30;
        while (this->stats.nTries() < this->max_trials) {
            this->classicMetropolis(true);
            if (this->stats.nTries() % this->stats.rAdjust() == 0) {
                if (this->stats.acceptRatio() > 0.5) {
                    this->rmax *= 1.05;
                } else {
                    this->rmax *= 0.95;
                }
                std::cout << "Max search - Iteration: " << this->stats.nTries() << " Accept ratio: " << this->stats.acceptRatio() << " Rmax: " << this->rmax << "  Energy: " << this->current_energy << std::endl;
                this->stats.reset();
            }
            if (this->current_energy > max_energy) {
                max_energy = this->current_energy;
            }
        }
        this->reset();
        return max_energy;
    }

    void reset() {
        this->rmax = 0.1;
        this->stats.globalReset();
        this->readXyzFile(this->config.xyzFile());
    }

    double minEnergy() {
        double min_energy = 1e30;
        while (this->stats.nTries() < this->max_trials) {
            this->classicMetropolis();
            if (this->stats.nTries() % this->stats.rAdjust() == 0) {
                if (this->stats.acceptRatio() > 0.5) {
                    this->rmax *= 1.05;
                } else {
                    this->rmax *= 0.95;
                }
                std::cout << "Min search - Iteration: " << this->stats.nTries() << " Accept ratio: " << this->stats.acceptRatio() << " Rmax: " << this->rmax << "  Energy: " << this->current_energy << std::endl;
                this->stats.reset();
            }
            if (this->current_energy < min_energy) {
                min_energy = this->current_energy;
            }
        }
        this->reset();
        return min_energy;
    }

    void adjustR() {
        if (this->stats.acceptRatio() > 0.5) {
            this->rmax *= 1.05;
        } else {
            this->rmax *= 0.95;
        }
    }

    void run() {
        ofstream output("output/trajectory.xyz");
        while (this->stats.nTries() < this->max_trials) {
            this->classicMetropolis();
            if (this->stats.nTries() % this->stats.rAdjust() == 0) {
                // this->adjustR();
                std::cout << "Iteration: " << this->stats.nTries() << " Accept ratio: " << this->stats.acceptRatio() << " Rmax: " << this->rmax << "  Energy: " << this->current_energy / this->n_atoms << std::endl;
                output << this->n_atoms << "\n\n";
                for (int i = 0; i < this->n_atoms; i ++) {
                    output << this->atom_types[i] << "  " << this->atom_positions[i][0] << "  " << this->atom_positions[i][1] << "  " << this->atom_positions[i][2] << endl;
                }
                cout.flush();
                output.flush();
                this->writeSummary();
                this->updateRdfs();
                // this->updateMsd();
                this->n_outputs ++;
                this->stats.reset();
                // this->current_energy = this->lennardJonesEnergy(this->atom_positions);
                this->writeXSF();
            }
            // this->updateEnergyDistro();
        }
        cout << "Min energy: " << this->min_energy << "  Max energy: " << this->max_energy << endl;
        output.close();
    }

    void runTargetedEnergy() {
        double min_energy = -1.0;
        double max_energy = 1.0;
        int bins = 100;
        int incr = this->max_trials / bins;
        double step = ( max_energy - min_energy ) / bins;
        while (this->stats.nTries() < this->max_trials) {
            if (this->stats.nTries() % incr == 0) {
                this->rmax = 0.1;
                this->readXyzFile(this->config.xyzFile());
            }
            int which = this->stats.nTries() / incr;
            double e_target = min_energy + which * step;
            this->targetedEnergyMetropolis(e_target);
            if (this->stats.nTries() % this->stats.rAdjust() == 0) {
                if (this->stats.acceptRatio() > 0.5) {
                    this->rmax *= 1.05;
                } else {
                    this->rmax *= 0.95;
                }
                std::cout << "Target E: " << e_target << " Iteration: " << this->stats.nTries() << " Accept ratio: " << this->stats.acceptRatio() << " Rmax: " << this->rmax << "  Energy: " << this->current_energy / this->n_atoms << std::endl;
                this->stats.reset();
            }
            this->updateEnergyDistro();
            this->writeSummary();
            this->updateRdfs();
        }
    }

    void targetedEnergyMetropolis(double etarget) {
        int atom_index = this->rng.randint(0, this->n_atoms);
        std::vector<std::vector<double>> trial_atom_positions = this->atom_positions;
        for (int i = 0; i < 3; i ++) {
            trial_atom_positions[atom_index][i] += (2.0 * this->rng.random() - 1) * this->rmax;
            pbc_wrap(trial_atom_positions[atom_index][i], i);
        }

        double next_energy = this->lennardJonesEnergy(trial_atom_positions);
        double delta_e = abs(next_energy - etarget) - abs(this->current_energy - etarget);
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
        this->energies.push_back(this->current_energy);
        this->avg_distance = this->avgDistance();
        this->writeXSF();
        this->stats.increment();
    }

    double avgDistance() {
        double dist_sum = 0;
        int counter = 0;
        for (int i = 0; i < this->n_atoms; i ++) {
            for (int j = 0; j < this->n_atoms; j ++) {
                if (i != j) {
                    dist_sum += this->r(this->atom_positions[i], this->atom_positions[j]);
                    counter ++;
                }
            }
        }
        return dist_sum / counter;
    }

    void classicMetropolis(bool flipped=false) {
        // int atom_index = this->rng.randint(0, this->n_atoms);
        // std::vector<std::vector<double>> trial_atom_positions = this->atom_positions;
        for (int a = 0; a < this->n_atoms; a++ ) {
            std::vector<std::vector<double>> trial_atom_positions = this->atom_positions;

            for (int i = 0; i < 3; i ++) {
                trial_atom_positions[a][i] += (2.0 * this->rng.random() - 1) * this->rmax;
            // pbc_wrap(trial_atom_positions[atom_index][i], i);
            }
        // }
        
            // vector<vector<double>> old_subset_positions, new_subset_positions;
            // vector<int> old_nns = this->nearestNeighbors(a, this->atom_positions);
            // vector<int> new_nns = this->nearestNeighbors(a, trial_atom_positions);
            // set<int> all_nns;
            // vector<int> vec_all_nns;
            
            // for (auto nn : old_nns) {
            //     all_nns.insert(nn);
            // }
            
            // for (auto nn : new_nns) {
            //     all_nns.insert(nn);
            // }

            // old_subset_positions.push_back(this->atom_positions[a]);
            // new_subset_positions.push_back(trial_atom_positions[a]);
            // vec_all_nns.push_back(a);

            // for (auto nn : all_nns) {
            //     vec_all_nns.push_back(nn);
            //     old_subset_positions.push_back(this->atom_positions[nn]);
            //     new_subset_positions.push_back(trial_atom_positions[nn]);
            // }
            
            // // double next_energy = this->aenetEnergy(trial_atom_positions);
            // // double prev_en = this->aenetEnergy(old_subset_positions, 1);
            // double next_en, prev_en;
            // thread t1([&] {next_en = this->aenetEnergy(new_subset_positions, 0);});
            // thread t2([&] {prev_en = this->aenetEnergy(old_subset_positions, 1);});
            // t1.join();
            // t2.join();
            // double next_energy = (next_en - prev_en) + this->current_energy;
            // cout << next_en << "\t" << prev_en << "\n";
            // double next_energy = (this->localLennardJonesEnergy(trial_atom_positions, a) - this->localLennardJonesEnergy(this->atom_positions, a)) + this->current_energy;
            // cout << "Smaller: " << next_energy << " " << " diff: " << next_en - prev_en  << " Current: " << this->current_energy << endl;
            // double full_next_en = this->aenetEnergy(trial_atom_positions, 0);
            // double full_prev_en = this->aenetEnergy(this->atom_positions, 1);
            // cout << "Bigger: " << full_next_en << " diff: " << full_next_en - full_prev_en << " Current: " << full_prev_en << "\n\n";
            // cout << next_energy << endl;
            // double next_energy = this->current_energy + ( this->aenetEnergy(new_subset_positions, vec_all_nns) - this->aenetEnergy(old_subset_positions, vec_all_nns) );
            double next_energy = this->aenetEnergy(trial_atom_positions);
            // cout << this->current_energy << "\t" << next_energy << "\t" << this->r(trial_atom_positions[0], trial_atom_positions[1]) << endl;

            double delta_e = next_energy - this->current_energy;
            double delta_e_beta = delta_e * this->beta;

            if (delta_e_beta < 75.0) {
                if (delta_e_beta < 0.0) {
                    this->current_energy = next_energy;
                    this->atom_positions = trial_atom_positions;
                    this->stats.accept();
                    // this->getForces();
                } else if (exp( - delta_e_beta ) > this->rng.random() ) {
                    this->current_energy = next_energy;
                    this->atom_positions = trial_atom_positions;
                    this->stats.accept();
                    // this->getForces();
                }
            }
    /*        if (this->current_energy < this->min_energy) {
                this->min_energy = this->current_energy;
            }

            if (this->current_energy > this->max_energy) {
                this->max_energy = this->current_energy;
            }*/
            // this->writeXSF();
            // this->energies.push_back(this->current_energy);
            // this->avg_distance = this->avgDistance();
            // if (this->avg_distance > 8.0) {
            //     cout << this->atom_positions[0][0] << " " << this->atom_positions[0][1] << this->atom_positions[0][2] << "\n";
            //     cout << this->atom_positions[1][0] << " " << this->atom_positions[1][1] << this->atom_positions[1][2] << "\n";
            //   cout << this->r(this->atom_positions[0], this->atom_positions[1]) << "\n";
            // }
            this->stats.increment();
            // std::cout << "Iteration: " << this->stats.nTries() << " Accept ratio: " << this->stats.acceptRatio() << " Rmax: " << this->rmax << "  Energy: " << this->current_energy / this->n_atoms << std::endl;
                // output << this->n_atoms << "\n\n";
        }
    }
    
    void readXyzFile(std::string filename) {
        std::ifstream xyzfile(filename.c_str());
        if (!xyzfile) {
            std::cout << "ERROR: We could not open your XYZ file. Please check the input file 'input.yaml' and verify your XYZ file. Exiting...\n";
            exit(-1);
        }
        this->atom_types.clear();
        this->atom_positions.clear();
        xyzfile >> this->n_atoms;
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
        // vector<int> all_indices;
        // for (int i = 0; i < this->n_atoms; i ++) {
        //     all_indices.push_back(i);
        // }
        // std::cout << "INFO: Read in XYZ file.\n";
        this->current_energy = this->aenetEnergy(this->atom_positions);
        // cout << this->current_energy << endl;

        // this->getForces();
        // this->current_energy = this->lennardJonesEnergy(this->atom_positions);
        this->original_positions = this->atom_positions;
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
        this->rng = Rng(13);
        this->max_trials = config.getMaxTrials();
        this->allocateRAMDisk();
        this->readXyzFile(config.xyzFile());
        this->summary_output.open("output/summary.dat");
        this->msd_output.open("output/msd.dat");
        this->initDistros();
    }

    ~System() {
        this->deallocateRAMDisk();
        this->summary_output.close();
        this->msd_output.close();
    }

};

#endif
