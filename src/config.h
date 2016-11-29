#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <string>
#include <iostream>
#include "yaml-cpp/yaml.h"

class Config {	
	YAML::Node config = YAML::LoadFile("input.yaml");
	YAML::Node cell_dimensions = config["cell_dimensions"];
	YAML::Node rdfs = config["rdfs"];

	double xlat = cell_dimensions["x"].as<double>();
	double ylat = cell_dimensions["y"].as<double>();
	double zlat = cell_dimensions["z"].as<double>();
	double rdf_max_dist = rdfs["max_dist"].as<double>();
	double rdf_incr = rdfs["increment"].as<double>();

	std::string xyz_file = config["xyz_file"].as<std::string>();
	int seed = config["seed"].as<int>();
	int max_trials = config["max_trials"].as<int>();

public:
	std::string xyzFile() {
		return this->xyz_file;
	}

	double xDim() {
		return this->xlat;
	}

	double yDim() {
		return this->ylat;
	}

	double zDim() {
		return this->zlat;
	}

	int getSeed() {
		return this->seed;
	}

	int getMaxTrials() {
		return this->max_trials;
	}

	double rdfMaxDist() {
		return this->rdf_max_dist;
	}

	double rdfIncrement() {
		return this->rdf_incr;
	}

	void printConfig() {
		std::cout << "Cell dimensions: " << "x = " << this->xlat << ", y = " << this->ylat <<  ", z = " << this->zlat << std::endl;
		std::cout << "XYZ file: " << this->xyz_file << std::endl;
		std::cout << "Seed: " << this->seed << std::endl;
	}

};

#endif