#pragma once
#ifndef Session_h
#define Session_h
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "file_io.h"
using namespace std;

class Session
{
public:
	int algo;
	string structure_file;
	string rules_file;
	string sim_type;
	string phase_init;
	string spin_init;
	string species_init;
	bool use_poscar;
	int numb_passes;
	int eq_passes;
	float start_temp;
	float end_temp;
	float temp_inc;
	float thermostat_2;
	float sro_target;
	int shape[3] = { 0,0,0 };
	vector<int> atom_numbs;
	Session(string input_file);
};
#endif
