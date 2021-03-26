#pragma once
#ifndef mc_h
#define mc_h
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <map>
#include "rule.h"
#include "cluster.h"
#include "sim_cell.h"
using namespace std;

void fillDistList(vector<float> &dist_list, vector<Rule> mc_rules);
void fillRuleList(vector<Rule> &list, string rule_file, int offset);
void runMetropolis3(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules);
void runMetropolis4(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules);
void runMetropolis5(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules);
void runMetropolis6(float passes, float temp1, float temp2, float temp_inc, float thermostat, int eq_passes, SimCell &sim_cell, vector<Rule> &mc_rules);
void runMetropolis7(float passes, float temp1, float temp2, float temp_inc, float sro_target, float thermostat, int eq_passes, SimCell& sim_cell, vector<Rule>& mc_rules);
void writeSuperCell(vector<int> &atom_species, vector<int>& atom_spins, SimCell &sim_cell);
float evalLattice(map<string, float> &rule_map_spin, map<string, float> &rule_map_chem, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list);
float evalLatticeSpin(map<string, float>& rule_map_spin, vector<int>& atom_spin, vector<int>& atom_species, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
float evalSiteEnergyAll(int site, map<string, float> &rule_map_spin, map<string, float> &rule_map_chem, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list);
float evalSiteEnergySpin(int site, map<string, float> &rule_map_spin, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list);
float calcMag2(int site, vector<int> &atom_spin,vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list);
float calcMag2Diff(int site, int old_spin, vector<int>& atom_spin, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list);
float calc_struct(int site, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list);
float calcSpecies(int site, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list);
float delSiteEnergySpin(int site, int old_spin, map<string, float> &rule_map_spin, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list);
float stagMag(int site, int spin, SimCell& sim_cell);
#endif