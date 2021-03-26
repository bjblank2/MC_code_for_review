#include "monte_carlo.h"
#include <cmath>
#include <random>
extern SimCell sim_cell;

// Fill list of distances from mc_rules
void fillDistList(vector<float> &dist_list, vector<Rule> mc_rules) {
	// Loop through all mc_rules
	for (int i = 0; i < mc_rules.size(); i++) {
		// for two atom terms
		if (mc_rules[i].GetLength() == 2) {
			if (find(dist_list.begin(), dist_list.end(), mc_rules[i].GetDists()[0]) == dist_list.end()) {
				dist_list.push_back(mc_rules[i].GetDists()[0]);
			}
		}
		// for three atom terms
		else if (mc_rules[i].GetLength() == 3) {
			if (find(dist_list.begin(), dist_list.end(), mc_rules[i].GetDists()[0]) == dist_list.end()) {
				dist_list.push_back(mc_rules[i].GetDists()[0]);
			}
			if (find(dist_list.begin(), dist_list.end(), mc_rules[i].GetDists()[1]) == dist_list.end()) {
				dist_list.push_back(mc_rules[i].GetDists()[1]);
			}
			if (find(dist_list.begin(), dist_list.end(), mc_rules[i].GetDists()[2]) == dist_list.end()) {
				dist_list.push_back(mc_rules[i].GetDists()[2]);
			}
		}
	}
}

// Fill list of mc_rules objects from rule file
void fillRuleList(vector<Rule> &list, string rule_file, int offset) {
	
	vector<float> distances;
	vector<int> spins;
	vector<int> species;
	float energy_contribution = 0;
	int rule_type = 0;
	int rule_length = 0;
	string phase = "";
	string rule_line;
	vector<string> rule_lines;
	ifstream rule_list;
	rule_list.open(rule_file);
	// Parce rule list txt file
	if (rule_list.is_open()) {
		while (getline(rule_list, rule_line))
		{
			rule_lines.push_back(rule_line);
		}
		rule_list.close();
		for (int i = 0; i < rule_lines.size(); i++) {
			list.push_back(Rule(rule_lines[i]));
		}
	}
	else cout << "Unable to open file";
}

// Evaluate the energy controbution from a singal atom, includeing spin chem rules.
float evalSiteEnergyAll(int site, map<string, float> &rule_map_spin, map<string, float> &rule_map_chem, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	// Initalize variables
	float Kb = 0.00008617333262;
	float uB = .000057883818012;
	float H = 0;
	string key;
	float enrg = 0;
	map<string, float>::iterator rule_itr;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;
	int neighbor_site2;
	int neighbor_species2;
	int neighbor_spin2;
	int dist_2_ind;
	float neighbor_dist2;
	float neighbor_dist3;
	enrg -= 3 * uB*H*site_spin; // add energy from external magnetic field (if there is one... right now its hard coded to 0)
	key = "_" + to_string(site_species) + ",0,"; // make key for one atom term
	rule_itr = rule_map_chem.find(key);
	enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // apply rule for one atom term
	// for 2 atom terms
	// loop through all neighbors in the neighbor list (list only includes atoms that have distances listed in the dist_list (the list of distances used in a rule))
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		neighbor_site1 = neighbor_index_list[site][i];
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ","; // make key string
		rule_itr = rule_map_chem.find(key);
		enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second : 0.0; // if the key was in the chem rule map, apply the rule, otherwise do nothing
		rule_itr = rule_map_spin.find(key);
		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * site_spin * neighbor_spin1 : 0.0; // if the key was in the spin rule map, apply the rule, otherwise do nothing
		// for 3 atom terms
		// loop through all neighbors in the neighbor list
		for (int j = 0; j < neighbor_index_list[site].size(); j++) {
			if (i != j) {
				neighbor_site2 = neighbor_index_list[site][j];
				neighbor_species2 = atom_species[neighbor_site2];
				neighbor_spin2 = atom_spin[neighbor_site2];
				vector<int>::iterator dist_2_itr = find(neighbor_index_list[neighbor_site1].begin(), neighbor_index_list[neighbor_site1].end(), neighbor_site2);
				if (dist_2_itr != neighbor_index_list[neighbor_site1].end()) {
					dist_2_ind = distance(neighbor_index_list[neighbor_site1].begin(), dist_2_itr);
					neighbor_dist2 = neighbor_dist_list[neighbor_site1][dist_2_ind];
					neighbor_dist3 = neighbor_dist_list[site][j];
					// make key for three atom term
					key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_species2);
					key +="," + to_string(neighbor_dist1) + "," + to_string(neighbor_dist2) + "," + to_string(neighbor_dist3) + ",";
					rule_itr = rule_map_chem.find(key);
					enrg += (rule_itr != rule_map_chem.end()) ? rule_itr->second / 6 : 0.0; // apply rule (if possible). Devide by six because of tripple count degeneracy of three atom terms
				}
			}
		}
	}
	return enrg;
}

// Evaluate the energy contribution from a singal atom using only spin rules. (Assuming the spin is all that changed)
float evalSiteEnergySpin(int site, map<string, float> &rule_map_spin, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	// Initalize variables and constants
	float Kb = 0.00008617333262;
	float uB = .000057883818012;
	float H = 0;
	string key;
	float enrg = 0;
	map<string, float>::iterator rule_itr;
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	enrg -= 3 * uB*H*site_spin; // apply magnetic field (if any, right now its hard coded to zero)
	// for 2 atom terms
	// Loop through all neighbors in the neighbor list (neighbor list only includes distances included in the dist_list (the list of distances used in a rule))
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		// gather info to build key
		neighbor_site1 = neighbor_index_list[site][i];
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ","; // create key string
		rule_itr = rule_map_spin.find(key);
  		enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * site_spin * neighbor_spin1 : 0.0; // if key is found in rule list, apply the rule, else do nothing.
	}
	return enrg;
}

float calc_struct(int site, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	int site_species = atom_species[site];
	int count = 0;
	if (site_species == 2) {
		for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
			if (neighbor_dist_list[site][i] == 0.5 or neighbor_dist_list[site][i] == 0.75) {
				if (atom_species[neighbor_index_list[site][i]] == 1) {
					count += 1;
				}
			}
		}
	}
	return count / 6.0;
}

// calculate the contribution to the spin-product order paramiter from a single site
float calcMag2(int site, vector<int> &atom_spin, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list){
	int site_spin = atom_spin[site];
	float mag_product = 0;
	// loop through all neighbors in neighbor list
	for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
		if (abs(neighbor_dist_list[site][i] - 0.707106769) < 0.001) {
		//if (neighbor_dist_list[site][i] == .5  or neighbor_dist_list[site][i] == 0.75) {
			mag_product += site_spin * atom_spin[neighbor_index_list[site][i]]; // if neighbor is directly above/below or directly adjacent... do the product
		}
	}
	return mag_product; // devide by 6 for normalization
}

// calculate the contribution to the spin-product order paramiter from a single site overload
float calcMag2Diff(int site, int old_spin, vector<int>& atom_spin, vector<vector<int>>& neighbor_index_list, vector<vector<float>>& neighbor_dist_list) {
	int site_spin = atom_spin[site];
	float mag_product = 0;
	float mag_product_old = 0;
	// loop through all neighbors in neighbor list
	for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
		if (abs(neighbor_dist_list[site][i] - 0.707106769) < 0.001) {
		//if (neighbor_dist_list[site][i] == .5 or neighbor_dist_list[site][i] == 0.75) {
			mag_product += 2 * (site_spin - old_spin) * atom_spin[neighbor_index_list[site][i]]; // if neighbor is directly above/below or directly adjacent... do the product
		}
	}
	return (mag_product); // devide by 6 for normalization
}

// claculate the stag mag order paramiter contribution from a single atom
float stagMag(int site, int spin, SimCell& sim_cell) {
	int mag = 0;
	int check = fmod(sim_cell.atom_list[site].pos[0] * sim_cell.atom_list[site].pos[1] * 2,2); // check to see what sublattice the site is on
	if (check == 0) { mag += spin; }
	else { mag += 0.0; }
	return mag;
}

float stagMagDiff(int site, int spin, int old_spin, SimCell& sim_cell) {
	int mag = 0;
	int check = fmod(sim_cell.atom_list[site].pos[0] * sim_cell.atom_list[site].pos[1] * 2, 2); // check to see what sublattice the site is on
	if (check == 0) { mag += (spin - old_spin); }
	else { mag += 0.0; }
	return mag;
}
// Curently unused!!
float calcSpecies(int site, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	int site_species = 0;
	int neighbor_species = 0;
	if (atom_species[site] == 1) { site_species = 1; }
	else if (atom_species[site] == 2) { site_species = -1; }
	float species_product = 0;
	for (int i = 0; i < neighbor_dist_list[site].size(); i++) {
		if (neighbor_dist_list[site][i] == .5 or neighbor_dist_list[site][i] == 0.75) {
			if (atom_species[neighbor_index_list[site][i]] == 1) { neighbor_species = 1; }
			else if (atom_species[neighbor_index_list[site][i]] == 2) { neighbor_species = -1; }
			species_product += site_species * neighbor_species;
		}
	}
	return species_product / 6;
}
// Curently unused!!
float delSiteEnergySpin(int site, int old_spin, map<string, float> &rule_map_spin, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	float Kb = 0.000086173324;
	float uB = .000057883818012;
	float H = 0;
	string key;
	float enrg = 0;
	map<string, float>::iterator rule_itr;
	int neighbor_site1;
	int neighbor_species1;
	int neighbor_spin1;
	float neighbor_dist1;
	int site_species = atom_species[site];
	int site_spin = atom_spin[site];
	//	site_energy -= Kb * temp * log(8)*(1 - pow(site_phase, 2));
	enrg -= 3 * uB*H*site_spin;
	// for 2 atom terms
	for (int i = 0; i < neighbor_index_list[site].size(); i++) {
		neighbor_site1 = neighbor_index_list[site][i];
		neighbor_species1 = atom_species[neighbor_site1];
		neighbor_spin1 = atom_spin[neighbor_site1];
		neighbor_dist1 = neighbor_dist_list[site][i];
		key = "_" + to_string(site_species) + "," + to_string(neighbor_species1) + "," + to_string(neighbor_dist1) + ",";
		rule_itr = rule_map_spin.find(key);
		if (site_spin == 0) {
			enrg -= (rule_itr != rule_map_spin.end()) ? rule_itr->second * old_spin * neighbor_spin1 : 0.0;
		}
		else {
			enrg += (rule_itr != rule_map_spin.end()) ? rule_itr->second * 2 * site_spin * neighbor_spin1 : 0.0;
		}
	}
	return enrg;
}
// Overloaded def of evalLattice that is only used for BEG. Don't need to look at it

// Evaluate the total energy of the simulation cell in it's current state
float evalLattice(map<string, float> &rule_map_spin, map<string, float> &rule_map_chem, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	float enrg = 0;
	// loop thruogh simulation cell and add everything up
	for (int site = 0; site < atom_species.size(); site++) {
		enrg += evalSiteEnergyAll(site, rule_map_spin, rule_map_chem, atom_spin, atom_species, neighbor_index_list, neighbor_dist_list);
	}
	return enrg;
}

// Evaluate the spin contribution to the total energy of the simulation cell
float evalLatticeSpin(map<string, float> &rule_map_spin, vector<int> &atom_spin, vector<int> &atom_species, vector<vector<int>> &neighbor_index_list, vector<vector<float>> &neighbor_dist_list) {
	float enrg = 0;
	// loop through everything and add it up
	for (int site = 0; site < atom_species.size(); site++) {
		enrg += evalSiteEnergySpin(site, rule_map_spin, atom_spin, atom_species, neighbor_index_list, neighbor_dist_list);
	}
	return enrg;// / atom_species.size() * 16 - 66.8295069760671;
}

// Run MC with ONLY spin flips
void runMetropolis3(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules) {
	// Pull some basic information out of the sim_cell object (an object mostly used for setup that contains inforamtion on the
	// position, spin, atomic species, and neighbor list for each atom in the simulation cell)
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, float> rule_map_spin;
	map <string, float> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into a map for spin and a map for chemestry rules to make lookup faster.
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	// Loop through all mc_rules
	for (int i = 0; i < mc_rules.size(); i++) {
		string rule_key = "_"; // add starting charicter to key string for parcing purposes.
		// If rule is single atom term...
		if (mc_rules[i].GetLength() == 1) {
			rule_key = "_" + to_string(mc_rules[i].GetSpecies()[0]) + ",0,"; // make key string
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));// assign key atom pair 
		}
		// If rule is two atom term...
		else if (mc_rules[i].GetLength() == 2) {
			vector<int> species = mc_rules[i].GetSpecies();
			float dist = mc_rules[i].GetDists()[0];
			// make degenerate keys for [0,1]->[1,0] type rules
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont())); // assign key atom pair if chem type
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont())); // assign key atom pair if spin type
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont())); // assign key atom pair if chem type
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont())); // assign key atom pair if spin type
			}
		}
		// if rule is three atom term...
		else if (mc_rules[i].GetLength() == 3) {
			vector<int> trip = mc_rules[i].GetSpecies();
			vector<float> dists = mc_rules[i].GetDists();
			// make keys and assign key value pairs for all degenerate three atom ruels, (always chem rules)
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
	}
	// Extract species list, spin list, neighbor index list and neighbor dist list from sim_cell object to make them more accesable latter.
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	// initalize variables 
	float Kb = 0.00008617333262;
	float e_total = 0;
	float e_flip = 0;
	float e_site = 0;
	float e_site_new = 0;
	float new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	int spin_avg_Mn = 0;
	int spin_total_Mn = 0;
	float spin_avg_Mn2 = 0;
	float spin_total_Mn2 = 0;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	// Initalize random number generator. This part might give your computer some problems if it is a Mac FYI !!!!!!!!!!
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// Calculate the total energy of the initalized simulation cell
	float initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	// Calculate the total energy contribution form magnetic ECI
	float initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms << "\n";
	// Enter main MC loop
	//Begin decreasing temperature
	ofstream Output;
	Output.open("OUTPUT.txt");
	Output << "Phase: " << sim_cell.phase_init << "\n";
	Output << "Composition: " << sim_cell.species_numbs[0] << ", " << sim_cell.species_numbs[1] << ", " << sim_cell.species_numbs[2] << "\n";
	Output << "MC passes: " << passes << "\n";
	Output << "Begin MC using ALGO_3\n";
	for (int temp = temp1; temp >= temp2; temp+=temp_inc) {
		e_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		spin_avg_Mn = 0;
		spin_avg_Mn2 = 0;
		flip_count = 0;
		flip_count2 = 0;
		// loop through the simulation cell specified number of times
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			spin_total_Mn = 0;
			spin_total_Mn2 = 0;
			// loop through every atom in the cell
			for (int site = 0; site < numb_atoms; site++) {
				// Flip a spin
				old_spin = atom_spin_list[site]; // get pre-flip spin
				e_site = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list); // get pre-flip energy
				// pick a new spin at random (don't want to waist time picking same spin that site already has)
				spin_same = true;
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) {
						new_spin = -1;
					}
					else if (spin_rand <= 0.6666666666666666) {
						new_spin = 0;
					}
					else {
						new_spin = 1;
					}
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin; // flip the spin
				e_site_new = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list); // find energy after the flip
				e_flip = (e_site_new - e_site); // find net change in energy
				// apply chose keep/throw out using Metropolis
				if (e_flip < 0) {
					flip_count2 += 1;
				}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_flip));
					if (keep_rand < keep_prob) {
						flip_count += 1;
					}
					else {
						atom_spin_list[site] = old_spin;
						e_flip = 0;
					}
				}
				initial_enrg += e_flip; // update the total energy
				e_total += initial_enrg / numb_atoms; // normalize
				current_spin = atom_spin_list[site]; 
				spin_total += current_spin; // get spin total
				if (atom_species_list[site] == 1) {
					spin_total_Mn += current_spin; // Mn only spin totla
					spin_total_Mn2 += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list); // calc spin product order param
					//spin_total_Mn2 += stagMag(site, atom_spin_list[site], sim_cell); // calc stag mag order param
				}
			}
			spin_avg += spin_total;
			spin_avg_Mn += spin_total_Mn;
			spin_avg_Mn2 += spin_total_Mn2;
			e_avg += e_total / passes;
		}
		// Output normaized avarages of tracked quantities
		Output << temp;
		Output << " , ";
		Output << e_avg / numb_atoms;
		Output << " , ";
		Output << spin_avg / passes / numb_atoms;
		Output << " , ";
		Output << spin_avg_Mn / passes / sim_cell.species_numbs[1];
		Output << " , ";  
		Output << spin_avg_Mn2 / passes / sim_cell.species_numbs[1];
		Output << " , ";
		Output << flip_count;
		Output << " , ";
		Output << flip_count2;
		Output << "\n";
		cout << temp << "\n";
	}
	writeSuperCell(atom_species_list, atom_spin_list, sim_cell); // write the filnal state of the simulation cell to txt file (for debugging and sanity checking mostly)
}

// This is not the function I am testing with. You can look through it if you want but its not as important 
void runMetropolis4(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules) {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, float> rule_map_spin;
	map <string, float> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < mc_rules.size(); i++) {
		string rule_key = "_";
		if (mc_rules[i].GetLength() == 1) {
			rule_key = "_" + to_string(mc_rules[i].GetSpecies()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
		else if (mc_rules[i].GetLength() == 2) {
			vector<int> species = mc_rules[i].GetSpecies();
			float dist = mc_rules[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
		}
		else if (mc_rules[i].GetLength() == 3) {
			vector<int> trip = mc_rules[i].GetSpecies();
			vector<float> dists = mc_rules[i].GetDists();
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	float Kb = .0000861733035;
	float e_total = 0;
	float e_flip = 0;
	float e_site = 0;
	float e_site_new = 0;
	float new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	int spin_avg_Mn = 0;
	int spin_total_Mn = 0;
	float spin_avg_Mn2 = 0;
	float spin_total_Mn2 = 0;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// Find ground state configuration
	int temp_inc_relax = -10;
	for (int temp = temp1; temp >= temp2; temp += temp_inc_relax) {
		cout << temp << "\n";
		if (temp <= 20) {
			temp_inc_relax = -1;
		}
		for (int i = 0; i < passes; i++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						e_site = evalSiteEnergyAll(site, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_site += evalSiteEnergyAll(rand_index, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						e_site_new = evalSiteEnergyAll(site, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_site_new += evalSiteEnergyAll(rand_index, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_flip = e_site_new - e_site;
						if (e_flip < 0) {}
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb*temp)*(e_flip));
							if (keep_rand < keep_prob) {}
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
							}
						}
					}
				}
				// Flip Spin
				//e_site_old = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				if (e_flip < 0) {}
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_flip));
					if (keep_rand < keep_prob) {}
					else { atom_spin_list[site] = old_spin; }
				}
			}
		}
	}
	// begin actual MC run //
	writeSuperCell(atom_species_list, atom_spin_list, sim_cell);
	float initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + -4.176806541934177 << ", " << initial_spin_cont / numb_atoms << "\n";
	for (int temp = temp1; temp >= temp2; temp += temp_inc) {
		e_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		spin_avg_Mn = 0;
		spin_avg_Mn2 = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			spin_total_Mn = 0;
			spin_total_Mn2 = 0;
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Spin
				//e_site_old = evalSiteEnergySpin(temp, site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				if (e_flip < 0) { flip_count += 1; }
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_flip));
					if (keep_rand < keep_prob) { flip_count2 += 1; }
					else {
						atom_spin_list[site] = old_spin;
						e_flip = 0;
					}
				}
				initial_enrg += e_flip;
				e_total += initial_enrg / numb_atoms;
				current_spin = atom_spin_list[site];
				spin_total += current_spin;
				if (atom_species_list[site] == 1) {
					spin_total_Mn += current_spin;
					spin_total_Mn2 += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list);
				}
			}
			spin_avg += spin_total;
			spin_avg_Mn += spin_total_Mn;
			spin_avg_Mn2 += spin_total_Mn2;
			e_avg += evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list) / (passes);        
			//e_avg += e_total / (passes);
		}
		cout << temp;
		cout << " , ";
		cout << e_avg / numb_atoms;
		cout << " , ";
		cout << spin_avg / passes / numb_atoms;
		cout << " , ";
		cout << spin_avg_Mn / passes / sim_cell.species_numbs[1];
		cout << " , ";
		cout << spin_avg_Mn2 / passes / sim_cell.species_numbs[1];
		cout << " , ";
		cout << flip_count;
		cout << " , ";
		cout << flip_count2;
		cout << "\n";
	}
}

// This is not the function I am testing with. You can look through it if you want but its not as important 
void runMetropolis5(float passes, float temp1, float temp2, float temp_inc, SimCell &sim_cell, vector<Rule> &mc_rules) {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, float> rule_map_spin;
	map <string, float> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < mc_rules.size(); i++) {
		string rule_key = "_";
		if (mc_rules[i].GetLength() == 1) {
			rule_key = "_" + to_string(mc_rules[i].GetSpecies()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
		else if (mc_rules[i].GetLength() == 2) {
			vector<int> species = mc_rules[i].GetSpecies();
			float dist = mc_rules[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
		}
		else if (mc_rules[i].GetLength() == 3) {
			vector<int> trip = mc_rules[i].GetSpecies();
			vector<float> dists = mc_rules[i].GetDists();
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	float Kb = .0000861733035;
	float e_total = 0;
	float e_flip = 0;
	float e_site = 0;
	float e_site_new = 0;
	float new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	int spin_avg_Mn = 0;
	int spin_total_Mn = 0;
	float spin_avg_Mn2 = 0;
	float spin_total_Mn2 = 0;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// begin actual MC run //
	ofstream Output;
	Output.open("OUTPUT.txt");
	Output << "Beginning MC run using ALGO_5\n";
	float initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + -4.918174088267291 << ", " << initial_spin_cont / numb_atoms << "\n";
	Output << initial_enrg / numb_atoms + -4.918174088267291 << ", " << initial_spin_cont / numb_atoms << "\n";

	for (int temp = temp1; temp >= temp2; temp += temp_inc) {
		cout << temp << "\n";
		e_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		spin_avg_Mn = 0;
		spin_avg_Mn2 = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			spin_total_Mn = 0;
			spin_total_Mn2 = 0;
			float order_enrgy_running = 0;
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Species
				if (atom_species_list[site] != 0) {
					int rand_index = site;
					while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
					int old_species_site = atom_species_list[site];
					int old_species_rand = atom_species_list[rand_index];
					if (old_species_site != old_species_rand) {
						e_site = evalSiteEnergyAll(site, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_site += evalSiteEnergyAll(rand_index, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						atom_species_list[site] = old_species_rand;
						atom_species_list[rand_index] = old_species_site;
						e_site_new = evalSiteEnergyAll(site, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_site_new += evalSiteEnergyAll(rand_index, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
						e_flip = e_site_new - e_site;
						if (e_flip < 0) {}
						else {
							keep_rand = unif(rng);
							keep_prob = exp(-1 / (Kb*temp)*(e_flip));
							if (keep_rand < keep_prob) {}
							else {
								atom_species_list[site] = old_species_site;
								atom_species_list[rand_index] = old_species_rand;
								e_flip = 0;
							}
						}
						initial_enrg += e_flip;
					}
				}
				// Flip Spin
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				if (e_flip < 0) { flip_count += 1; }
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb*temp)*(e_flip));
					if (keep_rand < keep_prob) { flip_count2 += 1; }
					else {
						atom_spin_list[site] = old_spin;
						e_flip = 0;
					}
				}
				initial_enrg += e_flip;
				e_total += initial_enrg;
				current_spin = atom_spin_list[site];
				spin_total += current_spin;
				if (atom_species_list[site] == 1) { spin_total_Mn += current_spin; }
				spin_total_Mn2 += calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list); 
			}
			spin_avg += spin_total;
			spin_avg_Mn += spin_total_Mn;
			spin_avg_Mn2 += spin_total_Mn2;
			e_avg = e_total / passes;
		}
		Output << temp;
		Output << " , ";
		Output << e_avg / numb_atoms; // -4.918174088267291;
		Output << " , ";
		Output << spin_avg / passes / numb_atoms;
		Output << " , ";
		Output << spin_avg_Mn / passes / sim_cell.species_numbs[1];
		Output << " , ";
		Output << spin_avg_Mn2 / passes / sim_cell.species_numbs[2];
		Output << " , ";
		Output << flip_count;
		Output << " , ";
		Output << flip_count2;
		Output << "\n";
	}
	writeSuperCell(atom_species_list, atom_spin_list, sim_cell);
	Output.close();
}

// Uses seperate thermostat for Spin and Chem
void runMetropolis6(float passes, float temp1, float temp2, float temp_inc, float thermostat, int eq_passes, SimCell& sim_cell, vector<Rule>& mc_rules) {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, float> rule_map_spin;
	map <string, float> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < mc_rules.size(); i++) {
		string rule_key = "_";
		if (mc_rules[i].GetLength() == 1) {
			rule_key = "_" + to_string(mc_rules[i].GetSpecies()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
		else if (mc_rules[i].GetLength() == 2) {
			vector<int> species = mc_rules[i].GetSpecies();
			float dist = mc_rules[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
		}
		else if (mc_rules[i].GetLength() == 3) {
			vector<int> trip = mc_rules[i].GetSpecies();
			vector<float> dists = mc_rules[i].GetDists();
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	float Kb = .0000861733035;
	float e_total = 0;
	float e_flip = 0;
	float e_site = 0;
	float e_site_new = 0;
	float new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	int current_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float spin_total = 0;
	float spin_avg2 = 0;
	float spin_total2 = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	int spin_avg_Mn = 0;
	int spin_total_Mn = 0;
	float spin_avg_Mn2 = 0;
	float spin_total_Mn2 = 0;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// begin actual MC run //
	ofstream Output;
	Output.open("OUTPUT.txt");
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: " << sim_cell.species_numbs[0] << ", " << sim_cell.species_numbs[1] << ", " << sim_cell.species_numbs[2] << "\n";
	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_6\n";
	float initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + -4.918174088267291 << ", " << initial_spin_cont / numb_atoms << "\n";
	Output << initial_enrg / numb_atoms + -4.918174088267291 << ", " << initial_spin_cont / numb_atoms << "\n";
	cout << thermostat << "\n";
	float sro;
	for (int i = 0; i < eq_passes; i++) {
		sro = 0;
		for (int site = 0; site < numb_atoms; site++) {
			// Flip Species
			if (atom_species_list[site] != 0) {
				int rand_index = site;
				while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
				int old_species_site = atom_species_list[site];
				int old_species_rand = atom_species_list[rand_index];
				if (old_species_site != old_species_rand) {
					e_site = evalSiteEnergyAll(site, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
					e_site += evalSiteEnergyAll(rand_index, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
					atom_species_list[site] = old_species_rand;
					atom_species_list[rand_index] = old_species_site;
					e_site_new = evalSiteEnergyAll(site, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
					e_site_new += evalSiteEnergyAll(rand_index, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
					e_flip = e_site_new - e_site;
					if (e_flip < 0) {}
					else {
						keep_rand = unif(rng);
						keep_prob = exp(-1 / (Kb * thermostat) * (e_flip));
						if (keep_rand < keep_prob) {}
						else {
							atom_species_list[site] = old_species_site;
							atom_species_list[rand_index] = old_species_rand;
						}
					}
				}
			}
			// Flip Spin
			old_spin = atom_spin_list[site];
			spin_same = true;
			e_site = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
			while (spin_same == true) {
				spin_rand = unif(rng);
				if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
				else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
				else { new_spin = 1; }
				if (new_spin != old_spin) { spin_same = false; }
			}
			atom_spin_list[site] = new_spin;
			e_site_new = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
			e_flip = e_site_new - e_site;
			if (e_flip < 0) { flip_count += 1; }
			else {
				keep_rand = unif(rng);
				keep_prob = exp(-1 / (Kb * thermostat) * (e_flip));
				if (keep_rand < keep_prob) { flip_count2 += 1; }
				else {
					atom_spin_list[site] = old_spin;
				}
			}
			sro += calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
		}
	}
	Output << "SRO: " << sro / sim_cell.species_numbs[1] << "\n";
	//Begin REAL MC
	initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	Output << initial_enrg / numb_atoms + -4.918174088267291 << ", " << initial_spin_cont / numb_atoms << "\n";
	for (int temp = temp1; temp >= temp2; temp += temp_inc) {
		cout << temp << "\n";
		e_avg = 0;
		spin_avg = 0;
		spin_avg2 = 0;
		spin_avg_Mn = 0;
		spin_avg_Mn2 = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int i = 0; i < passes; i++) {
			e_total = 0;
			spin_total = 0;
			spin_total2 = 0;
			spin_total_Mn = 0;
			spin_total_Mn2 = 0;
			float order_enrgy_running = 0;
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Spin
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				if (e_flip < 0) { flip_count += 1; }
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb * temp) * (e_flip));
					if (keep_rand < keep_prob) { flip_count2 += 1; }
					else {
						atom_spin_list[site] = old_spin;
						e_flip = 0;
					}
				}
				initial_enrg += e_flip;
				e_total += initial_enrg;
				current_spin = atom_spin_list[site];
				spin_total += current_spin;
				if (atom_species_list[site] == 1) {
					spin_total_Mn += current_spin;
					spin_total_Mn2 += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list);
				}
			}
			spin_avg += spin_total;
			spin_avg_Mn += spin_total_Mn;
			spin_avg_Mn2 += spin_total_Mn2;
			e_avg = e_total / passes;
			//e_avg += evalLattice(temp1, rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list) / (passes);
		}
		Output << temp;
		Output << " , ";
		Output << e_avg / numb_atoms; // -4.918174088267291;
		Output << " , ";
		Output << spin_avg / passes / numb_atoms;
		Output << " , ";
		Output << spin_avg_Mn / passes / sim_cell.species_numbs[1];
		Output << " , ";
		Output << spin_avg_Mn2 / passes / sim_cell.species_numbs[1];
		Output << " , ";
		Output << flip_count;
		Output << " , ";
		Output << flip_count2;
		Output << "\n";
	}
	writeSuperCell(atom_species_list, atom_spin_list, sim_cell);
	Output.close();
}

void runMetropolis7(float passes, float temp1, float temp2, float temp_inc, float sro_target, float thermostat, int eq_passes, SimCell& sim_cell, vector<Rule>& mc_rules) {
	int numb_atoms = sim_cell.numb_atoms;
	int numb_neighbors = sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell);
	map <string, float> rule_map_spin;
	map <string, float> rule_map_chem;
	vector<int> atom_species_list;
	vector<int> atom_spin_list;
	vector<vector<int>> neighbor_index_list(numb_atoms, vector<int>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	vector<vector<float>> neighbor_dist_list(numb_atoms, vector<float>(sim_cell.atom_list[0].getNumbNeighbors(0, sim_cell), 0));
	// Turn rule list into map for spin and map for chem
	// Redundent rules are added for [0,1] and [1,0] ect... so that comparison is quicker in mc steps
	for (int i = 0; i < mc_rules.size(); i++) {
		string rule_key = "_";
		if (mc_rules[i].GetLength() == 1) {
			rule_key = "_" + to_string(mc_rules[i].GetSpecies()[0]) + ",0,";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
		else if (mc_rules[i].GetLength() == 2) {
			vector<int> species = mc_rules[i].GetSpecies();
			float dist = mc_rules[i].GetDists()[0];
			rule_key = "_" + to_string(species[0]) + "," + to_string(species[1]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			rule_key = "_" + to_string(species[1]) + "," + to_string(species[0]) + "," + to_string(dist) + ",";
			if (mc_rules[i].GetType() == 0) {
				rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
			else if (mc_rules[i].GetType() == 1) {
				rule_map_spin.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			}
		}
		else if (mc_rules[i].GetLength() == 3) {
			vector<int> trip = mc_rules[i].GetSpecies();
			vector<float> dists = mc_rules[i].GetDists();
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(trip[2]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[1]) + "," + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[0]) + "," + to_string(trip[1]) + "," + to_string(dists[2]) + "," + to_string(dists[0]) + "," + to_string(dists[1]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
			rule_key = "_" + to_string(trip[2]) + "," + to_string(trip[1]) + "," + to_string(trip[0]) + "," + to_string(dists[1]) + "," + to_string(dists[0]) + "," + to_string(dists[2]) + ",";
			rule_map_chem.insert(pair<string, float>(rule_key, mc_rules[i].GetEnrgCont()));
		}
	}
	// make atom_list more acessable (and a map) for species and spin and neighbors
	for (int i = 0; i < sim_cell.numb_atoms; i++) {
		atom_species_list.push_back(sim_cell.atom_list[i].getSpecies());
		atom_spin_list.push_back(sim_cell.atom_list[i].getSpin());
		for (int j = 0; j < numb_neighbors; j++) {
			neighbor_index_list[i][j] = sim_cell.atom_list[i].getNeighborIndex(j, sim_cell);
			neighbor_dist_list[i][j] = sim_cell.atom_list[i].getNeighborDist(j, sim_cell);
		}
	}
	// Start MC stuff
	float Kb = .0000861733035;
	float e_flip = 0;
	float spin_flip = 0;
	float spin_flipMn = 0;
	float spin_flipMn2 = 0;
	float stag_flip = 0;
	float e_site = 0;
	float e_site_new = 0;
	float new_enrg = 0;
	float spin_rand = 0;
	float keep_rand = 0;
	int old_spin = 0;
	int new_spin = 0;
	bool spin_same;
	float e_avg = 0;
	float spin_avg = 0;
	float stag_avg = 0;
	float keep_prob = 0;
	int flip_count = 0;
	int flip_count2 = 0;
	float spin_avg_Mn = 0;
	float spin_avg_Mn2 = 0;
	map<int, int>::iterator atom_itr;
	map<int, int>::iterator spin_itr;
	std::mt19937_64 rng;
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	rng.seed(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	// begin actual MC run //
	ofstream Output;
	Output.open("OUTPUT.txt");
	Output << "Phase: " << sim_cell.phase_init;
	Output << "Composition: " << sim_cell.species_numbs[0] << ", " << sim_cell.species_numbs[1] << ", " << sim_cell.species_numbs[2] << "\n";
	Output << "MC passes: " << passes << "\n";
	Output << "Beginning MC EQ run using ALGO_7\n";
	Output << "EQ passes: " << eq_passes << ", EQ Temp: " << thermostat << "\n";
	float initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	float initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	cout << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms <<", " << initial_enrg << "\n";
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	cout << "SRO Target: " << sro_target << "\n";
	float sro_final = 0;
	float sro_initial = 0;
	float sro_site_new = 0;
	float sro_site_old = 0;
	float sro_site_flip = 0;
	float sro_flip = 0;
	for (int i = 0; i < numb_atoms; i++) {
		sro_initial += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	for (int i = 0; i < eq_passes; i++) {
		for (int site = 0; site < numb_atoms; site++) {
			// Flip Species
			if (atom_species_list[site] != 0) {
				int rand_index = site;
				while (rand_index == site or atom_species_list[rand_index] == 0) { rand_index = rand() % numb_atoms; }
				int old_species_site = atom_species_list[site];
				int old_species_rand = atom_species_list[rand_index];
				if (old_species_site != old_species_rand) {
					sro_site_old = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
					sro_site_old += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
					atom_species_list[site] = old_species_rand;
					atom_species_list[rand_index] = old_species_site;
					sro_site_new = calc_struct(site, atom_species_list, neighbor_index_list, neighbor_dist_list);
					sro_site_new += calc_struct(rand_index, atom_species_list, neighbor_index_list, neighbor_dist_list);
					sro_site_flip = 2 * (sro_site_new - sro_site_old) / sim_cell.species_numbs[2];
					sro_flip = abs((sro_site_flip + sro_initial) - sro_target) - abs(sro_initial - sro_target);
					if (sro_flip < 0) { sro_initial += sro_site_flip; }
					else {
						keep_rand = unif(rng);
						keep_prob = exp(-1 / (Kb * thermostat) * (sro_flip));
						if (keep_rand < keep_prob) { sro_initial += sro_site_flip; }
						else {
							atom_species_list[site] = old_species_site;
							atom_species_list[rand_index] = old_species_rand;
						}
					}
				}
			}
		}
	}
	for (int i = 0; i < numb_atoms; i++) {
		sro_final += calc_struct(i, atom_species_list, neighbor_index_list, neighbor_dist_list) / sim_cell.species_numbs[2];
	}
	cout << "SRO: " << sro_final << "\n";
	Output << "SRO: " << sro_final << "\n";
	//Begin REAL MC
	initial_enrg = evalLattice(rule_map_spin, rule_map_chem, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	initial_spin_cont = evalLatticeSpin(rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
	Output << initial_enrg / numb_atoms + 0 << ", " << initial_spin_cont / numb_atoms << "\n";
	float initial_spin = 0.0;
	float initial_spin_Mn = 0.0;
	float initial_spin_2_Mn = 0.0;
	float initial_stag_mag = 0;
	for (int site = 0; site < numb_atoms; site++) {
		initial_spin += atom_spin_list[site];
		initial_spin_2_Mn += calcMag2(site, atom_spin_list, neighbor_index_list, neighbor_dist_list);
		initial_stag_mag += stagMag(site, atom_spin_list[site], sim_cell);
		if (atom_species_list[site] == 1) {
			initial_spin_Mn += atom_spin_list[site];
		}
	}
	cout << "spin " << initial_spin/numb_atoms << " spin2 " << initial_spin_Mn/numb_atoms << " spin 2 prod " << initial_spin_2_Mn/numb_atoms <<"\n";
	for (int temp = temp1; temp >= temp2; temp += temp_inc) {
		e_avg = 0;
		spin_avg = 0;
		spin_avg_Mn = 0;
		spin_avg_Mn2 = 0;
		stag_avg = 0;
		flip_count = 0;
		flip_count2 = 0;
		for (int pass = 0; pass < passes; pass++) {
			for (int site = 0; site < numb_atoms; site++) {
				// Flip Spin
				old_spin = atom_spin_list[site];
				spin_same = true;
				e_site = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				while (spin_same == true) {
					spin_rand = unif(rng);
					if (spin_rand <= 0.3333333333333333) { new_spin = -1; }
					else if (spin_rand <= 0.6666666666666666) { new_spin = 0; }
					else { new_spin = 1; }
					if (new_spin != old_spin) { spin_same = false; }
				}
				atom_spin_list[site] = new_spin;
				e_site_new = evalSiteEnergySpin(site, rule_map_spin, atom_spin_list, atom_species_list, neighbor_index_list, neighbor_dist_list);
				e_flip = e_site_new - e_site;
				spin_flip = new_spin - old_spin;
				spin_flipMn2 = calcMag2Diff(site, old_spin, atom_spin_list, neighbor_index_list, neighbor_dist_list);
				stag_flip = stagMagDiff(site, new_spin, old_spin, sim_cell);
				if (e_flip < 0) { flip_count += 1; }
				else {
					keep_rand = unif(rng);
					keep_prob = exp(-1 / (Kb * temp) * (e_flip));
					if (keep_rand < keep_prob) { flip_count2 += 1; }
					else { atom_spin_list[site] = old_spin; e_flip = 0; spin_flip = 0; spin_flipMn2 = 0; stag_flip = 0; }
				}
				initial_enrg += e_flip; 
				initial_spin += spin_flip;
				initial_spin_2_Mn += spin_flipMn2;
				initial_stag_mag += stag_flip;
				if (atom_species_list[site] == 1) { initial_spin_Mn += spin_flip; }
				if (pass >= passes * .2) {
					e_avg += initial_enrg;
					spin_avg += initial_spin;
					spin_avg_Mn2 += initial_spin_2_Mn;
					stag_avg += initial_stag_mag;
					if (atom_species_list[site] == 1) {
						spin_avg_Mn += initial_spin_Mn;
					}
				}
			}
		}
		cout << temp << ": " << e_avg / (pow(numb_atoms,2) * passes * .8) - 0 << ": " << numb_atoms <<": " << initial_spin_2_Mn / numb_atoms << ": " << e_avg << "\n";
		Output << temp;
		Output << " , ";
		Output << e_avg / (pow(numb_atoms, 2) * passes  *.8) - 0;
		Output << " , ";
		Output << abs(spin_avg) / (pow(numb_atoms,2) * passes * .8);
		Output << " , ";
		Output << abs(spin_avg_Mn) / (pow(sim_cell.species_numbs[1],2) * passes * .8);
		Output << " , ";
		Output << spin_avg_Mn2 / (numb_atoms * sim_cell.species_numbs[1] * passes * .8);
		Output << " , ";
		Output << stag_avg / (numb_atoms * sim_cell.species_numbs[1] * passes * .8);
		Output << " , ";
		Output << flip_count;
		Output << " , ";
		Output << flip_count2;
		Output << "\n";
	}
	writeSuperCell(atom_species_list, atom_spin_list, sim_cell);
	Output.close();
}

// Write the current state of the simulation cell to a txt file. Mostly for debugging.
void writeSuperCell(vector<int> &atom_species, vector<int>& atom_spins, SimCell &sim_cell) {
	ofstream OUT_file;
	OUT_file.open("OUTCAR.txt");
	if (OUT_file.is_open()) {
		OUT_file << "Ni2Mn(1-X)In(X)\n" << 5.76 * sim_cell.sup_cell[0] <<"\n";
		OUT_file << sim_cell.unit_LC[0] << " 0 0\n";
		OUT_file << "0 " << sim_cell.unit_LC[1] << " 0\n" << "0 0 ";
		OUT_file << sim_cell.unit_LC[2] << "\n";
		OUT_file << "Ni Mn In\n";
		OUT_file << sim_cell.species_numbs[0] << " " << sim_cell.species_numbs[1] << " " << sim_cell.species_numbs[2] << "\n" << "Direct \n";
		for (int i = 0; i < sim_cell.numb_atoms; i++) {
			OUT_file << sim_cell.atom_list[i].pos[0] / (sim_cell.sup_cell[0] * sim_cell.unit_LC[0]) << " ";
			OUT_file << sim_cell.atom_list[i].pos[1] / (sim_cell.sup_cell[1] * sim_cell.unit_LC[1]) << " ";
			OUT_file << sim_cell.atom_list[i].pos[2] / (sim_cell.sup_cell[2] * sim_cell.unit_LC[2]) << " ";
			OUT_file << atom_species[i] << " " << atom_spins[i] << "\n";
		}
	}
	OUT_file.close();
}