#include <iostream>
#include <string>
#include "rule.h"
#include "monte_carlo.h"
#include "Session.h"
using namespace std;
// The MC code
// The paramiters for the MC run ( Composition, Simulation cell size, ect...) are all set in INPUT.txt
// These settings are read in and stored in the SimCell object. The SimCell object is then reffrenced to set up the actual MC run.
SimCell sim_cell;
int main(void) {
	Session ses("INPUT.txt"); // settings for MC run from input file
	vector<Rule> mc_rules; // Initalize the "Rules" that will be used to properly apply the ECI to each site in the simulation cell
	fillRuleList(mc_rules, ses.rules_file, 0); // Populate mc_rules. Rules are defined as follows, [species (0,1,2 for Ni,Mn,In)]:[NA if monomer,dist of dymer, AB,AC,BC if trimer]:no-spin(0) or spin(1): value
	cout << ses.shape[0] << ',' << ses.shape[1] << ',' << ses.shape[2] << '\n';
	cout << ses.atom_numbs[0] << ',' << ses.atom_numbs[1] << ',' << ses.atom_numbs[2] << '\n';
	vector<float> dist_list{ 0.0 }; // Initalize list of discances used in mc_rules
	fillDistList(dist_list, mc_rules); // Fill the dist_list 
	// Create the simulation cell object. Arguments (POSCAR_file, dist_list, shape, species numbs, cutoff (currently unused), sim_type (also unused), phase_init (aust/mart), spin_init (AFM/FM/RAND), species_init (Ordered/Random), bool use_poscar) 
	sim_cell.initSimCell(ses.structure_file, dist_list, ses.shape, ses.atom_numbs, 1, ses.sim_type, ses.phase_init, ses.spin_init, ses.species_init, ses.use_poscar); // Create and initalize the simulation cell
  	cout << "Testing" << '\n';
	cout << "beginning MC" << '\n';
	// All of the following MC functinos have the same argument format and represent diffrent implemetations of the MC algorithm: 
	//	(Number of MC moves on each site in the simulation cell per temperature increment, Initial temp, final temp, temp inc, sim_cell object, mc_rules object)
	if (ses.algo == 3) { runMetropolis3(ses.numb_passes, ses.start_temp, ses.end_temp, ses.temp_inc, sim_cell, mc_rules); }
	else if (ses.algo == 4) { runMetropolis4(ses.numb_passes, ses.start_temp, ses.end_temp, ses.temp_inc, sim_cell, mc_rules); }
	else if (ses.algo == 5) { runMetropolis5(ses.numb_passes, ses.start_temp, ses.end_temp, ses.temp_inc, sim_cell, mc_rules); }
	else if (ses.algo == 6) { runMetropolis6(ses.numb_passes, ses.start_temp, ses.end_temp, ses.temp_inc, ses.thermostat_2, ses.eq_passes, sim_cell, mc_rules); }
	else if (ses.algo == 7) { runMetropolis7(ses.numb_passes, ses.start_temp, ses.end_temp, ses.temp_inc, ses.sro_target, ses.thermostat_2, ses.eq_passes, sim_cell, mc_rules); }
	cout << "MC Finished...";
	int exit;
	std::cin >> exit;
}