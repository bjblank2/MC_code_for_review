#include "Session.h"

Session::Session(string input_file) {
	vector<string> input_lines;
	string input_settings;
	ifstream session_input;
	vector<string> setting;
	session_input.open(input_file);
	if (session_input.is_open()) {
		while (getline(session_input, input_settings))
		{
			input_lines.push_back(input_settings);
		}
		session_input.close();
	}
	else cout << "Unable to open input file";
	for (int i = 0; i < input_lines.size(); i++) {
		setting = split(input_lines[i], " ");
		if (setting[0] == "STRCTURE") { structure_file = setting[2]; }
		else if (setting[0] == "ALGO") { algo = stoi(setting[2]); }
		else if (setting[0] == "RULES_FILE") { rules_file = setting[2]; }
		else if (setting[0] == "SIM_TYPE") { sim_type = setting[2]; }
		else if (setting[0] == "PHASE_INIT") { phase_init = setting[2]; }
		else if (setting[0] == "SPIN_INIT") { spin_init = setting[2]; }
		else if (setting[0] == "SPECIES_INIT") { species_init = setting[2]; }
		else if (setting[0] == "NUMB_PASSES") { numb_passes = stoi(setting[2]); }
		else if (setting[0] == "START_TEMP") { start_temp = stof(setting[2]); }
		else if (setting[0] == "END_TEMP") { end_temp = stof(setting[2]); }
		else if (setting[0] == "TEMP_INC") { temp_inc = stof(setting[2]); }
		else if (setting[0] == "EQ_PASSES") { eq_passes = stoi(setting[2]); }
		else if (setting[0] == "THERMOSTAT_2") { thermostat_2 = stof(setting[2]); }
		else if (setting[0] == "SRO_TARGET") { sro_target = stof(setting[2]); }
		else if (setting[0] == "USE_POSCAR") { 
			if (setting[2] == "TRUE") { use_poscar = true; }
			else { use_poscar = false; }
		}
		else if (setting[0] == "ATOM_NUMBS") {
			atom_numbs.push_back(stoi(setting[2]));
			atom_numbs.push_back(stoi(setting[3]));
			atom_numbs.push_back(stoi(setting[4]));
		}
		else if (setting[0] == "SHAPE") {
			shape[0] = stoi(setting[2]);
			shape[1] = stoi(setting[3]);
			shape[2] = stoi(setting[4]);
		}

	}
}