#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <iostream>
#include "DIS.cpp"

int main(int argc, char const *argv[]) {
	DIS dis(318,
		{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Bottom, Flavor::Top},
		"CT18ANLO",
		20'000,
		Process {Process::Type::NeutrinoToLepton}
	);

	// dis.differential_cross_section({0.002}, {200}, "differential_cross_sections.csv");
	dis.differential_cross_section({0.002, 0.1, 0.2, 0.5}, {200, 500, 1'000, 10'000}, "differential_cross_sections.csv");

	// dis.y_max = 0.9;	
	// dis.cross_section({2e2, 5e2, 1e3, 2e3, 3e3, 5e3, 1e4, 2e4, 3e4}, "output2.csv");

	// dis.compute_all_structure_function({0.002, 0.1, 0.2, 0.5}, {200, 500, 1'000, 10'000}, "structure_functions.csv");
	// dis.compute_all_structure_function({0.002}, {200}, "structure_functions.csv");

	return 0;
}

#endif