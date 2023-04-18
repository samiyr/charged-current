#ifndef CROSS_SECTION_H
#define CROSS_SECTION_H

#include <iostream>
#include "DIS.cpp"

int main(int argc, char const *argv[]) {
	DIS dis(318,
		{Flavor::Up, Flavor::Down, Flavor::Charm, Flavor::Strange, Flavor::Top, Flavor::Bottom},
		"CT18NLO",
		100'000
	);
	dis.y_max = 0.9;
	// dis.xF3({1e-3, 1e-2, 1e-1, 0.2, 0.4, 0.5, 0.75, 0.9}, 3000);
	
	dis.cross_section({2e2, 5e2, 1e3, 2e3, 3e3, 5e3, 1e4, 2e4, 3e4}, "output2");

	return 0;
}

#endif