#ifndef ROW_H
#define ROW_H

#include <iostream>

struct Row {
	double x;
	double value;
	double error;
	double chi_squared;
	size_t iterations;

	friend std::ostream& operator<<(std::ostream& os, Row const & r) {
		return os << r.x << "," << r.value << "," << r.error << "," << r.chi_squared << "," << r.iterations;
	}
};

#endif