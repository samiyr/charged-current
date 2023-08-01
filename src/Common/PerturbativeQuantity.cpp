#ifndef PERTURBATIVE_QUANTITY_H
#define PERTURBATIVE_QUANTITY_H

#include <iostream>

struct PerturbativeQuantity {
	double lo;
	double nlo;

	friend std::ostream& operator<<(std::ostream &os, const PerturbativeQuantity &r) {
		return os << r.lo << ", " << r.nlo;
	}
};

const PerturbativeQuantity operator+(const PerturbativeQuantity lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs.lo + rhs.lo, lhs.nlo + rhs.nlo};
}
const PerturbativeQuantity operator-(const PerturbativeQuantity lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs.lo - rhs.lo, lhs.nlo - rhs.nlo};
}
const PerturbativeQuantity operator*(const PerturbativeQuantity lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs.lo * rhs.lo, lhs.nlo * rhs.nlo};
}
const PerturbativeQuantity operator/(const PerturbativeQuantity lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs.lo / rhs.lo, lhs.nlo / rhs.nlo};
}

const PerturbativeQuantity operator+(const double lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs + rhs.lo, lhs + rhs.nlo};
}
const PerturbativeQuantity operator-(const double lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs - rhs.lo, lhs - rhs.nlo};
}
const PerturbativeQuantity operator*(const double lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs * rhs.lo, lhs * rhs.nlo};
}
const PerturbativeQuantity operator/(const double lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs / rhs.lo, lhs / rhs.nlo};
}

const PerturbativeQuantity operator+(const PerturbativeQuantity lhs, const double& rhs) {
	return PerturbativeQuantity {lhs.lo + rhs, lhs.nlo + rhs};
}
const PerturbativeQuantity operator-(const PerturbativeQuantity lhs, const double& rhs) {
	return PerturbativeQuantity {lhs.lo - rhs, lhs.nlo - rhs};
}
const PerturbativeQuantity operator*(const PerturbativeQuantity lhs, const double& rhs) {
	return PerturbativeQuantity {lhs.lo * rhs, lhs.nlo * rhs};
}
const PerturbativeQuantity operator/(const PerturbativeQuantity lhs, const double& rhs) {
	return PerturbativeQuantity {lhs.lo / rhs, lhs.nlo / rhs};
}

#endif