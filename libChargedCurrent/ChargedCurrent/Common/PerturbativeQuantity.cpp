#ifndef PERTURBATIVE_QUANTITY_H
#define PERTURBATIVE_QUANTITY_H

#include <iostream>

// The perturbative order, from leading order (LO) to next-to-next-to leading order (NNLO).
enum class PerturbativeOrder {
	LO = 0, NLO = 1, NNLO = 2
};

// A quantity with LO, NLO and NNLO contributions. Supports basic arithmetic operations and printing. The NNLO contribution is by default zero.
struct PerturbativeQuantity {
	double lo;
	double nlo;
	double nnlo = 0.0;

	friend std::ostream& operator<<(std::ostream &os, const PerturbativeQuantity &r) {
		return os << r.lo << ", " << r.nlo << ", " << r.nnlo;
	}
};

static inline const PerturbativeQuantity operator+(const PerturbativeQuantity lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs.lo + rhs.lo, lhs.nlo + rhs.nlo, lhs.nnlo + rhs.nnlo};
}
static inline const PerturbativeQuantity operator-(const PerturbativeQuantity lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs.lo - rhs.lo, lhs.nlo - rhs.nlo, lhs.nnlo - rhs.nnlo};
}
static inline const PerturbativeQuantity operator*(const PerturbativeQuantity lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs.lo * rhs.lo, lhs.nlo * rhs.nlo, lhs.nnlo * rhs.nnlo};
}
static inline const PerturbativeQuantity operator/(const PerturbativeQuantity lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs.lo / rhs.lo, lhs.nlo / rhs.nlo, lhs.nnlo / rhs.nnlo};
}

static inline const PerturbativeQuantity operator+(const double lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs + rhs.lo, lhs + rhs.nlo, lhs + rhs.nnlo};
}
static inline const PerturbativeQuantity operator-(const double lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs - rhs.lo, lhs - rhs.nlo, lhs - rhs.nnlo};
}
static inline const PerturbativeQuantity operator*(const double lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs * rhs.lo, lhs * rhs.nlo, lhs * rhs.nnlo};
}
static inline const PerturbativeQuantity operator/(const double lhs, const PerturbativeQuantity& rhs) {
	return PerturbativeQuantity {lhs / rhs.lo, lhs / rhs.nlo, lhs / rhs.nnlo};
}

static inline const PerturbativeQuantity operator+(const PerturbativeQuantity lhs, const double& rhs) {
	return PerturbativeQuantity {lhs.lo + rhs, lhs.nlo + rhs, lhs.nnlo + rhs};
}
static inline const PerturbativeQuantity operator-(const PerturbativeQuantity lhs, const double& rhs) {
	return PerturbativeQuantity {lhs.lo - rhs, lhs.nlo - rhs, lhs.nnlo - rhs};
}
static inline const PerturbativeQuantity operator*(const PerturbativeQuantity lhs, const double& rhs) {
	return PerturbativeQuantity {lhs.lo * rhs, lhs.nlo * rhs, lhs.nnlo * rhs};
}
static inline const PerturbativeQuantity operator/(const PerturbativeQuantity lhs, const double& rhs) {
	return PerturbativeQuantity {lhs.lo / rhs, lhs.nlo / rhs, lhs.nnlo / rhs};
}

#endif