#ifndef PERTURBATIVE_RESULT_H
#define PERTURBATIVE_RESULT_H

struct PerturbativeResult {
	double lo;
	double nlo;
};

const PerturbativeResult operator+(const PerturbativeResult lhs, const PerturbativeResult& rhs) {
	return PerturbativeResult {lhs.lo + rhs.lo, lhs.nlo + rhs.nlo};
}
const PerturbativeResult operator-(const PerturbativeResult lhs, const PerturbativeResult& rhs) {
	return PerturbativeResult {lhs.lo - rhs.lo, lhs.nlo - rhs.nlo};
}
const PerturbativeResult operator*(const PerturbativeResult lhs, const PerturbativeResult& rhs) {
	return PerturbativeResult {lhs.lo * rhs.lo, lhs.nlo * rhs.nlo};
}
const PerturbativeResult operator/(const PerturbativeResult lhs, const PerturbativeResult& rhs) {
	return PerturbativeResult {lhs.lo / rhs.lo, lhs.nlo / rhs.nlo};
}

const PerturbativeResult operator+(const double lhs, const PerturbativeResult& rhs) {
	return PerturbativeResult {lhs + rhs.lo, lhs + rhs.nlo};
}
const PerturbativeResult operator-(const double lhs, const PerturbativeResult& rhs) {
	return PerturbativeResult {lhs - rhs.lo, lhs - rhs.nlo};
}
const PerturbativeResult operator*(const double lhs, const PerturbativeResult& rhs) {
	return PerturbativeResult {lhs * rhs.lo, lhs * rhs.nlo};
}
const PerturbativeResult operator/(const double lhs, const PerturbativeResult& rhs) {
	return PerturbativeResult {lhs / rhs.lo, lhs / rhs.nlo};
}

const PerturbativeResult operator+(const PerturbativeResult lhs, const double& rhs) {
	return PerturbativeResult {lhs.lo + rhs, lhs.nlo + rhs};
}
const PerturbativeResult operator-(const PerturbativeResult lhs, const double& rhs) {
	return PerturbativeResult {lhs.lo - rhs, lhs.nlo - rhs};
}
const PerturbativeResult operator*(const PerturbativeResult lhs, const double& rhs) {
	return PerturbativeResult {lhs.lo * rhs, lhs.nlo * rhs};
}
const PerturbativeResult operator/(const PerturbativeResult lhs, const double& rhs) {
	return PerturbativeResult {lhs.lo / rhs, lhs.nlo / rhs};
}

#endif