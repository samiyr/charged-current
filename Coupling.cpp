#ifndef COUPLING_H
#define COUPLING_H

#include <cmath>

struct Coupling {
	static double alpha_s(const double Q2, const size_t nf, const double lambda_QCD = Constants::lambda_QCD) {
		const double b0 = Coupling::beta_0(Q2, nf);
		const double b1 = Coupling::beta_1(Q2, nf);
		const double log = Coupling::L(Q2, lambda_QCD);
		const double prefactor = 1.0 / (b0 * log);
		const double term1 = (b1 * std::log(log)) / (b0 * b0 * log);
		return prefactor * (1.0 - term1);
	}

	private:
	static double beta_0(const double Q2, const size_t nf) {
		return (33.0 - 2.0 * nf) / (12 * M_PI);
	}
	static double beta_1(const double Q2, const size_t nf) {
		return (153.0 - 19.0 * nf) / (24 * M_PI * M_PI);
	}
	static double L(const double Q2, const double lambda_QCD) {
		return std::log(Q2 / (lambda_QCD * lambda_QCD));
	}
};


#endif