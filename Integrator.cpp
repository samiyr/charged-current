#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <gsl/gsl_monte_vegas.h>
#include <vector>
#include <cassert>
#include <iostream>

class Integrator {
	private:
	const gsl_rng_type *rng_type;
	gsl_rng *rng;
	public:
	double (* integrand)(double [], size_t, void *);
	std::vector<double> lower;
	std::vector<double> upper;
	void *params;

	const size_t points;
	const double max_chi_squared_deviation = 0.1;
	const double max_relative_error = 0.1;
	const int iter_max = 100;

	Integrator(double (* _integrand)(double [], size_t, void *), std::vector<double> _lower, std::vector<double> _upper, const size_t _points, void *_params = NULL) 
	: integrand(_integrand), lower(_lower), upper(_upper), points(_points), params(_params) {
		gsl_rng_env_setup();

		rng_type = gsl_rng_default;
		rng = gsl_rng_alloc(rng_type);
	}

	struct Result {
		double value;
		double error;
		double chi_squared;

		friend std::ostream& operator<<(std::ostream& os, Result const & r) {
			return os << r.value << " +/- " << r.error << " (" << r.chi_squared << ")";
		}
	};

	Result integrate() {
		assert(lower.size() == upper.size());
		assert(lower.size() > 0);

		const auto dim = lower.size();
		double integral, error;

		gsl_monte_function function;

		function.f = integrand;
		function.dim = dim;
		function.params = params;

		gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(dim);

		gsl_monte_vegas_integrate(&function, lower.data(), upper.data(), dim, points, rng, state, &integral, &error);

		int iteration = 0;
		bool iteration_limit_reached = false;

		double best_chi_squared = gsl_monte_vegas_chisq(state);
		double best_integral = integral;
		double best_error = error;

		while (true) {
			iteration++;
			if (iteration > iter_max) {
				iteration_limit_reached = true;
				break;
			}
			gsl_monte_vegas_integrate(&function, lower.data(), upper.data(), dim, points, rng, state, &integral, &error);

			const double chi_squared = gsl_monte_vegas_chisq(state);
			if (abs(chi_squared - 1.0) < abs(best_chi_squared - 1.0)) {
				best_chi_squared = chi_squared;
				best_integral = integral;
				best_error = error;
			}

			if (abs(gsl_monte_vegas_chisq(state) - 1.0) < max_chi_squared_deviation && abs(error / integral) < max_relative_error) {
				break;
			}
		}

		double final_integral = integral;
		double final_error = error;
		double final_chi_squared = gsl_monte_vegas_chisq(state);

		if (iteration_limit_reached) {
			final_integral = best_integral;
			final_error = best_error;
			final_chi_squared = best_chi_squared;
		}

		return Result {final_integral, final_error, final_chi_squared};
	}
};

#endif