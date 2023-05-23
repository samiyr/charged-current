#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <gsl/gsl_monte_vegas.h>
#include <vector>
#include <cassert>
#include <iostream>

template <typename Function>
struct ParametrizedFunctor {
	Function function;
	void *params;

	constexpr static double invoke(double input[], size_t dim, void *p) {
		ParametrizedFunctor *functor = static_cast<ParametrizedFunctor *>(p);
		return functor->function(input, dim, functor->params);
	}
};

template <typename Integrand>
class Integrator {
	private:
	const gsl_rng_type *rng_type;
	gsl_rng *rng;
	public:
	Integrand integrand;
	std::vector<double> lower;
	std::vector<double> upper;

	const size_t points;
	void *params;
	const double max_chi_squared_deviation;
	const double max_relative_error;
	const unsigned int iter_max;

	bool verbose = false;
	bool grid_warmup = true;

	Integrator(
		Integrand _integrand, 
		std::vector<double> _lower, 
		std::vector<double> _upper, 
		const size_t _points, 
		void *_params, 
		const double _max_chi_squared_deviation, 
		const double _max_relative_error,
		const unsigned int _iter_max
	) : integrand(_integrand), 
	lower(_lower), 
	upper(_upper), 
	points(_points), 
	params(_params),
	max_chi_squared_deviation(_max_chi_squared_deviation),
	max_relative_error(_max_relative_error),
	iter_max(_iter_max) {
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
		double integral = 0.0;
		double error = 0.0;

		gsl_monte_function function;

		ParametrizedFunctor<Integrand> functor { integrand, params };
		function.f = &ParametrizedFunctor<Integrand>::invoke;
		// function.f = integrand;
		function.dim = dim;
		function.params = &functor;
		// function.params = params;

		gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(dim);

		gsl_monte_vegas_params vegas_params;
		gsl_monte_vegas_params_get(state, &vegas_params);
		vegas_params.verbose = int(verbose) - 1;
		gsl_monte_vegas_params_set(state, &vegas_params);

		if (grid_warmup) {
			gsl_monte_vegas_integrate(&function, lower.data(), upper.data(), dim, std::max(points / 1000, size_t(10)), rng, state, &integral, &error);
		}

		int iteration = 0;
		bool iteration_limit_reached = false;
		double chi_squared = gsl_monte_vegas_chisq(state);

		double best_chi_squared = chi_squared;
		double best_integral = integral;
		double best_error = error;
		while (abs(chi_squared - 1.0) >= max_chi_squared_deviation || abs(error / integral) >= max_relative_error) {
			iteration++;
			if (iteration > iter_max) {
				iteration_limit_reached = true;
				break;
			}
			gsl_monte_vegas_integrate(&function, lower.data(), upper.data(), dim, points, rng, state, &integral, &error);

			chi_squared = gsl_monte_vegas_chisq(state);
			if (abs(chi_squared - 1.0) < abs(best_chi_squared - 1.0)) {
				best_chi_squared = chi_squared;
				best_integral = integral;
				best_error = error;
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

		gsl_monte_vegas_free(state);

		return Result {final_integral, final_error, final_chi_squared};
	}
};

#endif