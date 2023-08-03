#ifndef CUBA_INTEGRATOR_H
#define CUBA_INTEGRATOR_H

#include <vector>
#include <cassert>
#include <cuba.h>
#include <iostream>
#include <gsl/gsl_monte_vegas.h>
#include "Integration/IntegrationParameters.cpp"

enum class IntegrationMethod {
	CubaVegas, CubaSuave, CubaCuhre, GSLVegas
};

template <typename Integrand>
struct Integrator {
	template <typename Function>
	struct CubaFunctor {
		Function function;
		void *params;
		const double *lower;
		const double *upper;
		const double jacobian;

		static int invoke(const int *ndim, const double x[], [[maybe_unused]] const int *ncomp, double f[], void *userdata) {
			CubaFunctor *functor = static_cast<CubaFunctor *>(userdata);
			const std::size_t dim = std::size_t(*ndim);

			const double *lower_limits = functor->lower;
			const double *upper_limits = functor->upper;
			const double jacobian = functor->jacobian;

			double *input = new double[dim];

			for (std::size_t i = 0; i < dim; i++) {
				const double lower = lower_limits[i];
				const double upper = upper_limits[i];
				const double range = upper - lower;

				input[i] = lower + x[i] * range;
			}

			const double function_value = jacobian * functor->function(input, dim, functor->params);
			f[0] = function_value;

			return 0;
		}
	};

	template <typename Function>
	struct GSLFunctor {
		Function function;
		void *params;

		static double invoke(double input[], std::size_t dim, void *p) {
			GSLFunctor *functor = static_cast<GSLFunctor *>(p);
			return functor->function(input, dim, functor->params);
		}
	};

	struct Result {
		const double value;
		const double error;
		const double reliability;

		friend std::ostream &operator <<(std::ostream &os, const Result &r) {
			return os << r.value << " +- " << r.error << " (" << r.reliability << ")";
		}
	};

	const Integrand integrand;
	const std::vector<double> lower;
	const std::vector<double> upper;
	void *params;
	const IntegrationMethod method;

	IntegrationParameters::Cuba cuba;
	IntegrationParameters::GSL gsl;

	int verbose = false;

	Integrator(
		const Integrand integrand, 
		const std::vector<double> lower, 
		const std::vector<double> upper, 
		const IntegrationParameters parameters,
		void *params = nullptr,
		const IntegrationMethod method = IntegrationMethod::CubaSuave
	) : integrand(integrand), lower(lower), upper(upper), params(params), method(method), cuba(parameters.cuba), gsl(parameters.gsl) {
		if (method == IntegrationMethod::GSLVegas) {
			gsl_rng_env_setup();
			rng_type = gsl_rng_default;
			rng = gsl_rng_alloc(rng_type);
		}
	}

	Integrator(
		const Integrand integrand,
		const std::vector<double> lower,
		const std::vector<double> upper,
		void *params = nullptr,
		const IntegrationMethod method = IntegrationMethod::CubaSuave
	) : Integrator(integrand, lower, upper, IntegrationParameters(), params, method) { }

	~Integrator() {
		if (method == IntegrationMethod::GSLVegas) {
			gsl_rng_free(rng);
		}
	}

	Result integrate() {
		switch (method) {
		case IntegrationMethod::CubaVegas: return cuba_vegas_integrate();
		case IntegrationMethod::CubaSuave: return cuba_suave_integrate();
		case IntegrationMethod::CubaCuhre: return cuba_cuhre_integrate();
		case IntegrationMethod::GSLVegas: return gsl_vegas_integrate();
		}
	}

	private:
	const gsl_rng_type *rng_type;
	gsl_rng *rng = nullptr;
	gsl_monte_function gsl_function;

	static void set_cuba_cores(const int cores, const int max) {
		cubacores(&cores, &max);
	}
	static void disable_cuba_parallelization() {
		set_cuba_cores(0, 0);
	}

	int cuba_flags() const {
		int flag_value = 0;

		// Bits 0 and 1
		if (verbose < 0 || verbose > 3) {
			throw std::runtime_error("Cuba verbose flag must be an integer between 0 and 3, was given " + std::to_string(verbose) + " instead");
		}
		flag_value += verbose;
		// Bit 2
		flag_value += 4 * (!cuba.use_all_samples);
		// Bit 3
		flag_value += 8 * (!cuba.apply_smoothing);
		// Bit 4
		flag_value += 16 * (!cuba.delete_state_file_after_integration);
		// Bit 5
		flag_value += 32 * (cuba.reset_state);
		// Bits 8-31
		flag_value += 64 * cuba.random_number_generator_level;

		return flag_value;
	}

	double compute_jacobian(const std::vector<double> &lower_bounds, const std::vector<double> &upper_bounds) {
		double jacobian = 1.0;

		for (auto [lower, upper] : std::views::zip(lower_bounds, upper_bounds)) {
			jacobian *= upper - lower;
		}

		return jacobian;
	}

	Result cuba_vegas_integrate() {
		disable_cuba_parallelization();

		assert(lower.size() == upper.size());
		assert(lower.size() > 0);

		const int dim = Conversion::size_to_int(lower.size());

		double integral[] = {0.0};
		double error[] = {0.0};
		double chi_squared_probability[] = {0.0};

		CubaFunctor<Integrand> functor { integrand, params, lower.data(), upper.data(), compute_jacobian(lower, upper) };
		auto cuba_integrand = &CubaFunctor<Integrand>::invoke;

		int neval, fail;

		Vegas(
			dim, 1, cuba_integrand, &functor, 1, 
			cuba.maximum_relative_error, cuba.maximum_absolute_error,
			cuba_flags(), cuba.seed, cuba.minimum_evaluations, cuba.maximum_evaluations,
			cuba.vegas.starting_iterations, cuba.vegas.continuing_iterations, cuba.vegas.batch_size, cuba.vegas.grid_number,
			nullptr, nullptr, &neval, &fail,
			integral, error, chi_squared_probability
		);

		return Result { integral[0], error[0], chi_squared_probability[0] };
	}

	Result cuba_suave_integrate() {
		disable_cuba_parallelization();

		assert(lower.size() == upper.size());
		assert(lower.size() > 0);

		const int dim = Conversion::size_to_int(lower.size());

		double integral[] = {0.0};
		double error[] = {0.0};
		double chi_squared_probability[] = {0.0};

		CubaFunctor<Integrand> functor { integrand, params, lower.data(), upper.data(), compute_jacobian(lower, upper) };
		auto cuba_integrand = &CubaFunctor<Integrand>::invoke;

		int nregions, neval, fail;

		Suave(
			dim, 1, cuba_integrand, &functor, 1, 
			cuba.maximum_relative_error, cuba.maximum_absolute_error,
			cuba_flags(), cuba.seed, cuba.minimum_evaluations, cuba.maximum_evaluations,
			cuba.suave.new_subdivision_evaluations, cuba.suave.minimum_samples, cuba.suave.flatness,
			nullptr, nullptr, &nregions, &neval, &fail,
			integral, error, chi_squared_probability
		);

		return Result { integral[0], error[0], chi_squared_probability[0] };
	}

	Result cuba_cuhre_integrate() {
		disable_cuba_parallelization();

		assert(lower.size() == upper.size());
		assert(lower.size() > 0);

		const int dim = Conversion::size_to_int(lower.size());

		double integral[] = {0.0};
		double error[] = {0.0};
		double chi_squared_probability[] = {0.0};

		CubaFunctor<Integrand> functor { integrand, params, lower.data(), upper.data(), compute_jacobian(lower, upper) };
		auto cuba_integrand = &CubaFunctor<Integrand>::invoke;

		int nregions, neval, fail;

		Cuhre(
			dim, 1, cuba_integrand, &functor, 1,
			cuba.maximum_relative_error, cuba.maximum_absolute_error,
			cuba_flags(), cuba.minimum_evaluations, cuba.maximum_evaluations,
			cuba.cuhre.cubature_rule,
			nullptr, nullptr, &nregions, &neval, &fail,
			integral, error, chi_squared_probability
		);

		return Result { integral[0], error[0], chi_squared_probability[0] };
	}

	Result gsl_vegas_integrate() {
		assert(lower.size() == upper.size());
		assert(lower.size() > 0);

		const auto dim = lower.size();
		double integral = 0.0;
		double error = 0.0;

		GSLFunctor<Integrand> functor { integrand, params };
		gsl_function.f = &GSLFunctor<Integrand>::invoke;
		gsl_function.dim = dim;
		gsl_function.params = &functor;

		gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(dim);

		gsl_monte_vegas_params vegas_params;
		gsl_monte_vegas_params_get(state, &vegas_params);
		vegas_params.verbose = verbose - 1;
		gsl_monte_vegas_params_set(state, &vegas_params);

		std::vector<double> lower_copy = lower;
		std::vector<double> upper_copy = upper;

		double *lower_ptr = lower_copy.data();
		double *upper_ptr = upper_copy.data();

		if (gsl.grid_warmup) {
			gsl_monte_vegas_integrate(&gsl_function, lower_ptr, upper_ptr, dim, std::max(gsl.points / 1000, std::size_t(10)), rng, state, &integral, &error);
		}

		unsigned int iteration = 0;
		bool iteration_limit_reached = false;
		double chi_squared = gsl_monte_vegas_chisq(state);

		double best_chi_squared = chi_squared;
		double best_integral = integral;
		double best_error = error;

		// std::cout << "Warm-up integral: " << integral << " +- " << error << " (" << chi_squared << ")" << IO::endl;

		while (abs(chi_squared - 1.0) >= gsl.max_chi_squared_deviation || abs(error / integral) >= gsl.max_relative_error) {
			iteration++;
			if (iteration > gsl.iter_max) {
				iteration_limit_reached = true;
				break;
			}
			gsl_monte_vegas_integrate(&gsl_function, lower_ptr, upper_ptr, dim, gsl.points, rng, state, &integral, &error);

			chi_squared = gsl_monte_vegas_chisq(state);
			// std::cout << "Iteration " << iteration << " integral: " << integral << " +- " << error << " (" << chi_squared << ")" << IO::endl;
			if (abs(chi_squared - 1.0) < abs(best_chi_squared - 1.0)) {
				best_chi_squared = chi_squared;
				best_integral = integral;
				best_error = error;
			}
		}

		double final_integral = integral;
		double final_error = error;
		double final_chi_squared = gsl_monte_vegas_chisq(state);
		// std::cout << "Final iteration integral: " << final_integral << " +- " << final_error << " (" << final_chi_squared << ")" << IO::endl;

		if (iteration_limit_reached) {
			final_integral = best_integral;
			final_error = best_error;
			final_chi_squared = best_chi_squared;
		}

		// std::cout << "Final integral: " << final_integral << " +- " << final_error << " (" << final_chi_squared << ")" << IO::endl << IO::endl;
		gsl_monte_vegas_free(state);

		return Result {final_integral, final_error, final_chi_squared};
	}
};

#endif