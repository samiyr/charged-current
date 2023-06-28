#ifndef CUBA_INTEGRATOR_H
#define CUBA_INTEGRATOR_H

#include <vector>
#include <cassert>
#include <cuba.h>
#include <iostream>

enum class CubaMethod {
	Vegas, Suave, Cuhre
};

template <typename Integrand>
struct CubaIntegrator {
	template <typename Function>
	struct Functor {
		Function function;
		void *params;
		const double *lower;
		const double *upper;

		static int invoke(const int *ndim, const double x[], const int *ncomp, double f[], void *userdata) {
			Functor *functor = static_cast<Functor *>(userdata);
			const size_t dim = size_t(*ndim);

			const double *lower_limits = functor->lower;
			const double *upper_limits = functor->upper;

			const size_t size = dim * sizeof(double);
			double *input = new double[size];
			double jacobian = 1.0;
			for (size_t i = 0; i < dim; i++) {
				const double lower = lower_limits[i];
				const double upper = upper_limits[i];

				const double range = upper - lower;
				jacobian *= range;
				input[i] = lower + x[i] * range;
			}

			// #pragma clang diagnostic push
			// #pragma clang diagnostic ignored "-Wsign-conversion"
			const double function_value = jacobian * functor->function(input, dim, functor->params);
			// #pragma clang diagnostic pop
			f[0] = function_value;
			return 0;
		}
	};

	struct Result {
		const double value;
		const double error;
		const double chi_squared_probability;

		friend std::ostream &operator <<(std::ostream &os, const Result &r) {
			return os << r.value << " +- " << r.error << " (" << r.chi_squared_probability << ")";
		}
	};

	struct VegasParameters {
		int starting_iterations = 1000;
		int continuing_iterations = 1000;
		int batch_size = 1000;
		int grid_number = 0;
	};

	struct SuaveParameters {
		int new_subdivision_evaluations = 500;
		int minimum_samples = 20;
		double flatness = 100.0;
	};

	struct CuhreParameters {
		int cubature_rule = 0;
	};

	const Integrand integrand;
	const std::vector<double> lower;
	const std::vector<double> upper;
	void *params;
	const CubaMethod method;

	int flags = 0;
	int seed = 0;
	int minimum_evaluations = 0;
	int maximum_evaluations = 1'000'000;

	double maximum_relative_error = 1e-2;
	double maximum_absolute_error = 1e-12;

	VegasParameters vegas;
	SuaveParameters suave;
	CuhreParameters cuhre;

	// TODO: implement this
	bool verbose = false;

	CubaIntegrator(
		const Integrand integrand, 
		const std::vector<double> lower, 
		const std::vector<double> upper, 
		void *params = nullptr,
		const CubaMethod method = CubaMethod::Suave
	) : integrand(integrand), lower(lower), upper(upper), params(params), method(method) { }

	Result integrate() {
		switch (method) {
		case CubaMethod::Vegas: return vegas_integrate();
		case CubaMethod::Suave: return suave_integrate();
		case CubaMethod::Cuhre: return cuhre_integrate();
		}
	}

	private:
	static void set_cuba_cores(const int cores, const int max) {
		cubacores(&cores, &max);
	}
	static void disable_cuba_parallelization() {
		set_cuba_cores(0, 0);
	}

	Result vegas_integrate() {
		disable_cuba_parallelization();

		assert(lower.size() == upper.size());
		assert(lower.size() > 0);

		#pragma clang diagnostic push
		#pragma clang diagnostic ignored "-Wshorten-64-to-32"
		const int dim = int(lower.size());
		#pragma clang diagnostic pop

		double integral[] = {0.0};
		double error[] = {0.0};
		double chi_squared_probability[] = {0.0};

		Functor<Integrand> functor { integrand, params, lower.data(), upper.data() };
		auto cuba_integrand = &Functor<Integrand>::invoke;

		int neval, fail;

		Vegas(
			dim, 1, cuba_integrand, &functor, 1, 
			maximum_relative_error, maximum_absolute_error,
			flags, seed, minimum_evaluations, maximum_evaluations,
			vegas.starting_iterations, vegas.continuing_iterations, vegas.batch_size, vegas.grid_number,
			nullptr, nullptr, &neval, &fail,
			integral, error, chi_squared_probability
		);

		return Result { integral[0], error[0], chi_squared_probability[0] };
	}

	Result suave_integrate() {
		disable_cuba_parallelization();

		assert(lower.size() == upper.size());
		assert(lower.size() > 0);

		#pragma clang diagnostic push
		#pragma clang diagnostic ignored "-Wshorten-64-to-32"
		const int dim = int(lower.size());
		#pragma clang diagnostic pop

		double integral[] = {0.0};
		double error[] = {0.0};
		double chi_squared_probability[] = {0.0};

		Functor<Integrand> functor { integrand, params, lower.data(), upper.data() };
		auto cuba_integrand = &Functor<Integrand>::invoke;

		int nregions, neval, fail;

		Suave(
			dim, 1, cuba_integrand, &functor, 1, 
			maximum_relative_error, maximum_absolute_error,
			flags, seed, minimum_evaluations, maximum_evaluations,
			suave.new_subdivision_evaluations, suave.minimum_samples, suave.flatness,
			nullptr, nullptr, &nregions, &neval, &fail,
			integral, error, chi_squared_probability
		);

		return Result { integral[0], error[0], chi_squared_probability[0] };
	}

	Result cuhre_integrate() {
		disable_cuba_parallelization();

		assert(lower.size() == upper.size());
		assert(lower.size() > 0);

		#pragma clang diagnostic push
		#pragma clang diagnostic ignored "-Wshorten-64-to-32"
		const int dim = int(lower.size());
		#pragma clang diagnostic pop

		double integral[] = {0.0};
		double error[] = {0.0};
		double chi_squared_probability[] = {0.0};

		Functor<Integrand> functor { integrand, params, lower.data(), upper.data() };
		auto cuba_integrand = &Functor<Integrand>::invoke;

		int nregions, neval, fail;

		Cuhre(
			dim, 1, cuba_integrand, &functor, 1,
			maximum_relative_error, maximum_absolute_error,
			flags, minimum_evaluations, maximum_evaluations,
			cuhre.cubature_rule,
			nullptr, nullptr, &nregions, &neval, &fail,
			integral, error, chi_squared_probability
		);

		return Result { integral[0], error[0], chi_squared_probability[0] };
	}
};

#endif