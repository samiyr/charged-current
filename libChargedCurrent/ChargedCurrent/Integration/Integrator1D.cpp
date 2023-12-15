#ifndef INTEGRATOR_1D_H
#define INTEGRATOR_1D_H

#include <iostream>

#include <gsl/gsl_integration.h>

struct Integrator1D {
	/// The purpose of this functor is to enable passing arbitrary callable types to the integration routine.
	/// By default, both Cuba and GSL integration libraries take a function pointer as the integrand.
	/// Non-capturing lambdas can be implicitly converted to a function pointer, but capturing lambdas
	/// cannot be. To support capturing lambdas, this trick with a functor is used. The idea is to pass a functor
	/// instance as the arbitrary parameter void* to the integration routine. The actual integration parameters are then
	/// stored inside the functor, alongside the actual integrand function and potentially other data as well.
	/// At the integration callsite, instead of passing the actual integrand, the invoke method is passed instead.
	/// Inside the invoke method, the functor instance previously passed as the arbitrary parameter of type void*
	/// is obtained. Using the functor instance the actual integrand function stored in it can then be called.
	template <typename Function>
	struct GSL1DFunctor {
		Function function;
		void *params;

		static double invoke(double x, void *p) {
			GSL1DFunctor *functor = static_cast<GSL1DFunctor *>(p);
			return functor->function(x, functor->params);
		}
	};

	struct Result {
		const double value;
		const double error;
		const std::size_t nevals;
		const int return_code;

		friend std::ostream &operator <<(std::ostream &os, const Result &r) {
			return os << r.value << " +- " << r.error;
		}
	};

	void *params;
	const std::size_t interval_count;

	Integrator1D(
		void *params = nullptr,
		const std::size_t interval_count = 100
	) : params(params), interval_count(interval_count) {
		workspace = gsl_integration_cquad_workspace_alloc(interval_count);
	}

	~Integrator1D() {
		gsl_integration_cquad_workspace_free(workspace);
	}

	template <typename Integrand>
	Result integrate(const Integrand &integrand, const double lower, const double upper, const double epsabs, const double epsrel) const {
		double integral = 0.0;
		double error = 0.0;
		std::size_t nevals = 0;

		GSL1DFunctor<Integrand> functor { integrand, params };
		gsl_function function;

		function.function = &GSL1DFunctor<Integrand>::invoke;
		function.params = &functor;

		const int return_code = gsl_integration_cquad(&function, lower, upper, epsabs, epsrel, workspace, &integral, &error, &nevals);

		return Result { integral, error, nevals, return_code };
	}

	private:
	gsl_integration_cquad_workspace *workspace;
};

#endif