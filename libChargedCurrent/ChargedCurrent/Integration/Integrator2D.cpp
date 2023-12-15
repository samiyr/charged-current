#ifndef INTEGRATOR_2D_H
#define INTEGRATOR_2D_H

#include <iostream>

#include "Integration/Integrator1D.cpp"

struct Integrator2D {
	void *params;
	const std::size_t interval_count;

	Integrator2D(
		void *params = nullptr,
		const std::size_t interval_count = 10'000
	) : params(params), interval_count(interval_count) { }

	template <typename Integrand>
	Integrator1D::Result integrate(
		const Integrand &integrand, 
		const double lower1, const double upper1, const double lower2, const double upper2, 
		const double epsabs, const double epsrel
	) const {
		const Integrator1D x_integrator(params);
		const Integrator1D y_integrator(params);

		const auto x_integrand = [&](double x, void *) {
			const auto y_integrand = [&integrand, x](double y, void *) {
				return integrand(x, y);
			};
			return y_integrator.integrate(y_integrand, lower2, upper2, epsabs, epsrel).value;
		};
		
		return x_integrator.integrate(x_integrand, lower1, upper1, epsabs, epsrel);
	}
};

#endif