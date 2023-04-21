#ifndef SIDIS_COMPUTATION_H
#define SIDIS_COMPUTATION_H

#include "PDFInterface.cpp"
#include "FFInterface.cpp"
#include <vector>
#include "Flavor.cpp"
#include "PerturbativeResult.cpp"
#include "PDFFF.cpp"

class SIDISComputation {
	public:
	PDFInterface pdf1;
	FFInterface ff1;

	PDFInterface pdf2;
	FFInterface ff2;

	FlavorVector active_flavors;
	FlavorVector active_antiflavors;

	double sqrt_s;
	double s;

	const size_t points;
	const double max_chi_squared_deviation;
	const double max_relative_error;
	const unsigned int iter_max;

	const Process process;

	SIDISComputation (
		const double _sqrt_s, 
		const std::vector<FlavorType> _active_flavors, 
		const std::vector<FlavorType> _active_antiflavors, 
		const std::string _pdf_set_name,
		const std::string _ff_set_name,
		const size_t _points,
		const double _max_chi_squared_deviation, 
		const double _max_relative_error,
		const unsigned int _iter_max,
		const Process _process
	) : sqrt_s(_sqrt_s), 
	s(_sqrt_s * _sqrt_s), 
	active_flavors(_active_flavors), 
	active_antiflavors(_active_antiflavors), 
	pdf1(_pdf_set_name), 
	ff1(_ff_set_name),
	pdf2(_pdf_set_name), 
	ff2(_ff_set_name),
	points(_points), 
	max_chi_squared_deviation(_max_chi_squared_deviation), 
	max_relative_error(_max_relative_error),
	iter_max(_iter_max),
	process(_process) { }

	struct CommonParams {
		PDFInterface &pdf;
		FFInterface &ff;

		const std::vector<FlavorType> &active_flavors;
		const std::vector<FlavorType> &active_antiflavors;

		const std::vector<FlavorType> &active_upper_flavors;
		const std::vector<FlavorType> &active_lower_flavors;
		const std::vector<FlavorType> &active_upper_antiflavors;
		const std::vector<FlavorType> &active_lower_antiflavors;

		const double Q2;
		const double nlo_coefficient;
		const double s;

		const Process &process;
	};

	PerturbativeResult F2(const double x, const double z, const double Q2) {
		const std::vector<FlavorType> active_upper_flavors = upper_flavors(active_flavors);
		const std::vector<FlavorType> active_lower_flavors = lower_flavors(active_flavors);
		const std::vector<FlavorType> active_upper_antiflavors = upper_flavors(active_antiflavors);
		const std::vector<FlavorType> active_lower_antiflavors = lower_flavors(active_antiflavors);
	
		double alpha_s = pdf1.alpha_s(Q2);
		double nlo_coefficient = alpha_s / (2 * M_PI);

		pdf1.evaluate(x, Q2);
		ff1.evaluate(z, Q2);

		const double xq_zq = PDF_FF::xq_zq_sum(pdf1, ff1, active_upper_flavors, active_lower_flavors, process);
	}
};

#endif