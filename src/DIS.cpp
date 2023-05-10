#ifndef DIS_EVALUATION_H
#define DIS_EVALUATION_H

// #include "PDFInterface.cpp"
#include "Row.cpp"
#include <fstream>
#include <string>
#include "DISComputation.cpp"

template <typename PDFInterface>
struct DIS {
	const double sqrt_s;
	const double s;
	double y_max = 1.0;

	const std::vector<FlavorType> active_flavors;

	PDFInterface pdf;

	bool parallelize = true;
	int number_of_threads = 8;

	const size_t points;
	double max_chi_squared_deviation = 0.2;
	double max_relative_error = 1e-5;
	unsigned int iter_max = 10;

	const Process process;

	DIS (const double _sqrt_s, const std::vector<FlavorType> _active_flavors, const PDFInterface _pdf, const size_t _points, const Process _process)
	: sqrt_s(_sqrt_s), s(_sqrt_s * _sqrt_s), 
	active_flavors(_active_flavors), 
	pdf(_pdf),
	points(_points),
	process(_process)
	{ }

	PerturbativeResult compute_structure_function(StructureFunction F, const double x, const double Q2) {
		DISComputation dis(sqrt_s, active_flavors, pdf, points, max_chi_squared_deviation, max_relative_error, iter_max, process);
		return dis.structure_function(F, x, Q2);
	}
	PerturbativeResult F2(const double x, const double Q2) {
		return compute_structure_function(StructureFunction::F2, x, Q2);
	}
	PerturbativeResult FL(const double x, const double Q2) {
		return compute_structure_function(StructureFunction::FL, x, Q2);
	}
	PerturbativeResult xF3(const double x, const double Q2) {
		return compute_structure_function(StructureFunction::xF3, x, Q2);
	}

	// void compute_structure_function(StructureFunction F, const std::vector<double> x_bins, const double Q2, const std::string filename) {
	// 	const size_t step_count = x_bins.size();
	// 	std::vector<Row> result;
	// 	result.reserve(step_count);

	// 	int calculated_values = 0;

	// 	std::ofstream file(filename);

	// 	#pragma omp parallel for if(parallelize) num_threads(number_of_threads)
	// 	for (int i = 0; i < step_count; i++) {
	// 		DISComputation dis(sqrt_s, active_flavors, active_antiflavors, pdf.pdf_set_name, points, max_chi_squared_deviation, max_relative_error, iter_max, process);
	// 		const double x = x_bins[i];
	// 		const double value = dis.structure_function(F, x, Q2).nlo;
	// 		Row row {x, value, 0, 0, 0};

	// 		#pragma omp critical
	// 		{
	// 			result.push_back(row);
	// 			file << row << std::endl;
	// 			file.flush();

	// 			calculated_values++;
	// 			std::cout << "Calculated value " << calculated_values << " / " << step_count << " (index " << i << ", x " << x << ")" << std::endl;
	// 		}
	// 	}

	// 	file.close();
	// }
	void compute_all_structure_function(const std::vector<double> x_bins, const std::vector<double> Q2_bins, const std::string filename) {
		const size_t x_step_count = x_bins.size();
		const size_t Q2_step_count = Q2_bins.size();

		int calculated_values = 0;

		std::ofstream file(filename);

		#pragma omp parallel for if(parallelize) num_threads(number_of_threads) collapse(2)
		for (size_t i = 0; i < x_step_count; i++) {
			for (size_t j = 0; j < Q2_step_count; j++) {
				DISComputation dis(sqrt_s, active_flavors, pdf, points, max_chi_squared_deviation, max_relative_error, iter_max, process);
				const double x = x_bins[i];
				const double Q2 = Q2_bins[j];
				
				const PerturbativeResult value_F2 = dis.structure_function(StructureFunction::F2, x, Q2);
				const PerturbativeResult value_FL = dis.structure_function(StructureFunction::FL, x, Q2);
				const PerturbativeResult value_xF3 = dis.structure_function(StructureFunction::xF3, x, Q2);
				const PerturbativeResult value_F1 = (value_F2  - value_FL) / (2 * x);
				const PerturbativeResult value_F3 = value_xF3 / x;

				#pragma omp critical
				{
					file << x << ", " << Q2 << ", " << 
						value_F1.lo << ", " << value_F2.lo << ", " << value_F3.lo << ", " << value_xF3.lo << ", " << value_FL.lo << ", " <<
						value_F1.nlo << ", " << value_F2.nlo << ", " << value_F3.nlo << ", " << value_xF3.nlo << ", " << value_FL.nlo
						<< std::endl;
					file.flush();

					calculated_values++;
					std::cout << "Calculated value " << calculated_values << " / " << x_step_count * Q2_step_count << " (x = " << x << ", Q2 = " << Q2 << ")" << std::endl;
				}
			}
		}

		file.close();
	}
	// void F2(const std::vector<double> x_bins, const double Q2, const std::string filename) {
	// 	return compute_structure_function(StructureFunction::F2, x_bins, Q2, filename);
	// }
	// void FL(const std::vector<double> x_bins, const double Q2, const std::string filename) {
	// 	return compute_structure_function(StructureFunction::FL, x_bins, Q2, filename);
	// }
	// void xF3(const std::vector<double> x_bins, const double Q2, const std::string filename) {
	// 	return compute_structure_function(StructureFunction::xF3, x_bins, Q2, filename);
	// }

	// void cross_section(const std::vector<double> Q2_bins, const std::string filename) {
	// 	const size_t step_count = Q2_bins.size();
	// 	std::vector<Row> result;
	// 	result.reserve(step_count);

	// 	int calculated_values = 0;

	// 	std::ofstream file(filename);

	// 	#pragma omp parallel for if(parallelize) num_threads(number_of_threads)
	// 	for (int i = 0; i < step_count; i++) {
	// 		const double Q2 = Q2_bins[i];
			
	// 		DISComputation dis(sqrt_s, active_flavors, active_antiflavors, pdf.pdf_set_name, points, max_chi_squared_deviation, max_relative_error, iter_max, process);

	// 		const double x_min = Q2 / (s * y_max);
	// 		const double value = dis.x_integrated_cross_section(Q2, x_min);
	// 		Row row {Q2, value, 0, 0, 0};

	// 		#pragma omp critical
	// 		{
	// 			result.push_back(row);
	// 			file << row << std::endl;
	// 			file.flush();

	// 			calculated_values++;
	// 			std::cout << "Calculated value " << calculated_values << " / " << step_count << " (index " << i << ", Q2 " << Q2 << ")" << std::endl;
	// 		}
	// 	}

	// 	file.close();
	// }

	void differential_cross_section(const std::vector<double> x_bins, const std::vector<double> Q2_bins, const std::string filename) {
		const size_t x_step_count = x_bins.size();
		const size_t Q2_step_count = Q2_bins.size();

		int calculated_values = 0;

		std::ofstream file(filename);

		#pragma omp parallel for if(parallelize) num_threads(number_of_threads) collapse(2)
		for (size_t i = 0; i < x_step_count; i++) {
			for (size_t j = 0; j < Q2_step_count; j++) {
				DISComputation dis(sqrt_s, active_flavors, pdf, points, max_chi_squared_deviation, max_relative_error, iter_max, process);
				const double x = x_bins[i];
				const double Q2 = Q2_bins[j];
				
				const PerturbativeResult differential_cs = dis.differential_cross_section(x, Q2);

				#pragma omp critical
				{
					file << x << ", " << Q2 << ", " << differential_cs.lo << ", " << differential_cs.nlo << std::endl;
					file.flush();

					calculated_values++;
					std::cout << "Calculated value " << calculated_values << " / " << x_step_count * Q2_step_count << " (x = " << x << ", Q2 = " << Q2 << ")" << std::endl;
				}
			}
		}

		file.close();
	}
};

#endif