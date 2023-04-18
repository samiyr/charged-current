#ifndef DIS_EVALUATION_H
#define DIS_EVALUATION_H

#include "PDFEvaluation.cpp"
#include "Row.cpp"
#include <fstream>
#include <string>
#include "DISComputation.cpp"

struct DIS {
	PDFEvaluation pdf;

	const std::vector<FlavorType> active_flavors;
	const std::vector<FlavorType> active_antiflavors;

	bool parallelize = true;
	int number_of_threads = 8;

	double sqrt_s;
	double s;
	double y_max = 1.0;

	const size_t points;

	DIS (const double _sqrt_s, const std::vector<FlavorType> _active_flavors, const std::string pdf_set, const size_t _points)
	: sqrt_s(_sqrt_s), s(_sqrt_s * _sqrt_s), 
	active_flavors(_active_flavors), 
	active_antiflavors(conjugate_flavors(_active_flavors)), 
	pdf(pdf_set, 0),
	points(_points) 
	{ }

	double compute_structure_function(DISComputation::StructureFunction F, const double x, const double Q2) {
		DISComputation dis(sqrt_s, active_flavors, active_antiflavors, pdf.pdf_set_name, points);
		return dis.structure_function(F, x, Q2);
	}
	double F2(const double x, const double Q2) {
		return compute_structure_function(DISComputation::StructureFunction::F2, x, Q2);
	}
	double FL(const double x, const double Q2) {
		return compute_structure_function(DISComputation::StructureFunction::FL, x, Q2);
	}
	double xF3(const double x, const double Q2) {
		return compute_structure_function(DISComputation::StructureFunction::xF3, x, Q2);
	}

	void compute_structure_function(DISComputation::StructureFunction F, const std::vector<double> x_bins, const double Q2, const std::string filename) {
		const size_t step_count = x_bins.size();
		std::vector<Row> result;
		result.reserve(step_count);

		int calculated_values = 0;

		std::ofstream file(filename);

		#pragma omp parallel for if(parallelize) num_threads(number_of_threads)
		for (int i = 0; i < step_count; i++) {
			DISComputation dis(sqrt_s, active_flavors, active_antiflavors, pdf.pdf_set_name, points);
			const double x = x_bins[i];
			const double value = dis.structure_function(F, x, Q2);
			Row row {x, value, 0, 0, 0};

			#pragma omp critical
			{
				result.push_back(row);
				file << row << std::endl;
				file.flush();

				calculated_values++;
				std::cout << "Calculated value " << calculated_values << " / " << step_count << " (index " << i << ", x " << x << ")" << std::endl;
			}
		}

		file.close();
	}
	void F2(const std::vector<double> x_bins, const double Q2, const std::string filename) {
		return compute_structure_function(DISComputation::StructureFunction::F2, x_bins, Q2, filename);
	}
	void FL(const std::vector<double> x_bins, const double Q2, const std::string filename) {
		return compute_structure_function(DISComputation::StructureFunction::FL, x_bins, Q2, filename);
	}
	void xF3(const std::vector<double> x_bins, const double Q2, const std::string filename) {
		return compute_structure_function(DISComputation::StructureFunction::xF3, x_bins, Q2, filename);
	}

	void cross_section(const std::vector<double> Q2_bins, const std::string filename) {
		const size_t step_count = Q2_bins.size();
		std::vector<Row> result;
		result.reserve(step_count);

		int calculated_values = 0;

		std::ofstream file(filename);

		#pragma omp parallel for if(parallelize) num_threads(number_of_threads)
		for (int i = 0; i < step_count; i++) {
			const double Q2 = Q2_bins[i];
			
			DISComputation dis(sqrt_s, active_flavors, active_antiflavors, pdf.pdf_set_name, points);

			const double x_min = Q2 / (s * y_max);
			const double value = dis.x_integrated_cross_section(Q2, x_min);
			Row row {Q2, value, 0, 0, 0};

			#pragma omp critical
			{
				result.push_back(row);
				file << row << std::endl;
				file.flush();

				calculated_values++;
				std::cout << "Calculated value " << calculated_values << " / " << step_count << " (index " << i << ", Q2 " << Q2 << ")" << std::endl;
			}
		}

		file.close();
	}
};

#endif