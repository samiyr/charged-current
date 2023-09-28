#ifndef ANALYSIS_UTILITY_H
#define ANALYSIS_UTILITY_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <ranges>

#include "Analysis/Parameters.cpp"

#include "PDF/Interfaces/LHASetInterface.cpp"

struct UtilityAnalysis {
	const AnalysisParameters params;

	UtilityAnalysis(const AnalysisParameters params) noexcept : params(params) {}

	void pdf_values(
		const std::vector<double> x_bins, 
		const std::vector<double> Q2_bins, 
		const std::filesystem::path base_output, 
		const std::string comment = "") {
		
		LHASetInterface pdf(params.pdf_set);

		const std::size_t member_count = pdf.size();
		const std::size_t x_step_count = x_bins.size();
		const std::size_t Q2_step_count = Q2_bins.size();
		const std::size_t total_count = member_count * x_step_count * Q2_step_count;

		int calculated_values = 0;

		std::streamsize original_precision = std::cout.precision();

		for (typename decltype(pdf)::size_type member_index = 0; member_index < member_count; member_index++) {
			const auto &pdf_member = pdf[member_index];

			const std::string path_trail = IO::leading_zeroes(pdf_member.set_member_number, 4);
			std::filesystem::path full_filename = base_output.stem();
			full_filename /= path_trail;
			full_filename.replace_extension(base_output.extension());
			std::filesystem::path output = base_output;
			output.replace_filename(full_filename);

			IO::create_directory_tree(output);
			std::ofstream file(output);

			output_run_info(file, pdf, comment);
			file << "#pdf_member = " << member_index << IO::endl;

			file << "x,Q2,y,E,xtbar,xbbar,xcbar,xsbar,xubar,xdbar,xg,xd,xu,xs,xc,xb,xt" << IO::endl;

			for (std::size_t i = 0; i < x_step_count; i++) {
				for (std::size_t j = 0; j < Q2_step_count; j++) {
					const double x = x_bins[i];
					const double Q2 = Q2_bins[j];

					pdf_member.evaluate(x, Q2);					

					file << x << ", " << Q2 << ", " << Q2 << ", " << Q2;
					for (const double xf : pdf_member.get_flavor_values()) {
						file << ", " << xf;
					}
					file << IO::endl;

					file.flush();

					calculated_values++;

					const int base_precision = 5;
					const int Q2_precision = base_precision - static_cast<int>(Math::number_of_digits(static_cast<int>(Q2)));

					std::cout << std::fixed << std::setprecision(base_precision);
					std::cout << "[Utility] " << IO::leading_zeroes(calculated_values, Math::number_of_digits(total_count));
					std::cout << " / " << total_count;
					std::cout << " [" << "pdf member " << IO::leading_zeroes(member_index + 1, Math::number_of_digits(member_count));
					std::cout << " / " << member_count << "]";
					std::cout << " (x = " << x;
					std::cout << std::setprecision(Q2_precision) << ", Q2 = " << Q2 << ")";
					std::cout << "\r" << std::flush;
				}
			}
			file.close();
		}
		std::cout << std::setprecision(static_cast<int>(original_precision)) << IO::endl;
	}

	private:
	void output_run_info(std::ofstream &file, const auto &pdf, const std::string comment) {		
		file << "#pdf = " << pdf.set_name << " [" << typeid(pdf).name() << "]" << IO::endl;
		file << "#parallelize = false" << IO::endl;
		file << "#number_of_threads = 1" << IO::endl;
		if (!comment.empty()) {
			file << "#comment = " << comment << IO::endl;
		}
	}
};

#endif