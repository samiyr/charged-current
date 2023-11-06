#ifndef INTERPOLATING_FUNCTION_H
#define INTERPOLATING_FUNCTION_H

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "Utility/Math.cpp"
#include "Utility/Utility.cpp"

struct InterpolatingFunction {	
	const std::filesystem::path grid_path;

	InterpolatingFunction(const std::filesystem::path grid_path, const gsl_interp_type *interpolation_type = gsl_interp_linear, const std::string separator = "-----") : grid_path(grid_path) {
		std::cout << grid_path << IO::endl;
		std::ifstream grid_file(grid_path);

		std::size_t current_index = 0;

		const std::size_t line_count = static_cast<std::size_t>(std::count(std::istreambuf_iterator<char>(grid_file), std::istreambuf_iterator<char>(), '\n'));

		std::vector<double> grid_values(line_count);
		std::vector<double> x_points;

		bool in_grid_info = true;

		std::string line;

		while (std::getline(grid_file, line)) {
			std::cout << line << IO::endl;
			std::cout << current_index << IO::endl;
			current_index++;

			if (current_index == 1) { continue; }
			
			if (line == separator) {
				in_grid_info = false;
				continue;
			}

			if (in_grid_info) {
				std::istringstream string_stream(line);
				double value;

				while (string_stream >> value) {
					x_points.push_back(value);
				}
			} else {
				grid_values.push_back(std::stod(line));
			}
		}

		initialize_interpolation(x_points, grid_values, interpolation_type);
	}

	InterpolatingFunction(const std::vector<double> x_points, const std::vector<double> grid_values, const gsl_interp_type *interpolation_type = gsl_interp_linear) {
		initialize_interpolation(x_points, grid_values, interpolation_type);
	}

	~InterpolatingFunction() {
		gsl_spline_free(spline);
		gsl_interp_accel_free(x_accelerator);
	}

	double operator()(const double x) const {
		return gsl_spline_eval(spline, x, x_accelerator);
	}

	private:
	gsl_spline *spline;
	gsl_interp_accel *x_accelerator;

	void initialize_interpolation(const std::vector<double> x_points, const std::vector<double> grid_values, const gsl_interp_type *interpolation_type) {
		spline = gsl_spline_alloc(interpolation_type, x_points.size());
		x_accelerator = gsl_interp_accel_alloc();

		gsl_spline_init(spline, x_points.data(), grid_values.data(), x_points.size());
	}
};


#endif