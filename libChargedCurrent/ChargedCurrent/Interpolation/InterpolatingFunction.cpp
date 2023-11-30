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

	InterpolatingFunction(
		const std::filesystem::path grid_path, 
		const gsl_interp_type *interpolation_type = gsl_interp_akima, 
		const std::string separator = "-----") : grid_path(grid_path), interpolation_type(interpolation_type), separator(separator) { }

	// Initializes the interpolating function directly. An instance constructed using this method cannot be copied across threads.
	InterpolatingFunction(const std::vector<double> grid_points, const std::vector<double> grid_values, const gsl_interp_type *interpolation_type = gsl_interp_akima) {
		initialize_interpolation(grid_points, grid_values, interpolation_type);
		initialized = true;
	}

	~InterpolatingFunction() {
		if (initialized) {
			gsl_spline_free(spline);
			gsl_interp_accel_free(accelerator);
		}
	}

	double operator()(const double x) const {
		if (!initialized) { initialize(); }

		const double min = spline->interp->xmin;
		const double max = spline->interp->xmax;
		if (x < min) { return operator()(min); }
		if (x > max) { return operator()(max); }
		return gsl_spline_eval(spline, x, accelerator);
	}

	private:
	mutable bool initialized = false;
	mutable gsl_spline *spline;
	mutable gsl_interp_accel *accelerator;
	const gsl_interp_type *interpolation_type;
	mutable std::string separator;

	void initialize() const {
		std::ifstream grid_file(grid_path);
		if (!grid_file) {
			throw std::runtime_error("Grid file '" + grid_path.generic_string() + "' could not be opened.");
		}

		std::size_t current_index = 0;

		const std::size_t line_count = static_cast<std::size_t>(std::count(std::istreambuf_iterator<char>(grid_file), std::istreambuf_iterator<char>(), '\n'));
		grid_file.seekg(0);

		std::string line;

		std::vector<double> grid_values;
		grid_values.reserve(line_count);

		std::vector<double> grid_points;

		bool in_grid_info = true;

		while (std::getline(grid_file, line)) {
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
					grid_points.push_back(value);
				}
			} else {
				grid_values.push_back(std::stod(line));
			}
		}

		initialize_interpolation(grid_points, grid_values, interpolation_type);

		initialized = true;
	}

	void initialize_interpolation(const std::vector<double> grid_points, const std::vector<double> grid_values, const gsl_interp_type *interpolation_type) const {
		spline = gsl_spline_alloc(interpolation_type, grid_points.size());
		accelerator = gsl_interp_accel_alloc();

		gsl_spline_init(spline, grid_points.data(), grid_values.data(), grid_points.size());
	}
};


#endif