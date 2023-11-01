// #ifndef INTERPOLATING_FUNCTION_H
// #define INTERPOLATING_FUNCTION_H

// #include <filesystem>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <algorithm>
// #include <iterator>
// #include <tuple>
// #include <concepts>
// #include <ranges>
// #include <vector>
// #include <span>

// #include "Utility/Math.cpp"
// #include "Utility/Utility.cpp"

// struct InterpolatingFunction3D {
// 	using Vertex = std::tuple<std::size_t, std::size_t, std::size_t>;
// 	using Point = std::tuple<double, double, double>;
	
// 	const std::filesystem::path grid_path;

// 	InterpolatingFunction3D(const std::filesystem::path grid_path, const std::string separator = "-----") : grid_path(grid_path) {
// 		std::ifstream grid_file(grid_path);
// 		std::size_t current_index = 0;

// 		const std::size_t line_count = static_cast<std::size_t>(std::count(std::istreambuf_iterator<char>(grid_file), std::istreambuf_iterator<char>(), '\n'));

// 		std::vector<double> grid_values(line_count);

// 		bool in_grid_info = true;

// 		for (std::string line = ""; std::getline(grid_file, line);) {
// 			current_index++;

// 			if (current_index == 1) { continue; }
			
// 			if (line == separator) {
// 				in_grid_info = false;
// 				continue;
// 			}

// 			if (in_grid_info) {
// 				std::istringstream string_stream(line);
// 				double value;

// 				while (string_stream >> value) {
// 					switch (current_index) {
// 					case 2:
// 						x_points.push_back(value);
// 						break;
// 					case 3:
// 						y_points.push_back(value);
// 						break;
// 					case 4:
// 						z_points.push_back(value);
// 						break;
// 					default:
// 						break;
// 					}
// 				}
// 			} else {
// 				grid_values.push_back(std::stod(line));
// 			}
// 		}
// 	}

// 	double at(std::same_as<std::vector<double>::size_type> auto... indices) const {
// 		const std::array<std::vector<double>::size_type, sizeof...(indices)>index_array{indices...};
// 		return at(index_array);
// 	}

// 	double at(const std::vector<std::vector<double>::size_type> &index_array) const {
// 		std::vector<double>::size_type index = std::inner_product(std::begin(index_array), std::end(index_array), std::begin(dimension_conversions), std::vector<double>::size_type{0});
// 		return grid[index];
// 	}

// 	double operator()(std::same_as<double> auto... input) {
// 		const Point evaluation_point{input...};

// 		Vertex closest_point;

// 		for (const auto &[value, points] : std::views::zip(evaluation_point, grid_points)) {
// 			const auto closest_index = Collections::closest(points, value);
// 			closest_point.push_back(static_cast<Vertex::value_type>(closest_index));
// 		}

// 		return interpolate(closest_point, evaluation_point);
// 	}

// 	private:
// 	std::vector<double> x_points;
// 	std::vector<double> y_points;
// 	std::vector<double> z_points;
// 	Kokkos::mdspan<double, std::dynamic_extent, std::dynamic_extent, std::dynamic_extent> grid;

// 	// double interpolate(Vertex vertex, Point evaluation_point) {
// 	// 	const Vertex v000 = vertex;

// 	// 	const Vertex v001 = push_along_axes({0}, vertex);
// 	// 	const Vertex v010 = push_along_axes({1}, vertex);
// 	// 	const Vertex v100 = push_along_axes({2}, vertex);

// 	// 	const Vertex v011 = push_along_axes({0, 1}, vertex);
// 	// 	const Vertex v101 = push_along_axes({0, 2}, vertex);

// 	// 	const Vertex v110 = push_along_axes({1, 2}, vertex);

// 	// 	const Vertex v111 = push_along_axes({0, 1, 2}, vertex);

// 	// 	const Point g000 = grid_point(v000);
// 	// 	const Point g001 = grid_point(v001);
// 	// 	const Point g010 = grid_point(v010);
// 	// 	const Point g100 = grid_point(v100);
// 	// 	const Point g011 = grid_point(v011);
// 	// 	const Point g101 = grid_point(v101);
// 	// 	const Point g110 = grid_point(v110);
// 	// 	const Point g111 = grid_point(v111);

// 	// 	const double xd = (evaluation_point[0] - g000[0]) / (g001[0] - g000[0]);
// 	// 	const double yd = (evaluation_point[1] - g000[1]) / (g010[1] - g000[1]);
// 	// 	const double zd = (evaluation_point[2] - g000[2]) / (g100[2] - g000[2]);

// 	// 	const double c00 = linear_interpolation(at(v000), at(v100), xd);
// 	// 	const double c01 = linear_interpolation(at(v001), at(v101), xd);
// 	// 	const double c10 = linear_interpolation(at(v010), at(v110), xd);
// 	// 	const double c11 = linear_interpolation(at(v011), at(v111), xd);

// 	// 	const double c0 = linear_interpolation(c00, c10, yd);
// 	// 	const double c1 = linear_interpolation(c01, c11, yd);

// 	// 	const double c = linear_interpolation(c0, c1, zd);

// 	// 	return c;

// 	// 	// const double base_value = at(vertex);
// 	// 	// const Point lower_grid_point = grid_point(vertex);

// 	// 	// std::vector<double> linear_interpolations;

// 	// 	// for (std::size_t axis = 0; axis < vertex.size(); axis++) {
// 	// 	// 	const Vertex pushed_vertex = push_along_axis(axis, vertex);
// 	// 	// 	const double pushed_value = at(pushed_vertex);
// 	// 	// 	const Point upper_grid_point = grid_point(pushed_vertex);

// 	// 	// 	const double length = upper_grid_point[axis] - lower_grid_point[axis];
// 	// 	// 	const double relative_position = (evaluation_point[axis] - lower_grid_point[axis]) / length;
// 	// 	// 	const double interpolation = linear_interpolation(base_value, pushed_value, relative_position);

// 	// 	// 	linear_interpolations.push_back(interpolation);
// 	// 	// }
// 	// }

// 	// double linear_interpolation(const double x1, const double x2, const double t) {
// 	// 	return std::lerp(x1, x2, t);
// 	// }

// 	// std::vector<std::vector<double>::size_type> push_along_axis(const std::size_t axis, const Vertex &vertex) {
// 	// 	std::vector<std::vector<double>::size_type> output = vertex;
// 	// 	output[axis]++;
// 	// 	return output;
// 	// }
// 	// std::vector<std::vector<double>::size_type> push_along_axes(const std::vector<std::size_t> &axes, const Vertex &vertex) {
// 	// 	std::vector<std::vector<double>::size_type> output = vertex;
		
// 	// 	for (const std::size_t axis : axes) {
// 	// 		output[axis]++;
// 	// 	}

// 	// 	return output;
// 	// }

// 	// Point grid_point(const Vertex &vertex) {
// 	// 	Point output;

// 	// 	for (const auto &[hypercube_vertex, grid_point_list] : std::views::zip(vertex, grid_points)) {
// 	// 		output.push_back(grid_point_list[hypercube_vertex]);
// 	// 	}

// 	// 	return output;
// 	// }
// };


// #endif