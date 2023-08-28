// #ifndef FILE_PATH_H
// #define FILE_PATH_H

// #include <filesystem>
// #include <optional>

// struct FilePath {
// 	FilePath(
// 		const std::string filename, 
// 		const std::string extension = "", 
// 		const std::optional<std::filesystem::path> parent_directory = std::nullopt, 
// 		const bool create_parent = true
// 	) {
// 		if (parent_directory) {
// 			_path = *parent_directory / filename;
// 		} else {
// 			_path = filename;
// 		}
		
// 		_path.replace_extension(extension);

// 		// if (create_parent) {
// 		// 	create_directory_tree();
// 		// }
// 	}
	
// 	void create_directory_tree() const {
// 		const auto directory = path().parent_path();
// 		std::filesystem::create_directories(directory);
// 	}

// 	std::filesystem::path path() const { return _path; }

// 	private:
// 	std::filesystem::path _path;
// };

// #endif