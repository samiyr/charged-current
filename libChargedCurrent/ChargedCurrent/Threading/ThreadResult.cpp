#ifndef THREAD_RESULT_H
#define THREAD_RESULT_H

#include <fstream>
#include <iostream>
#include <filesystem>

#include "Utility/Utility.cpp"

#include "Threading/BS_thread_pool.hpp"
#include "Threading/shared_multi_future.cpp"

struct ThreadResult {
	shared_multi_future<std::string> future;
	std::filesystem::path output;
	std::string header;

	std::vector<std::string> get() {
		return future.get();
	}

	void write() {
		IO::create_directory_tree(output);
		std::ofstream file(output);

		file << header;

		for (const std::string &string : get()) {
			file << string;
		}
	}
};

struct ThreadResultCollection {
	std::vector<ThreadResult> results;

	void add(ThreadResult result) {
		results.push_back(result);
	}

	void add(std::vector<ThreadResult> &&vector) {
		results.insert(results.end(), vector.begin(), vector.end());
	}

	void write() {
		for (ThreadResult &result : results) {
			result.write();
		}
	}
};

#endif