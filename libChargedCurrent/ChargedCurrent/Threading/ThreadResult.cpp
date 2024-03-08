#ifndef THREAD_RESULT_H
#define THREAD_RESULT_H

#include <fstream>
#include <iostream>

#include "Threading/BS_thread_pool.hpp"

struct ThreadResult {
	BS::multi_future<std::string> &future;
	std::ofstream &file;

	std::vector<std::string> get() {
		return future.get();
	}

	void write() {
		for (const std::string &string : get()) {
			file << string;
		}

		file.close();
	}
};

struct ThreadResultCollection {
	std::vector<std::reference_wrapper<ThreadResult>> results;

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