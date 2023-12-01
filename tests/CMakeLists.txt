# Source files
set(SOURCE_FILES main.cpp)

set(TEST_PROJECT_NAME ${PROJECT_NAME}Tests)

enable_testing()

# Executable
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)
add_executable(${TEST_PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${TEST_PROJECT_NAME} PUBLIC libChargedCurrent)

# Thread sanitizer
SET(USE_THREAD_SANITIZER FALSE CACHE BOOL "Enable thread sanitizer")
if (${USE_THREAD_SANITIZER})
	add_compile_options(-fPIE -pie -g)
	add_compile_options(-fsanitize=thread -DLIBOMP_TSAN_SUPPORT=TRUE)
	add_link_options(-fsanitize=thread)
	message("Thread sanitizer enabled")
endif()
	
# Address sanitizer
SET(USE_ADDRESS_SANITIZER FALSE CACHE BOOL "Enable address sanitizer")
if (${USE_ADDRESS_SANITIZER})
	add_compile_options(-fPIE -pie -g)
	add_compile_options(-fsanitize=address)
	add_link_options(-fsanitize=address)
	message("Address sanitizer enabled")
endif()

# UB sanitizer
SET(USE_UB_SANITIZER FALSE CACHE BOOL "Enable UB sanitizer")
if (${USE_UB_SANITIZER})
	add_compile_options(-fPIE -pie -g)
	add_compile_options(-fsanitize=undefined)
	add_link_options(-fsanitize=undefined)
	message("Undefined behavior sanitizer enabled")
endif()

SET(USE_CLANG_TIDY FALSE CACHE BOOL "Enable clang-tidy")
if (${USE_CLANG_TIDY})
	set(CMAKE_CXX_CLANG_TIDY 
		clang-tidy;
		-format-style='file';
		-header-filter=${CMAKE_CURRENT_SOURCE_DIR};
	)
	message("clang-tidy enabled")
endif()

# Flags
target_compile_options(
	${TEST_PROJECT_NAME} PRIVATE 
	-Wall -Wconversion -Wextra -Wpedantic
	$<$<CONFIG:Release>:-O3 -DNDEBUG -ffast-math -march=native>
)

message("Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}")

# GoogleTest
find_package(GTest REQUIRED)
if (GTest_FOUND)
	target_include_directories(${TEST_PROJECT_NAME} PRIVATE ${GTEST_INCLUDE_DIRS})
	target_link_libraries(${TEST_PROJECT_NAME} PRIVATE GTest::gtest_main)
endif()
include(GoogleTest)
gtest_discover_tests(${TEST_PROJECT_NAME})

# Unset cached values
unset(USE_THREAD_SANITIZER CACHE)
unset(USE_ADDRESS_SANITIZER CACHE)
unset(USE_UB_SANITIZER CACHE)
unset(USE_CLANG_TIDY CACHE)