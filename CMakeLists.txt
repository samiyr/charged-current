# Basic definitions
set(PROJECT_NAME "ChargedCurrentDIS")
cmake_minimum_required(VERSION 3.26.3..3.26.4)
project(${PROJECT_NAME} VERSION 1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Load .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Source files
set(SOURCE_FILES src/ChargedCurrentDIS.cpp)

# Thread sanitizer
SET(USE_THREAD_SANITIZER FALSE CACHE BOOL "Enable thread sanitizer")
if (${USE_THREAD_SANITIZER})
	add_compile_options(-fPIE -pie -g)
	add_compile_options(-fsanitize=thread -DLIBOMP_TSAN_SUPPORT=TRUE)
	add_link_options(-fsanitize=thread)
endif()

# Address sanitizer
SET(USE_ADDRESS_SANITIZER FALSE CACHE BOOL "Enable address sanitizer")
if (${USE_ADDRESS_SANITIZER})
	add_compile_options(-fPIE -pie -g)
	add_compile_options(-fsanitize=address)
	add_link_options(-fsanitize=address)
endif()

# Executable
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# Flags
target_compile_options(${PROJECT_NAME} PRIVATE -ffast-math -Wall -Wconversion -Wextra -Wpedantic $<$<CONFIG:Release>:-O3> -g)
# target_compile_options(${PROJECT_NAME} PRIVATE -ffast-math -Wall -Wconversion -Wextra -Wpedantic $<$<CONFIG:Release>:-O3> $<$<CONFIG:Debug>:-g>)

# OpenMP
if (APPLE)
	set(CMAKE_C_COMPILER clang)
	set(CMAKE_CXX_COMPILER clang++)

	if (CMAKE_C_COMPILER_ID MATCHES "Clang\$")
		set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
		set(OpenMP_C_LIB_NAMES "omp")
		set(OpenMP_omp_LIBRARY omp)
	endif()

	if (CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
		set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
		set(OpenMP_CXX_LIB_NAMES "omp")
		set(OpenMP_omp_LIBRARY omp)
	endif()

	target_link_directories(${PROJECT_NAME} PRIVATE /opt/homebrew/lib /opt/homebrew/opt/llvm/lib)
endif()

find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

# GSL
find_package(GSL REQUIRED)
if (GSL_FOUND)
	target_link_libraries(${PROJECT_NAME} PRIVATE GSL::gsl GSL::gslcblas)
endif()

# LHAPDF
find_package(LHAPDF REQUIRED)
if (LHAPDF_FOUND)
	target_link_libraries(${PROJECT_NAME} PRIVATE LHAPDF)
endif()

# Pythia
find_package(Pythia8 REQUIRED)
if (PYTHIA8_FOUND)
	target_include_directories(${PROJECT_NAME} PRIVATE ${PYTHIA8_INCLUDE_DIRS})
	target_link_libraries(${PROJECT_NAME} PRIVATE ${PYTHIA8_LIBRARIES})
endif()

# Unset cached values
unset(USE_THREAD_SANITIZER CACHE)
unset(USE_ADDRESS_SANITIZER CACHE)