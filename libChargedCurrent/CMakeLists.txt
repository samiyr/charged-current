set(LIBRARY_NAME libChargedCurrent)
add_library(${LIBRARY_NAME} STATIC "")
target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ChargedCurrent)

add_subdirectory(ChargedCurrent/Analysis)
add_subdirectory(ChargedCurrent/Common)
add_subdirectory(ChargedCurrent/Decay)
add_subdirectory(ChargedCurrent/DIS)
add_subdirectory(ChargedCurrent/Integration)
add_subdirectory(ChargedCurrent/PDF)
add_subdirectory(ChargedCurrent/SIDIS)
add_subdirectory(ChargedCurrent/Utility)

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

# Link-time optimization
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)

if (ipo_supported)
    message(STATUS "IPO / LTO enabled")
    set_property(TARGET ${LIBRARY_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${ipo_error}>")
endif()

target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Flags
target_compile_options(
	${LIBRARY_NAME} PUBLIC 
	-Wall -Wconversion -Wextra -Wpedantic -g
	$<$<CONFIG:Release>:-O3 -DNDEBUG -ffast-math -march=native>
)

message("Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}")

# OpenMP
if (APPLE)
	target_link_directories(${LIBRARY_NAME} PUBLIC /opt/homebrew/lib /opt/homebrew/opt/llvm/lib)

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

endif()

SET(OMP_OVERRIDE FALSE CACHE BOOL "Override find_package(OpenMP)")
if(OMP_OVERRIDE)
	target_compile_options(${LIBRARY_NAME} PUBLIC -fopenmp)
else()
	find_package(OpenMP REQUIRED)
	if (OpenMP_CXX_FOUND)
		target_link_libraries(${LIBRARY_NAME} PUBLIC OpenMP::OpenMP_CXX)
	else()
		target_compile_options(${LIBRARY_NAME} PUBLIC -fopenmp)
	endif()
endif()

# GSL
find_package(GSL REQUIRED)
if (GSL_FOUND)
	target_link_libraries(${LIBRARY_NAME} PUBLIC GSL::gsl GSL::gslcblas)
endif()

# LHAPDF
find_package(LHAPDF REQUIRED)
if (LHAPDF_FOUND)
	target_link_libraries(${LIBRARY_NAME} PUBLIC LHAPDF)
endif()

# Pythia
find_package(Pythia8 REQUIRED)
if (PYTHIA8_FOUND)
	target_include_directories(${LIBRARY_NAME} PUBLIC ${PYTHIA8_INCLUDE_DIRS})
	target_link_libraries(${LIBRARY_NAME} PUBLIC ${PYTHIA8_LIBRARIES})
endif()

# Cuba
find_package(Cuba REQUIRED)
if (Cuba_FOUND)
	target_include_directories(${LIBRARY_NAME} PUBLIC ${Cuba_INCLUDE_DIRS})
	target_link_libraries(${LIBRARY_NAME} PUBLIC ${Cuba_LIBRARIES})
endif()
	
# Unset cached values
unset(USE_THREAD_SANITIZER CACHE)
unset(USE_ADDRESS_SANITIZER CACHE)
unset(USE_UB_SANITIZER CACHE)
unset(USE_CLANG_TIDY CACHE)
unset(OMP_OVERRIDE CACHE)