# Source files
set(SOURCE_FILES main.cpp)

# Executable
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC libChargedCurrent)

# Thread sanitizer
SET(USE_THREAD_SANITIZER FALSE CACHE BOOL "Enable thread sanitizer")
if (${USE_THREAD_SANITIZER})
	# add_compile_options(-fPIE -pie -g)
	# add_compile_options(-fsanitize=thread -DLIBOMP_TSAN_SUPPORT=TRUE)
	# add_link_options(-fsanitize=thread)
	target_compile_options(${PROJECT_NAME} PUBLIC -fPIE -g)
	target_compile_options(${PROJECT_NAME} PUBLIC -fsanitize=thread -DLIBOMP_TSAN_SUPPORT=TRUE)
	target_link_options(${PROJECT_NAME} PUBLIC -fsanitize=thread)
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
    set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${ipo_error}>")
endif()

# Flags
target_compile_options(
	${PROJECT_NAME} PRIVATE 
	-Wall -Wconversion -Wextra -Wpedantic -g
	$<$<CONFIG:Release>:-O3 -DNDEBUG -ffast-math -march=native>
)

message("Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}")

# Unset cached values
unset(USE_THREAD_SANITIZER CACHE)
unset(USE_ADDRESS_SANITIZER CACHE)
unset(USE_UB_SANITIZER CACHE)
unset(USE_CLANG_TIDY CACHE)