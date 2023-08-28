add_subdirectory(Coefficients)

target_sources(${LIBRARY_NAME} 
	PRIVATE 
		SIDIS.cpp
		SIDISComputation.cpp
		SIDISFunctions.cpp
)