add_subdirectory(Coefficients)

target_sources(${LIBRARY_NAME} 
	PRIVATE 
		DIS.cpp
		DISComputation.cpp
		DISFunctions.cpp
)