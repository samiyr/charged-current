# Basic definitions
set(PROJECT_NAME "ChargedCurrentDIS")
cmake_minimum_required(VERSION 3.26.3..3.27.2)
project(${PROJECT_NAME} VERSION 1.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# Load .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Load subdirectories
add_subdirectory(libChargedCurrent)
add_subdirectory(ChargedCurrent)
add_subdirectory(tests)