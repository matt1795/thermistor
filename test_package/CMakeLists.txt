cmake_minimum_required(VERSION 2.8.12)
project(ThermistorTest CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CTest)
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)

enable_testing()
conan_basic_setup()

file(GLOB SOURCE_FILES "src/*.cpp")

if(WITH_GMOCK)
    add_definitions("-DWITH_GMOCK")
endif()

add_executable(${PROJECT_NAME}
    src/ntc.cpp
    src/circuit.cpp)

target_link_libraries(${PROJECT_NAME} ${CONAN_LIBS})
target_include_directories(${PROJECT_NAME} PRIVATE include)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
