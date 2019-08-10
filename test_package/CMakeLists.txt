cmake_minimum_required(VERSION 2.8.12)
project(PackageTest CXX)

find_package(thermistor REQUIRED)
add_executable(main main.cpp)
target_link_libraries(main thermistor::thermistor)
