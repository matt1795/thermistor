cmake_minimum_required(VERSION 2.8.12)
project(ThermistorTest CXX)

include(CTest)
include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)

enable_testing()
conan_basic_setup()

file(GLOB SOURCE_FILES "src/*.cpp")

if(WITH_GMOCK)
    add_definitions("-DWITH_GMOCK")
endif()

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${CONAN_LIBS})
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
