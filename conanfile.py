from conans import ConanFile, CMake, tools


class ThermistorConan(ConanFile):
    name = "thermistor"
    version = "0.1"
    license = "MIT"
    author = "Matthew Knight <mgk1795@gmail.com>"
    description = "<Description of Thermistor here>"
    topics = ("thermistor", "embedded", "")
    exports_sources = "include/*"
    requires = "gcem/1.12.0@matt1795/testing"

    def package(self):
        self.copy("*", ".")

    def package_info(self):
        self.cpp_info.cxxflags.append("-std=c++17")
