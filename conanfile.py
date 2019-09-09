from conans import ConanFile, CMake, tools


class ThermistorConan(ConanFile):
    name = "thermistor"
    version = "1.0"
    license = "MIT"
    author = "Matthew Knight <mgk1795@gmail.com>"
    description = "Header-only, compile-time thermistor lookup table library"
    homepage = "https://github.com/matt1795/thermistor"
    url = "https://github.com/matt1795/thermistor"
    topics = ("thermistor", "embedded", "compile-time")
    exports_sources = "include/*"
    requires = "gcem/1.12.0@matt1795/stable"

    def package(self):
        self.copy("*", ".")

    def package_info(self):
        self.cpp_info.cxxflags.append("-std=c++17")
