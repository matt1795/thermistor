import os

from conans import ConanFile, CMake, tools

class ThermistorTestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"
    requires = "gtest/1.8.0@bincrafters/stable", "thermistor/1.0@matt1795/stable"

    def configure(self):
        self.options["gtest"].build_gmock = False

    def build(self):
        cmake = CMake(self)
        cmake.definitions['WITH_GMOCK'] = self.options['gtest'].build_gmock
        cmake.configure()
        cmake.build()

    def test(self):
        if not tools.cross_building(self.settings):
            os.chdir("bin")
            self.run(".%sThermistorTest" % os.sep)

