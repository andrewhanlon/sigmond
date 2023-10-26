#stolen from https://stackoverflow.com/questions/53397121/python-pip-packaging-how-to-move-built-files-to-install-directory
import os
import re
import sys
import sysconfig
import site
import platform
import subprocess
import pathlib
import shutil

from distutils.version import LooseVersion
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as build_ext_orig

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.join(os.path.abspath(sourcedir), "src", "sigmond", "source" )

class CMakeBuild(build_ext_orig):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            raise RuntimeError("Sorry, pyScannerBit doesn't work on Windows platforms. Please use Linux or OSX.")

        print(self.extensions)
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # extdir = os.path.join(extdir,"sigmond")
        print(extdir)
        # shutil.copyfile("/home/sarahski/latticeQCD/formerges/sigmond.pip/sigmond/libsigmond.so", extdir )

        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                       '-DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF',
                       '-Wno-dev',
                    #   '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=' + extdir,
                    #   '-DSCANNERBIT_STANDALONE=True',
                    #   '-DCMAKE_INSTALL_RPATH=$ORIGIN',
                    #   '-DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON',
                    #   '-DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON',
                    #   '-DCMAKE_INSTALL_PREFIX:PATH=' + extdir,
                     ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            # build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            # build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # # First cmake
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        # Main build
        subprocess.check_call(['cmake', '--build', "."] , cwd=self.build_temp)
        # subprocess.check_call(['make'] , cwd=self.build_temp)
        # # Install
        # subprocess.check_call(['cmake', '--build', '.', '--target', 'install'], cwd=self.build_temp)

setup(
    name='sigmond',
    version="0.0.0.dev1",
    author="Sarah Skinner",
    # Add yourself if you contribute to this package
    author_email="sarakski@andrew.cmu.edu",
    description='A python interface to the for the Sigmond analysis software.',
    long_description='',
    packages=['sigmond'],
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux"
    ],
    ext_modules=[CMakeExtension('sigmond/pysigmond/pysigmond')],
    python_requires='>=3.6',
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)


# [tool.mypy]
# files = "setup.py"
# python_version = "3.7"
# strict = true
# show_error_codes = true
# enable_error_code = ["ignore-without-code", "redundant-expr", "truthy-bool"]
# warn_unreachable = true

# [[tool.mypy.overrides]]
# module = ["ninja"]
# ignore_missing_imports = true


# [tool.pytest.ini_options]
# minversion = "6.0"
# addopts = ["-ra", "--showlocals", "--strict-markers", "--strict-config"]
# xfail_strict = true
# filterwarnings = [
#     "error",
#     "ignore:(ast.Str|Attribute s|ast.NameConstant|ast.Num) is deprecated:DeprecationWarning:_pytest",
# ]
# testpaths = ["tests"]

# [tool.cibuildwheel]
# test-command = "pytest {project}/tests"
# test-extras = ["test"]
# test-skip = ["*universal2:arm64"]
# # Setuptools bug causes collision between pypy and cpython artifacts
# before-build = "rm -rf {project}/build"

# [tool.ruff]
# extend-select = [
#   "B",    # flake8-bugbear
#   "B904",
#   "I",    # isort
#   "PGH",  # pygrep-hooks
#   "RUF",  # Ruff-specific
#   "UP",   # pyupgrade
# ]
# extend-ignore = [
#   "E501",   # Line too long
# ]
# target-version = "py37"