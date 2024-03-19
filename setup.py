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
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext as build_ext_orig

#https://stackoverflow.com/questions/47360113/compile-c-library-on-pip-install
        
class CMakeExtension(Extension):
    def __init__(self, name, sources=[]):
        Extension.__init__(self, name, sources=[])
        print(name, sources)
        self.sourcedir = os.path.join(os.path.abspath(''), "src", "sigmond", "source" )

class CMakeBuild(build_ext_orig):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        new_lib = os.path.basename(self.get_ext_fullpath(ext.name))
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        # extdir = os.path.join(os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name))),"sigmond")

        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable,
                       '-DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF',
                       '-Wno-dev',
                      '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=' + extdir,
                      '-DVERSION_INFO=0.0.1',
                    #   '-DCMAKE_INSTALL_RPATH=$ORIGIN',
                    #   '-DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON',
                    #   '-DCMAKE_INSTALL_RPATH_USE_LINK_PATH:BOOL=ON',
                    #   '-DCMAKE_INSTALL_PREFIX:PATH=' + extdir,
                    #   '-DLIBRARY_OUTPUT_NAME=sigmond.so',#+out_file,
                     ]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            cmake_args += ["-DCMAKE_CXX_FLAGS='-DDEFAULTENSFILE=\'\"nonesense\"\' -Wall /std:c++17 /Ox /EHa'"] 
            cmake_args += ["-DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE"] # needs lapack and hdf5 flags too
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            # build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            cmake_args += ["-DCMAKE_CXX_FLAGS='-DDEFAULTENSFILE=\'\"\"\' -Wall -std=c++17 -O3 -llapack -lhdf5'"]
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

setup(
    name='sigmond',
    version="0.0.0.dev1",
    author="Sarah Skinner",
    author_email="sarakski@andrew.cmu.edu",
    description='A python interface to the for the Sigmond analysis software.',
    long_description='',
    packages=['sigmond'],
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux"
    ],
    ext_modules=[CMakeExtension('sigmond',['src/sigmond/source/pysigmond/pysigmond.cc'])],
    python_requires='>=3.6',
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)