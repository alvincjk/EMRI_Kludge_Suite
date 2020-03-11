#!/usr/bin/env python

from __future__ import print_function

import sys, os, platform
from os.path import join as pjoin

from setuptools import setup
from setuptools import Extension
import distutils.sysconfig
from Cython.Distutils import build_ext

from Cython.Build import cythonize

import numpy

# maybe need rpath links to shared libraries on Linux
# if platform.system() == 'Linux':
#    linkArgs = ['-Wl,-R{}/lib'.format(...)]
# else:


def find_in_path(name, path):
    """Find a file in a search path"""

    # Adapted fom http://code.activestate.com/recipes/52224
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None


def locate_cuda():
    """Locate the CUDA environment on the system

    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDAHOME env variable. If not found,
    everything is based on finding 'nvcc' in the PATH.
    """

    # First check if the CUDAHOME env variable is in use
    if "CUDAHOME" in os.environ:
        home = os.environ["CUDAHOME"]
        nvcc = pjoin(home, "bin", "nvcc")
    else:
        # Otherwise, search the PATH for NVCC
        nvcc = find_in_path("nvcc", os.environ["PATH"])
        if nvcc is None:
            raise EnvironmentError(
                "The nvcc binary could not be "
                "located in your $PATH. Either add it to your path, "
                "or set $CUDAHOME"
            )
        home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {
        "home": home,
        "nvcc": nvcc,
        "include": pjoin(home, "include"),
        "lib64": pjoin(home, "lib64"),
    }
    for k, v in iter(cudaconfig.items()):
        if not os.path.exists(v):
            raise EnvironmentError(
                "The CUDA %s path could not be " "located in %s" % (k, v)
            )

    return cudaconfig


def customize_compiler_for_nvcc(self):
    """Inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.

    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on.
    """

    # Tell the compiler it can processes .cu
    self.src_extensions.append(".cu")

    # Save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # Now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == ".cu":
            # use the cuda for .cu files
            self.set_executable("compiler_so", CUDA["nvcc"])
            # use only a subset of the extra_postargs, which are 1-1
            # translated from the extra_compile_args in the Extension class
            postargs = extra_postargs["nvcc"]
        else:
            postargs = extra_postargs["gcc"]

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # Reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # Inject our redefined _compile method into the class
    self._compile = _compile


# Run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)


try:
    CUDA = locate_cuda()
    run_cuda = True

except OSError:
    run_cuda = False


# Obtain the numpy include directory. This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

lib_gsl_dir = "/opt/local/lib"
include_gsl_dir = "/opt/local/include"

if run_cuda:
    gpu_ext = Extension(
        "gpuAAK",
        sources=[
            "src/suite/gpuAAK.cc",
            "src/cuda/kernel.cu",
            "src/cuda/interpolate.cu",
            "src/cuda/AAK_manager.cu",
            "pygpuAAK/GPUAAK.pyx",
        ],
        library_dirs=[lib_gsl_dir, CUDA["lib64"], "/usr/local/lib", "./lib"],
        libraries=[
            "cudart",
            "cublas",
            "cusparse",
            "cufft",
            "gsl",
            "gslcblas",
            "KS",
            "IEKG",
            "LB",
            "NR",
            "RRGW",
            "GKG",
            "Circ",
        ],
        language="c++",
        runtime_library_dirs=[CUDA["lib64"]],
        # This syntax is specific to this build system
        # we're only going to use certain compiler args with nvcc
        # and not with gcc the implementation of this trick is in
        # customize_compiler()
        extra_compile_args={
            "gcc": ["-std=c99", "-Wno-unused-function"],  # , "-g"]
            "nvcc": [
                "-arch=sm_70",
                "--ptxas-options=-v",
                "-c",
                "--compiler-options",
                "'-fPIC'",
            ],  # ,"-G", "-g"] # for debugging
        },
        include_dirs=[
            numpy_include,
            include_gsl_dir,
            CUDA["include"],
            "./include",
            "./src",
            "./src/suite",
            "./src/cuda",
        ],
    )


linkArgs = []

cpu_ext = Extension(
    "AAKwrapper.AAKwrapper",
    ["AAKwrapper/AAKwrapper.pyx"],
    language="c++",
    include_dirs=["./include", numpy.get_include()],
    libraries=["gslcblas", "gsl", "KS", "IEKG", "LB", "NR", "RRGW", "GKG", "Circ",],
    library_dirs=["/usr/local/lib", "./lib"],
    extra_compile_args={"gcc": ["-Wno-unused-function"]},
    extra_link_args=linkArgs,
)

if run_cuda:
    extensions = [gpu_ext, cpu_ext]
    packages = ["AAKwrapper", "pygpuAAK"]
    package_dir = {"AAKwrapper": "AAKwrapper", "pygpuAAK": "pygpuAAK"}
    py_modules = ["pygpuAAK.pygpuAAK", "pygpuAAK.convert"]

else:
    extensions = [cpu_ext]
    packages = ["AAKwrapper"]
    package_dir = {"AAKwrapper": "AAKwrapper"}
    py_modules = []

setup(
    name="AAKwrapper",
    version="0.1",
    description="A Python wrapper for AAK",
    # author = 'Michele Vallisneri',
    # author_email = 'vallis@vallis.org',
    packages=packages,
    package_dir=package_dir,
    py_modules=py_modules,
    ext_modules=cythonize(extensions),
    cmdclass={"build_ext": custom_build_ext},
    zip_safe=False,
)
