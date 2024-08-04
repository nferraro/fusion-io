# Copyright 2013-2024 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class FusionIo(CMakePackage):
    """Fusion-IO"""

    homepage = ""
    git = "https://github.com/nferraro/fusion-io"

    maintainers("changliu")

    license("BSD-3-Clause")

    # We will use the scorec/core master branch as the 'nightly' version
    # of pumi in spack.  The master branch is more stable than the
    # scorec/core develop branch and we prefer not to expose spack users
    # to the added instability.
    version("master", submodules=True, branch="cmake_changliu")

    variant("python", default=True, description="Enable Python support")
    variant("trace", default=True, description="Build trace program")

    depends_on("mpi")
    depends_on("hdf5")
    depends_on("cmake@3:", type="build")

    extends("python", when="+python")

    def cmake_args(self):
        spec = self.spec

        args = [
            "-DCMAKE_C_COMPILER=%s" % spec["mpi"].mpicc,
            "-DCMAKE_CXX_COMPILER=%s" % spec["mpi"].mpicxx,
            "-DCMAKE_Fortran_COMPILER=%s" % spec["mpi"].mpifc,
            self.define_from_variant("ENABLE_TRACE", "trace"),
        ]

        if self.spec.variants["python"].value:
            args.extend(["-DENABLE_PYTHON:BOOL=ON"])
            args.append(self.define("PYTHON_MODULE_INSTALL_PATH", python_platlib))

        return args

