# Membrane Dynamics in 3D using Discrete Differential Geometry (Mem3DG)

[![Testing and release](https://github.com/RangamaniLabUCSD/Mem3DG/actions/workflows/ci.yaml/badge.svg?branch=main)](https://github.com/RangamaniLabUCSD/Mem3DG/actions/workflows/ci.yaml)
[![Zenodo](https://zenodo.org/badge/244037679.svg)](https://zenodo.org/doi/10.5281/zenodo.10359392)
[![PyPI](https://img.shields.io/pypi/v/pymem3dg)](https://pypi.org/project/pymem3dg/)
[![Conda-forge](https://anaconda.org/conda-forge/pymem3dg/badges/version.svg)](https://anaconda.org/conda-forge/pymem3dg)

Mem3DG is a flexible software package to model the membrane and its dynamics using unstructured meshes.
This work is currently under heavy development, please star this repository to follow along!

## Acknowledging your use of Mem3DG

Mem3DG is developed by [Cuncheng Zhu](https://github.com/cuzhucuncheng), [Christopher T. Lee](https://ctlee.github.io/), with contributions from others.
Development of Mem3DG is funded in part by AFOSR MURI FA9550-18-1-0051, and a Hartwell Foundation Postdoctoral Fellowship.

## Installation

```
git submodule update --init --recursive
mkdir build
cd build
cmake -DBUILD_PYDDG=ON -DWITH_NETCDF=ON -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release
```

Source released can also be obtained from [PyPi](https://pypi.org/project/pymem3dg/).


## Building pymem3dg

The python library can be built using pip:
```
pip install . -v
```
Options can be passed to scikit-build by exporting environment variables.
```
export SKBUILD_CONFIGURE_OPTIONS="-DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo"
```

## Temporary notes for setting up netcdf (especially on windows...)

1. Download `vcpkg` and follow the instructions to install
2. Install 32 or 64 bit version of `netcdf-c` and `netcdf-cxx4` libraries depending upon your configuration.

   `vcpkg install netcdf-c:x64-windows netcdf-cxx4:x64-windows eigen3:x64-windows`

   Remove the `:x64-windows` from the above string for the 32 bit libraries.

3. Configure the vcpkg CMake toolchain `vcpkg integrate install`
4. Copy and paste the `-DCMAKE_TOOLCHAIN_File="..."` string as an input into your CMake configuration.
5. Build as normal

The toolchain options can be passed through `setup.py` accordingly:

1.  python setup.py build -- -DCMAKE_TOOLCHAIN_FILE="C:/Users/Kieran/vcpkg/scripts/buildsystems/vcpkg.cmake" -G "Visual Studio 16 2019" -T host=x64 -A x64 -- /m:6

## Dependencies

* Meshes are represented using [Geometry-Central](https://geometry-central.net/).
* Optional trajectory output uses [NetCDF-cxx4](https://github.com/Unidata/netcdf-cxx4).
### Visualization
* Optional trajectory visualization in Python uses [Polyscope](https://polyscope.run/py/).
* Optional trajectory reading in Python uses [NetCDF4-python](https://github.com/Unidata/netcdf4-python).
