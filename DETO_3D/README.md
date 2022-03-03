# DETO3D

This code is currently in development. It is writen in C++ and requires the LAMMPS submodule. It is designed with the intention of giving the user as much freedom as possible to set up and run there own Topology Optimisation

## Dependancies

| Name   | Description                        | Version                          |
|--------|------------------------------------|----------------------------------|
| Git    | Distributed version control system | Most recent versions should work |
| CMake  | Cross-platform build tool          | 3.10 or above  (same req. as LAMMPS)                    |
| MPI    | Implementation of the MPI standard | Same requirements as LAMMPS      |

## Downloading DETO_3D

use:
```
$ cd <path>
$ git clone git@github.com:Connor-OS/DETO.git
$ cd DETO/DETO_3D
$ git submodule init
$ git submodule update
```
*path* should be the place where you want the source code to be downloaded.

## Building DETO_3D

We rely on [CMake](https://cmake.org) to provide cross-platform build support.

### Linux and MacOS

first we compile LAMMPS. from inside the DETO_3D folder, use:
```
$ cd lammps
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release -DWITH_JPEG=Off -DWITH_PNG=Off -DWITH_FFMPEG=Off -DBUILD_LIB=On -DBUILD_OMP=Off -DPKG_MISC=yes -DPKG_GRANULAR=yes ../cmake
$ cmake --build .
$ cd ../..
```
Then to install DETO_3D, use:
```
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
```
## Folder Layout

| Folder | Description                |
|--------|----------------------------|
| src    | Source files               |
| test   | Collection of test cases   |
| lammps | LAMMPS source folder       |


