# ![PICKSC Logo](http://exodus.physics.ucla.edu/~uclapic/repo_images/PICKSC-logo-OSHUN.png)

# OSHUN

OSHUN is a parallel Vlasov-Fokker-Planck plasma simulation code that employs an arbitrary-order spherical harmonic velocity-space decomposition.

This is the UCLA Plasma Simulation Group's official open-source repository for OSHUN. This repository houses a 1D C++ version and a 1D pure Python version (with optional C-modules for improved performance).  A 2D C++ version will be added shortly.

# Upon cloning the repository

If you clone this repository, we ask that you __please contact__ Ben Winjum (bwinjum@ucla.edu).  The development of OSHUN relies on grant funding for which we are required to report code usage, and we would greatly appreciate being able to maintain accurate user-number tallies.

Please also feel free to send an email to this address if you would like assistance using the code, if you have any questions related to the code base, or if you have suggestions for improvements.  We are in the process of establishing a user forum for discussions related to the use and development of OSHUN and will be more than happy to include you in such a forum if you desire.

# Installation for 1D C++ code

After cloning the repository, use a terminal (i.e. command line) to navigate to the OSHUN repository. Then go into the 1d_cpp subdirectory. From there:

```
./install.sh
```

Please note: The build may need to be configured for compilers and HDF5 libraries on your system. To set these configuration options, edit ```/1d_cpp/source/makefile```.

The build process creates an executable ```oshun-1d.e``` in the bin directory. Example inputdecks are available in the input directory.  The executable and inputdeck files can be copied and moved to any desired location. Output will be placed in the current directory.

There is a HTML manual in ```1d_cpp/manual.html```. There are also loaders and plotters to perform Braginskii-style transport tests using inputdecks provided in ```1d_cpp/inputdecks-for-example-1```

# Installation for 1D Python code

The 1D Python version of OSHUN uses the Cmake build system. Make sure this is installed on your system. After cloning the repository, use a terminal (i.e. command line) and navigate to the OSHUN repository. Then go into the 1d_python subdirectory. From there:

```
mkdir build
cd build
cmake ..
make
```

(Note: on windows, type ```nmake``` instead of make). (Also note: The build can be heavily configured - although its usually not necessary. To see the configuration options, play with the cmake a bit by replacing ```cmake ..``` with ```ccmake ..``` or ```cmake-gui ..```)

The build process creates a file ```oshun.zip``` in the build directory. This file consists of the contents of the ```/1d_python/source/python/oshun_modules``` directory along with the compiled results (i.e. shared libraries) of the C-code in the  ```/1d_python/source/c/oshun_modules``` directory. Let's run a sample simulation (calculation of the Spitzer-Harm heat conduction limit) in the ```/1d_python/examples``` directory. Still assuming we are in the build directory, type:
```
python oshun.zip ../examples/sample-spitzer.py
```
(With Windows replace all ```/``` with ```\```)

Output should immediately be in the ```output``` subdirectory. The ```oshun.zip``` file can be copied and moved to any desired location. Input decks can exist in any directory as well. Just make sure the fully-qualified names are passed to python (i.e. ```python /full/path/to/oshun.zip /full/path/to/inputdeck```). Output will be placed in the current directory.
