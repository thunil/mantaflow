# Mantaflow #

Mantaflow is an open-source framework targeted at fluid simulation research in Computer Graphics.
Its parallelized C++ solver core, python scene definition interface and plugin system allow for quickly prototyping and testing new algorithms.

In addition, it provides a toolbox of examples for deep learning experiments with fluids. E.g., it contains examples
how to build convolutional neural network setups in conjunction with the [tensorflow framework](https://www.tensorflow.org).

For more information on how to install, run and code with Mantaflow, please head over to our home page at
[http://mantaflow.com](http://mantaflow.com)

![mantaflow logo](resources/mantaflow-logo1.png)

## Building from source ##

This installation guide focusses on Ubuntu 24.04 as a distribution. The process will however look very similar with other distributions, the main differences being the package manager and library package names.

First, install a few pre-requisites:

    sudo apt install pt-get install cmake g++ git python3-dev qt5-qmake libqt5opengl5-dev libtbb-dev libopenvdb-dev

If you want to enable CUDA support, additionally get the latest toolkit from nVidia, and install the appropriate developer driver (be careful though, these driver tend to wreck X11 - get some installation instructions from the web if this is the first time you install CUDA on Linux)

Then, change to the directory to install the source code in, and obtain the current sources from the repository (or alternatively download and extract a source code package)

    git clone https://github.com/thunil/mantaflow.git

To build the project using CMake, set up a build directory and choose the build options (explanation):

    mkdir mantaflow/build
    cd mantaflow/build
    cmake .. -DGUI=ON -DOPENVDB=ON
    make -j8

To build Mantaflow with OpenVDB and Threading Building Blocks (TBB) support, you need to disable GUI. There is an issue with duplicated symbols between Qt and Threading Building Blocks (TBB):

    cmake .. -DTBB=ON -DOPENVDB=ON

That's it! You can now test mantaflow using an example scene

    ./manta ../scenes/simpleplume.py

## Run tests ##

Generate reference date and run tests:

    cd tools/tests
    MANTA_GEN_TEST_DATA=1 ./runTests.py ../../build/manta

Run tests:

    ./runTests.py ../../build/manta
