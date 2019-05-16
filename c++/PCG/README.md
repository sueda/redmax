# REDMAX: Efficient & Flexible Approach for Articulated Dynamics

> C++ implementation including Projected Block Jacobi Preconditioner
---
## Installation

The code was developed on Windows and has been tested on MacOS.

### Clone

- Clone this repo to your local machine using 
```shell
$ git clone https://github.com/sueda/redmax.git
$ cd redmax/c++/PCG/
$ mkdir build
```

### Dependencies

Required:

* [CMake](https://cmake.org/ "CMake")
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page "Eigen")
* [JsonCpp](https://github.com/open-source-parsers/jsoncpp "JsonCpp")
* [Intel® MKL](https://software.intel.com/en-us/mkl "MKL")

To run the simulation "online" with rendered output:

* [GLM](https://glm.g-truc.net/0.9.9/index.html "GLM")
* [GLFW](http://www.glfw.org/ "GLFW")
* [GLEW](http://glew.sourceforge.net/ "GLEW")

Optional:

PARDISO is strongly recommended if you want to reproduce some results in the paper.

* [PARDISO](https://www.pardiso-project.org/ "PARDISO")

### Setup

- Compile the code using CMake:

```shell
$ cd build
$ cmake [options] ..
$ make
```
---
## Usage

### Features
### Command Line Switches

---

## Contact

If you would like to contact us for anything regarding REDMAX feel free to email us. 

---
