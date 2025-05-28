# RRT* Path Finding Algorithm

This is an implementation of the RRT* algorithm for path planning. This code is designed to be run seperately from the main codebase and is intended for testing and development purposes only.

The code is written in C++ and uses the SFML library for visualization. The implementation includes a simple 2D environment with some tracks. 

## Running the Prototype

To run the prototype, you will need to meet the dependencies listed in the CMakeLists.txt file. The prototype can be run using the following commands:

(while in the prototype directory)
```bash
mkdir build
cd build
```

After that, run the following command to generate the Makefile:

```bash
cmake ..
```

Inside the build directory, run the following command to compile the project:

```bash
make
```

Finally, execute the compiled program:

```bash
./rrt
```

## The Executable

Running the executable will open a window displaying the environment including the track, the start point and the finish line rectangle. 

The RRT* algorithm will then begin to search for a path from the start point until it reaches the rectangle from the opposite side. The path tree will be displayed in real-time in white, and whenever a new path is found, it will be displayed in red. 

The algorithm will continuously run until a certain number of iterations after the first path is found. The algorithm will also frequently print out the current path length and the number of iterations to the terminal.