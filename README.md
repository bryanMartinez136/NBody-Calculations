# NBody-Calculations
N-Body Simulation
This is a CUDA-based N-Body simulation that calculates the motion of multiple particles under gravitational forces and friction. The simulation uses parallel processing to accelerate the computation.

Overview
The N-Body simulation is a computational model used to study the gravitational interactions between multiple particles or celestial bodies. In this implementation, we use CUDA to parallelize the calculations for increased performance.

Code Structure
The code is organized into several sections, each serving a specific purpose:

Constants and Definitions: The code starts by defining various constants, including the number of bodies, physical properties, and computational parameters.

Device Initialization: The init kernel initializes random states for CUDA device threads, which are used to generate random initial positions and velocities for the bodies.

Math Functions: The code defines several mathematical functions for vector operations and normalization.

Random Initialization of Bodies: The randoms kernel generates random initial positions and velocities for each body while ensuring they are appropriately normalized.

Force Calculation: The forces kernel computes gravitational forces between bodies and accumulates the forces on each body.

Update Positions and Velocities: The update kernel updates the positions and velocities of the bodies based on the computed forces.

Main Function: The main function coordinates the simulation by calling the kernels in a loop over a specified number of time steps. It also handles the initialization of bodies, printing output in PDB format, and measuring the execution time.

Usage
To use this N-Body simulation:

Compile the code, ensuring you have CUDA installed and properly configured on your system.

shell
Copy code
nvcc nbody.cu -o nbody
Run the executable with the desired number of time steps:

shell
Copy code
./nbody <number_of_timesteps>
Output
The simulation outputs the positions of particles at each time step in PDB format. These positions represent the simulated trajectory of the bodies.

License
This code is provided under an unspecified license. Please refer to the code's license or contact the author for more information.
