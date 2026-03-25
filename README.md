# Project-directory
My project directory for the course Advanced Scientific Programming with Python (3 hp)

Background: My Phd work is centered around control theory for multi-agent systems. It is standard practice in this field to provide simulation examples at the end of every paper to show that your algorithms that you have analyzed/derived work as intended.

Project description: This project focuses on organizing and version-controlling a simulation framework i have developed in Matlab for testing multi-agent formation control algorithms. The goal is to migrate a messy local library of scripts (which together simulate the movement of agents) into a GitHub repository to ensure research reproducibility and historical tracking.

My hope with this project is that I will learn how to manage my (continuously updated) code in structured way, and to present it to the public (reviewers/other researchers in my field) in an accessible and professional way.

# Multi-Agent Formation Control Simulation

## Project Description
This project is a Matlab-based simulation environment for multi-agent systems. It allows for the simulation of agents (unicycle models) that coordinate to track a shared centroid reference trajectory. The system uses real-block DFT matrices and low-frequency projections obtained from truncating a set of Fourier descriptors to specify and maintain formation shapes.

## Features
- **Trajectory Options**: Supports centroid tracking for ellipses, circles, lines, and lemniscates (infinity signs).
- **Control Architecture**: Includes a consensus algorithm and descriptor-based formation control.
- **Visualizations**: Built-in functions for plotting trajectories and animating simulations.

## File Structure
The repository includes the following main files:
- `main_simulator.m`: Main entry point for running a single simulation.
- `xdot_full'.m` Estimator and agent dynamics
- `rk4_step.m` and `simulate_rk4.m`: Runge-Kutta 4th order (RK4) integration functions.
- `init_....m`: Functions for setting inital conditions.
- `plot_....m` and `animate_simulation.m`: Animation and trajectory visualization functions.

## How to Run
1. Open Matlab and navigate to the project directory.
2. Run `main_simulator.m` to start a default simulation.
3. Adjust parameters in `init_params.m` or change trajectories in the `centroid_reference` section of the main script.
