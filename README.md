# Codes related to paper "Data-driven control of a magnetohydrodynamic flow"

This repository contains all code used for the data generation, learning, simulation, and experimental validation of the Koopman-based MPC strategy for magnetohydrodynamic (MHD) flow control as presented in the paper "Data-driven control of a magnetohydrodynamic flow" by Uchytil et al. ([arXiv:2507.12479](https://arxiv.org/abs/2507.12479)).

## Structure

- `trajectory-generation/` – Scripts to generate the simulated flow data used for training
- `learning/` – Koopman model training pipeline (symmetry augmentation, POD projection, EDMDc)
- `control-simulation/` – Simulated Koopman-MPC control of the MHD flow
- `control-experimental/` – Scripts for running the closed-loop control on the physical platform (non-executable without hardware)
- `figure-reproduction/` – Scripts for reproducing experiment figures (velocity and vorticity control)


## Further code
This repository depends on two external Julia packages:
- [`AlternatingMPC.jl`](https://github.com/aa4cc/AlternatingMPC.jl) – Alternating optimization solver for Koopman-MPC
- [`MHDSim.jl`](https://github.com/aa4cc/MHDSim.jl) – Simulation environment for the MHD flow


## Dataset
The `dataset/` folder (not included in this repository) contains all training trajectories, and experimental data used in the paper. It is not included in this repository due to its size but required for generating the figures and training the Koopman model. Please obtain it from Zenodo: [10.5281/zenodo.15782793](https://doi.org/10.5281/zenodo.15782793), and place in the root directory of this repository.


## Notes
Hardware-specific code (in `control-experimental/`) is included for documentation but cannot be run without access to the experimental platform.
