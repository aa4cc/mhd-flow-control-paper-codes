# Simulated Control Scripts
This folder contains the code used to run closed-loop control experiments on the magnetohydrodynamic (MHD) flow platform described in the paper.

It reproduces, in simulation, the control strategy demonstrated experimentally in the paper.

## Overview

The experiments are present in two scripts:
- `velocity-control.jl` implements the velocity control strategy using Koopman MPC.
- `vorticity-control.jl` implements the vorticity control strategy using Koopman MPC.

## Requirements
To run these scripts, you first need to install the required Julia packages. You can do this by running

```bash
julia --project=. -e '
using Pkg
Pkg.add(url="https://github.com/aa4cc/MHDSim.jl")
Pkg.add(url="https://github.com/aa4cc/AlternatingMPC.jl")
Pkg.instantiate()
'
```

The scripts also require access to the Koopman predictor, it is assumed to be present in `../learning/model/koopman_model.jld`. 
If it is not available, please, first, run the learning script to train the Koopman model.


## Running the Scripts

You can run the individual scripts as follows:

For the velocity control:
```bash
julia --project=. --threads=auto velocity-control.jl <experiment_name>
```
where `<experiment_name>` is the name of the experiment you want to run `two-vortices`, `jet`, `sides`, or `cross`.

For the vorticity control:
```bash
julia --project=. --threads=auto vorticity-control.jl 
```
The results will be saved in the `paraview` directory with the specified experiment name as `vtu`, and `pvd` files for visualization in ParaView.


