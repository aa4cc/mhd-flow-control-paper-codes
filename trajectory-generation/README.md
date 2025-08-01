# Training Trajectory Generation

This directory contains the code for generating training trajectories for the MHD flow control paper. The code is designed to create a set of trajectories that is used to train the Koopman operator model.

## Requirements
To run these scripts, you first need to install the required Julia packages. You can do this by running

```bash
julia --project=. -e '
using Pkg
Pkg.develop(url="https://github.com/aa4cc/MHDSim.jl")
Pkg.instantiate()
'
```


## You can run the trajectory generation script as follows:

```bash
julia --project=. --threads=auto generate_trajectories.jl  <K> <T>
```

where:
- `<K>` is the number of trajectories to generate.
- `<T>` is the length of each trajectory in seconds.


## Output

The script will generate trajectories based on the specified parameters and save them in the `trajectories` directory. Each trajectory is saved in a separate file named `traj-<i>.jld`, where `<i>` is the trajectory index. It contains the following data:

- `vec_field`: The velocity field data for the trajectory.
- `psi`: The scaling factors of the individual electromagnets.
- `phi`: The voltages applied to the electrodes.

The velocity field, and control inputs are sampled at fixed rate of `Δt = 0.5` seconds, thus contain `Int(T ÷ Δt) + 1` and `Int(T ÷ Δt)` samples, respectively.

Parameters for the trajectories, such as the number of electrodes, electromagnets, and the desired switch time, are defined in the script. 


