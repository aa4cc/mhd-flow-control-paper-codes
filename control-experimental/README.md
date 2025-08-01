# Experimental Control Scripts

This folder contains the code used to run closed-loop control experiments on the magnetohydrodynamic (MHD) flow platform described in the paper.

⚠️ These scripts require direct access to the experimental hardware and are **not executable** outside the original setup. They are provided for transparency and reproducibility.

## Structure

- `flow_control.jl`, `flow_control_vorticity.jl` – Main scripts for running velocity and vorticity control experiments
- `comm/` – Low-level communication routines for interacting with electrodes, and electromagnets
- `config/` – Configuration files (TOML and Julia) specifying experiment parameters, saving paths, and other settings
- `references/` – Reference velocity fields used as control objectives
- `models/` – Precomputed Koopman models used within the MPC controller
- Other `.jl` files – Supporting modules for state estimation (Kalman filter), camera service, PIV processing, control logic, and plotting