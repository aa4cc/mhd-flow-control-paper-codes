# Experiment Figure Reproduction

Reproduce the velocity and vorticity control experiment figures for the paper.

The scripts require the experimental data to be available at `../dataset/experiments/`.

## Setup

From this directory:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage

Run the following commands to generate the figures:

### Velocity figure
```bash
julia --project=. generate-experiment-figure-velocity.jl
```

### Vorticity figure
```bash
julia --project=. generate-experiment-figure-vorticity.jl
```

## Output
The generated figures will be saved in the `figures/` directory under the names
`velocity.pdf` and `vorticity.pdf`.