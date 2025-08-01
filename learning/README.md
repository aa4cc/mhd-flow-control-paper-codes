# Koopman Model Training

Trains the Koopman linear predictor used for MPC on the MHD flow.  
Steps performed:
1. Symmetry augmentation  
2. POD projection  
3. EDMDc identification  

The script requires trajectory data present under `../dataset/trajectories/`.

## Run

From this directory:

1. Activate and instantiate the environment
```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

2. Execute the training script:

```bash
julia --project=. train_koopman_model.jl
```

## Output

Creates `model/koopman_model.jld` with:
- `F`, `G`, `H`: Koopman model matrices  
- `Xi`: POD basis  
- `nd`: delay-embedding dimension  
