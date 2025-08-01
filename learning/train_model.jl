using JLD
using LinearAlgebra

const TRAJECTORY_COUNT = 7 # Number of original trajectories
const POD_DIMENSION = 512 # Number of POD modes
const ND = 4 # Delay embedding dimension
const NEL = 4 # Number of electrodes
const NMAG = 4 # Number of electromagnets

const TRAJECTORY_PATH = "../dataset/trajectories/"

include("src/symmetry.jl")
include("src/POD.jl")
include("src/EDMDc.jl")

mktempdir() do tmp

    println("Step 1: Running symmetry augmentation...")
    aug_scratch = joinpath(tmp, "augmented") # Path for augmented trajectories
    mkpath(aug_scratch)
    traj_aug_cnt = symmetrize_trajectories(TRAJECTORY_PATH, aug_scratch, TRAJECTORY_COUNT)

    println("Step 2: Running POD computation...")
    pod_scratch = joinpath(tmp, "pod") # Path for POD projected trajectories
    mkpath(pod_scratch)
    Ξ = compute_POD(aug_scratch, pod_scratch, traj_aug_cnt, POD_DIMENSION)

    println("Step 3: Running EDMDc computation...")
    F, G, H̃ = compute_EDMDc(pod_scratch, traj_aug_cnt, POD_DIMENSION, ND, NEL, NMAG)

    H = Ξ * H̃ # Output in velocity space

    println("Saving Koopman model...")
    mkpath("model")

    JLD.save("model/koopman_model.jld", "F", F, "G", G, "H", H, "Xi", Ξ, "nd", ND)

end





