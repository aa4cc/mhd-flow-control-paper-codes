#=
Main for measuring the vector field and/or control
=#
using GLMakie, JLD2, CodecZlib
using JLD
using Dates
include("VFPlot.jl")
include("PIV.jl")
include("Camera.jl")
include("preprocess.jl")
include("config/config.jl")
include("comm/elecComm.jl")
include("comm/coilComm.jl")
include("Controller.jl")


const Δt = 0.5

function cam_handler(config::Dict{<:Any})
    # Handles camera functions
    cam = nothing
    try
        cam = initCam(config["CAMERA"])
    catch err
        if isa(err, LoadError) || isa(err, PyCall.PyError)
            @error err
            @error "Camera probably not connected"
            return nothing, nothing, nothing
        else
            rethrow(err)
        end
    end

    function GRAB()
        rotl90(grab_as_uint(cam))
    end

    function CLEAR()
        clear_buffer(cam)
    end

    function END()
        closeCam(cam)
    end
    
    return GRAB, CLEAR, END
end

function PIV_handler(config::Dict{<:Any})
    # Handles calculating PIV
    PIV_data = nothing
    if config["PIV"]["PIV_TYPE"] == "DEFAULT"
        PIV_data = PIVData(config["CAMERA"]["Width"], config["CAMERA"]["Height"], (config["PIV"]["WINDOW_SIZE_R"], config["PIV"]["WINDOW_SIZE_C"]), sub_pixel_function=config["PIV"]["SUB_PIXEL_FUNCTION"])
    elseif config["PIV"]["PIV_TYPE"] == "OVERLAYED"
        PIV_data = OverlayedPIVData(config["CAMERA"]["Width"], config["CAMERA"]["Height"], (config["PIV"]["WINDOW_SIZE_R"], config["PIV"]["WINDOW_SIZE_C"]), (config["PIV"]["OVERLAY_R"], config["PIV"]["OVERLAY_C"]), sub_pixel_function=config["PIV"]["SUB_PIXEL_FUNCTION"])
    else
        @error "PIV Type:$(config["PIV_TYPE"]) is not yet implemented or is not a valid option"
    end
    fps = config["CAMERA"]["AcquisitionFrameRate"]
    vec_field = zeros(Float32, 2, PIV_data.rows, PIV_data.cols)
    prev_img = nothing
    function VF!(img::Array{<:Real})
        if(isnothing(prev_img))
            prev_img = img
        end
        compute_vector_field!(PIV_data, vec_field, prev_img, img; fps=fps)
        prev_img = img
    end
    return vec_field, VF!
end

function VIS_handler(config::Dict{<:Any}, vec_field_size::Tuple{Int, Int})
    #=
    Handles all visualization
    _____________________________________
    vec_field_size is (rows, cols)
    =#
    fig = Figure(size = (800, 400))
    ax_vf = GLMakie.Axis(fig[1, 1])
    ax_cam = GLMakie.Axis(fig[1, 3])
    update_vid = show_frame(zeros(Float32, config["CAMERA"]["Height"], config["CAMERA"]["Width"]), ax_cam, video=true)
    update_vf = plot_field(zeros(Float32, 2, vec_field_size...), ax_vf, fig, video=true)
    display(fig)
    function PLOT(img::Array{<:Real}, vec_field::Array{<:Real})
        update_vid(img)
        update_vf(vec_field)
        yield()
    end
    return PLOT, fig
end

function COM_handler(config::Dict{<:Any}, electrodes_order::Vector{<:Int})
    # handles communication with control boards
    elec_com = init_comm()
    coil_com = open_port(config["MEASURE"]["COIL_COMM_PATH"])
    function SET(electrodes::Vector{<:Int}, coils::Vector{<:Number})
        set_channels(elec_com, electrodes_order , electrodes) 
        set_coil(coil_com, 1, 4, coils[1])
        set_coil(coil_com, 2, 2, coils[2])
        set_coil(coil_com, 3, 1, coils[3])
        set_coil(coil_com, 4, 0, coils[4])
    end
    function END_COM()
        close_port(coil_com)
    end
    return SET, END_COM
end

const NUM_ELECTRODES = 4
const NUM_COILS = 4

function main()
    # --------------- INIT STUFF
    config = get_config("config/control_flow.toml", "config/control_flow_default.toml")
    GRAB, CLEAR, END = cam_handler(config)
    if isnothing(GRAB)
        return
    end
    vec_field, VF! = PIV_handler(config)
    PLOT, fig = VIS_handler(config, size(vec_field)[2:end])
    
    SET, END_COM = COM_handler(config,[1, 3, 6, 8])

    # Check if we are saving data and prepare storage
    save = config["MEASURE"]["SAVE"]
    num_samples = config["MEASURE"]["NUM_SAMPLES"]
    print("num_samples: $num_samples")
    plot = true
    if(save)
        plot = config["MEASURE"]["PLOT"]
        vec_field_save = zeros(Float32, (size(vec_field)..., num_samples))
        electrodes_save = zeros(Int, (NUM_ELECTRODES, num_samples))
        coils_save = zeros(Float32, (NUM_COILS, num_samples))
        video_save = zeros(UInt8, config["CAMERA"]["Height"], config["CAMERA"]["Width"], num_samples)
    end
    
    quit = false
    on(events(fig).window_open) do w_open
        quit = !(w_open)
    end
    frame_cnt = 0
    file_name = config["MEASURE"]["FILE_NAME"]
    # ------ Timing stuff
    t_prev = time()
    fps = config["CAMERA"]["AcquisitionFrameRate"]
    dt = 1 / fps
    dt_err = (config["MEASURE"]["DT_DEVIATION"] / 100) * dt
    sample = 0
    start_sending = false

    # ------ Control stuff
    t_control = nothing
    controller_steps = config["CONTROL"]["STEPS"]
    controller = setup_mpc()



    @info "save: $save | plot: $plot\nnum samples: $num_samples | file name: $file_name"

    # ------- coil and electrode vector
    electrodes = zeros(Int, (NUM_ELECTRODES))
    coils = zeros(Float32, (NUM_COILS))

    #img_uint = zeros(UInt8, config["CAMERA"]["Height"], config["CAMERA"]["Width"])
    #img = zeros(Float32, config["CAMERA"]["Height"], config["CAMERA"]["Width"])

    t = @elapsed while !quit

        img_uint = GRAB()
        img = convert.(Float32, img_uint) ./ 255
        #img = GRAB()
        equalize!(img)
        VF!(img)
        # ------------------ Remove outliers ------------
        vec_field .= validate(vec_field)
        if(plot)
            PLOT(img, vec_field)
        end
        frame_cnt = frame_cnt + 1
        # ------------- Wait for loop to have ~constant dt
        dt_elapsed = time() - t_prev
        #@show dt_elapsed
        t_prev = time()
        if(!start_sending && (dt_elapsed >= (dt - dt_err)) && (dt_elapsed <= (dt + dt_err)))
            @info "loop stable: dt=$dt_elapsed"
            start_sending = true
            continue
        end
        # ------------ Once stable start sending and measuring
        if(start_sending)
            # ----------- Place function that generates electrodes and coils signal here
            # ----------- To end the loop set quit = true ------------------------------
            # ----------- Vec field and control commands get automatically saved
            # electrodes, coils = some_function(vec_field)
            # ---------------------------------------------------------------------------

            t_now = time()

            if isnothing(t_control) || (t_now - t_control >= Δt)

                if !isnothing(t_control)
                    @info "Loop time:  $(t_now - t_control)"
                end

                t_control = time()

                electrodes, coils = update_controller!(controller, Float64.(vec_field))

                electrodes = clamp.(electrodes, 0, 10)
                coils = clamp.(coils, -1.0, 1.0)

                electrodes = round.(Int, 255 * electrodes / 10)

                coils = Float32.(clamp.(fi_inv.(coils) / 0.4330, -1.0, 1.0))

                SET(electrodes, coils)

                @show electrodes, coils

            end
            
        end

        if(start_sending && save)
            # each measuremnt gets saved here
            sample = sample + 1
            # check for overflow
            if(controller.k > controller_steps + 1)
                quit = true
                @info "max samples reached, quit!"
                break
            end
            vec_field_save[:, :, :, sample] .= vec_field
            electrodes_save[:, sample] .= copy(electrodes)
            coils_save[:, sample] .= copy(coils)
            video_save[:, :, sample] .= img_uint
        end
    end

    # turn off coils and electrodes
    SET(zeros(Int, (NUM_ELECTRODES)), zeros(Float32, (NUM_COILS)))
    sleep(0.2)
    # coils sometimes dont turn off, call twice just to be safe
    SET(zeros(Int, (NUM_ELECTRODES)), zeros(Float32, (NUM_COILS)))
    END_COM()

    if (save)
        # save the actual data, file_name can be set in config/control_flow.toml
        JLD2.save(file_name * "-$(now())" * ".jld2", Dict("vec_field" => vec_field_save[:, :, :, 1:min(sample, end)], "electrodes" => electrodes_save[:, 1:min(sample, end)], "coils" => coils_save[:, 1:min(sample, end)], "video" => video_save[:, :, 1:min(sample, end)]); compress=true)
    end

    
    if(config["MEASURE"]["CREATE_VIDEO"] && save)
        # plot / gif from the data can be made here
    end

    JLD2.save(file_name * "-controller-" * "$(now()).jld", "X", controller.X[:, 1:controller.k-1], "phi", controller.Φ[:, 1:controller.k-1], "psi", controller.Ψ[:, 1:controller.k-1])
    
    END()
    @info "Avarage FPS:$(1/(t/frame_cnt))"
end

#if abspath(PROGRAM_FILE) == @__FILE__
    main()
#end

