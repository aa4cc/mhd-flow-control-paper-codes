#using Conda
using PyCall
using GLMakie
pylon = pyimport("pypylon.pylon")

function initCam(config::Dict{<:Any})
    #=
    init camera paramters from config dictionary
    =#
    camera = pylon.InstantCamera(pylon.TlFactory.GetInstance().CreateFirstDevice())
    camera.Open()

    camera.MaxNumBuffer = config["MaxNumBuffer"]

    camera.AcquisitionFrameRateEnable.SetValue(config["AcquisitionFrameRateEnable"])
    camera.AcquisitionFrameRate.SetValue(config["AcquisitionFrameRate"])

    camera.Width.SetValue(config["Width"]) # MAX -> 1152
    camera.Height.SetValue(config["Height"])
    camera.OffsetX.SetValue(config["OffsetX"])
    camera.PixelFormat.SetValue(config["PixelFormat"])

    camera.ReverseX.Value = config["ReverseX"]
    camera.ReverseY.Value = config["ReverseY"]

    camera.ExposureAuto.SetValue(config["ExposureAuto"])
    camera.ExposureTime.SetValue(config["ExposureTime"]) 
    camera.GainAuto.SetValue(config["GainAuto"])
    camera.Gain.SetValue(config["Gain"])

    camera.StartGrabbing(pylon.GrabStrategy_OneByOne)
    return camera
end

function initCam()
    #=
    init camera paramters to default values
    =#
    camera = pylon.InstantCamera(pylon.TlFactory.GetInstance().CreateFirstDevice())
    camera.Open()

    camera.MaxNumBuffer = 50
    camera.AcquisitionFrameRateEnable.SetValue(true)
    camera.AcquisitionFrameRate.SetValue(30.0)

    camera.Width.SetValue(1024) # MAX -> 1152
    camera.Height.SetValue(1024)
    camera.OffsetX.SetValue(128)
    camera.PixelFormat.SetValue("Mono8") 

    camera.ExposureAuto.SetValue("Off")
    camera.ExposureTime.SetValue(10000.0) 
    camera.GainAuto.SetValue("Off")
    camera.Gain.SetValue(2.0)

    camera.StartGrabbing()
    return camera
end


function clear_buffer(camera)
    #=
    Clears camera buffer
    =#
    camera.StopGrabbing()
    camera.StartGrabbing(pylon.GrabStrategy_OneByOne)
end


function grab(camera)
    #=
    grabs single frame from camera and converts it to values between 0.0 - 1.0
    =#
    grabResult = camera.RetrieveResult(5000, pylon.TimeoutHandling_ThrowException)
    frame = Nothing
    if grabResult.GrabSucceeded()
        frame = convert(Array{Float32}, grabResult.Array) ./ 255
        grabResult.Release()
    else
        println("Error!")
    end
    return frame
end

function grab_as_uint(camera)
    #=
    Grabs frame and returns it as UInt8 Array
    =#
    grabResult = camera.RetrieveResult(5000, pylon.TimeoutHandling_ThrowException)
    frame = Nothing
    if grabResult.GrabSucceeded()
        frame = convert(Array{UInt8}, grabResult.Array)
        grabResult.Release()
    else
        println("Error!")
    end
    return frame
end

function grab!(camera, img)
    #=
    Grabs frame in place
    =#
    grabResult = camera.RetrieveResult(5000, pylon.TimeoutHandling_ThrowException)
    if grabResult.GrabSucceeded()
        #=
        for i = 1:size(grabResult.Array, 1) j = 1:size(grabResult.Array, 2)
            img[i, j] = grabResult.Array[i, j] / 255
        end
        =#
        img .= convert(Array{Float32}, grabResult.Array) ./ 255
        grabResult.Release()
        return true 
    else
        print("Error!")
        return false
    end
end


function show_frame(frame::Array{<:Real}; video::Bool = true)
    #=
    Shows frame, either single image or video, after calling update(), "yield()" needs to be called
    rotr90 used, because plotting in julia is weird, using this corresponds to indexing the img array as row, col like one would expect
    - this may not work its better to use the function with custom axis
    =#
    scene = GLMakie.Scene(camera=GLMakie.campixel!, size=size(frame))
    if(video)
        obs_img = GLMakie.Observable(frame)
        function update(frame::Array{<:Real})
            obs_img[] = rotr90(frame, 1)
        end
        GLMakie.image!(scene, rotr90(frame, 1)) # This is because plotting is weird
        display(scene)
        return update
    else
        GLMakie.image!(scene, rotr90(frame, 1))
        wait(display(scene))
    end
end

function show_frame(frame::Array{<:Real}, ax::T; video::Bool = true) where {T<:GLMakie.Axis}
    #=
    Shows frame, either single image or video, after calling update(), "yield()" needs to be called
    rotr90 used, because plotting in julia is weird, using this corresponds to indexing the img array as row, col like one would expect
    =#
    if(video)
        obs_img = GLMakie.Observable(frame)
        GLMakie.image!(ax, obs_img)
        function update(frame::Array{<:Real})
            obs_img[] = rotr90(frame, 1)
        end
        return update
    else
        GLMakie.image!(ax, rotr90(frame))
    end
end


function closeCam(camera)
    #=
    Closes communication with cam, use this, since this may cause problems when cam not closed
    =#
    camera.Close()
end

#=
--------- Simple Example on how to use ------------
=#

function main()
    record = false
    frames = 300
    frame = 0
    if(record)
        vid_save = zeros(Float32, (1024, 1024, frames))
    end
    cam = initCam()

    try
        img = grab(cam)

        fig = Figure(size = (600, 600))
        ax_cam = GLMakie.Axis(fig[1, 1])

        update_vid = show_frame(zeros(Float32, 1024, 1024), ax_cam, video=true)
        # display the GLMakie scene
        display(fig)
        
        run = true
        # until the window is open
        on(events(fig).window_open) do w_open
            run = w_open
        end
        t1 = time()
        while run
            # read from camera
            # rotl90 used just to view from the same angle as in real world setup
            img = rotl90(grab(cam))
    
            if ((time() - t1 > 10) && record)
                frame += 1
                if (frame >= frames)
                    @info "recording done"
                    break
                end
                @show frame
                vid_save[:, :, frame] = copy(img)
            end
            # update image in the window
            update_vid(img)
            yield()  # To update plot
            
        end
        if (record)
            JLD2.save("measurements/video_with_tracker_aruco.jld2", Dict("vid"=>vid_save); compress=true)
            @info "vid saved as measurements/raw_video.jld2"
        end
    finally
        # in any case, close cam after use
        @info "close cam"
        closeCam(cam)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    using Images
    using BenchmarkTools
    using JLD2, CodecZlib
    main()
end