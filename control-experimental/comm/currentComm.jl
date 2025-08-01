using PyCall

comm_path = "/"*relpath((@__FILE__)*"/../..","/")*"/comm/current_comm.py"
@pyinclude(comm_path)

function init_comm_current() :: PyObject
    return py"CurrentUSB"()
end

function toggle_led(comm::PyObject)
    comm.toggle_led()    
end

function set_led(comm::PyObject, value::Number)
    comm.set_led(value)
end

function set_channels(comm::PyObject, channel_indexes::Vector{<:Int}, channels_pwms::Vector{<:Number})
    comm.set_channels(channel_indexes, channels_pwms)
end

