using PyCall


comm_path = "/"*relpath((@__FILE__)*"/../..","/")*"/comm/electrodes_v2.py"
@pyinclude(comm_path)

function init_comm()
    return py"ElectrodeUSB"()
end

function toggle_led(comm::PyObject)
    comm.toggle_led()    
end

function set_led(comm::PyObject, value)
    comm.set_led(value)
end

function set_channels(comm::PyObject, channel_indexes::Vector{<:Int}, channels_pwms::Vector{<:Int})
    comm.set_channels(channel_indexes, channels_pwms)
end

function set_channels(comm::PyObject, channels_pwms::Vector{<:Int})
    comm.set_all_channels(channels_pwms)
end

function clear_channels(comm::PyObject)
    comm.clear_channels()
end

function measure_current(comm::PyObject, channel_indexes)
    return comm.measure_current(channel_indexes)
end

function measure_current(comm::PyObject)
    return comm.measure_all_currents()
end

function get_firmware_version(comm::PyObject)
    return comm.get_firmware_version()
end


if abspath(PROGRAM_FILE) == @__FILE__
    usb_comm = init_comm()
    @show typeof(usb_comm)
    toggle_led(usb_comm)
    sleep(1)
    toggle_led(usb_comm)
    sleep(1)
    set_channels(usb_comm, [12, 10, 7, 5], [255, 255, 255, 255])
    sleep(1)
    clear_channels(usb_comm)
end