using LibSerialPort


const tab= [0,49345,49537,320,49921,960,640,49729,50689,1728,1920,51009,1280,50625,
            50305,1088,52225,3264,3456,52545,3840,53185,52865,3648,2560,51905,52097,
            2880,51457,2496,2176,51265,55297,6336,6528,55617,6912,56257,55937,6720,
            7680,57025,57217,8000,56577,7616,7296,56385,5120,54465,54657,5440,55041,
            6080,5760,54849,53761,4800,4992,54081,4352,53697,53377,4160,61441,12480,
            12672,61761,13056,62401,62081,12864,13824,63169,63361,14144,62721,13760,
            13440,62529,15360,64705,64897,15680,65281,16320,16000,65089,64001,15040,
            15232,64321,14592,63937,63617,14400,10240,59585,59777,10560,60161,11200,
            10880,59969,60929,11968,12160,61249,11520,60865,60545,11328,58369,9408,
            9600,58689,9984,59329,59009,9792,8704,58049,58241,9024,57601,8640,8320,
            57409,40961,24768,24960,41281,25344,41921,41601,25152,26112,42689,42881,
            26432,42241,26048,25728,42049,27648,44225,44417,27968,44801,28608,28288,
            44609,43521,27328,27520,43841,26880,43457,43137,26688,30720,47297,47489,
            31040,47873,31680,31360,47681,48641,32448,32640,48961,32000,48577,48257,
            31808,46081,29888,30080,46401,30464,47041,46721,30272,29184,45761,45953,
            29504,45313,29120,28800,45121,20480,37057,37249,20800,37633,21440,21120,
            37441,38401,22208,22400,38721,21760,38337,38017,21568,39937,23744,23936,
            40257,24320,40897,40577,24128,23040,39617,39809,23360,39169,22976,22656,
            38977,34817,18624,18816,35137,19200,35777,35457,19008,19968,36545,36737,
            20288,36097,19904,19584,35905,17408,33985,34177,17728,34561,18368,18048,
            34369,33281,17088,17280,33601,16640,33217,32897,16448]


function open_port(port::String)
    sp = LibSerialPort.open(port, 230400)
    return sp
end

function close_port(port::LibSerialPort.SerialPort)
    close(port)
end 

function crc(data::Array{})
    buff = UInt8[data...]
    count = length(buff)
    dl::UInt8 = 0
    dh::UInt8 = 0
    for k = 1:count
        i = UInt16(((buff[k] ⊻ dl)&255))
        dl = (dh ⊻ (tab[i + 1] % 256))
        dh = tab[i + 1] >>> 8
    end
    return dl, dh
end

function send(port::LibSerialPort.SerialPort, adr::Int, cmd::Int, par)
    #=
    port: serial port,
    adr: coil address
    cmd: command number
    par: command parameters
    =#
    message = UInt8[1, 128 + adr, 128 + cmd, par..., 4]
    dl, dh = crc(message)
    #println(Int[[message..., dl, dh]...])
    write(port, UInt8[message..., dl, dh])
end 

function set_coils(port::LibSerialPort.SerialPort, adr::Int, duties::Vector{<:Number})
    dutys = [duties...]
    if length(dutys) == 1
        dutys = ones(4) .* dutys[1]
    end
    if length(dutys) != 4
        return
    end
    params = Array{UInt8}(undef, 3*4)
    for i=1:4
        dutys[i] = min(1.0, max(-1, dutys[i]))
        dutys[i] = dutys[i] * 2047
        if(dutys[i] < 0)
            dutys[i] = 4096 + dutys[i]
        end
        for (index, b) in enumerate(string(round(Int16, dutys[i]), base = 16, pad = 3))
            params[(i-1) * 3 + index] = UInt8(codepoint(b))
        end
    end
    #println(Int[params...])
    send(port, adr, 17, params)
end

function set_coil(port::LibSerialPort.SerialPort, adr::Int, coil::Int, duty)
    duty = min(1.0, max(-1.0, duty))
    coil = min(3, max(0, coil))
    coil = Int(codepoint(string(coil)[1]))

    duty *= 2047
    duty = round(Int16, duty)
    if(duty < 0)
        duty = -duty
        duty = 4096 - duty
    end
    params = Array{UInt8}(undef, 3)
    for (index, b) in enumerate(string(round(Int16, duty), base = 16, pad = 3))
        params[index] = UInt8(codepoint(b))
    end
    send(port, adr, 16, [coil, params...])
end

function send_response(port::LibSerialPort.SerialPort, adr::Int, cmd::Int, par::Array{})
    set_read_timeout(port, 1)
    if(bytesavailable(port) > 0)
        read(port)
    end
    send(port, adr, cmd, par)
    header = 0
    try
        header = read(port, UInt8)
    catch e
        if isa(e, LibSerialPort.Timeout)
            error("No reply received in timeout!")
        end
    end
    if header != 2
        set_read_timeout(port, 0.05)
        error("Invalid header -> Unparsable: ", Int(header))
    end

    adr_r = read(port, UInt8) - 128
    if(adr_r != adr)
        set_read_timeout(port, 0.05)
        error("Invalid adress -> Reply: ", adr_r + 128)
    end

    cmd_r = read(port, UInt8) - 128
    if(cmd_r != cmd)
        set_read_timeout(port, 0.05)
        error("Invalid command -> Reply: ", cmd_r + 128)
    end
    reply = zeros(UInt8, 64)
    p = 1
    d = read(port, UInt8)
    while(d != 4)
        reply[p] = d
        p += 1
        d = read(port, UInt8)
    end
    reply = reply[1:p-1]
    crc_A = [0, 0]
    crc_A[1], crc_A[2] = crc([2, (adr_r + 128), (cmd_r + 128), reply..., 4])
    crc_B = [0, 0]
    crc_B[1] = read(port, UInt8)
    crc_B[2] = read(port, UInt8)
    if(crc_A != crc_B)
        set_read_timeout(port, 0.05)
        error("Invalid CRC: ", crc_A, crc_B)
    end
    set_read_timeout(port, 0.05)
    return reply
end

function discover(port::LibSerialPort.SerialPort, rng::Int = 128)
    modules = []
    for m = 1:rng
        println("Looking for address: ", m)
        try
            reply = send_response(port, m, 15, [])
            append!(modules, m)
            println("Module found at address: ", m)
        catch e
            println(e)
        end
    end
    println(modules)
end

function test()
    port = open_port("/dev/tty.usbserial-1120")
    #discover(port, 4)
    #set_coil(port, 1, 1, 0.8)
    set_coils(port, 1, [0.8, 0.8, 0.8, 0.8])
    sleep(1)
    set_coils(port, 1, [0.0, 0.0, 0.0, 0.0])
    #set_coil(port, 1, 1, 0.0)
    close_port(port)
end

if abspath(PROGRAM_FILE) == @__FILE__
    test()
end