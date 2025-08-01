using TOML

function get_config(path::String, path_default::String)
    conf = TOML.parsefile(path)
    default = TOML.parsefile(path_default)
    check_types(conf, default)
    return conf
end

function check_types(conf::Dict{<:Any}, default::Dict{<:Any})
    for key in keys(default)
        if !haskey(conf, key)
            @warn "Key: $key missing in config file, replacing with default"
            conf[key] = default[key]
        else
            if(typeof(conf[key]) == typeof(default[key]))
                if(isa(conf[key], Dict))
                    check_types(conf[key], default[key])
                end
            elseif(isa(conf[key], Union{Int16, Int32, Int64}) && isa(default[key], Union{Float32, Float64}))
                @warn "Value of key: $key is Int should be Float... changing to Float64"
                conf[key] = Float64(conf[key])
            else
                @warn "Value of key: $key has different type than default: $(typeof(conf[key])) != $(typeof(default[key])) replacing with default"
                conf[key] = default[key]
            end
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    @show get_config("config.toml", "default_config.toml")
end