#=
Misc functions: output, time
=#

#########################################################################################

using Dates

#########################################################################################

function debug(message, is_debug_output) 
    if is_debug_output println(message) end
end

function debug_display(message, is_debug_output) 
    if is_debug_output display(message) end
end

function progress(message, is_debug_output) 
    if is_debug_output println(message) end 
end

#########################################################################################

function elapsed_time_string(time_ns)
    seconds = time_ns/1e9
    munutes = seconds/60
    hours = munutes/60
    return "Elapsed time: $(seconds)s = $(munutes)m = $(hours)h"
end

function get_easy_time()
    return Dates.format(now(), "HH:MM:SS") 
end