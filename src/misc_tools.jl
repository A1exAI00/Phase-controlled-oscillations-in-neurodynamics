module misc_tools
    using Statistics: mean


########################################################################
########################################################################

    
    function mean_tail(arr, tail_start)
        return mean(arr[Int(round(end*tail_start)):end])
    end
    
    
    function out_bounds(arr, bounds)
        return any(i->(i<bounds[1] || i>bounds[2]), arr)
    end


########################################################################
########################################################################
    
    
    function first_greater_then(seq, check)
        for (i,x) in enumerate(seq)
            if x≥check
                return x
            end
        end
        return NaN
    end
    
    
    function calc_maxes(seq, t_span)
        indexes_of_max = []
        for i in 1:length(seq)-2
            if seq[i] < seq[i+1] && seq[i+1] > seq[i+2]
                push!(indexes_of_max, i+1)
            end
        end
        return t_span[indexes_of_max]
    end
    
    
    function mesure_T(seq, t_span)
        times_of_max = calc_maxes(seq, t_span)
        T_arr = diff(times_of_max)
        return mean(T_arr)
    end


########################################################################
########################################################################

    
    function debug(message, is_debug_output) 
        if is_debug_output println(message) end
    end
    
    
    function debug_display(message, is_debug_output) 
        if is_debug_output display(message) end
    end
    
    
    function progress(message, is_debug_output) 
        if is_debug_output println(message) end 
    end


########################################################################
########################################################################

    
    function elapsed_time_string(time_ns)
        seconds = time_ns/1e9
        munutes = seconds/60
        hours = munutes/60
        return "Elapsed time: $(seconds)s = $(munutes)m = $(hours)h"
    end
    
    
########################################################################
########################################################################
    
    
    export calc_maxes, first_greater_then, mean_tail, out_bounds, mesure_T, debug, debug_display, progress, elapsed_time_string
end
