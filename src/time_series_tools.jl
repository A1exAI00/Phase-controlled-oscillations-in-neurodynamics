#=
A collection of function, that manipulate time series.
=#

using Statistics: mean

########################################################################
########################################################################

function mean_tail(arr, tail_start)
    return mean(arr[Int(round(end*tail_start)):end])
end

function out_bounds(arr, bounds)
    return any(i->(i<bounds[1] || i>bounds[2]), arr)
end

function first_greater_then(seq, check)
    for (i,x) in enumerate(seq)
        if x≥check
            return x
        end
    end
    return NaN
end

function times_of_max(seq, t_seq)
    indexes = indexes_of_max(seq)
    return t_seq[indexes]
end

function indexes_of_max(seq)
    indexes = []
    for i in 1:length(seq)-2
        if seq[i] < seq[i+1] && seq[i+1] > seq[i+2]
            push!(indexes, i+1)
        end
    end
    return indexes
end

function mesure_T(seq, t_span)
    times = times_of_max(seq, t_span)
    T_arr = diff(times)
    return mean(T_arr)
end
