module TimeAxisArrays
import IntervalSets: ClosedInterval
using AxisArrays
using ExtendedDates
using PrettyTables

export TimeAxisArray, .., diff, lag, lead

const Instant = ExtendedDates.Dates.Instant

struct TimeAxisArray
    data::AxisArray
end

function TimeAxisArray(data_arg::AbstractMatrix,
                       time::AbstractVector,
                       variables::AbstractVector{Symbol},
                       )
    return TimeAxisArray(AxisArray(data_arg, time=time, variables=variables)) 
end

get_axisarray(TA::TimeAxisArray) = getfield(TA, :data)

Base.broadcastable(TA::TimeAxisArray) = TA
struct TimeAxisArrayStyle <: Base.Broadcast.BroadcastStyle end

Base.Broadcast.BroadcastStyle(::Type{<:TimeAxisArray}) =
    TimeAxisArrayStyle()

Base.Broadcast.BroadcastStyle(::TimeAxisArrayStyle, ::Base.Broadcast.BroadcastStyle) = TimeAxisArrayStyle()
Base.Broadcast.BroadcastStyle(::Base.Broadcast.BroadcastStyle, ::TimeAxisArrayStyle) = TimeAxisArrayStyle()
Base.Broadcast.BroadcastStyle(::TimeAxisArrayStyle, ::TimeAxisArrayStyle) = TimeAxisArrayStyle()

import Base: getproperty, getindex, length, setindex!, show, size, diff

function getproperty(TA::TimeAxisArray, s::Symbol)
    return TimeAxisArray(getindex(get_axisarray(TA), :, [s]))
end

function format(p::D) where {D <: Instant}
    return ExtendedDates.Dates.format(p)
end

format(p) = string(p)

function Base.show(io::IO, TA::TimeAxisArray)
    aa = get_axisarray(TA)
    pretty_table(hcat(format.(aa.axes[1].val), aa.data), header=vcat("Period",collect(aa.axes[2])))
end

Base.size(TA::TimeAxisArray) = size(get_axisarray(TA))

Base.getindex(TA::TimeAxisArray, idx::Symbol...) = TimeAxisArray(view(get_axisarray(TA), :, [idx...]))
function Base.getindex(TA::TimeAxisArray, idx::D) where D <: Instant
    aa = get_axisarray(TA)
    return TimeAxisArray(reshape(view(aa, idx, :), 1, size(aa, 2)),
                         [idx],
                         aa.axes[2].val)
end
Base.getindex(TA::TimeAxisArray, idx::ClosedInterval{<: Instant}) = TimeAxisArray(view(get_axisarray(TA), idx, :))
Base.getindex(TA::TimeAxisArray, idx1::D, idx2::Symbol) where D <: Instant = TimeAxisArray(hcat(get_axisarray(TA)[idx1, idx2]), [idx1], [idx2])
Base.getindex(TA::TimeAxisArray, idx1, idx2...) = TimeAxisArray(view(get_axisarray(TA), idx1, [idx2...]))

function Base.setindex!(TA::TimeAxisArray, x, idx1, idx2)
    d = get_axisarray(TA)
    if !isa(idx2, Colon) && !checkvariables(d.axes[2].val, idx2)
        d = addvariables(d, idx2)
    end
    d = addperiods(d, idx1)
    setindex!(d, x, idx1, idx2)
    return TimeAxisArray(d)
end

function addperiods(AA::AxisArray, idx::AbstractVector)
    p1 = max(AA.axes[1][1], idx[1])
    p2 = min(AA.axes[1][end], idx[end])
    (p1 == idx[1]) && (p2 == idx[end]) && return AA 
    n, m = size(AA.data)
    d = hcat( repeat([missing], p1 - idx[1], m),
              AA.data,
              repeat([missing], idx[end] - p2, m))
    return AxisArray(d, time=collect(p1:p2), variable=AA.axes[2].val)
end

function addperiods(AA::AxisArray, idx)
    isa(idx, Colon) && return AA
    p1 = max(AA.axes[1][1], idx)
    p2 = min(AA.axes[1][end], idx)
    (p1 == p2) && return AA 
    n, m = size(AA.data)
    d = hcat( repeat([missing], p1 - idx, m),
              AA.data,
              repeat([missing], idx - p2, m)
              )
    return AxisArray(d, time=collect(p1:p2), variable=AA.axes[2].val)
end

function checkvariables(old, new)
    isa(new, Symbol) && return true
    for n in new
        if !(n in old)
            return false
        end
    end
    return true
end

function addvariables(AA::AxisArray, idx)
    pv = union(AA.axes[2].val, idx)
    n, m = size(AA.data)
    k = length(pv) - m
    (k == 0) && return(AA)
    d = vcat(AA.data,
             repeat([missing], n, m))
    return AxisArray(d, time=p1:p2, variable=AA.axes[2].val)
end

function Base.diff(TA::TimeAxisArray, k::Int64=1)
    d = get_axisarray(TA)
    n, m = size(d)
    abs(k) > n - 1 && return(TimeAxisArray(repeat([missing], n, m), d.axes[1].val, d.axes[2].val))
    d1 = d[time=k+1:n]
    return TimeAxisArray(AxisArray(vcat(repeat([missing], k, m),
                                        d1 - d[time=1:n-k]),
                                   time=d.axes[1].val,
                                   variables=d.axes[2].val))
end

function TimeAxisArray(f::Function, TA::TimeAxisArray, args...)
    aa = get_axisarray(TA)
    f1(x) = f(x, args...)
    return TimeAxisArray(map(f1, aa), aa.axes[1].val, aa.axes[2].val)
end

function lag(TA::TimeAxisArray, k::Int64=1)
    d = get_axisarray(TA)
    n, m = size(d)
    abs(k) > n - 1 && return(TimeAxisArray(repeat([missing], n, m), d.axes[1].val, d.axes[2].val))
    return TimeAxisArray(AxisArray(vcat(repeat([missing], k, m),
                                        d[time=1:n-k]),
                                   time=d.axes[1].val,
                                   variables=d.axes[2].val))
end

function lead(TA::TimeAxisArray, k::Int64=1)
    d = get_axisarray(TA)
    n, m = size(d)
    abs(k) > n - 1 && return(TimeAxisArray(repeat([missing], n, m), d.axes[1].val, d.axes[2].val))
    return TimeAxisArray(AxisArray(vcat(d[time=k+1:n],
                                        repeat([missing], k, m)),
                                   time=d.axes[1].val,
                                   variables=d.axes[2].val))
end


end # module TimeAxisArrays
