module TimeAxisArrays
import IntervalSets: ClosedInterval
using AxisArrays
using ElasticArrays
using ExtendedDates
using PrettyTables
using Tables

export TimeAxisArray, .., diff, lag, lead

const Instant = ExtendedDates.Dates.Instant

abstract type AbstractTimeAxisArray{T} <: AbstractArray{T, 2} end
    
struct TimeAxisArray{T}  <: AbstractTimeAxisArray{T}  
    data::AxisArray
end

mutable struct MTimeAxisArray{T} <: AbstractTimeAxisArray{T}  
    data::AxisArray
end

function TimeAxisArray(data_arg::AbstractMatrix{T},
                       time::AbstractVector,
                       variables::AbstractVector{Symbol},
                       ) where T
    return TimeAxisArray{T}(AxisArray(data_arg, time=time, variables=variables)) 
end

function add_column!(tdf::TimeAxisArray, label, data)
    append!(tdf.data.data, data)
    push!(tdf.data.axes[2].val, label)
tdf
end

Tables.istable(TA::AbstractTimeAxisArray) = true
Tables.schema(TA::AbstractTimeAxisArray) = Tables.Schema(vcat("Periods", TA.data.axes[2].val),
                                                         vcat(typeof(TA.data.axes[1]),
                                                              fill(eltype(typeof(TA.data)),
                                                                   length(TA.data.axes[2].val))))
Tables.columnaccess(::Type{<:AbstractTimeAxisArray}) = true
Tables.columns(TA::AbstractTimeAxisArray) = TA
# required AbstractColumns object methods
Tables.getcolumn(TA::AbstractTimeAxisArray, ::Type{T}, col::Int, nm::Symbol) where {T} = (col == 1) ? format.(TA.data.axes[1].val) : TA.data.data[:, col]
Tables.getcolumn(TA::AbstractTimeAxisArray, nm::Symbol) = (nm == :Periods) ? format.(TA.data.axes[1].val) : TA.data.data[:, nm]
Tables.getcolumn(TA::AbstractTimeAxisArray, i::Int) = (i == 1) ? format.(TA.data.axes[1].val) : TA.data.data[:, i]
function Tables.columnnames(TA::AbstractTimeAxisArray)
    @show TA.data
    return vcat(:Periods, TA.data.axes[2].val)    
end

# declare that any AbstractTimeAxisArray defines its own `rows` method
rowaccess(::Type{<:AbstractTimeAxisArray}) = true
# just return itself, which means AbstractTimeAxisArray must iterate `AbstractRow`-compatible objects
rows(TA::AbstractTimeAxisArray) = TA
# a custom row type; acts as a "view" into a row of an AbstractVecOrMat
struct MatrixRow <: Tables.AbstractRow
    row::Int
    source::AbstractTimeAxisArray
end
# the iteration interface, at a minimum, requires `eltype`, `length`, and `iterate`
# for `AbstractTimeAxisArray` `eltype`, we're going to provide a custom row type
Base.eltype(TA::AbstractTimeAxisArray{T}) where T = T
Base.length(TA::AbstractTimeAxisArray) = size(TA.data, 1)

Base.iterate(TA::AbstractTimeAxisArray, st=1) = st > length(TA) ? nothing : (MatrixRow(st, TA), st + 1)

# required `AbstractRow` interface methods (same as for `AbstractColumns` object before)
# but this time, on our custom row type
getcolumn(m::MatrixRow, ::Type, col::Int, nm::Symbol) =
    getfield(getfield(m, :source), :data)[getfield(m, :row), col]
getcolumn(m::MatrixRow, i::Int) =
    getfield(getfield(m, :source), :data)[getfield(m, :row), i]
getcolumn(m::MatrixRow, nm::Symbol) =
    getfield(getfield(m, :source), :data)[getfield(m, :row), getfield(getfield(m, :source), :lookup)[nm]]
columnnames(m::MatrixRow) = names(getfield(m, :source))

# Broadcasting
Base.broadcastable(TA::AbstractTimeAxisArray) = TA

struct AbstractTimeAxisArrayStyle <: Base.Broadcast.BroadcastStyle end

Base.Broadcast.BroadcastStyle(::Type{<:AbstractTimeAxisArray}) =
    AbstractTimeAxisArrayStyle()

Base.Broadcast.BroadcastStyle(::AbstractTimeAxisArrayStyle, ::Base.Broadcast.BroadcastStyle) = AbstractTimeAxisArrayStyle()
Base.Broadcast.BroadcastStyle(::Base.Broadcast.BroadcastStyle, ::AbstractTimeAxisArrayStyle) = AbstractTimeAxisArrayStyle()
Base.Broadcast.BroadcastStyle(::AbstractTimeAxisArrayStyle, ::AbstractTimeAxisArrayStyle) = AbstreactTimeAxisArrayStyle()

function copyto_widen!(res::AbstractVector{T}, bc::Base.Broadcast.Broadcasted,
                       pos, col) where T
    for i in pos:length(Base.axes(bc)[1])
        val = bc[CartesianIndex(i, col)]
        S = typeof(val)
        if S <: T || promote_type(S, T) <: T
            res[i] = val
        else
            newres = similar(Vector{promote_type(S, T)}, length(res))
            copyto!(newres, 1, res, 1, i-1)
            newres[i] = val
            return copyto_widen!(newres, bc, i + 1, col)
        end
    end
    return res
end

function getcolbc(bcf::Base.Broadcast.Broadcasted{Style}, colind) where {Style}
    # we assume that bcf is already flattened and unaliased
    newargs = map(bcf.args) do x
        Base.Broadcast.extrude(x isa AbstractTimeAxisArray ? x[:, colind] : x)
    end
    Base.Broadcast.Broadcasted{Style}(bcf.f, newargs, bcf.axes)
end

function Base.copy(bc::Base.Broadcast.Broadcasted{AbstractTimeAxisArrayStyle})
    ndim = length(Base.axes(bc))
    if ndim != 2
        throw(DimensionMismatch("cannot broadcast an AbstractTimeAxisArrayStyle into $ndim dimensions"))
    end

    
    bcf = Base.Broadcast.flatten(bc)
    first_ta = true
    local periods, continuous
    colnames = []
    for ta in bcf.args
        if ta isa AbstractTimeAxisArray
            if first_ta
                periods = getfield(ta.data, :axes)[1].val
                first_ta = false
            elseif getfield(ta.data, :axes)[1].val != periods
                error("TimeAxisArrays don't have the same periods")
            end
        end
    end
    colnames = unique!([ta.data.axes[2].val for ta in bcf.args if ta isa AbstractTimeAxisArray])
    if length(colnames) > 1 && any(length(c) != 1 for c in colnames)
        wrongnames = setdiff(union(colnames...), intersect(colnames...))
        if isempty(wrongnames)
            throw(ArgumentError("Column names in broadcasted AbstractTimeAxisArray " *
                                "must have the same order"))
        else
            msg = join(wrongnames, ", ", " and ")
            throw(ArgumentError("Column names in broadcasted AbstractTimeAxisArray must match. " *
                                "Non matching column names are $msg"))
        end
    end
    nrows = length(Base.axes(bcf)[1])
    ncols = length(colnames[1])
    ta = TimeAxisArray{Float64}(AxisArray(Matrix(undef, nrows, ncols), periods, colnames[1]))
    for i in Base.axes(bcf)[2]
        if nrows == 0
            col = Any[]
        else
            bcf′ = getcolbc(bcf, i)
            v1 = bcf′[CartesianIndex(1, i)]
            startcol = similar(Vector{typeof(v1)}, nrows)
            startcol[1] = v1
            col = copyto_widen!(startcol, bcf′, 2, i)
        end
        ta[:, colnames[1][i]] = col
    end
    return ta
end

function Base.copyto!(ta::AbstractTimeAxisArray, bc::Base.Broadcast.Broadcasted)
    copyto!(ta.data, bc)
end

import Base: getproperty, getindex, length, ndims, setindex!, show, size, diff

Base.ndims(::AbstractTimeAxisArray) = 2
Base.ndims(::Type{<:AbstractTimeAxisArray}) = 2

function getproperty(TA::TimeAxisArray{T}, s::Symbol) where T
    data = getfield(TA, :data)
    if s == :data
        return data
    elseif s == :axes
        return getfield(TA.data, :axes)
    end
    return TimeAxisArray{T}(getindex(data, :, [s]))
end

function format(p::D) where {D <: Instant}
    return ExtendedDates.Dates.format(p)
end

format(p) = string(p)

function Base.show(io::IO, TA::TimeAxisArray)
    aa = TA.data
    @show aa
    pretty_table(aa, nosubheader=true)
end

Base.size(TA::TimeAxisArray) = size(TA.data)

Base.getindex(TA::TimeAxisArray{T}, idx::Symbol...) where T = TimeAxisArray{T}(view(TA.data, :, [idx...]))
function Base.getindex(TA::TimeAxisArray{T}, idx::D) where {T, D <: Instant}
    aa = TA.data
    return TimeAxisArray(reshape(view(aa, idx, :), 1, size(aa, 2)),
                         [idx],
                         aa.axes[2].val)
end
Base.getindex(TA::AbstractTimeAxisArray{T}, idx::ClosedInterval{<: Instant}) where T = TimeAxisArray{T}(view(TA.data, idx, :))
Base.getindex(TA::AbstractTimeAxisArray{T}, idx1::D, idx2::Symbol) where {T, D <: Instant} = TimeAxisArray(hcat(TA.data[idx1, idx2]), [idx1], [idx2])
Base.getindex(TA::AbstractTimeAxisArray{T}, idx::CartesianIndex) where T = getindex(TA.data, idx)
function Base.getindex(TA::AbstractTimeAxisArray{T}, idx1...) where T
    @show idx1
    return TimeAxisArray{T}(view(TA.data, idx1...))
end

function Base.setindex!(TA::TimeAxisArray{T}, x, idx1, idx2) where T
    d = TA.data
    if !isa(idx2, Colon) && !checkvariables(d.axes[2].val, idx2)
        d = addvariables(d, idx2)
    end
    d = addperiods(d, idx1)
    setindex!(d, x, idx1, idx2)
    return TimeAxisArray{T}(d)
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

function Base.diff(TA::TimeAxisArray{T}, k::Int64=1) where T
    d = TA.data
    n, m = size(d)
    abs(k) > n - 1 && return(TimeAxisArray(repeat([missing], n, m), d.axes[1].val, d.axes[2].val))
    d1 = d[time=k+1:n]
    return TimeAxisArray{T}(AxisArray(vcat(repeat([missing], k, m),
                                           d1 - d[time=1:n-k]),
                                      time=d.axes[1].val,
                                      variables=d.axes[2].val))
end

function TimeAxisArray(f::Function, TA::TimeAxisArray{T}, args...) where T
    aa = TA.data
    f1(x) = f(x, args...)
    return TimeAxisArray{T}(map(f1, aa), aa.axes[1].val, aa.axes[2].val)
end

function lag(TA::TimeAxisArray{T}, k::Int64=1) where T
    d = TA.data
    n, m = size(d)
    abs(k) > n - 1 && return(TimeAxisArray{T}(repeat([missing], n, m), d.axes[1].val, d.axes[2].val))
    return TimeAxisArray{T}(AxisArray(vcat(repeat([missing], k, m),
                                        d[time=1:n-k]),
                                   time=d.axes[1].val,
                                   variables=d.axes[2].val))
end

function lead(TA::TimeAxisArray{T}, k::Int64=1) where T
    d = TA.data
    n, m = size(d)
    abs(k) > n - 1 && return(TimeAxisArray{T}(repeat([missing], n, m), d.axes[1].val, d.axes[2].val))
    return TimeAxisArray{T}(AxisArray(vcat(d[time=k+1:n],
                                        repeat([missing], k, m)),
                                   time=d.axes[1].val,
                                   variables=d.axes[2].val))
end


end # module TimeAxisArrays
