# ==================================================================================================== #
#                                           SVertex.jl                                                 #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#   SVertex (Sparse Vertex) type, that stores a minimal representation                                 #
# -------------------------------------------- TODO -------------------------------------------------- #
#   This is a more general formulation of the expand function. Not finished...                         #
# ==================================================================================================== #

import Base: show, display, @propagate_inbounds
import Base: axes, getindex, setindex!, similar, eachindex
import Base: eltype, length, size, iterate
# import Base: indices#, endof, start, next, done, iteratorsize, iteratoreltype


"""
    SVertex

Datastructure storing a two-particle (three-frequency) quantity in a minimal representation.


"""
struct SVertex{T} <: AbstractDict{T,1}
    len::Int            # length
    nBose::Int          # total number of bose frequencies
    nFermi::Int         # total number of fermi frequencies for each bose frequency
    nB::UInt32          # triple-to-int parameter
    nB2::UInt32         # triple-to-int parameter
    offset::UInt32      # triple-to-int parameter
    data::Array{T,1}    # independent data segments
    map::Dict{UInt32,UInt32}  # map from full range to data
    SVertex(map::Dict{UInt32,UInt32}, data::Array{T,1}, base::UInt32, offset::UInt32, nBose::Int, nFermi::Int) where {T} =
        new{T}(length(data), nBose, nFermi, base, base*base, offset, data, map)
end

# SVertex(redm::ReduceMap, data, base, offset) = SVertex(redm.map, data, UInt32(base), UInt32(offset), 0, 0)
# SVertex(redm::ReduceMap, data, base, offset, nBose, nFermi) = SVertex(redm.map, data, UInt32(base), UInt32(offset), nBose, nFermi)

similar(sv::SVertex, ::Type{T}, dims::Dims) where T = SVertex(sv.map, Array{T}(undef,sv.len), sv.nB, sv.offset, dims[1], dims[2])
similar(sv::SVertex, dims::Dims) = similar(sv, eltype(sv), dims)
similar(sv::SVertex, ::Type{T}) where T = similar(sv, T, size(sv))
similar(sv::SVertex) = similar(sv, eltype(sv), size(sv))

# ------------------ Indexing ------------------
IndexStyle(::SVertex) = IndexLinear()

axes(sv::SVertex) = keys(sv.map)

@propagate_inbounds getindex(sv::SVertex, i::Int, j::Int, k::Int) = 
    get(sv.data, sv.map[tripleToInt(i,j,k,sv.offset,sv.nB,sv.nB2)], missing)
@propagate_inbounds getindex(sv::SVertex, i::Number) = 
    get(sv.data, sv.map[UInt32(i)], missing)
@propagate_inbounds getindex(sv::SVertex, I) = [sv[i] for i in I]

@propagate_inbounds setindex!(sv::SVertex{T}, val::T, i::Int, j::Int, k::Int) where T = 
    sv.data[sv.map[tripleToInt(i,j,k,sv.offset,sv.nB,sv.nB2)]] = val
@propagate_inbounds setindex!(sv::SVertex{T}, val::T, i::Number) where T = 
    sv.data[sv.map[UInt32(i)]] = val
@propagate_inbounds function setindex!(sv::SVertex{T}, val::T, I) where T
    for i in I 
        sv.data[sv.map[UInt32(i)]] = val
    end
end

# ----------------- Iteration ------------------
@propagate_inbounds iterate(sv::SVertex) = iterate(sv.data)
@propagate_inbounds iterate(sv::SVertex, i) = iterate(sv.data, i)
eltype(sv::SVertex{T}) where T = T
length(sv::SVertex) = sv.len
size(sv::SVertex) = (sv.nBose, sv.nFermi, sv.nFermi)
size(sv::SVertex, ind::Int) = (sv.nBose, sv.nFermi, sv.nFermi)[ind...]

# -------------------- IO ----------------------
# Adapted from dict.jl
function Base.show(io::IO, ::MIME"text/plain", sv::SVertex{T}) where {T}
    recur_io = IOContext(io, :SHOWN_SET => sv,
                             :typeinfo => eltype(sv))
    n = 0
    limit = get(io, :limit, false)::Bool
    println(io, "SparseVertex{$T} with $(sv.len) elements. Encoded with base $(sv.nB), offset $(sv.offset)")
    for i in eachindex(sv.data)
        println(io, "  ", sv.data[i])
        limit && n >= 10 && (print(io, "…"); break)
    end
end

# -------------------- Auxilliary ----------------------
@fastmath @inline tripleToInt(i::UInt32, j::UInt32, k::UInt32, nB::UInt32, nB2::UInt32)::UInt32 = 
    UInt32(i*nB2 + j*nB + k)
@fastmath @inline tripleToInt(i::Int, j::Int, k::Int, offset::UInt32, nB::UInt32, nB2::UInt32)::UInt32 =
    UInt32((i+offset)*nB2 + (j+offset)*nB + (k+offset))

@fastmath @inline function intToTriple(::Type{T}, z::UInt32, offset, nB::UInt32) where {T<:Integer}
    r,k = divrem(z,nB)
    i,j = divrem(r,nB)
    return (convert(T,i)-offset,convert(T,j)-offset,convert(T,k)-offset)
end
intToTriple(z::UInt32) = intToTriple(Int64, z)

indices(sv::SVertex) = sort([intToTriple(Int, el, sv.offset, sv.nB) for el in keys(sv.map)])


function full(sv::SVertex{T}) where T
    sv.nBose*sv.nFermi*sv.nFermi != length(sv.map) && throw(ArgumentError("The total number of bosonic and fermionic frequencies must match the length of the map."))
    ind = indice(sv)
    return ind, full(sv, ind)
end


function full!(res::Array{T,3}, sv::SVertex{T}, indices) where T
    i = 1
    for nb in 1:sv.nBose
        for nf in 1:sv.nFermi
            for nfp in 1:sv.nFermi
                res[nb,nf,nfp] = sv[indices[i]...]
                i += 1
            end
        end
    end
end
