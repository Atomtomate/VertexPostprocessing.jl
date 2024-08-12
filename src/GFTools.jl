# ==================================================================================================== #
#                                           GFTools.jl                                                 #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#   Functions for computations involving  Green's functions.                                           #
# -------------------------------------------- TODO -------------------------------------------------- #
# ==================================================================================================== #


# ---------------- General Helpers -----------------
"""
    get_conjsymm_f(f::Array{Complex{Float64},1}, i::Int64)

Returns f(νn) for a 1-base array `f` with the property conj(f(-ν)) = f(ν) and ν being a fermionic Matsubara Frequecy.
`i` is the index of the Matsubara fequency `(2n+1)π/β`.
Example:
```
get_conjsymm_f(f, 0) # returns f[1]
get_conjsymm_f(f, 1) # returns f[2]
get_conjsymm_f(f,-1) # returns conj(f[1]), because of ν_{-1} = 2(-1) + 1 = -(2(0) + 1) = - ν_0
```
"""
@inline get_conjsymm_f(f::Array{Complex{Float64},1}, i::Int64) = (i < 0) ? conj(f[-i]) : f[i+1]


# ----------- Single Particle solution -------------
"""
    restore_1pt_GF(config_f::String, andpar_f::String; nFreq::Int=2000)

Reads config and anderson parameters and computes single particle Green's function.
Returns:
  - U::Float64
  - β::Float64
  - p::AndersonParameters
  - GImp::Vector{ComplexF64}
  - nden::Float64
"""
function restore_1pt_GF(config_f::String, andpar_f::String; nFreq::Int=2000)
    cfg         = TOML.parsefile(config_f)
    U           = cfg["parameters"]["U"]
    β           = cfg["parameters"]["beta"]
    ϵₖ, Vₖ, μ   = read_anderson_parameters(andpar_f)
    p  = AIMParams(ϵₖ, Vₖ)
    basis  = jED.Basis(length(ϵₖ) + 1);
    νnGrid = jED.OffsetVector([1im * (2*n+1)*π/β for n in 0:nFreq-1], 0:nFreq-1)
    model  = AIM(ϵₖ, Vₖ, μ, U)
    es     = Eigenspace(model, basis, verbose=false);
    GImp, nden = calc_GF_1(basis, es, νnGrid, β, with_density=true)
    indices = (-last(axes(GImp,1))-1):last(axes(GImp,1))
    GImp = OffsetVector(vcat(conj.(reverse(GImp.parent)), GImp.parent), indices)
    νnGrid =  jED.OffsetVector([1im * (2*n+1)*π/β for n in -nFreq:nFreq-1], -nFreq:nFreq-1)
    return U, β, p, νnGrid, GImp, μ, nden
end


# ------------------ Bare Bubble -------------------
"""
    compute_χ0(ω_range::AbstractArray{Int,1}, ν_range::AbstractArray{Int,1}, G::Array{Complex{Float64}, 1}, β::Float64; mode=:ph)

Computes bare susceptibility from Green's function `G` for particle-particle or particle-hole channel.
Warning! `G` must fulfill ``G^{-\\nu} = conj(G^{\\nu})``
"""
function compute_χ0(ω_range::AbstractArray{Int,1}, ν_range::AbstractArray{Int,1}, G::Array{Complex{Float64}, 1}, β::Float64; mode=:ph)
    !(mode in [:ph, :pp]) && error("unkown mode")
    χ0 = OffsetArray(Array{ComplexF64,2}(undef, length(ω_range), length(ν_range))) 
    for ωm in ω_range, νn in νVector_range
        χ0[ωm,νn] = (mode == :ph) ? -β*get_conjsymm_f(G, νn)*get_conjsymm_f(G, νn+ωm) : -β*get_conjsymm_f(G, νn)*get_conjsymm_f(G, ωm-νn-1)
    end
    return χ0
end

function compute_χ0(ω_range::AbstractArray{Int,1}, ν_range::AbstractArray{Int,1}, G::OffsetArray{Complex{Float64}, 1}, β::Float64; mode=:ph)
    !(mode in [:ph, :pp]) && error("unkown mode")
    χ0 = OffsetArray(Array{ComplexF64,2}(undef, length(ω_range), length(ν_range)), ω_range, ν_range) 
    for ωm in ω_range, νn in ν_range
        χ0[ωm,νn] = (mode == :ph) ? -β*G[νn]*G[νn+ωm] : -β*G[νn]*G[ωm-νn-1]
    end
    return χ0
end

# -------------- Irreducible Vertex ----------------
# ------------------ Full Vertex -------------------
"""
    computeF_ph(freqList::Vector, χm::Vector{T}, χd::Vector{T}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}, nBose::Int, nFermi::Int) where T

Computes `Fm` and `Fd` from `χm` and `χd`
"""
function computeF_ph(freqList::Vector, χm::Vector{T}, χd::Vector{T}, χ0::OffsetMatrix, nBose::Int, nFermi::Int) where T
    Fm = similar(χm)
    Fd = similar(χd)
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i]
        sub = ν == νp ? χ0[ω,ν] : 0.0
        Fm[i] = (-1.0/χ0[ω,ν])*(χm[i] - sub)*(1.0/χ0[ω,νp])
        Fd[i] = (-1.0/χ0[ω,ν])*(χd[i] - sub)*(1.0/χ0[ω,νp])
    end
    return reshape(Fm, 2*nFermi, 2*nFermi, 2*nBose+1), reshape(Fd, 2*nFermi, 2*nFermi, 2*nBose+1)
end

function computeF_pp(freqList::Vector, χ_s::Vector{T}, χ_t::Vector{T}, χ0::Dict{Tuple{Int,Int},Complex{Float64}}, nBose::Int, nFermi::Int) where T
    F_s = similar(χ_s)
    F_t = similar(χ_t)
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i]
        sub = ν == νp ? χ0[(ω,ν)] : 0.0
        F_s[i] = (-1.0/χ0[(ω,ν)])*(χ_s[i]+2*sub)*(1.0/χ0[(ω,νp)])
        F_t[i] = (-1.0/χ0[(ω,ν)])*(χ_t[i]-2*sub)*(1.0/χ0[(ω,νp)])
    end
    return reshape(F_s, 2*nFermi, 2*nFermi, 2*nBose+1), reshape(F_t, 2*nFermi, 2*nFermi, 2*nBose+1)
end

"""
    F_to_χ(freqList::Array, arr::Array{T}, gImp::Array{Complex{Float64}, 1}, β::Float64) where T <: Number
    F_to_χ!(res, freqList::Array, arr::Array{T}, gImp::Array{Complex{Float64}, 1}, β::Float64) where T <: Number

Subtracts disconnected term ``\\delta_{\\omega 0} G^\\nu G^{\\nu'}`` from ``F^{\\nu\\nu'\\omega}``
"""
function F_to_χ(freqList::Array, in::Array, gImp::OffsetArray{Complex{Float64}, 1}, β::Float64)
    res = similar(in)
    F_to_χ!(res, freqList, in, gImp, β)
    return res
end

function F_to_χ!(res::Array, freqList::Array, in::Array, gImp::OffsetArray{Complex{Float64}, 1}, β::Float64)
    @assert all(size(res) .== size(in))
    res[:] .= in[:]
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i]
        if ω == 0
            res[i] = in[i] - β*gImp[ν]*gImp[νp]
        end
    end
end


function computeχ0(ω_range::AbstractArray{Int,1}, ν_range::AbstractArray{Int,1}, gImp::OffsetArray, β::Float64; mode=:ph)
    !(mode in [:ph, :pp]) && error("unkown mode")
    χ0 = Dict{Tuple{Int,Int},Complex{Float64}}()
    for ω in ω_range, ν in ν_range
        χ0[(ω,ν)] = (mode == :ph) ? -β*gImp[ν]*gImp[ν+ω] : -β*gImp[ν]*gImp[ω-ν-1]
    end
    return χ0
end
function computeχ0(ω_range::AbstractArray{Int,1}, ν_range::AbstractArray{Int,1}, gImp::Array{Complex{Float64}, 1}, β::Float64; mode=:ph)
    !(mode in [:ph, :pp]) && error("unkown mode")
    χ0 = Dict{Tuple{Int,Int},Complex{Float64}}()
    for ω in ω_range, ν in ν_range
        χ0[(ω,ν)] = (mode == :ph) ? -β*get_conjsymm_f(gImp, ν)*get_conjsymm_f(gImp, ν+ω) : -β*get_conjsymm_f(gImp, ν)*get_conjsymm_f(gImp, ω-ν-1)
    end
    return χ0
end


# -------------- Irreducible Vertex ----------------
"""
    computeΓ_ph(freqList::Array{Tuple,3}, χ_r::Array{ComplexF64,3}, χ0::OffsetMatrix{ComplexF64}, nBose::Int64, nFermi::Int64)
    computeΓ_ph(freqList::Array{Tuple,1}, χ_r::Array{ComplexF64,1}, χ0::OffsetMatrix{ComplexF64}, nBose::Int64, nFermi::Int64)

freqList, χ_r are rank 3 Arrays of size [2*nBose + 1, 2*nFermi, 2*nFermi].
Linear arrays should be reshaped using [`reshape_lin_to_rank3`](@ref reshape_lin_to_rank3).
If rank `1` shapes are provided, the conversion will be done internally.
feeqList contains 3-Tuples with (ωn,νn,νpn), χ_r the corresponsing values of the susceptibilities at this frequencies.
Returns the irreducible vertex `Γ` in this channel.
"""
function computeΓ_ph(freqList::Array{Tuple{Int,Int,Int},3}, χ_r::Array{ComplexF64,3}, χ0::OffsetMatrix{ComplexF64}, nBose::Int64, nFermi::Int64)
    Γ_r = similar(χ_r)
    for ωi in axes(freqList,1)
        Γ_r[ωi,:,:] = inv(χ_r[ωi,:,:])
        for νi in axes(freqList, 2)
            for νpi in axes(freqList, 3)
                ωm, νn,νpn = freqList[ωi,νi,νpi] 
                if νn == νpn
                    Γ_r[ωi,νi,νi] -= 1.0/χ0[ωm,νn]
                end
            end
        end
    end
    return Γ_r
end

function computeΓ_ph(freqList::Array{Tuple{Int,Int,Int},1}, χ_r::Array{ComplexF64,1}, χ0::OffsetMatrix{ComplexF64}, nBose::Int64, nFermi::Int64)
    return computeΓ_ph( reshape_lin_to_rank3(freqList,nBose,nFermi), reshape_lin_to_rank3(χ_r,nBose,nFermi), χ0, nBose, nFermi)
end

function computeΓ_pp(freqList::Array, χs::Array{T,3}, χt::Array{T,3}, χ0::OffsetMatrix{ComplexF64}, shift, nBose::Int64, nFermi::Int64) where T
    Γs = fill!(Array{ComplexF64,3}(undef, 2*nFermi, 2*nFermi, 2*nBose+1), NaN)
    Γt = fill!(Array{ComplexF64,3}(undef, 2*nFermi, 2*nFermi, 2*nBose+1), NaN)
    # for (ωn,ω) in enumerate([0])
    
    for (ωi, ωn) in enumerate(-nBose:nBose)
        non_nan_slice = find_non_nan_matrix(χs[:,:,ωi], nFermi)
        if !isempty(non_nan_slice)
            # if any(isnan.(χs[non_nan_slice,non_nan_slice,ωi]))
            #     println("ERROR at $ωi/$ωn, NaN cut did not work. slice: $non_nan_slice")
            # end
            Γs[non_nan_slice,non_nan_slice,ωi] = 4 .* inv(χs[non_nan_slice,non_nan_slice,ωi])
            Γt[non_nan_slice,non_nan_slice,ωi] = 4 .* inv(χt[non_nan_slice,non_nan_slice,ωi])
            for (νi,νn) in enumerate((-nFermi:nFermi-1) .- shift*trunc(Int,ωn/2))
                Γs[νi,νi,ωi] += 2/χ0[ωn,νn]
                Γt[νi,νi,ωi] -= 2/χ0[ωn,νn]
            end
        end
    end
    return Γs, Γt
end
