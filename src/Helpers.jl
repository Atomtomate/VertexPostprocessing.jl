# ==================================================================================================== #
#                                            Helpers.jl                                                #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#  Generic helper functions                                                                            #
# -------------------------------------------- TODO -------------------------------------------------- #
# ==================================================================================================== #

"""
    reshape_lin_to_rank3(arr::Vector, nBose::Int, nFermi::Int)

Reshapes linear list to rank 3 with one bosonic and two fermionic Matsubara frequencies.
This uses a convention consistent through a set of Fortran codes.
"""
function reshape_lin_to_rank3(arr::Vector, nBose::Int, nFermi::Int)
    permutedims(reshape(arr, 2*nFermi, 2*nFermi, 2*nBose+1),[3,2,1])
end


"""
    Freq_to_OneToIndex(ωn::Int, νn::Int, νpn::Int, shift::Union{Bool,Int}, nBose::Int, nFermi::Int)

Transforms triple of Matsubara indices (one bosonic, two fermionic ones, see Georg Rohringer thesis for details)
into array indices starting from one. 
This also respects shifted grids, if `shift = True`.
`nFermi` and `nBose` are the maximum number of positive fermionic and bosonic indices.
"""
function Freq_to_OneToIndex(ωn::Int, νn::Int, νpn::Int, shift::Union{Bool,Int}, nBose::Int, nFermi::Int)
    ωn+nBose+1,νn+nFermi+1+trunc(Int, shift*ωn/2), νpn+nFermi+1+trunc(Int, shift*ωn/2)
end

function OneToIndex_to_Freq(ωi::Int, νi::Int, shift::Union{Bool,Int}, nBose::Int, nFermi::Int)
    ωn = ωi - nBose - 1
    νn = (νi - nFermi - 1) - shift * trunc(Int, ωn / 2)
    return ωn, νn
end

function νnGrid(ωm::Int, shift::Union{Bool,Int}, nFermi::Int)
    (-nFermi:(nFermi-1)) .- shift * trunc(Int, ωm/2)
end

"""
    find_non_nan_matrix(data::Matrix, nFermi::Int)

Indices for largest slice of matrix not containing any NaNs.
"""
function find_non_nan_matrix(data::Matrix, nFermi::Int)
    nan_list       = sort(map(x->x[1],filter(x->x[1] == x[2], findall(x-> !isnan(x), data))))
    nh = searchsorted(nan_list, nFermi)
    res = if length(nh) == 1
        t1 = diff(nan_list[nh[1]:end])
        t2 = -1 * diff(nan_list[nh[1]:-1:1])
        lim_up = findfirst(x->x != 1, t1)
        lim_lo = findfirst(x->x != 1, t2)
        lim_up = isnothing(lim_up) ? length(t1) : lim_up
        lim_lo = isnothing(lim_lo) ? length(t2) : lim_lo
        lim_non_nan = min(lim_up, lim_lo)
        (nFermi-lim_non_nan+1):(nFermi+lim_non_nan)
    else
        Int[]
    end
    return res
end
