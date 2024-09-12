# ==================================================================================================== #
#                                      ChannelTransforms.jl                                            #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#   Helper functions for transformations between deifferent spin and frequency  notations.             #
# -------------------------------------------- TODO -------------------------------------------------- #
# ==================================================================================================== #

"""
    uu_ud_TO_m_d(val_upup, val_updo)

Takes to objects with spin-indices `upup` and `updo` and returns
them in `m` and `d` notation, i.e. 
`f_d = (f_upup + f_updo)/1`
`f_m = (f_upup - f_updo)/1`
"""
function uu_ud_TO_m_d(val_uu, val_ud)
    val_m = val_uu .- val_ud
    val_d = val_uu .+ val_ud
    return val_m, val_d 
end

"""
    m_d_TO_uu_ud(val_upup, val_updo)

Takes to objects with spin-indices `upup` and `updo` and returns
them in `m` and `d` notation, i.e. 
`f_uu = (f_d + f_m)/2`
`f_ud = (f_d - f_m)/2`
"""
function m_d_TO_uu_ud(val_m, val_d)
    val_uu = (val_d .+ val_m)/2
    val_ud = (val_d .- val_m)/2
    return val_uu, val_ud
end

"""
    χph_to_χpp(freqList::Array, χph_upup::Array, χph_updo::Array, χ0::Dict{Tuple{Int,Int},Complex{Float64}}, shift, nBose::Int, nFermi::Int)

transforms susceptibilities in particle-hole notation for the `↑↑` and `↑↓` spin channels into particle-particle
notation in the singlet and triplet channels.
"""
function χph_to_χpp(freqList::Vector{Tuple{Int,Int,Int}}, χph_upup::Vector, χph_updo::Vector, χ0::OffsetMatrix, shift, nBose::Int, nFermi::Int)
    χpp_s = similar(χph_upup)
    χpp_t = similar(χph_upup)
    χpp_s = fill!(χpp_s, NaN)
    χpp_t = fill!(χpp_t, NaN)
    χph_upup_tmp = reshape_lin_to_rank3(χph_upup,nBose,nFermi)
    χph_updo_tmp = reshape_lin_to_rank3(χph_updo,nBose,nFermi)

    for i in 1:size(freqList,1)
        ωn, νn, νpn = freqList[i]
        ωi,νi,νpi = Freq_to_OneToIndex(ωn, νn, νpn, shift, nBose, nFermi)
        ωi_ph,νi_ph,νpi_ph = Freq_to_OneToIndex(ωn - νn - νpn - 1, νn, νpn, shift, nBose, nFermi)
        if !(any((ωi,νi,νpi) .< 1) || any((ωi_ph,νi_ph,νpi_ph) .< 1) || any((ωi,νi,νpi) .> (2*nBose,2*nFermi-1,2*nFermi-1)) || any((ωi_ph,νi_ph,νpi_ph) .> (2*nBose,2*nFermi-1,2*nFermi-1)))
            χpp_s[i] = - χ0[ωn,νn]*(νn==νpn) - χph_upup_tmp[ωi_ph,νi_ph,νpi_ph] + 2*χph_updo_tmp[ωi_ph,νi_ph,νpi_ph]
            χpp_t[i] = + χ0[ωn,νn]*(νn==νpn) + χph_upup_tmp[ωi_ph,νi_ph,νpi_ph]
        end
    end
    χpp_s, χpp_t
end
