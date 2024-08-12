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
    val_d = val_uu .+ val_ud
    val_m = val_uu .- val_ud
    return val_d, val_m
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
