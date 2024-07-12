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
`f_d = f_upup + f_updo`
`f_m = f_upup - f_updo`
"""
function uu_ud_TO_m_d(val_upup, val_updo)
    val_d = val_upup .+ val_updo
    val_m = val_upup .- val_updo
    return val_d, val_m
end
