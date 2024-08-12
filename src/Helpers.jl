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

