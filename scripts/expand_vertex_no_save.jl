# ==================================================================================================== #
#                                       calc_quantities.jl                                             #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#   Internal script for expansion of 2-particle Green's function.                                      #
#   See `expand_vertex.jl` and `combine_vertices.jl`                                                   #
# ==================================================================================================== #

println("Expanding Vertex")
fname = legacy_mode ? "vert_chi" : "2_part_gf_red"
nBose, nFermi, shift, freqList, TwoPartGF_upup, TwoPartGF_updo = expand_2PtGF_CSV(freqListFile, joinpath(dataDir, fname))
println("Done expanding!")
