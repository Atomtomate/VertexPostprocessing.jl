# ==================================================================================================== #
#                                   interpolate_sparse_vertex.jl                                       #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe, Marvin Leusch                                                     #
# ----------------------------------------- Description ---------------------------------------------- #
#   Produces an input file for LadderDAG.jl by interpolating a two particle Green's functions          #
#   that was calculated on a sparse (i.e. not all Matsubara frequencies) grid.                         #
#   ARGS[1]: Path tp frequency list file, produced by EquivalenceClassesCosntructor.jl                 #
#   ARGS[2]: dataPath. It is epxected that this directory contains 3 files:                            #
#     - 2_part_gf_red (two-particle GF in frequency notation according to G.Rohringer PhD Thesis)      #
#     - config.toml only field `U` and `beta` in the section `[parameters]` need to be set.            #
#     - hubb.andpar: contains                                                                          #
#         - bath levels starting from line 10                                                          #
#         - hoppings starting from line 10+NBathsites+1 (one comment line above between levels and hoppings)
#         - chemical potential mu in the last line                                                     #
#   ARGS[3]: Output path+name                                                                          #
# -------------------------------------------- TODO -------------------------------------------------- #
# ==================================================================================================== #

using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using VertexPostprocessing
using JLD2

freqListFile = joinpath(ARGS[1],"freqList.jld2")
dataDir      = ARGS[2]
outFileName  = ARGS[3]
legacy_mode  = false
fname = "2_part_gf_red"

println("Expanding Vertex 1")
nBose, nFermi, shift, freqList, TwoPartGF_upup, TwoPartGF_updo = expand_2PtGF_CSV(freqListFile, joinpath(dataDir, fname))
println("Done expanding Vertex 1!")


include(joinpath(@__DIR__,"calc_quantities.jl"))

println("Storing results in $outFileName")

if false
    jldopen(outFileName, "w") do f
        f["Γch"] = -Γd      # Convention for lDGA is Γ → - Γ
        f["Γsp"] = -Γm
        f["Φpp_s"] = Φs
        f["Φpp_t"] = Φt
        f["χDMFTch"] = reshape_lin_to_rank3(χd,nBose,nFermi)
        f["χDMFTsp"] = reshape_lin_to_rank3(χm,nBose,nFermi)
        f["χ_ch_asympt"] = χ_d_asympt
        f["χ_sp_asympt"] = χ_m_asympt
        f["χ_pp_asympt"] = χ_pp_asympt
        f["gImp"] = GImp.parent
        f["g0"] = G0W.parent
        f["ϵₖ"] = p.ϵₖ
        f["Vₖ"] = p.Vₖ
        f["μ"] = μ
        f["U"] = U
        f["β"] = β
        f["nden"] = nden
        f["E_kin_DMFT"] = E_kin_DMFT
        f["E_pot_DMFT"] = E_pot_DMFT
        f["grid_shift"] = shift
        f["grid_nBose"] = nBose
        f["grid_nFermi"] = nFermi
    end;
end;