# ==================================================================================================== #
#                                        expand_vertex.jl                                              #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#   Produces an input file for LadderDAG.jl.                                                           #
#   This script takes 2+1 argumentes:                                                                  #
#   ARGS[1]: Path tp frequency list file, produced by EquivalenceClassesCosntructor.jl                 #
#   ARGS[2]: dataPath. It is epxected that this directory contains 3 files:                            #
#     - 2_part_gf_red (two-particle GF in frequency notation according to G.Rohringer PhD Thesis)      #
#     - config.toml only field `U` and `beta` in the section `[parameters]` need to be set.            #
#     - hubb.andpar: contains                                                                          #
#         - bath levels starting from line 10                                                          #
#         - hoppings starting from line 10+NBathsites+1 (one comment line above between levels and hoppings)
#         - chemical potential mu in the last line                                                     #
#   ARGS[3]: legacy_mode (optional). This reads a vert_chi instead of a 2_part_gf_red                  #
# -------------------------------------------- TODO -------------------------------------------------- #
# ==================================================================================================== #

using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using VertexPostprocessing
using JLD2

if length(ARGS) < 3
    println("ERROR: expected 3 (+1) command line arguments. Please provide\n
                (1) path to directory of the frequency file generated from EquivalencyClassesConstructor.jl\n
                (2) path to DataDir, containing `config.toml`, `hubb.andpar` and `chi_asympt`.\n
                (3) Bool for legacy_mode (if True, vert_chi will be read, otherwise 2_part_gf_red). In legacy mode the G*G contribution has already been subtracted.\n
                (4) (optional) filename of the output file. Default is DMFT2_out.jld2")
    exit(1)
end

freqListFile = joinpath(ARGS[1],"freqList.jld2")
dataDir      = ARGS[2]
legacy_mode  = parse(Bool, ARGS[3])
ofname = length(ARGS) > 3 ? ARGS[4] : "DMFT2_out.jld2"

include(joinpath(@__DIR__,"expand_vertex_no_save.jl"))
include(joinpath(@__DIR__,"calc_quantities.jl"))

println("Storing results in $ofname")
# #-1.0 .*  Γr
jldopen(joinpath(dataDir, ofname), "w") do f
    f["Γch"] = -Γd      # Convention for lDGA is Γ → - Γ
    f["Γsp"] = -Γm
    f["Γs"] = Γs
    f["Γt"] = Γt
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

