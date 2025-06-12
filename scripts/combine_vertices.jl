# ==================================================================================================== #
#                                       combine_TwoPartGF.jl                                           #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#   Produces an input file for LadderDAG.jl by combining two two particle Green's functions..          #
#   WARNING: For now the resulting vertex has just the frequency list of the larger of the two!        #
#   This script takes 5 argumentes:                                                                    #
#   ARGS[1]: Path tp frequency list file, produced by EquivalenceClassesCosntructor.jl                 #
#   ARGS[2]: dataPath. It is epxected that this directory contains 3 files:                            #
#     - 2_part_gf_red (two-particle GF in frequency notation according to G.Rohringer PhD Thesis)      #
#     - config.toml only field `U` and `beta` in the section `[parameters]` need to be set.            #
#     - hubb.andpar: contains                                                                          #
#         - bath levels starting from line 10                                                          #
#         - hoppings starting from line 10+NBathsites+1 (one comment line above between levels and hoppings)
#         - chemical potential mu in the last line                                                     #
#   ARGS[3]: The same as ARGS[1], but for additional vertex                                            #
#   ARGS[4]: The same as ARGS[2], for for additional vertex                                            #
#   ARGS[5]: Output path+name                                                                          #
# -------------------------------------------- TODO -------------------------------------------------- #
#   Add capability of combining distinct (non-overlapping) frequency lists!                            #
# ==================================================================================================== #

using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using VertexPostprocessing
using JLD2

freqListFile = joinpath(ARGS[1],"freqList.jld2")
dataDir      = ARGS[2]
legacy_mode  = false
fname = "2_part_gf_red"

println("Expanding Vertex 1")
nBose_1, nFermi_1, shift_1, freqList_1, TwoPartGF_upup_1, TwoPartGF_updo_1 = expand_2PtGF_CSV(freqListFile, joinpath(dataDir, fname))
println("Done expanding Vertex 1!")

freqListFile = joinpath(ARGS[3],"freqList.jld2")
dataDir      = ARGS[4]
outFileName  = ARGS[5]

println("Expanding Vertex 2")
nBose_2, nFermi_2, shift_2, freqList_2, TwoPartGF_upup_2, TwoPartGF_updo_2 = expand_2PtGF_CSV(freqListFile, joinpath(dataDir, fname))
println("Done expanding Vertex 2!")


TwoPartGF_upup, TwoPartGF_updo = combine_TwoPartGF(
        freqList_1, freqList_2, 
        TwoPartGF_upup_1, TwoPartGF_upup_2,
        TwoPartGF_updo_1, TwoPartGF_updo_2
    );
println("Done combining!")
prio_on_1 = length(freqList_1) > length(freqList_2)
nBose = prio_on_1 ? nBose_1 : nBose_2 
nFermi = prio_on_1 ? nFermi_1 : nFermi_2 
shift = prio_on_1 ? shift_1 : shift_2 
freqList = prio_on_1 ? freqList_1 : freqList_2 

include(joinpath(@__DIR__,"calc_quantities.jl"))

println("Storing results in $outFileName")
# #-1.0 .*  Γr
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

