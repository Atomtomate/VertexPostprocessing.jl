# ==================================================================================================== #
#                                        expand_vertex.jl                                              #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#   Produces an input file for LadderDAG.jl.                                                           #
#   This script takes 2 argumentes:                                                                    #
#   ARGS[1]: Path tp frequency list file, produced by EquivalenceClassesCosntructor.jl                 #
#   ARGS[2]: dataPath. It is epxected that this directory contains 3 files:                            #
#     - 2_part_gf_red (two-particle GF in frequency notation according to G.Rohringer PhD Thesis)      #
#     - config.toml only field `U` and `beta` in the section `[parameters]` need to be set.            #
#     - hubb.andpar: contains                                                                          #
#         - bath levels starting from line 10                                                          #
#         - hoppings starting from line 10+NBathsites+1 (one comment line above between levels and hoppings)
#         - chemical potential mu in the last line                                                     #
# -------------------------------------------- TODO -------------------------------------------------- #
# ==================================================================================================== #

using VertexPostprocessing

if length(ARGS) != 2
    println("ERROR: expected 2 command line arguments. The comments in this script contain instructions.")
    exit(1)
end

freqListFile = ARGS[1]
dataDir      = ARGS[2]


println("Expanding Vertex")
nBose, nFermi, shift, freqList, TwoPartGF_upup, TwoPartGF_updo = expand_2PtGF_CSV(freqListFile, joinpath(datap_s0,"vert_chi"))
println("Done expanding!")

println("Calculating single particle Green's function")
U, β, p, νnGrid, GImp, μ, nden = restore_1pt_GF(joinpath(dataDir, "config.toml"), joinpath(dataDir, "hubb.andpar"); nFreq=2000)
E_kin_DMFT = calc_EKin_DMFT(νnGrid[0:end], p.ϵₖ, p.Vₖ, GImp[0:end], nden, U, β, μ)
E_pot_DMFT = calc_EPot_DMFT(νnGrid[0:end], p.ϵₖ, p.Vₖ, GImp[0:end], nden, U, β, μ)

println("Calculating derived quantities")
χ_upup = F_to_χ(freqList, TwoPartGF_upup, GImp, β)
χ_updo = F_to_χ(freqList, TwoPartGF_updo, GImp, β)
χm, χd = uu_ud_TO_m_d(χ_upup, χ_updo)

χ0_full = compute_χ0(-nBose:nBose, -(nFermi+2*nBose):(nFermi+2*nBose)-1, GImp, β)
Fm, Fd = computeF_ph(freqList, χm, χd, χ0_full, nBose, nFermi)
Γm = -1.0 .* computeΓ_ph(freqList, χm, χ0_full, nBose, nFermi)
Γd = -1.0 .* computeΓ_ph(freqList, χd, χ0_full, nBose, nFermi)
println("Done with ph channel!")
res = isfile(joinpath(dataDir, "chi_asympt")) ? read_chi_asympt(joinpath(dataDir, "chi_asympt")) : error("chi_asympt not found!")
χ_d_asympt, χ_m_asympt, χ_pp_asympt = res

