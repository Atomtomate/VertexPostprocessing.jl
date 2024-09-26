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

using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using VertexPostprocessing
using JLD2

if length(ARGS) != 3
    println("ERROR: expected 3 command line arguments. Please provide 
                (1) path to frequency file generated from EquivalencyClassesConstructor.jl 
                (2) path to DataDir, containing `config.toml`, `hubb.andpar` and `chi_asympt`.
                (3) Bool for legacy_mode (if True, vert_chi will be read, otherwise 2_part_gf_red). In legacy mode the G*G contribution has already been subtracted.")
    exit(1)
end

freqListFile = ARGS[1]
dataDir      = ARGS[2]
legacy_mode  = parse(Bool, ARGS[3])



println("Expanding Vertex")
fname = legacy_mode ? "vert_chi" : "2_part_gf_red"
nBose, nFermi, shift, freqList, TwoPartGF_upup, TwoPartGF_updo = expand_2PtGF_CSV(freqListFile, joinpath(dataDir, fname))
println("Done expanding!")

println("Calculating single particle Green's function")
U, β, p, νnGrid, G0W, GImp, μ, nden = restore_1pt_GF(joinpath(dataDir, "config.toml"), joinpath(dataDir, "hubb.andpar"); nFreq=2000)
E_kin_DMFT = calc_EKin_DMFT(νnGrid[0:end], p.ϵₖ, p.Vₖ, GImp[0:end], nden, U, β, μ)
E_pot_DMFT = calc_EPot_DMFT(νnGrid[0:end], p.ϵₖ, p.Vₖ, GImp[0:end], nden, U, β, μ)

println("Calculating derived quantities. Rank 3 quantities have the index convention [ωm, νn, νpn]")
χ_upup, χ_updo = if legacy_mode
    TwoPartGF_upup, TwoPartGF_updo
else
    G2_to_χ(freqList, TwoPartGF_upup, GImp, β), G2_to_χ(freqList, TwoPartGF_updo, GImp, β)
end
χm, χd = uu_ud_TO_m_d(χ_upup, χ_updo)

χ0_full = compute_χ0(-nBose:nBose, -(nFermi+2*nBose):(nFermi+2*nBose)-1, GImp, β)
Fm, Fd = computeF_ph(freqList, χm, χd, χ0_full, nBose, nFermi)
Γm = computeΓ_ph(freqList, χm, χ0_full, nBose, nFermi)
Γd = computeΓ_ph(freqList, χd, χ0_full, nBose, nFermi)
println("Done with ph channel!")
res = isfile(joinpath(dataDir, "chi_asympt")) ? read_chi_asympt(joinpath(dataDir, "chi_asympt")) : error("chi_asympt not found!")
_, χ_d_asympt, χ_m_asympt, χ_pp_asympt = res


# χ0_pp_full   = compute_χ0(-nBose:nBose, -(nFermi+2*nBose):(nFermi+2*nBose)-1, GImp, β; mode=:pp)
# χpp_s, χpp_t = χph_to_χpp(freqList, χ_upup, χ_updo, χ0_pp_full, shift, nBose, nFermi)
# Fs, Ft       = computeF_pp(freqList, χpp_s, χpp_t, χ0_pp_full, nBose, nFermi)
# Γs, Γt       = computeΓ_pp(freqList, χpp_s, χpp_t, χ0_pp_full, nBose, nFermi)
# Φs = Fs .- Γs
# Φt = Ft .- Γt
println("Done with pp channel!")
println("Storing results in DMFT_out.jld2")
# #-1.0 .*  Γr
jldopen(joinpath(dataDir,"DMFT2_out.jld2"), "w") do f
    f["Γch"] = -Γd      # Convention for lDGA is Γ → - Γ
    f["Γsp"] = -Γm
    # f["Φpp_s"] = Φs
    # f["Φpp_t"] = Φt
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
