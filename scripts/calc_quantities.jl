# ==================================================================================================== #
#                                       calc_quantities.jl                                             #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#   Internal script for computation of quantities derived from the 2-particle GF.                      #
#   See `expand_vertex.jl`                                                                             #
# ==================================================================================================== #

println("Calculating single particle Green's function")
U, β, p, νnGrid, G0W, GImp, ΣImp, μ, nden = restore_1pt_GF(joinpath(dataDir, "config.toml"), joinpath(dataDir, "hubb.andpar"); nFreq=2000)
gLoc = nothing
println("mu = $μ, U= $U, beta = $β, nden = $nden")
println("andpar: \n $p")
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
Fm = F_from_χ(χm, GImp, shift, nBose, nFermi,β, :m)
Fd = F_from_χ(χd, GImp, shift, nBose, nFermi,β, :d)
Γm = computeΓ_ph(freqList, χm, χ0_full, nBose, nFermi)
Γd = computeΓ_ph(freqList, χd, χ0_full, nBose, nFermi)
println("Done with ph channel!")
res = isfile(joinpath(dataDir, "chi_asympt")) ? read_chi_asympt(joinpath(dataDir, "chi_asympt")) : error("chi_asympt not found!")
_, χ_d_asympt, χ_m_asympt, χ_pp_asympt = res


χ0_pp_full   = compute_χ0(-nBose:nBose, -(nFermi+2*nBose):(nFermi+2*nBose)-1, GImp, β; mode=:pp)
χpp_s, χpp_t = χph_to_χpp(freqList, χ_upup, χ_updo, χ0_pp_full, shift, nBose, nFermi)
Fs, Ft       = computeF_pp(freqList, χpp_s, χpp_t, χ0_pp_full, nBose, nFermi)
Γs, Γt       = computeΓ_pp(freqList, χpp_s, χpp_t, χ0_pp_full, nBose, nFermi)
Φs = Fs .- Γs
Φt = Ft .- Γt
println("Done with pp channel!")

ΣLoc_m, ΣLoc_d = calc_local_EoM(Fm, Fd, GImp, U, β, nden, nBose, nFermi, shift)
loc_EoM = 0.5 .* (ΣLoc_m .+ ΣLoc_d)
loc_EoM_diff_abs = abs.(loc_EoM .- ΣImp[axes(ΣLoc_m,1)]) 
loc_EoM_diff_rel = round.(100 * loc_EoM_diff_abs ./ (abs.(loc_EoM) .+ abs.(ΣImp[axes(ΣLoc_m,1)])), digits=4)
fitCheck = if !isnothing(gLoc)
    abs.(gLoc .- G0W)
else
    VertexPostprocessing.OffsetVector(repeat([NaN], 6), 0:5)
end
sum_vk, min_eps_diff, min_vk, min_eps = andpar_check_values(p.ϵₖ, p.Vₖ)

println(repeat("=",80))
println(" ϵₖ = $(lpad.(round.(p.ϵₖ, digits=4),9)...)")
println(" Vₖ = $(lpad.(round.(p.Vₖ, digits=4),9)...)")
println("   1. min(|Vₖ|)           = ", min_vk)
println("   2. ∑Vₗ^2               = ", sum_vk)
println("   3. min(|ϵₖ|)           = ", min_eps)
println("   4. min(|ϵₖ - ϵₗ|)      = ", min_eps_diff)
println("   5. |GLoc - G0W|[ν]     = ", fitCheck[0:5])
println("   6. |ΣImp - Σ_2part|[ν] = ", loc_EoM_diff_abs[0:5])
println("   7. 6 (relative, %)     = ", loc_EoM_diff_rel[0:5])
println(repeat("=",80))
