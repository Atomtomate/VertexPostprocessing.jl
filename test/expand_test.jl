datap_fs0 = joinpath(@__DIR__,"test_data/data_s0_full")
nBose_fs0, nFermi_fs0, shift_fs0, freqList_fs0, TwoPartGF_upup_fs0, TwoPartGF_updo_fs0 = VertexPostprocessing.expand_2PtGF_CSV(joinpath(datap_fs0,"grid_full_b5_f5_s0.jld2"),joinpath(datap_fs0,"vert_chi"))
U, β, p, νnGrid, GImp, μ, nden = restore_1pt_GF(joinpath(datap_fs0,"config.toml"), joinpath(datap_fs0, "hubb.andpar"); nFreq=6*(nFermi_fs0+nBose_fs0))

datap_s0 = joinpath(@__DIR__,"test_data/data_s0")
nBose_s0, nFermi_s0, shift_s0, freqList_s0, TwoPartGF_upup_s0, TwoPartGF_updo_s0 = VertexPostprocessing.expand_2PtGF_CSV(joinpath(datap_s0,"grid_b5_f5_s0.jld2"),joinpath(datap_s0,"vert_chi"))

datap_s1 = joinpath(@__DIR__,"test_data/data_s1")
nBose_s1, nFermi_s1, shift_s1, freqList_s1, TwoPartGF_upup_s1, TwoPartGF_updo_s1 = VertexPostprocessing.expand_2PtGF_CSV(joinpath(datap_s1,"grid_b5_f5_s1.jld2"),joinpath(datap_s1,"vert_chi"))

# test for less than 0.1% discrepancy between exact and reduced
diff_F = abs.(TwoPartGF_updo_s0 .- TwoPartGF_updo_fs0) ./ abs.(TwoPartGF_updo_s0 + TwoPartGF_updo_fs0)
@test all(diff_F .< 0.001)
χ0_ph = compute_χ0(-4*nBose_fs0:4*nBose_fs0, -4*nFermi_fs0:4*nFermi_fs0, GImp, β; mode=:ph)
χ0_ph_test = [-β*GImp[0]*GImp[0], -β*GImp[1]*GImp[1], -β*GImp[-1]*GImp[0]]
@test all([χ0_ph[0,0], χ0_ph[0,1], χ0_ph[1,-1]] .≈ χ0_ph_test)

χ_upup_fs0 = F_to_χ(freqList_fs0, TwoPartGF_upup_fs0, GImp, β)
χ_updo_fs0 = F_to_χ(freqList_fs0, TwoPartGF_updo_fs0, GImp, β)
χ_d_fs0, χ_m_fs0 = uu_ud_TO_m_d(χ_upup_fs0, χ_updo_fs0) 
Γd_ph_fs0 = computeΓ_ph(freqList_fs0, χ_d_fs0, χ0_ph, nBose_fs0, nFermi_fs0);
Γm_ph_fs0 = computeΓ_ph(freqList_fs0, χ_m_fs0, χ0_ph, nBose_fs0, nFermi_fs0);

χ_upup_s0 = F_to_χ(freqList_s0, TwoPartGF_upup_s0, GImp, β)
χ_updo_s0 = F_to_χ(freqList_s0, TwoPartGF_updo_s0, GImp, β)
χ_d_s0, χ_m_s0 = uu_ud_TO_m_d(χ_upup_s0, χ_updo_s0) 
Γd_ph_s0 = computeΓ_ph(freqList_s0, χ_d_s0, χ0_ph, nBose_s0, nFermi_s0);
Γm_ph_s0 = computeΓ_ph(freqList_s0, χ_m_s0, χ0_ph, nBose_s0, nFermi_s0);

χ_upup_s1 = F_to_χ(freqList_s1, TwoPartGF_upup_s1, GImp, β)
χ_updo_s1 = F_to_χ(freqList_s1, TwoPartGF_updo_s1, GImp, β)
χ_d_s1, χ_m_s1 = uu_ud_TO_m_d(χ_upup_s1, χ_updo_s1) 
Γd_ph_s1 = computeΓ_ph(freqList_s1, χ_d_s1, χ0_ph, nBose_s1, nFermi_s1);
Γm_ph_s1 = computeΓ_ph(freqList_s1, χ_m_s1, χ0_ph, nBose_s1, nFermi_s1);

diff_chi = abs.((χ_d_s0 .- χ_d_fs0) ./ abs.(χ_d_s0 .+ χ_d_fs0))
@test all(diff_chi .< 0.001)
