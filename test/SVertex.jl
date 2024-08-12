function expand_test(TwoPartGF_upup, TwoPartGF_updo, freqList_map, freqList, parents, ops, nBose, nFermi)
    off(f) = f[1]+nBose+1,f[2]+nFermi+1,f[3]+nFermi+1
    Fup_full = -1 .* ones(Complex{Float64}, length(freqList))
    Fdo_full = -1 .* ones(Complex{Float64}, length(freqList))
    for (k,v) in freqList_map
        Fup_full[k] = TwoPartGF_upup[v]
        Fdo_full[k] = TwoPartGF_updo[v]
    end
    return Fup_full, Fdo_full
end

