@testset "freq_to_index" begin
    tt = falses(length(freqList))
    for (i,el) in enumerate(freqList)
        ωi_t,νi_t,νpi_t =  VertexPostprocessing.Freq_to_OneToIndex(el[1], el[2], el[3], shift, nBose, nFermi)
        tt[i] = all(el .== freqList_t[ωi_t,νi_t,νpi_t])
    end
    @test all(tt)
end