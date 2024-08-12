
function identityMap(nB_min, nB_max, nF_min, nF_max, shift, offset, base)
    freqList = [(i,j,k) for i in (nB_min:nB_max) for j in (nF_min:nF_max) .- trunc(Int64,shift*i/2) for k in (nF_min:nF_max) .- trunc(Int64,shift*i/2)]
    fullMap = Dict(tripleToInt(freqList[i][1],freqList[i][2],freqList[i][3],offset,base,base*base) => UInt32(i) for i in 1:length(freqList))
    return fullMap
end

function cut_margin(freqList, arr, Ncut_ω::Int, Ncut_ν::Int, Nν_max::Int, Nω_max::Int, shift::Int)
    freqList_cut = similar(freqList, 0)
    arr_cut      = similar(arr, 0)
    for i in 1:size(freqList,1)
        ω, ν, νp = freqList[i]
        # i < 10 && println(abs(ω), " ?>? ", Nω_max-Ncut_ω,", ", abs(ν+shift*ω/2), " ?>? ", Nν_max-Ncut_ν, ", ", abs(νp+shift*ω/2), " ?>? ", Nν_max-Ncut_ν)
        if abs(ω) <= Nω_max-Ncut_ω && ν+shift*ω/2 >=  -(Nν_max-Ncut_ν) && ν+shift*ω/2 < (Nν_max-Ncut_ν) &&  νp+shift*ω/2 >=  -(Nν_max-Ncut_ν) && νp+shift*ω/2 < (Nν_max-Ncut_ν) 
            push!(freqList_cut, freqList[i])
            push!(arr_cut, arr[i])
            # i < 100 && print("pushing: ",freqList[i])
        end
    end
    return freqList_cut, arr_cut
end

function Freq_to_OneToIndex(ωn::Int, νn::Int, νpn::Int, shift::Union{Bool,Int}, nBose::Int, nFermi::Int)
    (ωn+nBose+1,νn+nFermi+1+trunc(Int, shift*ωn/2), νpn+nFermi+1+trunc(Int, shift*ωn/2))
end


"""
    find_non_nan_matrix(data::Matrix, nFermi::Int)

Indices for slice of matrix not containing any NaNs.
"""
function find_non_nan_matrix(data::Matrix, nFermi::Int)
    nan_list       = sort(map(x->x[1],filter(x->x[1] == x[2], findall(x-> !isnan(x), data))))
    nh = searchsorted(nan_list, nFermi)
    res = if length(nh) == 1
        t1 = diff(nan_list[nh[1]:end])
        t2 = -1 * diff(nan_list[nh[1]:-1:1])
        lim_up = findfirst(x->x != 1, t1)
        lim_lo = findfirst(x->x != 1, t2)
        lim_up = isnothing(lim_up) ? length(t1) : lim_up
        lim_lo = isnothing(lim_lo) ? length(t2) : lim_lo
        lim_non_nan = min(lim_up, lim_lo)
        (nFermi-lim_non_nan):(nFermi+lim_non_nan)
    else
        []
    end
    return res
end

# ==================== GF Stuff ====================
# ------------------ Bare Bubble -------------------



# ------------------ Full Vertex -------------------
function FUpDo_from_χDMFT(χupdo, GImp, freqFile, β)
    FUpDo = similar(χupdo)
    jldopen(freqFile, "r") do freqFile
        freqList = freqFile["freqList"]
        n_iω = freqFile["nBose"]
        n_iν = freqFile["nFermi"]
        shift = freqFile["shift"]
        for f in freqList
            i = f[1] + n_iω+1
            j = f[2] + n_iν+1 + trunc(Int, shift*f[1]/2)
            k = f[3] + n_iν+1 + trunc(Int, shift*f[1]/2)
            FUpDo[i,j,k] = χupdo[i,j,k]/(β^2 * get_symm_f(GImp,f[2]) * get_symm_f(GImp,f[1]+f[2])
                               * get_symm_f(GImp,f[3]) * get_symm_f(GImp,f[1]+f[3]))
        end
    end
    return FUpDo
end




function χph_to_χpp(freqList::Array, χph_upup::Array, χph_updo::Array, χ0::Dict{Tuple{Int,Int},Complex{Float64}}, shift, nBose::Int, nFermi::Int)
    χpp_s = fill!(Array{eltype(χph_upup),3}(undef, 2*nFermi, 2*nFermi, 2*nBose+1), NaN)
    χpp_t = fill!(Array{eltype(χph_upup),3}(undef, 2*nFermi, 2*nFermi, 2*nBose+1), NaN)
    old_shape = size(χph_upup)
    χph_upup = reshape(χph_upup, 2*nFermi, 2*nFermi, 2*nBose+1)
    χph_updo = reshape(χph_updo, 2*nFermi, 2*nFermi, 2*nBose+1)

    for i in 1:size(freqList,1)
        ωn, νn, νpn = freqList[i]
        ωi,νi,νpi = Freq_to_OneToIndex(ωn, νn, νpn, shift, nBose, nFermi)
        ωi_ph,νi_ph,νpi_ph = Freq_to_OneToIndex(ωn - νn - νpn - 1, νn, νpn, shift, nBose, nFermi)
        if !(any((ωi,νi,νpi) .< 1) || any((ωi_ph,νi_ph,νpi_ph) .< 1) || any((ωi,νi,νpi) .> (2*nBose,2*nFermi-1,2*nFermi-1)) || any((ωi_ph,νi_ph,νpi_ph) .> (2*nBose,2*nFermi-1,2*nFermi-1)))
            χpp_s[νi,νpi,ωi] = - χ0[(ωn,νn)]*(νn==νpn) - χph_upup[νi_ph,νpi_ph,ωi_ph] + 2*χph_updo[νi_ph,νpi_ph,ωi_ph]
            χpp_t[νi,νpi,ωi] = + χ0[(ωn,νn)]*(νn==νpn) + χph_upup[νi_ph,νpi_ph,ωi_ph]
        end
    end
    χph_upup = reshape(χph_upup, old_shape)
    χph_updo = reshape(χph_updo, old_shape)
    χpp_s, χpp_t
end

function calc_E_ED(iνₙ, ϵₖ, Vₖ, GImp, n, U, β, μ; full=false)
    E_kin = 0.0
    E_pot = 0.0
    vk = sum(Vₖ .^ 2)
    Σ_hartree = n * U/2
    E_pot_tail = (U^2)/2 * n * (1-n/2) - Σ_hartree*(Σ_hartree-μ)
    E_kin_tail = vk

    for n in 1:length(GImp)
        Δ_n = sum((Vₖ .^ 2) ./ (iνₙ[n] .- ϵₖ))
        Σ_n = iνₙ[n] .- Δ_n .- 1.0 ./ GImp[n] .+ μ
        E_kin += 2*real(GImp[n] * Δ_n - E_kin_tail/(iνₙ[n]^2))
        E_pot += 2*real(GImp[n] * Σ_n - E_pot_tail/(iνₙ[n]^2))
    end
    E_kin = E_kin .* (2/β) - (β/2) .* E_kin_tail
    E_pot = E_pot .* (1/β) .+ 0.5*Σ_hartree .- (β/4) .* E_pot_tail
    return E_kin, E_pot
end
