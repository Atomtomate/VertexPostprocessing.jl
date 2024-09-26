# ==================================================================================================== #
#                                          IO_legacy.jl                                                #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#   Read/Writeoperations that keep compatibility with a set of private Fortran77 codes                 #
# -------------------------------------------- TODO -------------------------------------------------- #
# ==================================================================================================== #


# ========================================== Dependencies ============================================
"""
    read_vertex_csv(fname::String)

Reads two particle function from a human readable file with 7 columns of `Float64`, separated in a way,
that `readdlm` will recognize:
```
f1::Float64 f2::Float64 f3::Float64 val1_re::Float64 val1_im::Float64  val2_re::Float64 val2_im::Float64
```

Returns `3` arrays:
```
[f1, f2, f3], [val1_re .+ 1im .* val1_im], [val2_re .+ 1im .* val2_im]
```
"""
function read_two_pt_function(fname::String)
    inp = readdlm(fname, Float64)
    freq = Float64.(inp[:,1:3])
    Fupup = inp[:,4] + inp[:,5]*1im
    Fupdo = inp[:,6] + inp[:,7]*1im
    return freq, Fupup, Fupdo
end

"""
read_chi_asympt(fname::String)

Reads file with `7` columns, similar to [`read_two_pt_function`](@ref read_two_pt_function).
The first colums are the frequencies, the remaining correspond to `3` one-frequency objects. 
"""
# function read_chi_asympt(fname::String)
#     inp = readdlm(fname, Float64)
#     freq = Float64.(inp[:,1])
#     χ_upup = inp[:,4] + inp[:,5]*1im
#     χ_updo = inp[:,4] - inp[:,5]*1im
#     χpp = inp[:,6] + inp[:,7]*1im
#     return freq, χ_upup, χ_updo, χpp
# end
function read_chi_asympt(fname::String)
    inp = readdlm(fname, Float64)
    freq = Float64.(inp[:,1])
    χ_upup = inp[:,2] + inp[:,3]*1im
    χ_updo = inp[:,4] + inp[:,5]*1im
    χpp = inp[:,6] + inp[:,7]*1im
    χm, χd = uu_ud_TO_m_d(χ_upup, χ_updo)
    return freq, χd, χm, χpp
end

# function readGImp(filename; only_positive=false)

# function readGImp(filename; only_positive=false)
#     GFs = readdlm(filename)[:,1:end]


#     tmp = GFs[:,2] .+ GFs[:,3].*1im
#     tmpiνₙ = GFs[:,1] .* 1im
#     if only_positive
#         GImp = tmp
#         iνₙ  = tmpiνₙ
#     else
#         NH = size(tmp,1)
#         N  = 2*NH
#         GImp = zeros(Complex{Float64}, N)
#         iνₙ  = zeros(Complex{Float64}, N)
#         GImp[1:(NH)] = reverse(conj.(tmp[1:NH]))
#         GImp[(NH+1):N] = tmp[1:NH]
#         iνₙ[1:(NH)] = conj.(reverse(tmpiνₙ[1:(NH)]))
#         iνₙ[(NH+1):N] = tmpiνₙ[1:NH]
#     end
#     return iνₙ, GImp
# end


# function read_gm_wim(nFreq, filename; storedInverse, storeFull=false)
#     GFs = readdlm(filename)[:,2:end]

#     if size(GFs, 1)*(1 + 1*storeFull) < nFreq
#         throw(ArgumentError("nFermFreq in simulation parameters too large! Got $(size(GFString, 1)) lines but need data for $(nFreq) frequencies."))
#     end
#     tmp = GFs[:,1] .+ GFs[:,2].*1im

#     if storedInverse
#         tmp = 1 ./ tmp
#     end
#     
#     GF = Array{Complex{Float64}}(undef, nFreq)
#     if storeFull
#         NH = Int(nFreq/2)
#         GF[1:(NH-1)] = reverse(conj.(tmp[2:NH]))
#         GF[NH:nFreq] = tmp[1:NH]
#     else
#         GF = tmp[1:nFreq]
#     end
#     return GF
# end

# function read_densimp(file)
#     content = open(file) do f
#         readlines(f)
#     end
#     return parse(Float64,content[1])
# end

# function read_hubb_dat(file)
#     content = open(file) do f
#         readlines(f)
#     end
#     U = parse(Float64,split(replace(content[2], "d0" => ""), r"(\s|,)+")[1])
#     β = parse(Float64,split(replace(content[4], "d0" => ""), r"(\s|,)+")[1])
#     n = parse(Float64,split(replace(content[8], "d0" => ""), r"(\s|,)+")[2])
#     return U, β, n
# end

"""
    read_anderson_parameters(file)

Reads Anderson parameters from legacy Fortran human readable file
"""
function read_anderson_parameters(file)
    content = open(file) do f
        readlines(f)
    end
    
    in_epsk = false
    in_tpar = false
    ϵₖ = []
    Vₖ = []
    μ = 0
    for line in content
        if "Eps(k)" == strip(line)
            in_epsk = true
            continue
        elseif "tpar(k)" == strip(line)
            in_epsk = false
            in_tpar = true
            continue
        end
        
        if in_epsk
            push!(ϵₖ, parse(Float64, line))
        elseif in_tpar
            # skip last line, which is mu
            if length(Vₖ) < length(ϵₖ)
                push!(Vₖ, parse(Float64, line))
            else
                if occursin("#", line)
                    μ = parse(Float64, line[1:(findfirst("#", line))[1] - 1])
                else
                    μ = parse(Float64, line)
                end
            end
        end
    end
    return convert(Array{Float64,1}, ϵₖ), convert(Array{Float64,1}, Vₖ), μ
end

# function write_vert_chi(freqList::Array, ver::Array{Complex{Float64},1}, verdo::Array{Complex{Float64},1}, dirname::String, nBose::Int, nFermi::Int)
#     open(dirname * "/vert_chi", "w") do f
#         for (i,freq) in enumerate(freqList)
#             @printf(f, "  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f\n", 
#                     float(freq[1]), float(freq[2]), float(freq[3]),
#                 real(ver[i]), imag(ver[i]),
#                 real(verdo[i]), imag(verdo[i]))
#         end
#     end
# end

# function write_fort_dir(prefix::String, freqList::Array, arr_ch::Array{Complex{Float64},3}, arr_sp::Array{Complex{Float64},3}, dirname::String, nBose::Int, nFermi::Int)
#     mkpath(dirname)
#     for ωn in 1:size(arr_ch,1)
#         freqSegment = (ωn-1)*nFermi*nFermi+1:(ωn+0)*nFermi*nFermi
#         freq_sub_grid = freqList[freqSegment]
#         open(dirname * "/" * prefix * lpad(ωn-1, 3, "0"), "w") do f
#             for i in 1:nFermi
#                 for j in 1:nFermi
#                     @printf(f, "  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f\n", 
#                         float(freq_sub_grid[i][1]), float(freq_sub_grid[(i-1)*(nFermi)+j][2]), 
#                         float(freq_sub_grid[(i-1)*(nFermi)+j][3]),
#                         real(arr_ch[ωn, i, j]), imag(arr_ch[ωn, i, j]),
#                         real(arr_sp[ωn, i, j]), imag(arr_sp[ωn, i, j]))
#                 end
#             end
#         end
#     end
# end

# function write_fort_dir(prefix::String, freqList::Array, arr_ch::Array{Complex{Float64},1}, arr_sp::Array{Complex{Float64},1}, dirname::String, nBose::Int, nFermi::Int)
#     mkpath(dirname)
#     nF2 = nFermi*nFermi
#     for ωn in 1:nBose
#         freqSegment = (ωn-1)*nF2+1:(ωn+0)*nF2
#         freq_sub_grid = freqList[freqSegment]
#         open(dirname * "/" * prefix * lpad(ωn-1, 3, "0"), "w") do f
#             for i in 1:nF2
#                 ind = (ωn-1)*nF2+i
#                 @printf(f, "  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f  %18.10f\n", 
#                         float(freq_sub_grid[i][1]), float(freq_sub_grid[i][2]), float(freq_sub_grid[i][3]),
#                         real(arr_ch[ind]), imag(arr_ch[ind]),
#                         real(arr_sp[ind]), imag(arr_sp[ind]))
#             end
#         end
#     end
# end
