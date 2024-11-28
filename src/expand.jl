# ==================================================================================================== #
#                                           SVertex.jl                                                 #
# ---------------------------------------------------------------------------------------------------- #
#   Authors         : Julian Stobbe                                                                    #
# ----------------------------------------- Description ---------------------------------------------- #
#   SVertex (Sparse Vertex) type, that stores a minimal representation 
# -------------------------------------------- TODO -------------------------------------------------- #
# ==================================================================================================== #

# -------------------- Auxilliary ----------------------

# ----------------- Public Interface -------------------
"""
    expand_2pt_GF(TwoPartGF_upup::Vector, TwoPartGF_updo::Vector, transform!::Function, 

Expands two particle Green's function according to the relationships in `transform!`. 
Used by [`expand_2PtGF_CSV`](@ref expand_2PtGF_CSV)
"""
function expand_2pt_GF(TwoPartGF_upup::Vector, TwoPartGF_updo::Vector, transform!::Function, 
                freqList_map, freqList, parents, ops, nBose, nFermi)
    off(f) = f[1]+nBose+1,f[2]+nFermi+1,f[3]+nFermi+1
    Fup_full = Array{eltype(TwoPartGF_upup)}(undef, length(freqList))
    Fdo_full = Array{eltype(TwoPartGF_updo)}(undef, length(freqList))
    done = falses(length(freqList))
    for (k,v) in freqList_map
        Fup_full[k] = TwoPartGF_upup[v]
        Fdo_full[k] = TwoPartGF_updo[v]
        done[k] = true
    end
    open = Stack{eltype(parents)}()
    for i in 1:length(freqList)
        next = i
        while !done[next]
            push!(open, next)
            next = parents[next]
        end
        while length(open) > 0 
            prev = next
            next = pop!(open)
            done[next] = true
            transform!(Fup_full, Fdo_full, prev, next, ops)
        end
    end
    return Fup_full, Fdo_full
end

"""
    combine_TwoPartGF(
        freqList_1::Vector, freqList_2::Vector, 
        TwoPartGF_upup_1::Vector, TwoPartGF_upup_2::Vector,
        TwoPartGF_updo_1::Vector, TwoPartGF_updo_2::Vector
    )

Combines two 2-particle Green's functions.
The smaller frequency grid must be fully contained in the larger one.
"""
function combine_TwoPartGF(
        freqList_1::Vector, freqList_2::Vector, 
        TwoPartGF_upup_1::Vector, TwoPartGF_upup_2::Vector,
        TwoPartGF_updo_1::Vector, TwoPartGF_updo_2::Vector
    )
    prio_on_1 = length(freqList_1) > length(freqList_2)
    freqList = prio_on_1 ? freqList_1 : freqList_2
    freqList_t = deepcopy(freqList)
    ii = sortperm(freqList)
    freqListSorted = freqList[ii]
    freqList_small = !prio_on_1 ? freqList_1 : freqList_2
    TwoPartGF_upup = prio_on_1 ? TwoPartGF_upup_1 : TwoPartGF_upup_2
    TwoPartGF_updo = prio_on_1 ? TwoPartGF_updo_1 : TwoPartGF_updo_2
    TwoPartGF_upup_small = !prio_on_1 ? TwoPartGF_upup_1 : TwoPartGF_upup_2
    TwoPartGF_updo_small = !prio_on_1 ? TwoPartGF_updo_1 : TwoPartGF_updo_2

    indices = zeros(Int, length(freqList_small))
    for (i,f) in enumerate(freqList_small)
        ind = searchsortedfirst(freqListSorted, f) #findfirst(x-> all(x .== f), freqList)
        if isnothing(ind)
            error("Key from smaller frequency list not found in larger one! Non-overlapping lists not implemented yet!")
        else
            TwoPartGF_upup[ii[ind]] = TwoPartGF_upup[i]
            TwoPartGF_updo[ii[ind]] = TwoPartGF_updo[i]
            #Debug: freqList_t[ii[ind]] = f 
        end
    end
    return TwoPartGF_upup, TwoPartGF_updo
end
