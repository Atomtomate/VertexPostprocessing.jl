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
