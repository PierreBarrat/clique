type clique_type
    nodes::Array{Int64,1} # nodes present in the clique
    sample::Array{Int64,2}
    J::Array{Float64,2}
    map::Array{Int64,1} # map from clique position to full system position
    q::Int64
    L::Int64
end

#
#
function CreateZeroClique(L::Int64, q::Int64, M::Int64; Lmap = L)
    clique = clique_type(zeros(Int64,L), zeros(Int64,M,L), zeros(Float64,L*q,L*q), zeros(Int64,L), q, L)
    return clique
end


#
# Adds node k to clique_o, modifying clique_n
# Information on sample and couplings is taken from clique_ref, which should represent the whole system
function AddNode!(clique_n::clique_type, clique_o::clique_type, clique_ref::clique_type, k::Int64)
    Ln = clique_o.L + 1
    q = clique_o.q
    clique_n.L = Ln
    clique_n.q = q

    # List of nodes
    clique_n.nodes[1:(Ln-1)] = copy(clique_o.nodes)
    clique_n.nodes[Ln] = k
    sort!(clique_n.nodes)

    # Dealing with map
    pos_k = findin(clique_n.nodes, k)[1]
    clique_n.map = copy(clique_o.map)
    clique_n.map[k] = pos_k
    for l in (k+1):size(clique_n.map)[1]
        if clique_n.map[l] !=0 
            clique_n.map[l] += 1
        end
    end

    # Updating sample and J
    clique_n.sample = clique_ref.sample[:,clique_n.nodes]
    for pos_i in 1:Ln
        for pos_j in (pos_i+1):Ln
            i = clique_n.nodes[pos_i]
            j = clique_n.nodes[pos_j]
            clique_n.J[(pos_i-1)*q + (1:q), (pos_j-1)*q + (1:q)] = clique_ref.J[(i-1)*q + (1:q), (j-1)*q + (1:q)]
            clique_n.J[(pos_j-1)*q + (1:q), (pos_i-1)*q + (1:q)] = clique_ref.J[(j-1)*q + (1:q), (i-1)*q + (1:q)]
        end
    end
end 



"""
    RemoveNode!(clique_n::clique_type, clique_o::clique_type, k::Int64)

Removes node `k` from `clique_o` and stores result in `clique_n`. Node `k` corresponds to the total system mapping ; the node that is removed in practice is `clique_o.map[k]`.

# Arguments
-`clique_n::clique_type` : should have the correct dimensions.
-`clique_o::clique_type` : will not be modified.
-`k::Int64` : node to remove in the total system mapping.
"""
function RemoveNode!(clique_n::clique_type, clique_o::clique_type, k::Int64)
    Ln = clique_o.L - 1
    qn = clique_o.q
    Mn = size(clique_o.sample)[1]

    if Mn!=size(clique_n.sample)[1] || Ln != clique_n.L || qn != clique_n.q
        println("WARNING - CliqueTools.jl in RemoveNode!: clique_n does not seem to have the correct size.")
        clique_n.q = qn
        clique_n.L = Ln
    end

    kc = clique_o.map[k]

    clique_n.nodes = cat(1, clique_o.nodes[1:(kc-1)], clique_o.nodes[(kc+1):end])
    clique_n.map = copy(clique_o.map)
    clique_n.map[k] = 0
    clique_n.map[(k+1):end] = max(clique_n.map[(k+1):end]-1,0)
    clique_n.sample = clique_o.sample[:,cat(1,1:(clique_o.map[k]-1),(clique_o.map[k]+1):clique_o.L)]

    for i in 1:Ln
        for j in (i+1):Ln
            ic = clique_o.map[clique_n.nodes[i]]
            jc = clique_o.map[clique_n.nodes[j]]
            if ic == 0 || jc == 0
                prinln("problem with mapping")
            end
            clique_n.J[(i-1)*qn + (1:qn), (j-1)*qn + (1:qn)] = clique_o.J[(ic-1)*qn + (1:qn), (jc-1)*qn + (1:qn)]
            clique_n.J[(j-1)*qn + (1:qn), (i-1)*qn + (1:qn)] = clique_o.J[(jc-1)*qn + (1:qn), (ic-1)*qn + (1:qn)]
        end
    end
end


"""
    RemoveNode(clique_o::clique_type, k::Int64)

Removes node `k` from input `clique_o`. Output is a new `clique_type` object with reduced dimensions.
"""
function RemoveNode(clique_o::clique_type, k::Int64)

    (M,L) = size(clique_o.sample)
    Ln = L - 1
    qn = clique_o.q
    clique_n = CreateZeroClique(Ln,qn,M, Lmap = size(clique_o.map)[1])
    
    kc = clique_o.map[k] # Position of node to remove in old clique
    
    clique_n.nodes = cat(1, clique_o.nodes[1:(kc-1)], clique_o.nodes[(kc+1):end])
    clique_n.map = copy(clique_o.map)
    clique_n.map[k] = 0
    clique_n.map[(k+1):end] = max(clique_n.map[(k+1):end]-1,0)
    clique_n.sample = clique_o.sample[:,cat(1,1:(clique_o.map[k]-1),(clique_o.map[k]+1):clique_o.L)]
    
    for i in 1:Ln
        for j in (i+1):Ln
            ic = clique_o.map[clique_n.nodes[i]]
            jc = clique_o.map[clique_n.nodes[j]]
            if ic == 0 || jc == 0
                prinln("problem with mapping")
            end
            clique_n.J[(i-1)*qn + (1:qn), (j-1)*qn + (1:qn)] = clique_o.J[(ic-1)*qn + (1:qn), (jc-1)*qn + (1:qn)]
            clique_n.J[(j-1)*qn + (1:qn), (i-1)*qn + (1:qn)] = clique_o.J[(jc-1)*qn + (1:qn), (ic-1)*qn + (1:qn)]
        end
    end

    return clique_n
end

#
#
function CliqueCopy(clique::clique_type)

    clique_out = CreateZeroClique(clique.L,clique.q,size(clique.sample)[1])
    clique_out.nodes = copy(clique.nodes)
    clique_out.sample = copy(clique.sample)
    clique_out.J = copy(clique.J)
    clique_out.map = copy(clique.map)

    return clique_out
end


"""
    CliqueCorr!(i::Int64, j::Int64, clique::clique_type)

Computes frobenius norm of correlation between `i` and `j` in `clique`.
"""
function CliqueCorr(i::Int64, j::Int64, clique::clique_type)

    (f_1,f_2) = ReturnPairFreqs(i,j,clique.sample, q=clique.q)
    q = clique.q

    return f_2[0*q+(1:q), 1*q+(1:q)]
end
    # if exclude_gap
    #     Fcij = sqrt(sum((f_2[0*q + (2:q), 1*q + (2:q)] - f_1[0*q + (2:q)] * f_1[1*q + (2:q)]').^2))
    # else
    #     Fcij = sqrt(sum((f_2[0*q + (1:q), 1*q + (1:q)] - f_1[0*q + (1:q)] * f_1[1*q + (1:q)]').^2))
    # end

"""
    CliqueDistribution!(clique::clique_type)

Computes single site and pairwise frequencies for `clique`, using MCSwap sampling.
"""
function CliqueDistribution!(clique::clique_type)

    (clique.sample, n_it) = MCSwap(clique.sample, clique.J, clique.q, verbose = true)
    
    (f_1,f_2) = ReturnFreqs(clique.sample, q=clique.q)

    return (f_1,f_2)
end

# """
#     CliquePairDistribution!(i::Int64, j::Int64, clique::clique_type)
# Computes single site and pairwise frequencies for sites `i` and `j` in `clique`, using MCSwap sampling.
# """
# function CliquePairDistribution!(i,j,clique::clique_type)

#     (sample_n, n_it) = MCSwap(clique.sample, clique.J, clique.q, verbose = true)
#     clique.sample = sample_n
#     (f_1,f_2) = ReturnPairFreqs(i,j,clique.sample, q=clique.q)

#     return (f_1,f_2)
# end

"""
    TwoNodesFromFull(i::Int64, j::Int64, clique_ref::clique_type)

Remove all nodes from `clique_ref`, except for `i` and `j`. Return clique object.
"""
function TwoNodesFromFull(i::Int64, j::Int64, clique_ref::clique_type)
    M = size(clique_ref.sample)[1]
    L = clique_ref.L
    q = clique_ref.q

    clique_out = CliqueCopy(clique_ref)
    for k in 1:L
        if k!=i && k!=j
            clique_out = RemoveNode(clique_out, k)
        end
    end

    return clique_out
end