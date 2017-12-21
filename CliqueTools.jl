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





"""
    RemoveNode!(clique_n::clique_type, clique_o::clique_type, k::Int64)

Removes node `k` from `clique_o` and stores result in `clique_n`. Node `k` corresponds to the total system mapping ; the node that is removed in practice is `clique_o.map[k]`.

# Arguments
-`clique_n::clique_type` : should have the correct dimensions.
-`clique_o::clique_type` : will not be modified.
-`k::Int64` : node to remove, expressed in the total system mapping.
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
#    clique_n.sample = clique_o.sample[:,cat(1,1:(clique_o.map[k]-1),(clique_o.map[k]+1):clique_o.L)]
    for i in 1:(clique_o.map[k]-1)
        for m = 1:Mn
            clique_n.sample[m,i] = clique_o.sample[m,i]
        end
    end
    for i in (clique_o.map[k]+1):clique_o.L
        for m = 1:Mn
            clique_n.sample[m,i-1] = clique_o.sample[m,i]
        end
    end

    for i in 1:Ln
        for j in (i+1):Ln
            ic = clique_o.map[clique_n.nodes[i]]
            jc = clique_o.map[clique_n.nodes[j]]
            if ic == 0 || jc == 0
                prinln("problem with mapping")
            end
            # clique_n.J[(i-1)*qn + (1:qn), (j-1)*qn + (1:qn)] .= clique_o.J[(ic-1)*qn + (1:qn), (jc-1)*qn + (1:qn)]
            # clique_n.J[(j-1)*qn + (1:qn), (i-1)*qn + (1:qn)] .= clique_o.J[(jc-1)*qn + (1:qn), (ic-1)*qn + (1:qn)]

            ### Less memory
            for b = 1:qn
                for a = 1:qn
                    @inbounds clique_n.J[(i-1)*qn + a, (j-1)*qn + b] = clique_o.J[(ic-1)*qn + a, (jc-1)*qn + b]
                    @inbounds clique_n.J[(j-1)*qn + a, (i-1)*qn + b] = clique_o.J[(jc-1)*qn + a, (ic-1)*qn + b]
                end
            end
            ###
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

"""
"""
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

Returns pairwise frequencies for columns `i` and `j` of `clique.sample`. Output is a `qxq` matrix.
"""
function PairFreq(i::Int64, j::Int64, clique::clique_type)

    (f_1,f_2) = ReturnPairFreqs(i,j,clique.sample, q=clique.q)
    q = clique.q

    return f_2[0*q+(1:q), 1*q+(1:q)]
end


"""
   ReturnPairFreqs(i::Int64, j::Int64, msa::Array{Int64,2}; q=findmax(msa)[1], reweighting=0)

Return single point and pairwise frequencies for columns `i` and `j` in input `msa`. `i` is indexed as 1, and `j` as 2.
"""
function ReturnPairFreqs(i::Int64, j::Int64, msa::Array{Int64,2}; q::Int64=findmax(msa)[1], reweighting=0)

    (M,L) = size(msa)
    f_2 = zeros(Float64, 2*q, 2*q)
    f_1 = zeros(Float64, 2*q)
    # println("sample[1,:] = ", msa[1,:])
    # println("i = ",i," -- j = ",j)
    for m in 1:M
        f_1[0*q+msa[m,i]] += 1
        f_1[1*q+msa[m,j]] += 1

        f_2[0*q+msa[m,i], 0*q+msa[m,i]] += 1
        f_2[1*q+msa[m,j], 1*q+msa[m,j]] += 1
        f_2[0*q+msa[m,i], 1*q+msa[m,j]] += 1
        f_2[1*q+msa[m,j], 0*q+msa[m,i]] += 1
    end

    f_2 = f_2./M
    f_1 = f_1./M

    return (f_1,f_2)

end

"""
   ReturnPairFreqs!(f_2::Array{Float64,2},i::Int64, j::Int64, msa::Array{Int64,2}; q=findmax(msa)[1], reweighting=0)

Return single point and pairwise frequencies for columns `i` and `j` in input `msa`. Output is a `qxq` matrix.
"""
function ReturnPairFreqs!(f_2::Array{Float64,2}, i::Int64, j::Int64, msa::Array{Int64,2}, q::Int64, reweighting=0)

    (M) = size(msa)[1]
    eps::Float64 = 1/M
    # println("sample[1,:] = ", msa[1,:])
    # println("i = ",i," -- j = ",j)
    for a = 1:q
        for b = 1:q
            f_2[a,b] = 0
        end
    end
    for m in 1:M
        f_2[msa[m,i], msa[m,j]] += eps
    end

    # copy!(f_2,f_2/M)
    # f_2

end

"""
"""
function SampleFromClique!(i::Int64, j::Int64, clique::clique_type, n_it::Int64)
    (M,L) = size(clique.sample)
    tau::Int64 = 0
    # Parameters
    n_itp = n_it
    df::Int64 = 30 # M*L/df accepted swaps define one iteration
    min_n_acc::Int64 = floor(M*L/df)
    #Allocating space
    freq_out = zeros(Float64,clique.q,clique.q)
    freq_temp = zeros(Float64,clique.q,clique.q)

    # Estimating reasonnable sampling time
    tau = EstimateTau!(i, j, clique.sample, clique.J, clique.q, min_n_acc) # tau --> M*L/df accepted swaps
    tau_pair::Int64 = floor(df*tau/L/2) # tau_pair --> M/2 accepted swaps in those columns
    
    # Equilibrating for eqt tau
    eqt::Int64 = floor(df/5) # M*L/5 accepted swaps define eq. time
    DoSwap!(clique.sample,clique.J,eqt*tau,clique.q)

    # Sampling n_it times
    for it in 1:n_it
        DoSwap!(clique.sample,clique.J,tau,clique.q)
        ### No Swapping Pairs
        # ReturnPairFreqs!(freq_temp,i,j,clique.sample,clique.q)
        # for a = 1:clique.q
        #     for b = 1:clique.q
        #         freq_out[a,b] = freq_out[a,b] + freq_temp[a,b] / n_it 
        #     end
        # end
        ###

        ### Swapping Pairs
        for itp in 1:n_itp
            DoSwapPair!(clique.sample,clique.J,tau_pair,i,j,clique.q)
            ReturnPairFreqs!(freq_temp,i,j,clique.sample,clique.q)
            for a = 1:clique.q
                for b = 1:clique.q
                    freq_out[a,b] = freq_out[a,b] + freq_temp[a,b] / n_it / n_itp
                end
            end
        end
        ###
    end
    
    return freq_out
end