include("../Misc/AlignmentStat.jl")
include("../MCMC/MCMCSwap.jl")
include("../Clique/CliqueTools.jl")


function FindBestCliques(i::Int64, j::Int64, sample::Array{Int64,2}, J::Array{Float64,2}, q::Int64 ; parallel = false, file::String = "")

    (M,L) = size(sample)
    clique_tot = clique_type(collect(1:L), sample, J, collect(1:L), q, L)

    # Initialize log and freqs files
    if file!=""
        f = open(@sprintf("log_%s",file),"w")
        close(f)
        f = open(@sprintf("freqs_%s",file),"w")
        close(f)
    end

    (score_array, removed_nodes) = CliqueFromFull(clique_tot, i, j, parallel = parallel, file = file)

    return (score_array, removed_nodes)
end

"""
"""
function WriteLog(file::String,k,L,removed_node, score, freq_curr)
    if file!=""
        open(@sprintf("log_%s",file),"a") do f
            write(f,@sprintf("It. %d out of %d\n",k,L-2)) 
            write(f,@sprintf("Removed node %d -- New score = %f\n\n", removed_node,score) )
        end
        open(@sprintf("freqs_%s",file),"a") do f
            writedlm(f,freq_curr,' ')
        end
    end
end

"""
    CliqueFromFull(clique_ref::clique_type, i::Int64, j::Int64)

Iteratively removes nodes from `clique_ref` keeping the correlation between `i` and `j` maximal. Outputs values of correlations and removed nodes in format `(corr_array::Array{Float64,1}, removed_nodes::Array{Int64,1})`. 
"""
function CliqueFromFull(clique_ref::clique_type, i::Int64, j::Int64 ; parallel::Bool = false , file::String="")
    L = size(clique_ref.sample)[2]
    score_array = zeros(Float64, L-1)
    removed_nodes = zeros(Int64, L-2)
    t_clique = clique_ref

    MCSwap!(t_clique.sample, t_clique.J, t_clique.q, verbose = false, it_max = 30, random_control = false)
    freq_init = CliqueCorr(t_clique.map[i], t_clique.map[j], t_clique)
    score_array[1] = ScoresFromFreqs(freq_init ; f0 = freq_init, exclude_gap = true)[1]

    println("Initial score = ", score_array[1])
    WriteLog(file, 0, L, 0, score_array[1],freq_init)
    for k in 1:(L-2)
        println("It. ", k, " out of ",L-2)
        if parallel
            (t_clique, removed_nodes[k], score_array[k+1],freq_curr) = RemoveOptNodePar(t_clique, i, j, freq_init)
        else
            (t_clique, removed_nodes[k], score_array[k+1],freq_curr) = RemoveOptNode(t_clique, i, j, freq_init)
        end
        println("Removed node ", removed_nodes[k], "-- New score = ", score_array[k+1])
        println(" ")
        WriteLog(file,k,L,removed_nodes[k], score_array[k+1], freq_curr)
    end 

    return (score_array, removed_nodes)
end

"""
    RemoveOptNode(clique::clique_type, i::Int64, j::Int64) 
"""
function RemoveOptNode(clique::clique_type, i::Int64, j::Int64, freq_init = Array{Float64,2})
    (M,L) = size(clique.sample)
    t_clique = CreateZeroClique(L-1, clique.q, M)
    score_max::Float64 = -1
    opt_node::Int64 = 0
    q::Int64 = clique.q
    mapi::Int64 = clique.map[i]
    mapj::Int64 = clique.map[j]

    freq_array = Array(Float64,q * L, q)
    for k in 1:L
        if k!=mapi && k!=mapj
            RemoveNode!(t_clique,clique, clique.nodes[k])
            MCSwap!(t_clique.sample, t_clique.J, t_clique.q, verbose = false, it_max = 15, random_control = false)
            freq_array[(k-1)*q+(1:q),:] = CliqueCorr(t_clique.map[i], t_clique.map[j], t_clique)
        end
    end
    node_scores = ScoresFromFreqs(freq_array ; exclude_gap=true, f0 = freq_init)
    node_scores[mapi] = findmin(node_scores)[1]-1
    node_scores[mapj] = node_scores[mapi]

    (score_max, opt_node) = findmax(node_scores)
    while opt_node == mapi || opt_node == mapj
        opt_node = rand(1:L);
    end
    RemoveNode!(t_clique, clique, clique.nodes[opt_node])
    MCSwap!(t_clique.sample, t_clique.J, t_clique.q, verbose = false, it_max = 15, random_control = false)
    return (t_clique, clique.nodes[opt_node], score_max, freq_array[(opt_node-1)*q+(1:q),:])
end


"""
    RemoveOptNodePar(clique::clique_type, i::Int64, j::Int64 )
"""
function RemoveOptNodePar(clique::clique_type, i::Int64, j::Int64, freq_init = Array{Float64,2})
    (M,L) = size(clique.sample)
    score_max::Float64 = -1
    opt_node::Int64 = 0
    q::Int64 = clique.q
    mapi::Int64 = clique.map[i]
    mapj::Int64 = clique.map[j]

    freq_array = SharedArray(Float64,q * L, q)
    @sync @parallel for k in 1:L
        if k!=mapi && k!=mapj
            t_clique = RemoveNode(clique, clique.nodes[k])
            MCSwap!(t_clique.sample, t_clique.J, t_clique.q, verbose = false, it_max = 15, random_control = false)
            freq_array[(k-1)*q+(1:q),:] = CliqueCorr(t_clique.map[i], t_clique.map[j], t_clique)
        end
    end
    node_scores = ScoresFromFreqs(freq_array ; exclude_gap=true, f0 = freq_init)
    node_scores[mapi] = findmin(node_scores)[1]-1
    node_scores[mapj] = node_scores[mapi]

    (score_max, opt_node) = findmax(node_scores)
    while opt_node == mapi || opt_node == mapj
        opt_node = rand(1:L);
    end
    t_clique = RemoveNode(clique, clique.nodes[opt_node])
    MCSwap!(t_clique.sample, t_clique.J, t_clique.q, verbose = false, it_max = 15, random_control = false)
    return (t_clique, clique.nodes[opt_node], score_max, freq_array[(opt_node-1)*q+(1:q),:])
end

"""
"""
function ScoresFromFreqs(freq_array::Union{Array{Float64,2},SharedArray{Float64,2}} ; exclude_gap::Bool = true, f0::Array{Float64,2}=zeros(Float64,size(freq_array)[2],size(freq_array)[2]) )

    q::Int64 = size(freq_array)[2]
    L::Int64 = size(freq_array)[1]/q
    offset::Int64 = 1;
    if exclude_gap
        offset = 2;
    end

    node_scores = zeros(Float64,L);
    for k in 1:L
        node_scores[k] = -sqrt(sum( (freq_array[(k-1)*q+(offset:q),(offset:q)] - f0[offset:q,offset:q]).^2 ));
    end
    return node_scores
end 


""" 
    CliqueFromTwo(clique_i::clique_type, clique_ref::clique_type)

Iteratively add nodes to `clique_i` in order to maximize the correlation between initial nodes of `clique_i` at each step. Outputs values of correlations and added nodes in format `(corr_array::Array{Float64,1}, added_nodes::Array{Int64,1}).

# Arguments
- `clique_i::clique_type` should contain two nodes only.
- `clique_ref::clique_type` should contain a full system.
"""
function CliqueFromTwo(clique_i::clique_type, clique_ref::clique_type)

    (M,L_ref) = size(clique_ref.sample)
    q = clique_ref.q
    L = size(clique_i.nodes)[1]
    i = clique_i.nodes[1]
    j = clique_i.nodes[2]

    added_nodes = zeros(Int64, L_ref - L)
    corr_array = zeros(Float64, L_ref - L + 1)
    corr_array[1] = CliqueCorr!(clique_i.map[i], clique_i.map[j], clique_i)
    
    out_nodes = symdiff(clique_i.nodes, clique_ref.nodes)
    curr_L = L
    curr_corr = 0.0
    curr_clique = CliqueCopy(clique_i)
    pos = 1
    while curr_L < L_ref
        println("Length ",curr_L," out of ", L_ref)
        t_clique = CreateZeroClique(curr_L + 1, q, M)
        new_node = 0
        new_corr = 0
        new_sample = zeros(Int64, M, curr_L+1)
        for k in out_nodes
            print("k = ",k," -- ")
            AddNode!(t_clique, curr_clique, clique_ref, k)
            curr_corr = CliqueCorr!(t_clique.map[i], t_clique.map[j], t_clique)
            println("corr = ", curr_corr)
            if curr_corr > new_corr
                new_corr = curr_corr
                new_node = k
                new_sample = t_clique.sample
            end
        end
        curr_L += 1
        AddNode!(t_clique, curr_clique, clique_ref, new_node)
        curr_clique = CliqueCopy(t_clique)
        curr_clique.sample = new_sample
        corr_array[pos+1] = new_corr
        added_nodes[pos] = new_node
        out_nodes = symdiff(out_nodes, [new_node])
        pos +=1
    end

    return (corr_array, added_nodes)

end












