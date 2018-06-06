include("./CliqueMCMCSwap.jl")
include("./CliqueTools.jl")


"""
   CliqueMain(i::Int64, j::Int64, sample::Array{Int64,2}, J::Array{Float64,2}, q::Int64 ; file::String = "")

Main script to find cliques. Target pair is (`i`,`j`). `file` is a string used to identify output log and frequencies.
"""
function CliqueMain(i::Int64, j::Int64, sample::Array{Int64,2}, J::Array{Float64,2}, q::Int64 ; file::String = "")

    (M,L) = size(sample)
    clique_tot = clique_type(collect(1:L), sample, J, collect(1:L), q, L)

    # Initialize log and freqs files
    f = open(@sprintf("log_%s.txt",file),"w")
    close(f)
    f = open(@sprintf("freqs_%s.txt",file),"w")
    close(f)

    # Main calculation
    (score_array, removed_nodes) = CliqueFromFull(clique_tot, i, j, file = file)

    # Writing output
    writedlm(@sprintf("score_array_%s.txt",file), score_array, ' ')
    writedlm(@sprintf("removed_nodes_%s.txt",file), removed_nodes, ' ')

    #
    return (score_array, removed_nodes)
end

"""
"""
function WriteLog(file::String,k,L,removed_node, score, freq_curr)
    open(@sprintf("log_%s.txt",file),"a") do f
        write(f,@sprintf("It. %d out of %d\n",k,L-2)) 
        write(f,@sprintf("Removed node %d -- New score = %f\n\n", removed_node,score) )
    end
    open(@sprintf("freqs_%s.txt",file),"a") do f
        write(f,@sprintf("It. %d out of %d\n",k,L-2)) 
        writedlm(f,freq_curr,' ')
        write(f,@sprintf("\n"))
    end
end

"""
   CliqueFromFull(clique_ref::clique_type, i::Int64, j::Int64)

Iteratively removes nodes from `clique_ref` keeping the correlation between `i` and `j` maximal. Outputs values of correlations and removed nodes in format `(corr_array::Array{Float64,1}, removed_nodes::Array{Int64,1})`. 
"""
function CliqueFromFull(clique_ref::clique_type, i::Int64, j::Int64 ; file::String="")
    L = size(clique_ref.sample)[2]
    n_it::Int64 = 10
    score_array = zeros(Float64, L-1)
    removed_nodes = zeros(Int64, L-2)
    t_clique = clique_ref

    # Computing original score -- accounting for finite size and sampling time
    # f0 = PairFreq(t_clique.map[i], t_clique.map[j], t_clique)
    # MCSwap!(t_clique.sample, t_clique.J, t_clique.q, verbose = false, it_max = 30, random_control = false)
    # freq_init = PairFreq(t_clique.map[i], t_clique.map[j], t_clique)
    f0 = SampleFromClique!(t_clique.map[i], t_clique.map[j], t_clique,n_it)
    freq_init = SampleFromClique!(t_clique.map[i], t_clique.map[j], t_clique,n_it)
    score_array[1] = ScoresFromFreqs(freq_init, f0)[1]
    println("Initial score = ", score_array[1])
    WriteLog(file, 0, L, 0, score_array[1],freq_init)

    # Removing nodes iteratively - it is the iteration number
    for it in 1:(L-2)
        println("It. ", it, " out of ",L-2)
        # Computing
        (t_clique, removed_nodes[it], score_array[it+1],freq_curr) = RemoveOptNode(t_clique, i, j, freq_init)
        println("Removed node ", removed_nodes[it], "-- New score = ", score_array[it+1])
        println("Remaining sites: ", size(t_clique.sample)[2])
        println(" ")
        # Log
        WriteLog(file,it,L,removed_nodes[it], score_array[it+1], freq_curr)
    end 

    return (score_array, removed_nodes)
end

"""
   RemoveOptNode(clique::clique_type, i::Int64, j::Int64)

Attempts to remove one node from `clique`. `freq_init` is an array containing pairwise correlation between `i` and `j` in the original sample. 
"""
function RemoveOptNode(clique::clique_type, i::Int64, j::Int64, freq_init::Array{Float64,2})
    (M,L) = size(clique.sample)
    #Conveniencies
    q::Int64 = clique.q
    mapi::Int64 = clique.map[i]
    mapj::Int64 = clique.map[j]
    # Allocating space
    t_clique = CreateZeroClique(L-1, clique.q, M)
    freq_array = zeros(Float64,q * L, q)
    # Parameters
    n_it::Int64 = 10; # Number of swapping procedures to compute frequencies for one temptative node
    #Initialisation
    score_max::Float64 = -1
    opt_node::Int64 = 0

    # Computing pairwise frequencies for all temptative nodes k 
    for k in 1:L
        if k!=mapi && k!=mapj
            @printf("Node %d/%d...     \r",k,L)
            RemoveNode!(t_clique,clique, clique.nodes[k])
            freq_array[(k-1)*q+(1:q),:] = SampleFromClique!(t_clique.map[i], t_clique.map[j], t_clique,n_it)
        end
    end
    
    # Finding best node to remove
    node_scores = ScoresFromFreqs(freq_array ,freq_init)
    node_scores[mapi] = findmin(node_scores)[1]-1 # Setting scores of i and j to low values 
    node_scores[mapj] = node_scores[mapi]

    (score_max, opt_node) = findmax(node_scores)
    while opt_node == mapi || opt_node == mapj
        opt_node = rand(1:L);
        println("opt node is mapi or mapj"); println("opt node is mapi or mapj"); println("opt node is mapi or mapj")
    end

    # Removing node and equilibrating sample
    RemoveNode!(t_clique, clique, clique.nodes[opt_node])
    SampleFromClique!(t_clique.map[i], t_clique.map[j], t_clique,1)

    #
    return (t_clique, clique.nodes[opt_node], score_max, freq_array[(opt_node-1)*q+(1:q),:])
end


"""
    ScoresFromFreqs(freq_array::Array{Float64,2}, f0::Array{Float64,2}; exclude_gap::Bool = false)

Computes score for each temptative node. `freq_array` is (`L`*`q` x `q`) and contains measures pairwise frequencies for each removed node. 
"""
function ScoresFromFreqs(freq_array::Array{Float64,2}, f0::Array{Float64,2}; exclude_gap::Bool = true, scoretype = "corr")

    q::Int64 = size(freq_array)[2]
    L::Int64 = size(freq_array)[1]/q
    offset::Int64 = 1;
    if exclude_gap
        offset = 2;
    end

    node_scores = zeros(Float64,L);
    for k in 1:L
        if scoretype == "corr"
            # Minus frobenius norm of difference of correlation matrices. (I want to maximize scores).
            node_scores[k] = -sqrt(sum( (freq_array[(k-1)*q+(offset:q),(offset:q)] - f0[offset:q,offset:q]).^2 ));
        elseif scoretype == "MI"
            node_scores[k] = ComputeMI(freq_array[(k-1)*q+(offset:q),(offset:q)])# - ComputeMI(f0[offset:q,offset:q])
        elseif scoretype == "DKL"
            computeDKL(f0[offset:q,offset:q], freq_array[(k-1)*q+(offset:q),(offset:q)]))
        end
    end
    return node_scores
end 

"""
    ComputeMI(f)

Compute and return mutual information correpsonding to pairwise frequencies in `f`.
"""
function ComputeMI(f)
    q = size(f)[1]
    pc = 1e-5
    f = (1-pc)*f + pc/q/q
    m1 = sum(f,1)
    m2 = sum(f,2)
    return sum(f.*log.(f./(m2*m1)))
end













