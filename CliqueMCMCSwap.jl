include("./CliqueTools.jl")

"""
"""
function MCSwap!(sample::Array{Int64,2}, J::Array{Float64,2}, q::Int64 ; it_max::Int64 = 50    , it_min::Int64 = 5 , verbose::Bool=true , random_control::Bool = false)

    (M,L) = size(sample)

    tau = convert(Int64,floor(M*L/10))

    # Stop conditions
    eps = 0.01
    acc_rate_store = zeros(Float64,it_max)
    E = zeros(Float64, it_max+1)
    
    t::Int64 = 1
    for t in 1:it_max
        (E[t+1], acc_rate) = DoSwap!(sample, J, tau, q)
        E[t+1] = E[t+1] + E[t]
        acc_rate_store[t] = acc_rate
    end

end

"""
    function EstimateTau!(i::Int64, j::Int64, sample::Array{Int64,2}, J::Array{Float64,2}, q::Int64, min_n_acc::Int64)

Estimates reasonnable number of swaps as number of trials to achieve `min_n_acc` accepted moves.
"""
function EstimateTau!(i::Int64, j::Int64, sample::Array{Int64,2}, J::Array{Float64,2}, q::Int64, min_n_acc::Int64)

    (M,L) = size(sample)
    tau::Int64 = 0
    it_max = convert(Int64,floor(M))*10;
    n_acc = zeros(Int64, it_max+1)
    n_acc[1] = 0
    # min_n_acc = M*L/10;

    step = L;
    t::Int64 = 1
    # samplew = copy(sample);
    for t in 1:it_max
        n_acc[t+1] = DoSwap!(sample,J,step,q)[2]
        n_acc[t+1] = n_acc[t+1] + n_acc[t]
        if n_acc[t+1] > min_n_acc
            break
        end
    end
    tau = convert(Int64,t * L)
    return tau
end

"""
"""
function DoSwap!(sample::Array{Int64,2}, J::Array{Float64,2}, tau::Int64, q::Int64)

    (M,L) = size(sample)
    n_accept::Int64 = 0
    dE_tot::Float64 = 0
    dE::Float64 = 0
    m1::Int64 = 0
    m2::Int64 = 0
    i::Int64 = 0
    ia::Int64 = 0
    ib::Int64 = 0

    for t in 1:tau
        dE = 0
        i = rand(1:L)
        m1 = rand(1:M)
        m2 = rand(1:M)

        # If move does not change anything, accept it. 
        if(sample[m1,i]==sample[m2,i])
            n_accept+=1
        else
            ia = (i-1)*q + sample[m1,i]
            ib = (i-1)*q + sample[m2,i]
            for j in 1:L
                @fastmath @inbounds dE += J[(j-1)*q + sample[m1,j], ia] - J[(j-1)*q + sample[m2,j], ia] # Contribution from moving a from m1 to m2
                @fastmath @inbounds dE += J[(j-1)*q + sample[m2,j], ib] - J[(j-1)*q + sample[m1,j], ib] # Contribution from moving b from m2 to m1
            end
            if (dE < 0 || exp(-dE) > rand(Float64))
                # println("   --- OK ---")
                ia = sample[m1,i]
                sample[m1,i] = sample[m2,i]
                sample[m2,i] = ia
                n_accept+=1
                dE_tot += dE
            end
        end
    end

    return (dE_tot, n_accept)
end

"""
function DoSwapPair!(sample::Array{Int64,2}, J::Array{Float64,2}, tau::Int64, ir::Int64, jr::Int64, q::Int64)

    Same as `DoSwap!`, but only swaps in columns `ir` and `jr`. `sample` is changed in the process.
"""
function DoSwapPair!(sample::Array{Int64,2}, J::Array{Float64,2}, tau::Int64, ir::Int64, jr::Int64, q::Int64)

    (M,L) = size(sample)
    n_accept::Int64 = 0
    dE_tot::Float64 = 0
    dE::Float64 = 0
    m1::Int64 = 0
    m2::Int64 = 0
    i::Int64 = 0
    ia::Int64 = 0
    ib::Int64 = 0
    cols::Array{Int64,1} = [ir,jr]

    for t in 1:tau
    	dE = 0
        i = cols[rand(1:2)]
        m1 = rand(1:M)
        m2 = rand(1:M)

        # If move does not change anything, accept it. 
        if(sample[m1,i]==sample[m2,i])
            n_accept+=1
        else
            ia = (i-1)*q + sample[m1,i]
            ib = (i-1)*q + sample[m2,i]
            for j in 1:L
                @fastmath @inbounds dE += J[(j-1)*q + sample[m1,j], ia] - J[(j-1)*q + sample[m2,j], ia] # Contribution from moving a from m1 to m2
                @fastmath @inbounds dE += J[(j-1)*q + sample[m2,j], ib] - J[(j-1)*q + sample[m1,j], ib] # Contribution from moving b from m2 to m1
            end
            if (dE < 0 || exp(-dE) > rand(Float64))
                # println("   --- OK ---")
                ia = sample[m1,i]
                sample[m1,i] = sample[m2,i]
                sample[m2,i] = ia
                n_accept+=1
                dE_tot += dE
            end
        end
    end

    return (dE_tot, n_accept)
end


# function MCSwap(sample::Array{Int64,2}, J::Array{Float64,2}, q::Int64 ; it_max::Int64 = 50    , it_min::Int64 = 5 , verbose::Bool=true , random_control::Bool = false)

#     (M,L) = size(sample)
#     sample_n = deepcopy(sample)

#     tau = floor(Int64,M*L/5)

#     # Stop conditions
#     eps = 0.01
#     acc_rate_store = zeros(Float64,it_max)
#     E = zeros(Float64, it_max+1)
#     if random_control
#         sample_rand = ShuffleColumns(sample)
#         hamming_distance = zeros(Float64,it_max)
#         hamming_distance_rand = zeros(Float64,it_max)
#     end
#     flag_max = 5
#     flag = 0
    

#     if verbose
#         print("Starting MCSwap for ",it_max,"x",tau,"iterations ...")
#     end
    
#     t::Int64 = 1
#     for t in 1:it_max
#         (E[t+1], acc_rate) = DoSwap!(sample_n, J, tau, q)
#         E[t+1] = E[t+1] + E[t]
#         acc_rate_store[t] = acc_rate
        
#         if random_control
#             DoSwap!(sample_rand, J, tau, q)
#             hamming_distance[t] = sum(sample.!=sample_n) / L / M
#             hamming_distance_rand[t] = sum(sample.!=sample_rand) / L / M
#         end
#     end
#     if(t==it_max)
#         if verbose
#             print("MCSwap : maximum number of iterations reached (",it_max,")")
#         end
#     end
#     if verbose
#         @printf("MCSwap done - number of iterations %d (with %d swaps per iteration)\n", t, tau)
#     end
#     # return(sample_n, t, hamming_distance, hamming_distance_rand, E, sample_rand)
#     return (sample_n,t)

# end
