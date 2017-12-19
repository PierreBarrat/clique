"""
   ReturnFreqs(msa::Array{Int64,2}; q=findmax(msa)[1], reweighting=0)

Return single point and pairwise frequencies for input `msa`. 

# Warning

Reweighting is not yet implemented
"""
function ReturnFreqs(msa::Array{Int64,2}; q=findmax(msa)[1], reweighting=0)

    (M,L) = size(msa)
    f_2 = zeros(Float64,L*q, L*q)
    f_1 = zeros(Float64,L*q)

    for m in 1:M
        for j in 1:L
            f_1[(j-1)*q+msa[m,j]] += 1
            f_2[(j-1)*q+msa[m,j], (j-1)*q+msa[m,j]] += 1
            for i in (j+1):L
                f_2[(i-1)*q+msa[m,i], (j-1)*q+msa[m,j]] += 1
                f_2[(j-1)*q+msa[m,j], (i-1)*q+msa[m,i]] += 1
            end
        end
    end

    f_2 = f_2./M
    f_1 = f_1./M

    return (f_1,f_2)

end

"""
   ReturnPairFreqs(i::Int64, j::Int64, msa::Array{Int64,2}; q=findmax(msa)[1], reweighting=0)

Return single point and pairwise frequencies for columns `i` and `j` in input `msa`. `i` is indexed as 1, and `j` as 2.

# Warning

Reweighting is not yet implemented
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
    Fapc(A::Array{Float64,2}, q::Int64 ; APC::Bool = false)

Compute frobenius norm, with optional APC correction, for blocks of size `q` in input matrix `A`.
"""
function Fapc(A::Array{Float64,2}, q::Int64 ; APC::Bool = false)

    (L,L) = size(A)
    try
        L = convert(Int64,L/q)
    catch
        println("Size of input matrix is ",size(A))
        error("Input matrix A does not have the correct size.")
    end
    if !(A==A')
        error("Input matrix A is not symmetric")
    end

    S = zeros(Float64, L,L)
    F = zeros(Float64, convert(Int64, L*(L-1)/2))
    for i in 0:L-1
        for j in (i+1):L-1
            S[i+1,j+1] = vecnorm(A[i*q + (1:q), j*q + (1:q)])
            S[j+1,i+1] = S[i+1,j+1]
        end
    end

    pos = 1
    for i in 1:L
        for j in (i+1):L
            F[pos] = S[i,j] - (mean(S[:,j]) * mean(S[:,i]) / mean(S)) * APC
            pos += 1
        end
    end

    return F

end

"""
    ComputeWeights(msa::String ; theta::Float64 = 0.2)

Compute weights of each sequence in file `msa` using reweighting threshold `theta` (default 0.2).
"""
function ComputeWeights(msa::String; theta::Float64 = 0.2, output::String = @sprintf("%s_weights",msa))
    Y = readdlm(msa);
    if findmin(Y[2:end,:])[1]==0
        Y = Y[2:end,:];
        Y = Y +1;
    end
    Y = convert(Array{Int64,2},Y);
    ComputeWeights(Y,theta=theta,output=output)
end

function ComputeWeights(Y::Array{Int64,2};theta::Float64=0.2,output::String="weights.txt")
    Y = Y';
    M::Int64 = size(Y,2)
    L::Int64 = size(Y,1)
    h::Float64 = L * (1-theta)
    println(h)
    weights = ones(Float64, M,1);
    d::Int64=0
    for m = 1:M
        for l = (m+1):M
            d = 0
            for i = 1:L
                d += Y[i,m]==Y[i,l];
            end
            if d > h
                weights[m]+=1
                weights[l]+=1
            end
        end
    end

    writedlm(output, 1./weights, ' ')
end




# Full matrix of pairwise hamming distances between configurations
function pdist(alignment_1::Array{Int64,2}, alignment_2::Array{Int64,2} = alignment_1)

    (M1,N1) = size(alignment_1)
    (M2,N2) = size(alignment_2)

    if N1!=N2
        error("Both alignments should have the same number of columns")
        exit(1)
    end

    if N1<128
        pdist_array = zeros(Int8,M1,M2)
        temp = convert(Int8,0);
    else
        pdist_array = zeros(Int16,M1,M2)
        temp = convert(Int16,0);
    end
    
    for m in 1:M1
        for n in 1:M2
            for a = 1:N1
                temp += (alignment_1[m,a]!=alignment_2[n,a])
            end
            pdist_array[m,n] = temp;
            temp = 0;
        end
    end

    return pdist_array
end


# Smallest distances from alignment to itself, or from alignment 1 to alignment 2
function PdistMin(alignment_1::Array{Int64,2}, K::Int64 = 1; alignment_2 = alignment_1)

    (M1,N1) = size(alignment_1)
    (M2,N2) = size(alignment_2)

    if N1!=N2
        error("Both alignments should have the same number of columns")
        exit(1)
    end

    pdist_array = fill(N1+1, K, M1)
    index_array = zeros(Int64,K,M1)
    temp = 0;
    for m in 1:M1
        for n in 1:M2
            for a = 1:N1
                temp += (alignment_1[m,a]!=alignment_2[n,a])
            end
            k = K;
            while k>0 && temp < pdist_array[k,m]
                if k==K
                    pdist_array[k,m] = temp
                    index_array[k,m] = n
                else
                    pdist_array[k+1,m] = pdist_array[k,m]
                    index_array[k+1,m] = index_array[k,m]
                    pdist_array[k,m] = temp
                    index_array[k,m] = n
                end
                k-=1
            end
            temp = 0
        end
    end

    return (pdist_array, index_array)

end 



