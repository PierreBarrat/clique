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





