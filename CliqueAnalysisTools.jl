type PwFreq
   	f2::Array{Float64}
    q::Int64
end

type ScoreType
	MIdiff::Float64
	Cdiff::Float64
end

"""
"""
function ReadFreqs(ff::String, q::Int64)

	f2str = Array{PwFreq,1}(0)
	f2t = zeros(Float64,q,q)
	a = 1
	f = open(ff,"r")
	count = 1
	for l in eachline(f)
		if l!="" && l[1:2]=="It"
			a = 1
			f2t .= 0
		elseif l!=""
			map!(x->parse(x), view(f2t,a,:), split(l,' '))
			a += 1
		elseif l==""
		 	push!(f2str, PwFreq(copy(f2t),q))
		end
	end
	return f2str
end

"""
	ScoreFromFreqs(f2_0::Array{Float64,2}, f2::Array{Float64,2}, q::Int64)

Using initial frequencies `f2_0` and frequencies from clique `f2`, should compute every score imaginable. 
Currently
1. Difference between MI
2. Fro. norm of difference between off-diagonal correlations
"""
function ScoreFromFreqs(f2_0::Array{Float64,2}, f2::Array{Float64,2}, q::Int64)
	out = ScoreType(0.,0.)

	# Computing marginals and off-diag correlations
	f1_0_j = sum(f2_0,1)
	f1_0_i = sum(f2_0,2)
	f1_j = sum(f2,1)
	f1_i = sum(f2,2)
	C_0 = f2_0 - f1_0_i*f1_0_j
	C = f2 - f1_i*f1_j

	# Difference between MI
	out.MIdiff = ComputeMI(f2_0,q) - ComputeMI(f2,q)

	# Fro. norm of difference of C -- off diag. part
	out.Cdiff = vecnorm(C - C_0)

	return out
end

"""
	ComputeMI(f2::Array{Float64,2}, q::Int64)	
"""
function ComputeMI(f2::Array{Float64,2}, q::Int64)
	f1_i = sum(f2,2)
	f1_j = sum(f2,1)
	eps = 0.0001
	MI = 0.
	for a in 1:q
		for b in 1:q
			MI += f2[a,b]* log(f2[a,b]/f1_i[a]/f1_j[b]+eps)
		end
	end
	return MI
end

"""
    FittingQuality(f2_1::Array{Float64,2}, f1_1::Array{Float64,1}, f2_2::Array{Float64,2}, f1_2::Array{Float64,1}, q::Int64; withdiag=false)
Compares frequencies `(f1_1,f2_1)` to `(f1_2,f2_2)`. Output is 
1. Pearson correlation between connected correlations. (default: without diagonal elements)
2. Frobenius norm of the difference between connected correlations.
3. Same as 1. for magnetizations.
4. Same as 2. for magnetizations.
"""
function FittingQuality(f2_1::Array{Float64,2}, f1_1::Array{Float64,1}, f2_2::Array{Float64,2}, f1_2::Array{Float64,1}, q::Int64; withdiag=false)
    L::Int64=size(f2_1)[1]/q
    C1 = f2_1 - f1_1*f1_1'
    C2 = f2_2 - f1_2*f1_2'
    id = zeros(Bool,L*q,L*q)
    for i = 1:L
        for j = (i+1):L
            id[(j-1)*q+(1:q),(i-1)*q+(1:q)]=true
        end
    end
    cc = cor(C1[id],C2[id])
    froc = vecnorm(C1[id] - C2[id])
    cm = cor(f1_1,f1_2)
    from = vecnorm(f1_1-f1_2)
    return (cc,froc,cm,from)
end


"""
"""
function InEllipse(clique, dmat::Array{Float64,2}, i::Int64, j::Int64 ; alpha=1.5, threshold=-1)
	threshold<0?(threshold = alpha * dmat[i,j]):threshold=threshold
	# true?threshold=1;threshold=2
	return map(x-> ((dmat[i,x]>0&&dmat[j,x]>0)?(dmat[i,x] + dmat[j,x] < threshold):0), clique)

end