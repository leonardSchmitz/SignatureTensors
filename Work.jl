@testset "Truncated Tensor Algebra Tests :p2id" begin

    @testset "Constructor TTA" begin
        d = 6        # path dimension
        k = 5        # truncation level
        T = TruncatedTensorAlgebra(QQ, d, k, sequence_type=:p2id)
        @test T == TruncatedTensorAlgebra(QQ, d, k, sequence_type=:p2id)
        @test sequence_type(T) == :p2id
        @test base_dimension(T) == d
        @test base_algebra(T) == QQ
        @test truncation_level(T) == k
    end
   
function axis_core_3tensor_QQ_p2id(_d)
    # La misma lógica que tu axis_core_3tensor_QQ, pero aquí podrías cambiar coeficientes si p2id difiere
    C = zeros(QQ, _d, _d, _d)
    for al in 1:_d
        for be in 1:_d
            for ga in 1:_d
                if al == be && be == ga
                    C[al, be, ga] = QQ(1,6)
                end
                if (al < be && be == ga) || (al == be && be < ga)
                    C[al, be, ga] = QQ(1,2)
                end
                if al < be && be < ga
                    C[al, be, ga] = one(QQ)
                end
            end
        end
    end
    return C
end

 function membrane_signature_QQ(m::Int, n::Int)
    Caxis_m = axis_core_3tensor_QQ_p2id(m)
    Caxis_n = axis_core_3tensor_QQ_p2id(n)

    d_m = size(Caxis_m, 1)
    d_n = size(Caxis_n, 1)

    d = d_m * d_n

    M = zeros(QQ, d, d, d)

    for i in 1:d_m, j in 1:d_m, k in 1:d_m
        for a in 1:d_n, b in 1:d_n, c in 1:d_n

            a2 = (i-1)*d_n + a
            b2 = (j-1)*d_n + b
            c2 = (k-1)*d_n + c

            M[a2, b2, c2] += Caxis_m[i,j,k] * Caxis_n[a,b,c]
        end
    end

    return M
end


 @testset "Axis constructor in TTA for QQ :p2id" begin
     d = 6
     T = TruncatedTensorAlgebra(QQ, d, 4, sequence_type=:p2id)
     m=2
     n=3
    shape = (m, n)
     # Aquí usamos la función que genera el eje tipo p2id
     Caxis_d = sig(T, :axis, shape=shape)
    
     @test parent(Caxis_d) == T
     @test zero(T) + zero(T) == zero(T)
     @test Caxis_d == sig(T, :axis, shape=shape, algorithm=:Chen)
     @test Caxis_d == sig(T, :axis, shape=shape, algorithm=:AFS19)
    
     for i in 2:m-1
         @test Caxis_d[i] == one(QQ)
         @test Caxis_d[i, (i-1)] == zero(QQ)
         @test Caxis_d[i, i] == QQ(1,2)
         @test Caxis_d[i, i+1] == one(QQ)
     end
      
     @test membrane_signature_QQ(m, n) == Caxis_d[:,:,:]
    
          @test zero(T) + Caxis_d == Caxis_d
          @test one(T) * Caxis_d == Caxis_d
          @test Caxis_d * one(T) == Caxis_d
          @test inv(Caxis_d) * Caxis_d == one(T)
          @test Caxis_d * inv(Caxis_d) == one(T)
          @test inv(inv(Caxis_d)) == Caxis_d
      #    @test exp(log(Caxis_d)) == Caxis_d
          @test Caxis_d^3 == Caxis_d * Caxis_d * Caxis_d
      #    @test exp(log(Caxis_d) + log(Caxis_d)) == Caxis_d^2
      end




   end

    @testset "Axis constructor in TTA for polynomial rings" begin
        d = 6;
        R, a = polynomial_ring(QQ, :a => (1:d, 1:d));
        T = TruncatedTensorAlgebra(R,d,5, sequence_type=:p2id);
        m=2
        n=3
        shape = (m, n)
        Caxis_d = sig(T,:axis,shape=shape); 
        @test parent(Caxis_d) == T
        @test zero(T) + zero(T) == zero(T)
        @test Caxis_d == sig(T,:axis,algorithm=:Chen,shape=shape)
        @test Caxis_d == sig(T,:axis,algorithm=:AFS19,shape=shape) 
         
        @test zero(T) + Caxis_d == Caxis_d
        @test one(T)*Caxis_d == Caxis_d
        @test Caxis_d*one(T) == Caxis_d
        @test inv(Caxis_d)*Caxis_d == one(T)
        @test Caxis_d*inv(Caxis_d) == one(T)
        @test inv(inv(Caxis_d)) == Caxis_d
        #@test exp(log(Caxis_d)) == Caxis_d
        @test Caxis_d^3 == Caxis_d*Caxis_d*Caxis_d
        #@test exp(log(Caxis_d) + log(Caxis_d)) == Caxis_d^2
    end

    @testset "Pwln constructor in TTA" begin
        d = 6; k = 4; number_tests = 5;
        T = TruncatedTensorAlgebra(QQ,d,k,:p2id);
        ms = rand((1:7),number_tests); # maximal 7 number of segments
        As = [QQ.(rand((-20:20),d,m)) for m in ms];
        for A in As
          m = size(A)[2]  # number of segments in A
          S = sig(T,:pwln,coef=A); 
          @test S == sig(T,:pwbln,coef=A,algorithm=:Chen); 
          @test S == sig(T,:pwbln,coef=A,algorithm=:congruence); 
          C = sig(TruncatedTensorAlgebra(QQ,m,k,:p2id),:axis,shape); # axis tensor in dim m
          @test S == A*C
          @test zero(T) + S == S
          @test one(T)*S == S
          @test S*one(T) == S
          @test inv(S)*S == one(T)
          @test S*inv(S) == one(T)
          @test inv(inv(S)) == S
          @test exp(log(S)) == S
          @test S^3 == S*S*S
          @test exp(log(S) + log(S)) == S^2
        end 
    end    
    
    @testset "Converting base algebra via matrix tensor congruence" begin
        d = 6; m = 5; k = 4;
        R, a = polynomial_ring(QQ, :a => (1:d, 1:m));
        TmQQ = TruncatedTensorAlgebra(QQ,m,k);
        TdR = TruncatedTensorAlgebra(R,d,k);
        TmR = TruncatedTensorAlgebra(R,m,k);
        C = sig(TmQQ,:axis);
        S = sig(TdR,:pwln,coef=a);
        @test S == a*C
        A = QQ.(rand((-20:20),d,m));
        @test A*sig(TmR,:axis) ==  sig(TdR,:pwln,coef=A);
    end 




#typeof(membrane)=Vector{Matrix{Float64}}



d = 6
T = TruncatedTensorAlgebra(QQ, d, 4, sequence_type=:p2id)
m=2
n=3
shape = (m, n)
# Aquí usamos la función que genera el eje tipo p2id
Caxis_d = sig(T, :axis, shape=shape)
    
@test parent(Caxis_d) == T
@test zero(T) + zero(T) == zero(T)
@test Caxis_d == sig(T, :axis, shape=shape, algorithm=:Chen)
@test Caxis_d == sig(T, :axis, shape=shape, algorithm=:AFS19)
    
coef=[   3  -2  -2  1   2   5;
 -1  -3   4  3   5   0;
  2  -3   1  0   3   3;
  2   4  -4  0  -1  -4;
 -3   5  -1  1   4  -3;
 2   3  1  3   4  -5;

 ]

#shape


pwbln=sig(T, :pwbln,coef=coef,shape=shape)
pwbln2=sig(T, :pwbln,coef=coef,shape=shape,algorithm=:LS)

pwbln==pwbln2


d2=size(coef, 1)


membrane = Array{Int64}(undef, m, n, d2)

for di in 1:d2
    cont = 0
    for i in 1:m
        for j in 1:n
            cont += 1
            membrane[i, j, di] = coef[di, cont]
        end
    end
end

pwblnm=sig(T, :pwbln,coef=membrane,shape=shape)
pwblnm2=sig(T, :pwbln,coef=membrane,shape=shape,algorithm=:LS)

pwblnm==pwblnm2
pwbln==pwbln2
pwbln==pwblnm
pwbln==pwblnm2



d2=size(coef, 1)

membrane = Vector{Matrix{Float64}}(undef, d2)
   
for di in 1:d2
    M = Matrix{Float64}(undef, m, n)
    cont = 0
    for i in 1:m
        for j in 1:n
            cont += 1
            M[i,j] = coef[di, cont]
        end
    end
    membrane[di] = M
end

membrane
result = signature_levels(membrane, 4)

pwbln[:]
pwbln2[:]
result[2]

pwbln[:,:]
pwbln2[:,:]
result[3]


sig_lott(membrane, [1,1,1])

 
pwbln[:,:,1]==pwbln2[:,:,1]
pwbln2[:,:,1]
result[4][:,:,1]

pwbln==pwbln2

membrane






pwbln==pwbln2


pwbln[:,:,3]

pwbln2[:,:,3]




#sig_pwbln_p2id_Congruence(T,coef,shape=)


sig(T, :pwbln,coef=coef,shape=shape,algorithm=:LS)==sig(T, :pwbln,coef=coef,shape=shape)
sig(T, :pwbln,coef=coef,shape=shape,algorithm=:LS)





m=shape[1]
n=shape[2]
R=RealField()

T = TruncatedTensorAlgebra(R, 6, 4, sequence_type=:p2id)


membrane = [rand(m,n)*2 .- 1 for _ in 1:d]
membrane_tensor = cat(membrane...; dims=3)

coef = transpose(reshape(membrane_tensor,m*n,d))

coef2= coef*0

coef3=[  3  -2.0  -2.0  1.0   2   5;
 -1  -3   4  3   5   0;
  2  -3   1  0   3   3;
  2   4  -4  0  -1  -4;
 -3   5  -1  1   4  -3]

membrane_tensor = cat(membrane...; dims=3)

for di in 1:d
    cont=0
    for i in 1:m
        for j in 1:n
            cont=cont+1
            coef2[di,cont]=membrane_tensor[i,j,di]
        end
    end
end


pwbln=sig(T, :pwbln,coef=coef2,shape=shape)

Caxis_d2=sig(T, :axis, shape=(n,m))

Caxis_d = sig(T, :axis, shape=shape)
matrix=[1.0 1.0;
    1.0 1.0]
for i in 1:2
    for j in 1:2
        word=[i,j]
        matrix[i,j]=sig_lott(membrane, word)
    end
end

lotter=sig_lott(membrane, word)

pwbln[:,:]
matrix
Caxis_d[:,:]

coef
sig(T, :pwbln,coef=coef,shape=shape,algorithm=:LS)==sig(T, :pwbln,coef=coef,shape=shape)

sig(T, :pwbln,coef=coef,shape=shape,algorithm=:LS)


d=5
m=2
n=3

T = TruncatedTensorAlgebra(QQ, d, 3, sequence_type=:p2id);
A = QQ.(rand(-1*20:20, d, m,n))
R,a = polynomial_ring(QQ,:a=>(1:d,1:m,1:n));

AC = sig(T, :poly, shape=shape,coef=A); 
aC = sig(T, :poly, shape=shape,coef=a); 

T = TruncatedTensorAlgebra(QQ, d, 3, sequence_type=:p2id);
A = QQ.(rand(-1*20:20, d, m,n))
R,a = polynomial_ring(QQ,:a=>(1:d,1:m,1:n));
AC = sig(T, :poly, shape=shape,coef=A); 
aC = sig(T, :poly, shape=shape,coef=a); 













julia> using Oscar

julia> d = 4
4

julia> m = 2 
2

julia> n = 22
22

julia> n = 2
2


T = TruncatedTensorAlgebra(QQ, d, 3, sequence_type=:p2id);

A = QQ.(rand(-1*20:20, d, m,n))

R,a = polynomial_ring(QQ,:a=>(1:d,1:m,1:n));

AC = sig(T, :poly, shape=shape,coef=A); 
ERROR: BoundsError: attempt to access 3-element Vector{Array{QQFieldElem}} at index [4]
Stacktrace:
 [1] throw_boundserror(A::Vector{Array{QQFieldElem}}, I::Tuple{Int64})
   @ Base ./essentials.jl:15
 [2] getindex(A::Vector{Array{QQFieldElem}}, i::Int64)
   @ Base ./essentials.jl:919
 [3] applyMatrixToTTA(A::Matrix{QQFieldElem}, X::TruncatedTensorAlgebraElem{QQField, QQFieldElem})
   @ SignatureTensors ~/Documents/Berlin_projects/project_signature-tensors-in-OSCAR/Gabriel-project-signature-tensors/signature-tensors-in-OSCAR/src/TruncatedTensorAlgebra.jl:1967
 [4] sig2parPoly(T::TruncatedTensorAlgebra{QQField}, A::Array{QQFieldElem, 3})
   @ SignatureTensors ~/Documents/Berlin_projects/project_signature-tensors-in-OSCAR/Gabriel-project-signature-tensors/signature-tensors-in-OSCAR/src/TruncatedTensorAlgebra.jl:2216
 [5] sig(T::TruncatedTensorAlgebra{…}, path_type::Symbol; coef::Array{…}, shape::Function, composition::Vector{…}, regularity::Int64, algorithm::Symbol)
   @ SignatureTensors ~/Documents/Berlin_projects/project_signature-tensors-in-OSCAR/Gabriel-project-signature-tensors/signature-tensors-in-OSCAR/src/TruncatedTensorAlgebra.jl:453
 [6] top-level scope
   @ REPL[12]:1
Some type information was truncated. Use `show(err)` to see complete types.

julia> T
TruncatedTensorAlgebra{QQField}(Rational field, 4, 3, :p2id)






using Oscar

function polynomial_path_iis_signature(A::AbstractMatrix, k::Int; R)

    d, m = size(A)

    # Crear PolynomialRing en t sobre R
    T, tvec = polynomial_ring(R, [:t])
    t = tvec[1]   # tomar la variable única
   
    # Construir polinomios X_i(t) en T
    X = Vector{typeof(t)}(undef, d)
    for i in 1:d
        p = zero(T)
        for j in 1:m
            p += A[i,j] * t^j
        end
        X[i] = p
    end

    # Derivadas DX_i(t)
    DX = [derivative(X[i], t) for i in 1:d]
    
    # Vector de arrays para cada nivel
    results = Vector{Array{typeof(zero(R))}}(undef, k)
    
    # Nivel 1
    results[1] = [evaluate(X[i], [one(R)]) for i in 1:d]

    # Niveles 2..k
    for l in 2:k
        words = Iterators.product(ntuple(_->1:d, l)...)
        results[l] = Vector{typeof(zero(R))}(undef, d^l)
        idx = 1
        for w in words
            P = X[w[1]]
            for r in 2:l
                Q = P * DX[w[r]]
                P = integration(Q, t)
            end
            results[l][idx] = evaluate(P, [one(R)])
            idx += 1
        end
    end

    return results
end

function polynomial_path_iis_signature(A::AbstractMatrix, k::Int; R)

    d, m = size(A)

    T, tvec = polynomial_ring(R, [:t])
    t = tvec[1]

    X = Vector{typeof(t)}(undef, d)
    for i in 1:d
        p = zero(T)
        for j in 1:m
            p += A[i,j] * t^j
        end
        X[i] = p
    end

    DX = [derivative(X[i], t) for i in 1:d]

    E = typeof(one(R))

    res = Array{E}(undef, ntuple(_->d, k)...)

    for idx in CartesianIndices(res)

        w = Tuple(idx)

        P = X[w[1]]
        for r in 2:k
            Q = P * DX[w[r]]
            P = integration(Q, t)
        end

        res[idx] = evaluate(P, [one(R)])
    end

    return res
end


function polynomial_path_iis_signature(A::AbstractMatrix, k::Int; R)

    d, m = size(A)

    T, tvec = polynomial_ring(R, [:t])
    t = tvec[1]

    X = Vector{typeof(t)}(undef, d)
    for i in 1:d
        p = zero(T)
        for j in 1:m
            p += A[i,j] * t^j
        end
        X[i] = p
    end

    DX = [derivative(X[i], t) for i in 1:d]

    results = Vector{Array{typeof(one(R))}}(undef, k+1)
    results[1]= fill(one(R), ())    # array 0-dimensional

    for l in 1:k

        res = Array{typeof(one(R))}(undef, ntuple(_->d, l)...)

        for idx in CartesianIndices(res)

            w = Tuple(idx)

            P = X[w[1]]
            for r in 2:l
                Q = P * DX[w[r]]
                P = integration(Q, t)
            end

            res[idx] = evaluate(P, [one(R)])
        end

        results[l+1] = res
    end

    return results
end



A = [
1 0 0 0 0; 
0 1 0 0 0;
0 0 1 0 0;
0 0 0 1 0;
0 0 0 0 1 ]

k = 3
R3=typeof(QQ)




sig2 = polynomial_path_iis_signature(A, k; R=QQ)

T = TruncatedTensorAlgebra(QQ,size(A,1),k)

S = SignatureTensors.sig(T,:poly,coef=QQ.(A)); 

S2 = SignatureTensors.sig(T,:poly,coef=QQ.(A),algorithm=:ARS26); 


d = 4; k = 3; 
R=QQ
T = TruncatedTensorAlgebra(QQ,d,k)
m = 4; # maximal 7 number of segments
A = rand((-20:20),d,m)
QQ.(A)

S = sig(T,:poly,coef=QQ.(A)); 

S3 = SignatureTensors.sig(T,:poly,coef=QQ.(A)); 

S4 = SignatureTensors.sig(T,:poly,coef=QQ.(A),algorithm=:ARS26); 

S5 = SignatureTensors.sig(T,:poly,coef=QQ.(A),algorithm=:ARS26); 

sig3=polynomial_path_iis_signature(A, k; R=QQ)

R_new,_=polynomial_ring(R, [:t])
T, t = Oscar.PolynomialRing(R, "t")



using Oscar, SignatureTensors;         
d = 2; k = 3;                          # d dimension, k truncation level
T = TruncatedTensorAlgebra(QQ,d,k);   
sig(T,:axis) 
sig(T,:poly,coef=QQ.([1 2; 3 4]))

sig(T,:pwln,coef=QQ.([1 2; 3 4]))
sig(T,:pwmon,composition=[1,1])




d = 2; k = 4; m = 4; A = QQ.([ 6  -2  6   -10; 7  -4  10  -4])
S = sig(TruncatedTensorAlgebra(QQ,d,k),:pwln,coef=A);     # in L_{d,<=k,m}
C = sig(TruncatedTensorAlgebra(QQ,m,k),:axis);             # core tensor in T_mk
R, a = polynomial_ring(QQ, :a => (1:d, 1:m)); I = ideal(R,vec(S-a*C));
dim(I), degree(I)
(0, 4)



T = TruncatedTensorAlgebra(QQ,50,3);       # d=50, k=3
S = sig(T,:pwlin,coef=QQrandGL(d));        # QQrandGL(d) random invetible matrix 
@time A = learn(S,algorithm=:Sch25);       # solves (1) with fixed core tensor C 

8.293488 seconds (92.58 M allocations: 3.686 GiB, 21.73% gc time)

QQrandGL(d)




d = 2; k = 4; m1 = 1;m2 = 1; K,v = rational_function_field(QQ, :v => (1:d,1:d));
T = TruncatedTensorAlgebra(K,d,k);
s1 = sig(T,:segment,coef=v[:,1]); s2 = sig(T,:segment,coef=v[:,2]); 
b = bary([s1,s2]); 
R,a = polynomial_ring(K,:a => (1:d,1:3));           # Ansatz B24(1,1) = 3
y = sig(T,:pwln,coef=a); 
I = ideal(R,vec(y-b));
dim(I), degree(I)
(0,2)



d = 3; k = 3; m = 2; n = 2;
R, a = polynomial_ring(QQ, :a => (1:d, 1:m,1:n));
T = TruncatedTensorAlgebra(R,d,k,sequence_type=:p2id);   
S = sig(T,:poly,coef=a);
S[1,1,1]







function sig_rough_word_matrix(path::AbstractMatrix, word::Vector{Int})
    d, m = size(path)
    k = length(word)

    E = eltype(path)
    sig = fill(zero(E), k+1)
    sig[1] = one(E)

    for t in 2:m
        for i in k:-1:1
            delta = path[word[i], t] - path[word[i], t-1]
            sig[i+1] += delta * sig[i]
        end
    end

    return sig[end]
end



function sigRoughPath_TA(
    T::TruncatedTensorAlgebra{R},
    path::AbstractMatrix
) where R

   
    if T.sequence_type != :iis
        error("sigRoughPath_TA only defined for sequence_type = :iis")
    end

    d = base_dimension(T)
    k = truncation_level(T)
    Rbase = base_algebra(T)

    d_path, m = size(path)

    if d_path != d
        error("Path dimension does not match tensor algebra base dimension")
    end

   
    E = typeof(one(Rbase))

   
    elem_out = Vector{Array{E}}(undef, k+1)

   
    elem_out[1] = fill(one(E), ())

     function compute_level(r)
        tensor_dims = ntuple(_ -> d, r)
        Tlevel = fill(zero(E), tensor_dims...)

        function loop_word(word::Vector{Int}, level::Int)
            if level > r
                val = sig_rough_word_matrix(path, word)

                Tlevel[Tuple(word)...] = try
                    convert(E, val)
                catch
                    val
                end
            else
                for i in 1:d
                    word[level] = i
                    loop_word(word, level + 1)
                end
            end
        end

        loop_word(fill(0, r), 1)
        return Tlevel
    end

    for j in 1:k
        elem_out[j+1] = compute_level(j)
    end

    return TruncatedTensorAlgebraElem(T, elem_out)
end



A=As[1]

d = 6; k = 4; number_tests = 5;
T = TruncatedTensorAlgebra(QQ,d,k);
ms = rand((1:7),number_tests); # maximal 7 number of segments
As = [QQ.(rand((-20:20),d,m)) for m in ms];
A=As[1]

d,k,A

for A in As
        m = size(A)[2]  # number of segments in A
        S = sig(T,:pwln,coef=A); 
        @test S == sig(T,:pwln,coef=A,algorithm=:Chen); 
        @test S == sig(T,:pwln,coef=A,algorithm=:congruence);
        @test S == sig(T,:pwln,coef=A,algorithm=:LS); 
        C = sig(TruncatedTensorAlgebra(QQ,m,k),:axis); # axis tensor in dim m
        @test S == A*C
        @test zero(T) + S == S
        @test one(T)*S == S
        @test S*one(T) == S
        @test inv(S)*S == one(T)
        @test S*inv(S) == one(T)
        @test inv(inv(S)) == S
        @test exp(log(S)) == S
        @test S^3 == S*S*S
        @test exp(log(S) + log(S)) == S^2
end 