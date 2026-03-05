using BenchmarkTools
using Statistics 


# change default for `seconds` to 2.5
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 400

BenchmarkTools.DEFAULT_PARAMETERS.samples = 100


function benchmark_signature(ds::Vector{Int}, ms::Vector{Int}, k::Int; 
                             signature_path::Symbol=:pwln,seq_type::Symbol=:iis, num_samples::Int=100)


    results = Dict{Tuple{Int,Int}, Float64}()

    for d in ds
        for m in ms
            T = TruncatedTensorAlgebra(QQ, d, k,seq_type)  

            benchset = @benchmarkable sig(T, signature_path, coef=QQ.(rand(-20:20, d, m))) seconds=400 samples=num_samples 
            bench=run(benchset)
            # times in seconds
            results[(d,m)] = mean(bench).time*1000 
        end
    end

    # Median matrix
    median_matrix = zeros(length(ds), length(ms))
    for ((d, m), times) in results
        i = findfirst(==(d), ds)
        j = findfirst(==(m), ms)
        median_matrix[i, j] = median(times)
    end

    return median_matrix
end

# -------------------------------
# Example
# -------------------------------
#ds = [10, 20, 30, 40, 50, 60]
#ms = [10, 20, 30, 40, 50, 60]
k = 3

ds=[2,3]
ms=[2,3]
#ds 
#ms

median_matrix = benchmark_signature(ds, ms, 3, signature_path=:pwln, seq_type=:iis, num_samples=100)
println(median_matrix)
