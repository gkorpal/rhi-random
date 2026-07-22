using Oscar

"""
    fileParser(filename)

Reads in txt file, identifies prime p, and identifies all the degree forms and stores their coefficient matrices.

"""
function fileParser(filename::String)
    # Define the polynomial ring in 5 variables over ZZ.
    global R, (x1, x2, x3, x4, x5) = polynomial_ring(ZZ, 5)
    PolyElem = typeof(x1^2)

    prime = nothing
    types_data = Dict{Integer,NTuple{2,PolyElem}}()

    current_type = 0
    current_q = nothing
    current_deg = nothing

    # regex patterns.
    prime_re = r"^\s*p\s*=\s*(\d+)"
    type_re = r"^\s*Type\s+(\d+)"
    q_re = r"^\s*q\(ExE,θ\)\s*=\s*(.+)"
    deg_re = r"^\s*deg\(ExE,θ\)\s*=\s*(.+)"

    for line in eachline(filename)
        # Check for the prime line.
        m = match(prime_re, line)
        if m !== nothing
            prime = parse(BigInt, m.captures[1])
            continue
        end

        m = match(type_re, line)
        if m !== nothing
            if current_type != 0
                types_data[current_type] = (current_q, current_deg)
            end
            current_type = parse(BigInt, m.captures[1])
            current_q = nothing
            current_deg = nothing
            continue
        end

        if current_type != 0
            m = match(q_re, line)
            if m !== nothing
                q_str = strip(m.captures[1])
                current_q = eval(Meta.parse(q_str))
                continue
            end

            m = match(deg_re, line)
            if m !== nothing
                deg_str = strip(m.captures[1])
                current_deg = eval(Meta.parse(deg_str))
                continue
            end
        end
    end

    if current_type != 0
        types_data[current_type] = (current_q, current_deg)
    end

    return prime, types_data
end


"""
    minDeg(p)

Retrieve all the degree forms from the file "RHI1_p.txt" and compute minimum vector for the lattice corresponding to the gram matrix. Finally, returns the frequency distribution of all minimum vectors.

"""
function minDeg(p::Integer)

    monos = [
        [2, 0, 0, 0, 0],
        [1, 1, 0, 0, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 0, 0, 1],
        [0, 2, 0, 0, 0],
        [0, 1, 1, 0, 0],
        [0, 1, 0, 1, 0],
        [0, 1, 0, 0, 1],
        [0, 0, 2, 0, 0],
        [0, 0, 1, 1, 0],
        [0, 0, 1, 0, 1],
        [0, 0, 0, 2, 0],
        [0, 0, 0, 1, 1],
        [0, 0, 0, 0, 2],
    ]

    filename = "RHI1_$(p).txt"
    prime, forms = fileParser(filename)
    @assert prime == p "File mismatch"

    NDict = Dict{Integer,Integer}()

    for tp in keys(forms)
        _, q = forms[tp]
        coeffs = [coeff(q, mono) for mono in monos] #15 elements

        A = ZZ[
            2*coeffs[1] coeffs[2] coeffs[3] coeffs[4];
            coeffs[2] 2*coeffs[6] coeffs[7] coeffs[8];
            coeffs[3] coeffs[7] 2*coeffs[10] coeffs[11];
            coeffs[4] coeffs[8] coeffs[11] 2*coeffs[13]
        ]

        L = integer_lattice(gram = A .÷ 2)
        N = Int(minimum(L))
        NDict[N] = get(NDict, N, 0) + 1
    end
    return NDict
end


function getDeg(p::Integer)
    output_filename = "deg_$(p).txt"
    open(output_filename, "w") do io
        println(io, "p = ", p)
        try
            NDict = minDeg(p)
            println(io, "max(min deg) = ", maximum(keys(NDict)))
            println(io, "minimum degree frequency distribution")
            for (n, freq) in sort(collect(NDict), by = x -> x[2], rev = true)
                println(io, n, " => ", freq)
            end
            println(io, "")   
        catch e 
            println(io, "skipping $p since file doesn't exist.\n")
        end
    end
    return nothing
end


# uncomment to run the analysis
for p in 11:12:660
    if is_probable_prime(p)
        getDeg(p)
    end
end