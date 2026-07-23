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

    # regex patterns.
    prime_re = r"^\s*p\s*=\s*(\d+)"
    # Plain "q = ..." lines as in Gen5_p.txt, e.g.:
    #   q = x1^2 + 4*x2^2 + 4*x2*x3 + 12*x3^2 + ...
    # Each such line is its own form; there is no separate "Type"/"deg" line,
    # so the degree form is recovered as (q - x1^2)/4.
    q_re = r"^\s*q\s*=\s*(.+)$"

    for line in eachline(filename)
        # Check for the prime line.
        m = match(prime_re, line)
        if m !== nothing
            prime = parse(BigInt, m.captures[1])
            continue
        end

        # Each "q = ..." line is a complete, standalone form.
        m = match(q_re, line)
        if m !== nothing
            q_str = strip(m.captures[1])
            q_poly = eval(Meta.parse(q_str))
            deg_poly = divexact(q_poly - x1^2, 4)
            current_type += 1
            types_data[current_type] = (q_poly, deg_poly)
            continue
        end
    end

    return prime, types_data
end


"""
    minDeg(p)

Retrieve all the degree forms from the file "Gen5_p.txt" and compute minimum vector for the lattice corresponding to the gram matrix. Finally, returns the frequency distribution of all minimum vectors.

"""
function minDeg(p::Integer)

    # The degree forms produced by fileParser live in x2, x3, x4, x5
    # (x1 is absent, since deg = (q - x1^2)/4). These monomials pick out
    # the coefficients of the resulting rank-4 form in that basis.
    monos = [
        [0, 2, 0, 0, 0],  # x2^2
        [0, 1, 1, 0, 0],  # x2*x3
        [0, 1, 0, 1, 0],  # x2*x4
        [0, 1, 0, 0, 1],  # x2*x5
        [0, 0, 2, 0, 0],  # x3^2
        [0, 0, 1, 1, 0],  # x3*x4
        [0, 0, 1, 0, 1],  # x3*x5
        [0, 0, 0, 2, 0],  # x4^2
        [0, 0, 0, 1, 1],  # x4*x5
        [0, 0, 0, 0, 2],  # x5^2
    ]

    filename = "Gen5_$(p).txt"
    prime, forms = fileParser(filename)
    @assert prime == p "File mismatch"

    NDict = Dict{Integer,Integer}()

    for tp in keys(forms)
        _, deg = forms[tp]
        coeffs = [coeff(deg, mono) for mono in monos] #10 elements

        A = ZZ[
            2*coeffs[1] coeffs[2] coeffs[3] coeffs[4];
            coeffs[2] 2*coeffs[5] coeffs[6] coeffs[7];
            coeffs[3] coeffs[6] 2*coeffs[8] coeffs[9];
            coeffs[4] coeffs[7] coeffs[9] 2*coeffs[10]
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
            if e isa SystemError || e isa Base.IOError
                println(io, "skipping $p since file doesn't exist.\n")
            else
                println(io, "skipping $p due to error: ", sprint(showerror, e), "\n")
            end
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