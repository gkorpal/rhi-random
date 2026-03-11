# To use with Revise for automatic code reloading:
# julia> using Revise
# julia> includet("forms.jl")
# Now any changes to this file will be automatically tracked and reloaded

using Oscar


"""
    fileReader(filename)

Reads polarization files of the format "polarizations_p<prime>_FINAL.txt" where each line contains rhi_param values.

# Format
The file should have:
- Header comment: # p = <prime>
- Optional comment lines starting with #
- Data lines: rhi_param = [u, v, w, x, y, z]

# Examples
```jldoctest
julia> prime, polarizations = fileReader("polarizations_p23_FINAL.txt");

julia> prime
11

julia> length(polarizations)
5

julia> polarizations[1]
6-element Vector{Int64}:
   8
  28
   8
  -3
  -5
  -6
```
"""
function fileReader(filename::String)
    prime = nothing
    polarizations = Vector{Vector{Int}}()

    # Accept either "# p = 11" or "p = 11" headers.
    prime_re = r"^\s*#?\s*p\s*=\s*(\d+)"
    data_re = r"^\s*\[\d+\]\s+(.+)$"

    for line in eachline(filename)
        raw = strip(line)
        if isempty(raw)
            continue
        end

        m = match(prime_re, raw)
        if m !== nothing
            prime = parse(Int, m.captures[1])
            continue
        end

        if startswith(raw, "#")
            continue
        end

        data = match(data_re, raw)
        if data !== nothing
            raw = data.captures[1]
        end

        nums = split(raw)
        if length(nums) >= 6
            try
                values = [parse(Int, nums[i]) for i in 1:6]
                push!(polarizations, values)
            catch
                continue
            end
        end
    end

    return prime, polarizations
end



"""
    classNumbers(p)

Returns the class numbers h, (h(h+1)/2) and H. 

# Examples
```jldoctest
julia> @time classNumbers(19)
  0.000005 seconds
7
```

# References
[1] K. Hashimoto and T. Ibukiyama, On class numbers of positive definite binary quaternion Hermitian forms, J. Fac. Sci. Univ. Tokyo Sect. IA Math. **27** (1980), no. 3, 549--601; MR0603952
[2] T. Ibukiyama, T. Katsura and F. Oort, Supersingular curves of genus two and class numbers, Compositio Math. **57** (1986), no. 2, 127--152; MR0827350
[3] T. Katsura and F. Oort, Supersingular abelian varieties of dimension two or three and class numbers, *in Algebraic geometry, Sendai, 1985*, 253--281, Adv. Stud. Pure Math., 10, North-Holland, Amsterdam, ; MR0946242
[4] T. Ibukiyama and T. Katsura, On the field of definition of superspecial polarized abelian varieties and type numbers, Compositio Math. **91** (1994), no. 1, 37--46; MR1273924
[5] T. Ibukiyama, Principal polarizations of supersingular abelian surfaces, J. Math. Soc. Japan **72** (2020), no 4, 1161--1180; MR4165927
[6] T. Ibukiyama, Supersingular abelian varieties and quaternion hermitian lattices, in *Theory and Applications of Supersingular Curves and Supersingular Abelian Varieties*, 17--37, RIMS Kokyuroku Bessatsu, B90, Res. Inst. Math. Sci. (RIMS), Kyoto, 2022; MR4521511
"""
function classNumbers(p::Int)
    # Handle small p quickly H - (h(h+1))/2
    if p == 2 || p == 3
        return [1, 1, 0]
    elseif p == 5
        return [2, 1, 1]
    end

    # Compute H
    # Precompute Jacobi values and remainders
    a = (1 - jacobi_symbol(-1, p))
    b = (1 - jacobi_symbol(-2, p))
    c = (1 - jacobi_symbol(-3, p))
    r5 = p % 5
    d = (r5 == 4) ? (4//5) : 0
    H =
        ((p - 1) * (p + 12) * (p + 23))//2880 +
        (a * (2p + 13))//96 +
        (c * (p + 11))//36 +
        (b//8) +
        ((a * c)//12) +
        d

    # Compute h
    # Use a small lookup table for offset in h
    # Only 1, 5, 7, 11 matter; default 0 otherwise
    r12 = p % 12
    offset = if r12 == 1
        0
    elseif r12 == 5 || r12 == 7
        1
    elseif r12 == 11
        2
    else
        0
    end
    h = ((p - r12)//12 + offset)

    # Return result
    return Int(h), Int((h * (h + 1)) ÷ 2), Int(H)
end

"""
    RHI1(p, param)

Given Bp = (-1, -p | Q) and polarization param = [u0, v0, w0, x0, y0, z0] computes the coefficient matrix of the 5-ary refined Humbert invariant which represent 1. 


# Examples
```jldoctest
julia> @time RHI1(131, [5, 11, 4, -3, 0, -1])
```
"""
function RHI1(p::Int, param::Vector{Int})
    u0, v0, w0, x0, y0, z0 = param

    A = zero_matrix(QQ, 6, 6)

    # Fill diagonal
    A[1, 1] = 2 * v0^2
    A[2, 2] = 2 * u0^2
    A[3, 3] = 8 * w0^2 + 8 * w0 * z0 + 2 * z0^2 + 8
    A[4, 4] = 8 * x0^2 + 8 * x0 * y0 + 2 * y0^2 + 8
    A[5, 5] = 2 * x0^2 + 2 * x0 * y0 * p + 2 * x0 * y0 + QQ(1, 2) * y0^2 * p^2 + y0^2 * p + QQ(1, 2) * y0^2 + 2 * p + 2
    A[6, 6] = 2 * w0^2 - 2 * w0 * z0 * p + 2 * w0 * z0 + QQ(1, 2) * z0^2 * p^2 - z0^2 * p + QQ(1, 2) * z0^2 + 2 * p + 2

    # Fill upper triangular off-diagonals
    A[1, 2] = 2 * u0 * v0 - 4
    A[1, 3] = -4 * v0 * w0 - 2 * v0 * z0
    A[1, 4] = -4 * v0 * x0 - 2 * v0 * y0
    A[1, 5] = -2 * v0 * x0 - v0 * y0 * p - v0 * y0
    A[1, 6] = -2 * v0 * w0 + v0 * z0 * p - v0 * z0

    A[2, 3] = -4 * u0 * w0 - 2 * u0 * z0
    A[2, 4] = -4 * u0 * x0 - 2 * u0 * y0
    A[2, 5] = -2 * u0 * x0 - u0 * y0 * p - u0 * y0
    A[2, 6] = -2 * u0 * w0 + u0 * z0 * p - u0 * z0

    A[3, 4] = 8 * w0 * x0 + 4 * w0 * y0 + 4 * x0 * z0 + 2 * y0 * z0
    A[3, 5] = 4 * w0 * x0 + 2 * w0 * y0 * p + 2 * w0 * y0 + 2 * x0 * z0 + y0 * z0 * p + y0 * z0
    A[3, 6] = 4 * w0^2 - 2 * w0 * z0 * p + 4 * w0 * z0 - z0^2 * p + z0^2 + 4

    A[4, 5] = 4 * x0^2 + 2 * x0 * y0 * p + 4 * x0 * y0 + y0^2 * p + y0^2 + 4
    A[4, 6] = 4 * w0 * x0 + 2 * w0 * y0 - 2 * x0 * z0 * p + 2 * x0 * z0 - y0 * z0 * p + y0 * z0

    A[5, 6] = 2 * w0 * x0 + w0 * y0 * p + w0 * y0 - x0 * z0 * p + x0 * z0 - QQ(1, 2) * y0 * z0 * p^2 + QQ(1, 2) * y0 * z0

    # Mirror the upper triangular part to the lower triangular part
    for i = 2:6
        for j = 1:(i-1)
            A[i, j] = A[j, i]
        end
    end

    # Check A is positive semidefinite
    V = quadratic_space(QQ, A)
    D = diagonal(V)
    if !all(>=(0), D)
        println("Not semi-pd")
        return 0
    end
    #return A
    #println("det(A) = ", factor(ZZ(det(A))))

    # Transform A to reduce radical
    #T = Hecke._complete_to_basis(kernel(A; side = :left))
    #AA = T * A * transpose(T)

    #println(A)
    Azz = map_entries(x -> ZZ(x//2), A)
    AA = lll_gram(Azz)

    # Check rank, symmetry, and last row
    if rank(AA) != 5 || !is_symmetric(AA) || any(!iszero, AA[6, :])
        println("non sym")
        return 0
    end

    # Work with top-left 5×5 submatrix
    B = @view AA[1:5, 1:5]
    
    #if det(B) == (2^9) * p^2
        L = integer_lattice(; gram = B)
        if is_positive_definite(L) && minimum(L) == 1
            C = B .* 2
            #println("det(C) = ", factor(ZZ(det(C))))
            return C
        end
    #end
    return 0
end


"""
    RHI2(p, param)

Given Bp = (-1, -p | Q) and polarization param = [u0, v0, w0, x0, y0, z0] computes the coefficient matrix of the 5-ary refined Humbert invariant which represent 1. 


# Examples
```jldoctest
julia> @time RHI2(131, [5, 11, 4, -3, 0, -1])
```
"""
function RHI2(p::Int, param::Vector{Int})
    u0, v0, w0, x0, y0, z0 = param

    A = zero_matrix(QQ, 6, 6)

    # Fill diagonal
    A[1, 1] = 2 * v0^2
    A[2, 2] = 2 * u0^2
    A[3, 3] = 8 * w0^2 + 8 * w0 * z0 + 2 * z0^2 + 8
    A[4, 4] = 8 * x0^2 + 8 * x0 * y0 + 2 * y0^2 + 8
    A[5, 5] = 2 * x0^2 + 2 * x0 * y0 * p + 2 * x0 * y0 + QQ(1, 2) * y0^2 * p^2 + y0^2 * p + QQ(1, 2) * y0^2 + 2 * p + 2
    A[6, 6] = 2 * w0^2 - 2 * w0 * z0 * p + 2 * w0 * z0 + QQ(1, 2) * z0^2 * p^2 - z0^2 * p + QQ(1, 2) * z0^2 + 2 * p + 2

    # Fill upper triangular off-diagonals
    A[1, 2] = 2 * u0 * v0 - 4
    A[1, 3] = -4 * v0 * w0 - 2 * v0 * z0
    A[1, 4] = -4 * v0 * x0 - 2 * v0 * y0
    A[1, 5] = -2 * v0 * x0 - v0 * y0 * p - v0 * y0
    A[1, 6] = -2 * v0 * w0 + v0 * z0 * p - v0 * z0

    A[2, 3] = -4 * u0 * w0 - 2 * u0 * z0
    A[2, 4] = -4 * u0 * x0 - 2 * u0 * y0
    A[2, 5] = -2 * u0 * x0 - u0 * y0 * p - u0 * y0
    A[2, 6] = -2 * u0 * w0 + u0 * z0 * p - u0 * z0

    A[3, 4] = 8 * w0 * x0 + 4 * w0 * y0 + 4 * x0 * z0 + 2 * y0 * z0
    A[3, 5] = 4 * w0 * x0 + 2 * w0 * y0 * p + 2 * w0 * y0 + 2 * x0 * z0 + y0 * z0 * p + y0 * z0
    A[3, 6] = 4 * w0^2 - 2 * w0 * z0 * p + 4 * w0 * z0 - z0^2 * p + z0^2 + 4

    A[4, 5] = 4 * x0^2 + 2 * x0 * y0 * p + 4 * x0 * y0 + y0^2 * p + y0^2 + 4
    A[4, 6] = 4 * w0 * x0 + 2 * w0 * y0 - 2 * x0 * z0 * p + 2 * x0 * z0 - y0 * z0 * p + y0 * z0

    A[5, 6] = 2 * w0 * x0 + w0 * y0 * p + w0 * y0 - x0 * z0 * p + x0 * z0 - QQ(1, 2) * y0 * z0 * p^2 + QQ(1, 2) * y0 * z0

    # Mirror the upper triangular part to the lower triangular part
    for i = 2:6
        for j = 1:(i-1)
            A[i, j] = A[j, i]
        end
    end

    # Check A is positive semidefinite
    V = quadratic_space(QQ, A)
    D = diagonal(V)
    if !all(>=(0), D)
        println("Not semi-pd")
        return 0
    end
    #return A
    #println("det(A) = ", factor(ZZ(det(A))))

    # Transform A to reduce radical
    #T = Hecke._complete_to_basis(kernel(A; side = :left))
    #AA = T * A * transpose(T)

    #println(A)
    Azz = map_entries(x -> ZZ(x//2), A)
    AA = lll_gram(Azz)

    # Check rank, symmetry, and last row
    if rank(AA) != 5 || !is_symmetric(AA) || any(!iszero, AA[6, :])
        println("non sym")
        return 0
    end

    # Work with top-left 5×5 submatrix
    B = @view AA[1:5, 1:5]
    
    #if det(B) == (2^9) * p^2
        L = integer_lattice(; gram = B)
        if is_positive_definite(L) && minimum(L) > 1
            C = B .* 2
            return C
        end
    #end
    return 0
end


"""
    polyForm(M)

Given the coefficient matrix M of a quadratic form, this function computes the polynomial form using f = 1/2 * (X^t * M * X).

# Examples
```jldoctest
julia> M = RHI1(19, [14,31,7,2,4,2]);

julia> @time polyForm(M)
  0.000096 seconds (194 allocations: 7.609 KiB)
x1^2 + 12*x2^2 + 8*x2*x3 + 8*x2*x4 + 8*x2*x5 + 16*x3^2 - 8*x3*x4 + 16*x3*x5 + 24*x4^2 - 16*x4*x5 + 32*x5^2

julia> N = degForm(M)
[6    2    2    2]
[2    8   -2    4]
[2   -2   12   -4]
[2    4   -4   16]

julia> @time polyForm(N)
  0.000075 seconds (146 allocations: 5.695 KiB)
3*x1^2 + 2*x1*x2 + 2*x1*x3 + 2*x1*x4 + 4*x2^2 - 2*x2*x3 + 4*x2*x4 + 6*x3^2 - 4*x3*x4 + 8*x4^2
```
"""
function polyForm(M::ZZMatrix)
    n = number_of_rows(M)
    R, x = polynomial_ring(ZZ, n)
    f = R(0)
    # Use the formula: f = 1/2 * Σᵢ M[i,i]*x[i]^2 + Σ₍ᵢ<ⱼ₎ M[i,j]*x[i]*x[j]
    for i = 1:n
        f += (M[i, i] ÷ 2) * x[i]^2
        for j = (i+1):n
            f += M[i, j] * x[i] * x[j]
        end
    end
    return f
end


"""
    degForm(M)

Given the coefficient matrix M of RHI, it computes the coefficient matrix of the corresponding 4-ary quadratic form called degree form. 

Note that, RHI = X^2 + 4*deg.

# Examples
```jldoctest
julia> M = RHI(19, [14,31,7,2,4,2]);

julia> @time degForm(M)
  0.000013 seconds (64 allocations: 1.438 KiB)
[6    2    2    2]
[2    8   -2    4]
[2   -2   12   -4]
[2    4   -4   16]
```
"""
function degForm(M::ZZMatrix)
    # Check directly without allocating an intermediate matrix.
    if M[1, 1] == 2 && M[1, 2] == 0 && M[1, 3] == 0 && M[1, 4] == 0 && M[1, 5] == 0
        N = @view M[2:5, 2:5]  # Use a view to avoid copying the submatrix.
        return N .÷ 4
    end
    return nothing
end


"""
    allRHI1(p)

Given a prime p congruent to 11 mod 12

# Examples
```jldoctest
julia> @time allRHI1(11)
working with prime 11
saved data for prime 11
  0.668610 seconds (5.93 M allocations: 329.830 MiB, 20.30% gc time)
```
"""
function allRHI1(p::Int)
    start_time = time()
    println("working with prime ", p)

    # Open file once for writing.
    filename = "./RHI1_$(p).txt"

    file = open(filename, "w")
    try
        # Write global header immediately.
        println(file, "p = ", p, "\n")

        h, h2, H = classNumbers(p)
        println(file, "h(h+1)/2 = ", h2)
        println(file, "H = ", H, "\n")

        idx = 0   # Counting unique forms.
        total = 0 # Total RHIs computed.

        unique_forms = Vector{ZZMatrix}()
        pol_count = Dict{Int,Int}()

        prime, params = fileReader("polarizations_p$(p).txt")
        count = length(params)

        if prime == p
            # Create a temporary IOBuffer
            #m_buffer = IOBuffer()
            for param in params
                cmA = RHI1(p, param)  # ZZMatrix coefficient matrix.
                if cmA != 0
                    total += 1
                    is_unique = true
                    for (k, cmB) in enumerate(unique_forms)
                        if cmA == cmB
                            is_unique = false
                            pol_count[k] = get(pol_count, k, 1) + 1
                            break
                        end
                    end
                    if is_unique
                        LA = integer_lattice(gram = cmA .÷ 2)
                        for (k, cmB) in enumerate(unique_forms)
                            LB = integer_lattice(gram = cmB .÷ 2)
                            if is_isometric(LA, LB)
                                is_unique = false
                                pol_count[k] = get(pol_count, k, 1) + 1
                                break
                            end
                        end
                    end

                    if is_unique
                        idx += 1
                        println(file, "Type ", idx)
                        #println("Found no. ", idx)
                        s4 = param[4] >= 0 ? "+" : "-"
                        s5 = param[5] >= 0 ? "+" : "-"
                        s6 = param[6] >= 0 ? "+" : "-"
                        a4 = abs(param[4])
                        a5 = abs(param[5])
                        a6 = abs(param[6])

                        s4b = (-param[4]) >= 0 ? "+" : "-"
                        s5b = (-param[5]) >= 0 ? "+" : "-"
                        s6b = (-param[6]) >= 0 ? "+" : "-"
                        a4b = abs(-param[4])
                        a5b = abs(-param[5])
                        a6b = abs(-param[6])

                        println(
                            file,
                            "θ = [",
                            param[1],
                            "  ",
                            param[3],
                            s4,
                            a4,
                            "β₁",
                            s5,
                            a5,
                            "β₂",
                            s6,
                            a6,
                            "β₃]",
                        )
                        println(
                            file,
                            "    [",
                            param[3],
                            s4b,
                            a4b,
                            "β₁",
                            s5b,
                            a5b,
                            "β₂",
                            s6b,
                            a6b,
                            "β₃  ",
                            param[2],
                            "]",
                        )

                        push!(unique_forms, cmA)

                        q = polyForm(cmA)
                        println(file, "q(ExE,θ) = ", q)

                        cmC = degForm(cmA)
                        deg = polyForm(cmC)
                        println(file, "deg(ExE,θ) = ", deg, "\n")

                    end
                end
            end
        end                
        # Write the chunk for m to the file and flush.
        #write(file, String(take!(m_buffer)))
        #flush(file)

        println(file, "total polarizations checked: ", count)
        println(file, "total RHI's computed: ", total, "\n")
        println(file, "polarization leading to same type: ", pol_count, "\n")

        end_time = time()
        elapsed_time = end_time - start_time
        hours = floor(elapsed_time / 3600)
        minutes = floor((elapsed_time % 3600) / 60)
        seconds = round(elapsed_time % 60)
        println(file, "Total run time: ", hours, " hrs ", minutes, " min ", seconds, " sec")
    finally
        close(file)
    end
    println("saved data for prime ", p)
    return nothing
end


"""
    allRHI2(p)

Given a prime p congruent to 11 mod 12

# Examples
```jldoctest
julia> @time allRHI2(11)
working with prime 11
saved data for prime 11
  0.668610 seconds (5.93 M allocations: 329.830 MiB, 20.30% gc time)
```
"""
function allRHI2(p::Int)
    start_time = time()
    println("working with prime ", p)

    # Open file once for writing.
    filename = "./RHI2_$(p).txt"

    file = open(filename, "w")
    try
        # Write global header immediately.
        println(file, "p = ", p, "\n")

        h, h2, H = classNumbers(p)
        println(file, "h(h+1)/2 = ", h2)
        println(file, "H = ", H, "\n")

        idx = 0   # Counting unique forms.
        total = 0 # Total RHIs computed.

        unique_forms = Vector{ZZMatrix}()
        pol_count = Dict{Int,Int}()

        prime, params = fileReader("p=$(p)_polarizations.txt")
        count = length(params)

        if prime == p
            # Create a temporary IOBuffer
            #m_buffer = IOBuffer()
            for param in params
                cmA = RHI2(p, param)  # ZZMatrix coefficient matrix.
                if cmA != 0
                    total += 1
                    is_unique = true
                    for (k, cmB) in enumerate(unique_forms)
                        if cmA == cmB
                            is_unique = false
                            pol_count[k] = get(pol_count, k, 1) + 1
                            break
                        end
                    end
                    if is_unique
                        LA = integer_lattice(gram = cmA .÷ 2)
                        for (k, cmB) in enumerate(unique_forms)
                            LB = integer_lattice(gram = cmB .÷ 2)
                            if is_isometric(LA, LB)
                                is_unique = false
                                pol_count[k] = get(pol_count, k, 1) + 1
                                break
                            end
                        end
                    end

                    if is_unique
                        idx += 1
                        println(file, "Type ", idx)
                        #println("Found no. ", idx)
                        s4 = param[4] >= 0 ? "+" : "-"
                        s5 = param[5] >= 0 ? "+" : "-"
                        s6 = param[6] >= 0 ? "+" : "-"
                        a4 = abs(param[4])
                        a5 = abs(param[5])
                        a6 = abs(param[6])

                        s4b = (-param[4]) >= 0 ? "+" : "-"
                        s5b = (-param[5]) >= 0 ? "+" : "-"
                        s6b = (-param[6]) >= 0 ? "+" : "-"
                        a4b = abs(-param[4])
                        a5b = abs(-param[5])
                        a6b = abs(-param[6])

                        println(
                            file,
                            "θ = [",
                            param[1],
                            "  ",
                            param[3],
                            s4,
                            a4,
                            "β₁",
                            s5,
                            a5,
                            "β₂",
                            s6,
                            a6,
                            "β₃]",
                        )
                        println(
                            file,
                            "    [",
                            param[3],
                            s4b,
                            a4b,
                            "β₁",
                            s5b,
                            a5b,
                            "β₂",
                            s6b,
                            a6b,
                            "β₃  ",
                            param[2],
                            "]",
                        )

                        push!(unique_forms, cmA)

                        q = polyForm(cmA)
                        println(file, "q(ExE,θ) = ", q, "\n")

                        #cmC = degForm(cmA)
                        #deg = polyForm(cmC)
                        #println(file, "deg(ExE,θ) = ", deg, "\n")

                    end
                end
            end
        end                
        # Write the chunk for m to the file and flush.
        #write(file, String(take!(m_buffer)))
        #flush(file)

        println(file, "total polarizations checked: ", count)
        println(file, "total RHI's computed: ", total, "\n")
        println(file, "polarization leading to same type: ", pol_count, "\n")

        end_time = time()
        elapsed_time = end_time - start_time
        hours = floor(elapsed_time / 3600)
        minutes = floor((elapsed_time % 3600) / 60)
        seconds = round(elapsed_time % 60)
        println(file, "Total run time: ", hours, " hrs ", minutes, " min ", seconds, " sec")
    finally
        close(file)
    end
    println("saved data for prime ", p)
    return nothing
end


# Main execution: Run allRHI for all primes p between 11 and 690 where p ≡ 11 (mod 12)
# for p in 11:12:690
#     if is_prime(p)
#         allRHI(p)
#     end
# end


"""
    fileParser(filename)

Reads in txt file, identifies prime p, and identifies all the degree forms and stores their coefficient matrices.

# Examples
```jldoctest 
julia> prime, forms = fileParser("RHI_19_38_57.txt");

julia> types = sort(collect(keys(forms)));

julia> length(types)
8

julia> q_form, deg_form = forms[6];

julia> println("q(ExE,θ) = ", q_form)
q(E1xE2,θ) = x1^2 + 8*x2^2 + 8*x2*x4 + 12*x3^2 - 8*x3*x4 + 16*x4^2 + 76*x5^2

julia> println("deg(ExE,θ) = ", deg_form)
deg(E1xE2,θ) = 2*x1^2 + 2*x1*x3 + 3*x2^2 - 2*x2*x3 + 4*x3^2 + 19*x4^2
```
"""
function fileParser(filename::String)
    # Define the polynomial ring in 5 variables over ZZ.
    global R, (x1, x2, x3, x4, x5) = polynomial_ring(ZZ, 5)
    PolyElem = typeof(x1^2)

    prime = nothing
    types_data = Dict{Int,NTuple{2,PolyElem}}()

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
            prime = parse(Int, m.captures[1])
            continue
        end

        m = match(type_re, line)
        if m !== nothing
            if current_type != 0
                types_data[current_type] = (current_q, current_deg)
            end
            current_type = parse(Int, m.captures[1])
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
    minDeg(p,ell,a,b)

Retrieve all the degree forms from the file "RHI_p_a_b.txt" and compute minimum vector for the lattice corresponding to the gram matrix. Finally, returns the frequency distribution of all minimum vectors.

# Examples
```jldoctest 
julia> @time minDeg(19, 1, 2*19, 3*19)
  0.607081 seconds (1.49 M allocations: 75.458 MiB, 96.98% compilation time: 15% of which was recompilation)
Dict{Int64, Int64} with 4 entries:
  4 => 1
  2 => 2
  3 => 2
  1 => 3

julia> @time minDeg(19, 1, 2*19, 3*19)
  0.006532 seconds (14.66 k allocations: 759.359 KiB)
Dict{Int64, Int64} with 4 entries:
  4 => 1
  2 => 2
  3 => 2
  1 => 3
```
"""
function minDeg(p::Int)

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

    NDict = Dict{Int,Int}()

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

"""
    primeList(N, n)

Computes a list of n pseudoprimes bigger than N and congruent to 11 mod 12.
    
We take N > 11 so the moduli space of the supersingular surface has more than one element.
"""
function primeList(N::Int, n::Int)
    primes = Vector{Int}()
    p = N
    # Adjust p so that it is congruent to 11 mod 12.
    rem = p % 12
    if rem != 11
        p += (rem - 11) % 12
    end
    while length(primes) < n
        if is_probable_prime(p)
            push!(primes, p)
        end
        p += 12  # Only check numbers congruent to 11 mod 12.
    end
    return primes
end


function getDeg(p::Int, n::Int)
    primes = primeList(p, n)
    last_p = primes[end]
    output_filename = "deg_$(p)_$(last_p).txt"
    open(output_filename, "w") do io
        for p in primes
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
    end
    return nothing
end

# Representatives of 5-ary genus rep 1
function Genus5(p::Int)
    A = ZZ[2 0 0 0 0 ; 0 8 0 0 4; 0 0 8 4 0; 0 0 4 2*(p+1) 0; 0 4 0 0 2*(p+1)]
    LA = integer_lattice(; gram = A .÷2)
    Llist = genus_representatives(LA)
    count = 0
    output_filename = "Gen5_$(p).txt"
    open(output_filename, "w") do io
        println(io, "p = ", p, "\n")
        for L in Llist
            if minimum(L) == 1
                M = gram_matrix(L)
                q = polyForm(ZZ.(M .*2))
                println(io, "q = ", q, "\n")
                count += 1
            end
        end
        println(io, "Total Gen5 forms: ", length(Llist))
        println(io, "Total Gen5 forms rep 1: ", count)
    end
    return nothing
end

# Representatives of 4-ary norm genus
function Genus4(p::Int)
    A = ZZ[2 0 0 1; 0 2 1 0; 0 1 (p+1)//2 0; 1 0 0 (p+1)//2]
    LA = integer_lattice(; gram = A .÷2)
    Llist = genus_representatives(LA)
    output_filename = "Gen4_$(p).txt"
    open(output_filename, "w") do io
        println(io, "p = ", p, "\n")
        for L in Llist
            M = gram_matrix(L)
            q = polyForm(ZZ.(M .*2))
            println(io, "q = ", q, "\n")
        end
        println(io, "Total Gen4 forms: ", length(Llist))
    end
    return nothing
end