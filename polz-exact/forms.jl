using Oscar

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

function degForm(M::ZZMatrix)
    # Check directly without allocating an intermediate matrix.
    if M[1, 1] == 2 && M[1, 2] == 0 && M[1, 3] == 0 && M[1, 4] == 0 && M[1, 5] == 0
        N = @view M[2:5, 2:5]  # Use a view to avoid copying the submatrix.
        return N .÷ 4
    end
    return nothing
end

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


#allRHI1(p)


