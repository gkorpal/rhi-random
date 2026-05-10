# PrimitivePow4kThesis.jl
#
# The code for primitive representation of 4^k by a fixed
# positive definite integral quinary quadratic form Q_A(x)=x' A x.
#
# The order of operations in the algorithm:
#   1. Build the lattice L from A.
#   2. Check primitive 2-adic representability, cached by target residue modulo 2^N.
#   3. Check the odd bad local conditions q | det(A), q odd.
#   4. For the requested finite range of k's, do exact enumaration for global search.
#
# LLL/Cholesky are only used for speed.  Every found witness is verified back in
# the original coordinates and against the original Gram matrix A.

module PrimitivePow4Thesis

using Oscar
using LinearAlgebra

const ZZ = Oscar.ZZ

export Options, Problem,
       build_problem, decide, scan, summary, print_rows,
       qvalue, primitive_gcd, verify_witness

# --------------
# Data helpers
# --------------

Base.@kwdef struct Options
    two_adic_depth::Int = 6              # local test modulo 2^depth
    exact_cutoff::Int = 31              # 4^31=2^62 fits Int search
    assume_after_exact_cutoff::Bool = false
    max_nodes::Int = 300_000
    lll_delta::Float64 = 0.99
    lll_eta::Float64 = 0.501
    precision_bits::Int = 256
    tolerance::Float64 = 1e-20
end

struct Problem
    A::Matrix{Int}                       # original Gram matrix
    L::ZZLat                             # integral lattice with Gram A
    bad_primes::Vector{Int}              # odd primes where local checks are done
    options::Options

    # Local caches
    two_adic_method::String
    jordan_scales::Vector{Int}
    local2_by_target_residue::Dict{Int,Bool}
    odd_local_ok::Bool

    # Exact-search caches
    U::Matrix{Int}                       # x = U*z
    G::Matrix{Int}                       # G = U' A U
    R::Matrix{BigFloat}                  # Cholesky factor of G
end

pow4(k::Integer) = BigInt(1) << (2 * Int(k))

# ---------
# Functions
# ---------

function qvalue(A::AbstractMatrix{<:Integer}, x::AbstractVector{<:Integer})
    n = length(x)
    size(A) == (n, n) || error("dimension mismatch between A and x")
    s = BigInt(0)
    @inbounds for i in 1:n, j in 1:n
        s += BigInt(A[i, j]) * BigInt(x[i]) * BigInt(x[j])
    end
    return s
end

function primitive_gcd(x::AbstractVector{<:Integer})
    g = BigInt(0)
    for a in x
        g = gcd(g, abs(BigInt(a)))
        g == 1 && return BigInt(1)
    end
    return g
end

isprimitive(x::AbstractVector{<:Integer}) = primitive_gcd(x) == 1

function verify_witness(A::AbstractMatrix{<:Integer}, x, k::Integer)
    target = pow4(k)
    return (value = qvalue(A, x), target = target,
            primitive_gcd = primitive_gcd(x),
            ok = qvalue(A, x) == target && isprimitive(x))
end

function _check_quinary_form(Ain::AbstractMatrix{<:Integer})
    A = Matrix{Int}(Ain)
    size(A) == (5, 5) || error("expected a 5x5 Gram matrix")
    A == transpose(A) || error("Gram matrix must be symmetric")
    L = integer_lattice(; gram = matrix(ZZ, A))
    is_positive_definite(L) || error("Gram matrix must be positive definite")
    return A, L
end

function _odd_primes_from_det(A::Matrix{Int})
    d = abs(det(matrix(ZZ, A)))
    ps = Int[]
    for (p, _) in factor(d)
        q = Int(p)
        q > 2 && push!(ps, q)
    end
    return ps
end

# ------------------
# Local checks at 2
# ------------------

# There are only finitely many target residues modulo 2^N.
# This is only a finite-modulus cache, not a mathematical classification.
function _target_residue_mod_2N(k::Int, depth::Int)
    M = BigInt(1) << depth
    return Int(mod(pow4(k), M))
end

function _rat_mod_pow2(x, M::Int)
    if x isa QQFieldElem
        num = Int(ZZ(numerator(x)))
        den = Int(ZZ(denominator(x)))
        isodd(den) || error("even denominator in 2-adic Gram entry: $x")
        return mod(num * invmod(mod(den, M), M), M)
    else
        return mod(Int(ZZ(x)), M)
    end
end

function _primitive_residues_mod_2N(G, depth::Int)
    M = 1 << depth
    r = nrows(G)
    r == 0 && return Set{Int}()

    Gmod = Matrix{Int}(undef, r, r)
    for i in 1:r, j in 1:r
        Gmod[i, j] = _rat_mod_pow2(G[i, j], M)
    end

    residues = Set{Int}()
    ranges = ntuple(_ -> 0:(M - 1), r)
    for xs in Iterators.product(ranges...)
        all_even = true
        for t in xs
            if isodd(t)
                all_even = false
                break
            end
        end
        all_even && continue

        q = 0
        @inbounds for i in 1:r, j in 1:r
            q += xs[i] * Gmod[i, j] * xs[j]
        end
        push!(residues, mod(q, M))
    end
    return residues
end

function _local2_by_jordan(L::ZZLat, k::Int, depth::Int)
    M = 1 << depth
    target = Int(mod(pow4(k), M))
    _, blocks, scales = jordan_decomposition(L, ZZ(2))

    # Combine primitive residues block by block.  A vector in the whole lattice
    # is primitive as soon as one Jordan block contributes a primitive vector.
    nonprimitive = Set([0])
    primitive = Set{Int}()

    for G in blocks
        block_prim = _primitive_residues_mod_2N(G, depth)
        next_nonprimitive = Set{Int}()
        next_primitive = Set{Int}()

        for a in nonprimitive
            push!(next_nonprimitive, a)
            for b in block_prim
                push!(next_primitive, mod(a + b, M))
            end
        end
        for a in primitive
            push!(next_primitive, a)
            for b in block_prim
                push!(next_primitive, mod(a + b, M))
            end
        end

        nonprimitive = next_nonprimitive
        primitive = next_primitive
    end

    return (ok = target in primitive,
            method = "jordan_decomposition(L, ZZ(2)) modulo 2^$depth",
            scales = [Int(scales[i]) for i in 1:length(scales)])
end

function _local2_direct_fallback(A::Matrix{Int}, k::Int, depth::Int)
    # Fallback only.  Keep it small, because direct 64^5 enumeration is too big.
    d = min(depth, 4)
    M = 1 << d
    target = Int(mod(pow4(k), M))
    A2 = mod.(A, M)

    for x1 in 0:M-1, x2 in 0:M-1, x3 in 0:M-1, x4 in 0:M-1, x5 in 0:M-1
        iseven(x1) && iseven(x2) && iseven(x3) && iseven(x4) && iseven(x5) && continue
        x = (x1, x2, x3, x4, x5)
        q = 0
        @inbounds for i in 1:5, j in 1:5
            q += x[i] * A2[i, j] * x[j]
        end
        mod(q, M) == target && return (ok = true,
            method = "direct primitive search modulo 2^$d", scales = Int[])
    end
    return (ok = false, method = "direct primitive search modulo 2^$d", scales = Int[])
end

function _local2(A::Matrix{Int}, L::ZZLat, k::Int, depth::Int)
    try
        return _local2_by_jordan(L, k, depth)
    catch
        return _local2_direct_fallback(A, k, depth)
    end
end

function _build_local2_cache(A::Matrix{Int}, L::ZZLat, depth::Int)
    # It is enough to test k = 0, ..., ceil(depth/2): after that
    # 4^k is 0 modulo 2^depth.  The cache key is the residue itself.
    last_k_needed = cld(depth, 2)
    cache = Dict{Int,Bool}()
    method = ""
    scales = Int[]

    for k in 0:last_k_needed
        residue = _target_residue_mod_2N(k, depth)
        r = _local2(A, L, k, depth)
        cache[residue] = r.ok
        isempty(method) && (method = r.method)
        isempty(scales) && (scales = r.scales)
    end
    return cache, method, scales
end

# ---------------
# Odd primes dividing discriminant
# --------------

function _odd_local_unit_test(L::ZZLat, bad_primes::Vector{Int})
    # Since 4^k is a q-adic unit for odd q, this check is independent of k.
    unit_lattice = integer_lattice(; gram = matrix(ZZ, 1, 1, [ZZ(1)]))
    for q in bad_primes
        q > 2 || continue
        represents(genus(L, q), genus(unit_lattice, q)) || return false
    end
    return true
end

# -----------------------------------------------------------------------------
# LLL and exact global search
# -----------------------------------------------------------------------------

function _gram_schmidt_columns(B::Matrix{Float64})
    n = size(B, 2)
    Bstar = zeros(Float64, size(B, 1), n)
    mu = zeros(Float64, n, n)
    sqnorm = zeros(Float64, n)
    for i in 1:n
        v = copy(B[:, i])
        for j in 1:i-1
            mu[i, j] = dot(B[:, i], Bstar[:, j]) / sqnorm[j]
            v .-= mu[i, j] .* Bstar[:, j]
        end
        Bstar[:, i] .= v
        sqnorm[i] = dot(v, v)
    end
    return mu, sqnorm
end

function _lll_reduce_form(A::Matrix{Int}, opt::Options)
    n = size(A, 1)
    B = Matrix(transpose(cholesky(Symmetric(Float64.(A))).U))
    U = Matrix{Int}(I, n, n)

    k = 2
    while k <= n
        mu, _ = _gram_schmidt_columns(B)
        for j in k-1:-1:1
            c = round(Int, mu[k, j])
            if abs(mu[k, j]) > opt.lll_eta && c != 0
                B[:, k] .-= c .* B[:, j]
                U[:, k] .-= c .* U[:, j]
            end
        end

        mu, sqnorm = _gram_schmidt_columns(B)
        if sqnorm[k] >= (opt.lll_delta - mu[k, k-1]^2) * sqnorm[k-1]
            k += 1
        else
            B[:, [k-1, k]] = B[:, [k, k-1]]
            U[:, [k-1, k]] = U[:, [k, k-1]]
            k = max(k - 1, 2)
        end
    end

    G = transpose(U) * A * U
    R = setprecision(BigFloat, opt.precision_bits) do
        Matrix(cholesky(Symmetric(BigFloat.(G))).U)
    end
    return U, G, R
end

function _exact_search(P::Problem, k::Int)
    target = pow4(k)
    target <= BigInt(typemax(Int)) || return (
        decision = :UNKNOWN_TARGET_TOO_LARGE_FOR_INT_SEARCH,
        witness = nothing, value = nothing, gcd = nothing, nodes = 0)

    n = size(P.A, 1)
    z = zeros(Int, n)
    nodes = Ref(0)
    found = Ref{Union{Nothing,Vector{Int}}}(nothing)

    result = setprecision(BigFloat, P.options.precision_bits) do
        R = P.R
        T = BigFloat(target)
        eps = BigFloat(P.options.tolerance)

        function search(i::Int, partial::BigFloat)
            nodes[] += 1
            nodes[] > P.options.max_nodes && return :TOO_MANY_NODES
            partial > T + eps && return :CONTINUE

            if i == 0
                if qvalue(P.G, z) == target && isprimitive(z)
                    x = P.U * z
                    if isprimitive(x) && qvalue(P.A, x) == target
                        found[] = copy(x)
                        return :FOUND
                    end
                end
                return :CONTINUE
            end

            tail = BigFloat(0)
            @inbounds for j in i+1:n
                tail += R[i, j] * BigFloat(z[j])
            end

            remaining = T - partial
            remaining < -eps && return :CONTINUE
            radius = sqrt(max(remaining, BigFloat(0))) / abs(R[i, i])
            center = -tail / R[i, i]
            lo = ceil(Int, center - radius - eps)
            hi = floor(Int, center + radius + eps)

            for zi in lo:hi
                z[i] = zi
                term = R[i, i] * BigFloat(zi) + tail
                r = search(i - 1, partial + term * term)
                r == :FOUND && return :FOUND
                r == :TOO_MANY_NODES && return :TOO_MANY_NODES
            end
            z[i] = 0
            return :CONTINUE
        end

        search(n, BigFloat(0))
    end

    if result == :FOUND
        x = found[]
        return (decision = :YES_EXACT, witness = x,
                value = qvalue(P.A, x), gcd = primitive_gcd(x), nodes = nodes[])
    elseif result == :TOO_MANY_NODES
        return (decision = :UNKNOWN_TOO_MANY_NODES, witness = nothing,
                value = nothing, gcd = nothing, nodes = nodes[])
    else
        return (decision = :NO_EXACT, witness = nothing,
                value = nothing, gcd = nothing, nodes = nodes[])
    end
end

# -----------------------------
# Main function to check the primitive representation of 4^k
# --------------------------------

function build_problem(Ain::AbstractMatrix{<:Integer};
                       bad_primes::Union{Nothing,Vector{Int}} = nothing,
                       options::Options = Options())
    A, L = _check_quinary_form(Ain)
    bad = bad_primes === nothing ? _odd_primes_from_det(A) : copy(bad_primes)
    bad = sort(unique([q for q in bad if q > 2]))

    local2_cache, local2_method, scales =
        _build_local2_cache(A, L, options.two_adic_depth)
    odd_ok = _odd_local_unit_test(L, bad)
    U, G, R = _lll_reduce_form(A, options)

    return Problem(A, L, bad, options,
                   local2_method, scales, local2_cache, odd_ok,
                   U, G, R)
end

function decide(P::Problem, k::Integer)
    kk = Int(k)
    kk >= 0 || error("k must be nonnegative")
    target = pow4(kk)
    residue2 = _target_residue_mod_2N(kk, P.options.two_adic_depth)
    local2 = get(P.local2_by_target_residue, residue2, false)

    base = (k = kk, target = target,
            local2 = local2, local_odd = P.odd_local_ok,
            two_adic_target_residue = residue2,
            two_adic_modulus = BigInt(1) << P.options.two_adic_depth,
            two_adic_method = P.two_adic_method,
            jordan_scales = P.jordan_scales,
            bad_primes = P.bad_primes)

    if !local2
        return merge(base, (decision = :NO_LOCAL_AT_2,
            witness = nothing, value = nothing, gcd = nothing, nodes = 0,
            method = "local obstruction"))
    end

    if !P.odd_local_ok
        return merge(base, (decision = :NO_LOCAL_AT_ODD_BAD_PRIME,
            witness = nothing, value = nothing, gcd = nothing, nodes = 0,
            method = "local obstruction"))
    end

    if kk > P.options.exact_cutoff
        decision = P.options.assume_after_exact_cutoff ?
            :YES_ASSUMED_AFTER_EXACT_CUTOFF : :UNKNOWN_ABOVE_EXACT_CUTOFF
        return merge(base, (decision = decision,
            witness = nothing, value = nothing, gcd = nothing, nodes = 0,
            method = "beyond exact cutoff"))
    end

    ex = _exact_search(P, kk)
    return merge(base, (decision = ex.decision,
        witness = ex.witness, value = ex.value, gcd = ex.gcd,
        nodes = ex.nodes, method = "LLL + BigFloat exact ellipsoid search"))
end

scan(P::Problem, ks) = [decide(P, k) for k in ks]

function summary(rows)
    yes_exact = Int[]
    yes_assumed = Int[]
    no = Int[]
    unknown = Int[]

    for r in rows
        if r.decision == :YES_EXACT
            push!(yes_exact, r.k)
        elseif r.decision == :YES_ASSUMED_AFTER_EXACT_CUTOFF
            push!(yes_assumed, r.k)
        elseif r.decision in (:NO_EXACT, :NO_LOCAL_AT_2, :NO_LOCAL_AT_ODD_BAD_PRIME)
            push!(no, r.k)
        else
            push!(unknown, r.k)
        end
    end

    return (yes_exact_ks = yes_exact,
            yes_assumed_ks = yes_assumed,
            no_ks = no,
            unknown_ks = unknown,
            count_yes_exact = length(yes_exact),
            count_yes_assumed = length(yes_assumed),
            count_no = length(no),
            count_unknown = length(unknown))
end

function print_rows(rows)
    for r in rows
        println("k=", r.k,
                " target=", r.target,
                " local2=", r.local2,
                " odd=", r.local_odd,
                " decision=", r.decision,
                " nodes=", r.nodes,
                " witness=", r.witness)
    end
    return nothing
end

end # module PrimitivePow4kThesis
