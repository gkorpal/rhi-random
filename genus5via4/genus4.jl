using Oscar

# ============================================================================
#  Fast computation of quinary-genus representatives with minimum 1
#
#  q_5(t_0,...,t_4) = t_0^2 + 4( t_1^2 + t_1 t_4 + t_2^2 + t_2 t_3
#                                 + (p+1)/4 t_3^2 + (p+1)/4 t_4^2 )
#
#  MATHEMATICAL BACKGROUND (why this is fast)
#  --------------------------------------------------------------------------
#  Write  ω = (1+√-p)/2.  Since p ≡ 3 (mod 4), O_K = Z[ω] is the maximal
#  order of K = Q(√-p), of discriminant  -p, and its norm form is exactly
#
#         b_0(x,y) = N(x + yω) = x^2 + xy + (p+1)/4 y^2.
#
#  So q_5 decomposes as an ORTHOGONAL SUM
#
#         q_5  =  <1>  ⊥  4 b_0(t_1,t_4)  ⊥  4 b_0(t_2,t_3)
#              =  <1>  ⊥  Q_4         where   Q_4 := 4 b_0 ⊥ 4 b_0.
#
#  TRICK (rigorous, unconditional).
#  If v is a vector of norm 1 in an integral lattice L, then v is
#  automatically primitive and splits off orthogonally:
#      w ↦ w - B(v,w) v   maps L → v^⊥ ∩ L integrally  (since B(v,w) ∈ Z),
#  so   L ≅ <1> ⊥ (v^⊥ ∩ L).
#  Hence
#      { L ∈ gen(q_5) : min(L) = 1 }  ↔  gen(Q_4)          (bijectively!)
#      via                      L  =  <1> ⊥ L_4,   L_4 ∈ gen(Q_4).
#  We therefore never need to touch rank 5, and never need to filter by
#  minimum: EVERY class of gen(Q_4) contributes exactly one class of
#  gen(q_5) with minimum 1, and there are no others.
# ============================================================================


"""
    polyForm(M::ZZMatrix)

Given the *doubled* Gram matrix M of a quadratic form (i.e. M[i,i] = 2*Q(e_i),
M[i,j] = 2*B(e_i,e_j) for i<j... see convention below), return the
polynomial  f = 1/2 * X^T M X  as an element of Z[x_1,...,x_n].

Convention used throughout this file: if G = gram_matrix(L) is the *true*
Gram matrix of a lattice (G[i,i] = Q(e_i), no factor of 2), then the
associated polynomial is obtained via polyForm(2*G).
"""
function polyForm(M::ZZMatrix)
    n = number_of_rows(M)
    R, x = polynomial_ring(ZZ, n)
    f = R(0)
    for i = 1:n
        f += (M[i, i] ÷ 2) * x[i]^2
        for j = (i+1):n
            f += M[i, j] * x[i] * x[j]
        end
    end
    return f
end


"""
    doubled_binary_block(a::Int, b::Int, c::Int, scale::Int) -> ZZMatrix

Doubled Gram matrix (see `polyForm` convention) of the quadratic form
scale * (a x^2 + b x y + c y^2), i.e. of the binary form (a,b,c) scaled by
`scale`.  We use scale = 4 throughout, to match the "4 b_0" normalisation
appearing in q_5.
"""
function doubled_binary_block(a::Int, b::Int, c::Int, scale::Int)
    return ZZ[2*scale*a  scale*b ; scale*b  2*scale*c]
end


"""
    quaternary_form_doubled(p::Int) -> ZZMatrix

The doubled Gram matrix of Q_4 = 4 b_0 ⊥ 4 b_0 (the "obvious" member of
its own genus, obtained from the principal form (1,1,(p+1)/4) taken
twice). This is the lattice whose genus representatives we enumerate.
"""
function quaternary_form_doubled(p::Int)
    b0 = (1, 1, div(p + 1, 4))
    B = doubled_binary_block(b0..., 4)
    return block_diagonal_matrix([B, B])
end


"""
    Genus5(p::Int; verbose::Bool=true)

Compute representatives of the genus of q_5 with minimum 1, for a prime
p ≡ 11 (mod 12). Writes them (as quinary polynomials) to Gen5_\$(p).txt.

Via Trick: only a RANK-4 genus (that of Q_4) is ever computed/enumerated
-- no rank-5 neighbour search, and no minimum()/shortest-vector filtering
is ever required.
"""
function Genus5(p::Int; verbose::Bool=true)
    verbose && println("p = $p")

    L4 = integer_lattice(; gram = quaternary_form_doubled(p) .÷ 2)
    reps4 = genus_representatives(L4)

    output_filename = "Gen5_$(p).txt"
    open(output_filename, "w") do io
        println(io, "p = ", p)
        println(io, "method: full rank-4 neighbour search\n")
        for L in reps4
            M = gram_matrix(L)                       # true (undoubled) 4x4 Gram
            Mfull = block_diagonal_matrix([ZZ[2;;], ZZ.(M .* 2)])  # prepend <1>, doubled
            q = polyForm(Mfull)
            println(io, "q = ", q, "\n")
        end
        println(io, "Total Gen5 forms rep 1: ", length(reps4))
    end
    verbose && println("  -> wrote $(length(reps4)) forms to $output_filename\n")
    return nothing
end


# ----------------------------------------------------------------------------
# Driver: same range as the original script.
# ----------------------------------------------------------------------------
# for p in 11:12:660
#     if is_probable_prime(p)
#         Genus5(p; verbose=true)
#     end
# end

# ---------------------------------------------------------------------------
# Command-line interface
# Usage:
#     julia genus.jl 503
# ---------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) != 1
        println("Usage: julia genus.jl <prime>")
        exit(1)
    end

    p = parse(Int, ARGS[1])

    if p % 12 != 11
        error("Expected a prime congruent to 11 mod 12.")
    end

    Genus5(p; verbose=true)
end
