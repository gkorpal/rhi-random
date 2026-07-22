import Pkg
try
    import Oscar
catch
    Pkg.instantiate()
    import Oscar
end

module SuperspecialPolarizations_ThesisOscar

using LinearAlgebra
using Dates
import Oscar
# Hecke is provided inside Oscar; alias it so existing code can use `Hecke`.
const Hecke = Oscar.Hecke

export QElt, Hermitian2x2,
       one_quaternion, qi, qj, qk,
       order_basis_O, algebra_basis_E, O2_basis8,
       qconj, qmul, trd, nrd,
       order_coordinates, has_integral_order_coordinates,
       is_principal, normalize_hermitian,
       enumerate_principal_candidates,
       canonical_packet_forms,
       classNumbers, targetH, check_p_11mod12,
       clear_caches!, write_polarizations_file,
       certify_basis_conventions,
       reduce_packet_by_lll, hecke_packet_isometric, hecke_f1_isometric,
       classify_prime_lll_until_target, run_batch_lll_until_target
# ---------------------------------------------------------------
# Fixed quaternion/order conventions
# ----------------------------------
# B = (-1,-p)/Q, with i^2=-1, j^2=-p, k=ij=-ji.
# QElt(A0,A1,A2,A3) stores (A0 + A1*i + A2*j + A3*k)/2.
# Fixed maximal order basis:
#   O = Z< 1, i, (i+j)/2, (1+k)/2 >.
# Packet basis on O^2:
#   first 4 vectors on first coordinate, next 4 on second coordinate.
# Fixed algebra basis:
#   1, i, j, k.
# canonical_packet_forms(H,p) returns [F1,F2,F3,F4] in the fixed O-basis.
# -----------------------------------------------------------------------

struct QElt
    A0::Int
    A1::Int
    A2::Int
    A3::Int
end

Base.show(io::IO, q::QElt) = print(io, "(", q.A0, ", ", q.A1, ", ", q.A2, ", ", q.A3, ")")
Base.:(==)(x::QElt, y::QElt) = (x.A0 == y.A0 && x.A1 == y.A1 && x.A2 == y.A2 && x.A3 == y.A3)
Base.hash(x::QElt, h::UInt) = hash((x.A0, x.A1, x.A2, x.A3), h)
Base.:+(x::QElt, y::QElt) = QElt(x.A0 + y.A0, x.A1 + y.A1, x.A2 + y.A2, x.A3 + y.A3)
Base.:-(x::QElt, y::QElt) = QElt(x.A0 - y.A0, x.A1 - y.A1, x.A2 - y.A2, x.A3 - y.A3)
Base.:-(x::QElt) = QElt(-x.A0, -x.A1, -x.A2, -x.A3)

const ZEROQ = QElt(0, 0, 0, 0)

one_quaternion() = QElt(2, 0, 0, 0)
qi() = QElt(0, 2, 0, 0)
qj() = QElt(0, 0, 2, 0)
qk() = QElt(0, 0, 0, 2)

qconj(x::QElt) = QElt(x.A0, -x.A1, -x.A2, -x.A3)
trd(x::QElt) = x.A0
_qkey(x::QElt) = (x.A0, x.A1, x.A2, x.A3)

function _halfdiv(n::Int)
    iseven(n) || error("nonintegral quaternion coordinate: numerator=$n")
    return n ÷ 2
end

function nrd(x::QElt, p::Int)
    num = x.A0*x.A0 + x.A1*x.A1 + p*(x.A2*x.A2 + x.A3*x.A3)
    num % 4 == 0 || error("nonintegral reduced norm numerator=$num")
    return num ÷ 4
end

function qmul(x::QElt, y::QElt, p::Int)
    a0,a1,a2,a3 = x.A0, x.A1, x.A2, x.A3
    b0,b1,b2,b3 = y.A0, y.A1, y.A2, y.A3
    return QElt(
        _halfdiv(a0*b0 - a1*b1 - p*a2*b2 - p*a3*b3),
        _halfdiv(a0*b1 + a1*b0 + p*(a2*b3 - a3*b2)),
        _halfdiv(a0*b2 + a2*b0 - a1*b3 + a3*b1),
        _halfdiv(a0*b3 + a3*b0 + a1*b2 - a2*b1),
    )
end

function order_basis_O()
    return QElt[
        QElt(2,0,0,0),       # 1
        QElt(0,2,0,0),       # i
        QElt(0,1,1,0),       # (i+j)/2
        QElt(1,0,0,1),       # (1+k)/2
    ]
end

function algebra_basis_E()
    return QElt[
        QElt(2,0,0,0),       # 1
        QElt(0,2,0,0),       # i
        QElt(0,0,2,0),       # j
        QElt(0,0,0,2),       # k
    ]
end

function O2_basis8()
    B = order_basis_O()
    V = NTuple{2,QElt}[]
    for b in B
        push!(V, (b, ZEROQ))
    end
    for b in B
        push!(V, (ZEROQ, b))
    end
    return V
end

has_integral_order_coordinates(q::QElt) = iseven(q.A0 - q.A3) && iseven(q.A1 - q.A2)

function order_coordinates(q::QElt)
    has_integral_order_coordinates(q) || error("quaternion does not lie in the fixed maximal order: $q")
    z = q.A3
    y = q.A2
    w = (q.A0 - z) ÷ 2
    x = (q.A1 - y) ÷ 2
    return (w, x, y, z)
end

function certify_basis_conventions(p::Int=11)
    B = order_basis_O()
    for a in B, b in B
        has_integral_order_coordinates(qmul(a,b,p)) || return false
    end
    i, h, g = qi(), QElt(0,1,1,0), qj()
    qmul(i,i,p) == QElt(-2,0,0,0) || return false
    qmul(g,g,p) == QElt(-2p,0,0,0) || return false
    return true
end

# ------------------
# Hermitian forms
# ------------------

struct Hermitian2x2
    u::Int
    a::QElt
    v::Int
end

Base.show(io::IO, H::Hermitian2x2) = print(io, "Hermitian2x2(u=", H.u, ", a=", H.a, ", v=", H.v, ")")
Base.:(==)(X::Hermitian2x2, Y::Hermitian2x2) = (X.u == Y.u && X.v == Y.v && X.a == Y.a)
Base.hash(H::Hermitian2x2, h::UInt) = hash((H.u, _qkey(H.a), H.v), h)

function is_principal(H::Hermitian2x2, p::Int)
    return H.u > 0 && H.v > 0 && H.u*H.v - nrd(H.a,p) == 1
end

function normalize_hermitian(H::Hermitian2x2)
    if H.u > H.v
        return Hermitian2x2(H.v, qconj(H.a), H.u)
    elseif H.u == H.v
        ac = qconj(H.a)
        return _qkey(ac) < _qkey(H.a) ? Hermitian2x2(H.u, ac, H.v) : H
    else
        return H
    end
end

# -------------------------------------
# Candidate enumeration
# -------------------------------------
"""

    First candidate enumaration has been done by using Lemma 13 in 'On class numbers of positive definite binary
    quaternion hermitian forms', Ki-ichiro Hashimoto and Tomoyoshi Ibukiyama. This enables us to eliminate easy cases.

"""


const _SHELL_CACHE = Dict{Tuple{Int,Int},Vector{QElt}}()
const _PACKET_CACHE = Dict{Tuple{Int,Hermitian2x2},Vector{Matrix{Int}}}()

function clear_caches!()
    empty!(_SHELL_CACHE)
    empty!(_PACKET_CACHE)
    return nothing
end

function _units_norm1(p::Int)
    out = QElt[]
    for A0 in -2:2, A1 in -2:2, A2 in -2:2, A3 in -2:2
        if A0*A0 + A1*A1 + p*(A2*A2 + A3*A3) == 4 &&
           iseven(A0 - A3) && iseven(A1 - A2)
            push!(out, QElt(A0,A1,A2,A3))
        end
    end
    return sort!(unique(out), by=_qkey)
end

function _lemma3_residue_key(a::QElt, v::Int, p::Int)
    units = _units_norm1(p)
    keys = NTuple{4,Int}[]
    for u1 in units, u2 in units
        x = qmul(qconj(u1), qmul(a, u2, p), p)
        push!(keys, (mod(x.A0, 2v), mod(x.A1, 2v), mod(x.A2, 2v), mod(x.A3, 2v)))
        xc = qconj(x)
        push!(keys, (mod(xc.A0, 2v), mod(xc.A1, 2v), mod(xc.A2, 2v), mod(xc.A3, 2v)))
    end
    return minimum(keys)
end

function _order_elements_exact_norm(p::Int, N::Int)
    key = (p,N)
    haskey(_SHELL_CACHE, key) && return _SHELL_CACHE[key]

    vals = QElt[]
    if N == 0
        vals = QElt[ZEROQ]
    elseif N > 0
        B = Int(ceil(2sqrt(N))) + 2
        for w in -B:B, x in -B:B, y in -B:B, z in -B:B
            q = QElt(2w + z, 2x + y, y, z)
            nrd(q,p) == N && push!(vals, q)
        end
        vals = sort!(unique(vals), by=_qkey)
    end

    _SHELL_CACHE[key] = vals
    return vals
end

function _candidate_from_uv(p::Int, u::Int, v::Int)
    N = u*v - 1
    by_residue = Dict{Any,QElt}()
    for a in _order_elements_exact_norm(p, N)
        k = _lemma3_residue_key(a, v, p)
        get!(by_residue, k, a)
    end

    out = Hermitian2x2[]
    for a in values(by_residue)
        H = normalize_hermitian(Hermitian2x2(u,a,v))
        is_principal(H,p) && push!(out, H)
    end
    return sort!(unique(out), by = H -> (H.u, H.v, _qkey(H.a)))
end

function enumerate_principal_candidates(p::Int; Amax::Int)
    Amax >= 1 || return Hermitian2x2[]
    reps = Dict{Tuple{Int,Int,NTuple{4,Int}},Hermitian2x2}()
    for v in 1:Amax, u in 1:v
        for H in _candidate_from_uv(p,u,v)
            reps[(H.u,H.v,_qkey(H.a))] = H
        end
    end
    return sort!(collect(values(reps)), by = H -> (H.u, H.v, _qkey(H.a)))
end

# -----------------------------------------------------------------------
# Trace forms (auxiliary forms): F_a(x,y) = Trd(h_H(a*x,y)), a = 1,i,j,k
# -----------------------------------------------------------------------

"""
    The trace form and auxiliary form definitions are taken from the PhD thesis of Sarah Chisholm, University of Calgary, Section 3.8.1 and Lemma 97.
    In Lemma 97,  the result cites Lemma 6.2 of Lattice methods for algebraic modular forms on
    classical groups by  Greenberg and Voight.
"""

function _scale(q::QElt, c::Int)
    return QElt(c*q.A0, c*q.A1, c*q.A2, c*q.A3)
end

_leftmul_pair(a::QElt, x::NTuple{2,QElt}, p::Int) = (qmul(a,x[1],p), qmul(a,x[2],p))

function _hermitian_value(H::Hermitian2x2, x::NTuple{2,QElt}, y::NTuple{2,QElt}, p::Int)
    x1,x2 = x
    y1,y2 = y
    return _scale(qmul(x1, qconj(y1), p), H.u) +
           qmul(x1, qmul(H.a, qconj(y2), p), p) +
           qmul(x2, qmul(qconj(H.a), qconj(y1), p), p) +
           _scale(qmul(x2, qconj(y2), p), H.v)
end

function canonical_packet_forms(H::Hermitian2x2, p::Int)
    key = (p,H)
    haskey(_PACKET_CACHE, key) && return _PACKET_CACHE[key]
    is_principal(H,p) || error("canonical_packet_forms expects a principal positive Hermitian form; got $H")

    B = O2_basis8()
    A = algebra_basis_E()
    F = Matrix{Int}[]
    for a in A
        G = zeros(Int, 8, 8)
        for r in 1:8, s in 1:8
            G[r,s] = trd(_hermitian_value(H, _leftmul_pair(a, B[r], p), B[s], p))
        end
        push!(F, G)
    end
    _PACKET_CACHE[key] = F
    return F
end

# ============================================================
# Hashimoto--Ibukiyama class numbers (the number of principal polarizations)
# ============================================================

function check_p_11mod12(p::Int)
    Oscar.is_prime(p) || error("p must be prime; got p=$p")
    p % 12 == 11 || error("p must be 11 mod 12; got p=$p")
    return true
end

"""
    The class numbers are taken from 'Supersingular curves of genus two and class numbers' by Tomoyoshi Ibukiyama, Toshiyuki Katsura, and Frans Oort.
 
    classNumbers(p) returns the class numbers h, h(h+1)/2 and H.
    H is the number of all polarizations, and h(h+1)/2 is the number of reducible polarizations.
"""

function classNumbers(p::Integer)
    # Handle small p quickly.
    if p == 2 || p == 3
        return Int(1), Int(1), Int(0)
    elseif p == 5
        return Int(2), Int(1), Int(1)
    end

    # Compute H.
    a = 1 - Oscar.jacobi_symbol(-1, p)
    b = 1 - Oscar.jacobi_symbol(-2, p)
    c = 1 - Oscar.jacobi_symbol(-3, p)
    d = (p % 5 == 4) ? (4//5) : 0//1

    H =
        ((p - 1) * (p + 12) * (p + 23))//2880 +
        (a * (2p + 13))//96 +
        (c * (p + 11))//36 +
        (b//8) +
        ((a * c)//12) +
        d

    # Compute h.
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

    return Int(h), Int((h * (h + 1)) ÷ 2), Int(H)
end
targetH(p::Integer) = classNumbers(p)[3]

# ----------------------
# Output
# ----------------------

function _write_polarization_row(io, k::Int, H::Hermitian2x2, p::Int)
    w,x,y,z = order_coordinates(H.a)
    println(io, "[$k]\t$(H.u)\t$(H.v)\t$w\t$x\t$y\t$z\tdet1=$(is_principal(H,p))")
end

function write_polarizations_file(path::AbstractString, reps::Vector{Hermitian2x2}, p::Int)
    open(path, "w") do io
        println(io, "p=$p")
        println(io, "# format: u v w x y z where a = w + x*i + y*(i+j)/2 + z*(1+k)/2")
        for (k,H) in enumerate(reps)
            _write_polarization_row(io, k, H, p)
        end
    end
    return path
end



# ============================================================
# LLL-reduced packet isometry and dynamic batch classification
# ============================================================
# The functions below keep the original Hermitian2x2 representatives for
# output, but use a unimodular Z^8-basis change of the whole packet
# [F1,F2,F3,F4] before calling Hecke's Plesken--Souvignier routines.
# They also use exact F1-isometry as a safe gate:
#     packet isometry => F1-isometry.

_now() = Dates.format(Dates.now(), "HH:MM:SS")

# -------------------------------
# Integer matrix / Hecke helpers
# -------------------------------

function _zzmat(A::Matrix{Int})
    return Hecke.matrix(
        Hecke.ZZ,
        [[Hecke.ZZ(A[i,j]) for j in 1:size(A,2)] for i in 1:size(A,1)]
    )
end

_zzpacket(F::Vector{Matrix{Int}}) = Hecke.ZZMatrix[_zzmat(A) for A in F]
_mat_int(M) = Matrix{Int}([Int(M[i,j]) for i in 1:Hecke.nrows(M), j in 1:Hecke.ncols(M)])
_mat_big(M) = Matrix{BigInt}([BigInt(M[i,j]) for i in 1:Hecke.nrows(M), j in 1:Hecke.ncols(M)])
_big(A::Matrix{Int}) = BigInt.(A)
_diagmax(A::Matrix{Int}) = maximum(A[i,i] for i in 1:min(size(A,1), size(A,2)))
_det_big(A::Matrix{Int}) = det(BigInt.(A))

function _identity_big(n::Int)
    Ibig = zeros(BigInt, n, n)
    for i in 1:n
        Ibig[i,i] = 1
    end
    return Ibig
end

function _inv_unimodular_bigint(T::Matrix{BigInt})
    d = det(T)
    @assert d == 1 || d == -1 "matrix is not unimodular; det=$d"
    n = size(T, 1)
    R = Matrix{Rational{BigInt}}(T) \ Matrix{Rational{BigInt}}(_identity_big(n))
    @assert all(denominator(x) == 1 for x in R) "inverse is not integral"
    return Matrix{BigInt}([numerator(R[i,j]) for i in 1:n, j in 1:n])
end

_inv_unimodular_int(P::Matrix{Int}) = Int.(_inv_unimodular_bigint(BigInt.(P)))

# -----------------------------------------
# LLL transport of the whole four-form packet
# -----------------------------------------

function _lll_gram_transform(G::Matrix{Int})
    Gzz = _zzmat(G)
    if isdefined(Hecke, :lll_gram_with_transform)
        A, P = Hecke.lll_gram_with_transform(Gzz)
        return _mat_int(A), _mat_int(P), :Hecke_lll_gram_with_transform
    elseif isdefined(Oscar, :lll_gram_with_transform)
        A, P = Oscar.lll_gram_with_transform(Gzz)
        return _mat_int(A), _mat_int(P), :Oscar_lll_gram_with_transform
    else
        error("No lll_gram_with_transform found in Hecke or Oscar")
    end
end

function _detect_transport_convention(G::Matrix{Int}, Gred::Matrix{Int}, P::Matrix{Int})
    if P * G * transpose(P) == Gred
        return :row_PGPt, P
    end
    if transpose(P) * G * P == Gred
        return :col_PtGP, P
    end
    if _det_big(P) == 1 || _det_big(P) == -1
        Pinv = _inv_unimodular_int(P)
        if Pinv * G * transpose(Pinv) == Gred
            return :row_invP_G_invPt, Pinv
        end
        if transpose(Pinv) * G * Pinv == Gred
            return :col_invPt_G_invP, Pinv
        end
    end
    error("Could not detect LLL transport convention")
end

function _transport_packet(F::Vector{Matrix{Int}}, P::Matrix{Int}, conv::Symbol)
    length(F) == 4 || error("expected a 4-form packet")
    if conv == :row_PGPt || conv == :row_invP_G_invPt
        return [P * F[k] * transpose(P) for k in 1:4]
    elseif conv == :col_PtGP || conv == :col_invPt_G_invP
        return [transpose(P) * F[k] * P for k in 1:4]
    else
        error("unknown transport convention $conv")
    end
end

function reduce_packet_by_lll(F::Vector{Matrix{Int}})
    length(F) == 4 || error("expected a 4-form packet")
    old_bound = _diagmax(F[1])
    Gred, Praw, source = _lll_gram_transform(F[1])
    detP = _det_big(Praw)
    @assert detP == 1 || detP == -1 "LLL transform is not unimodular; det=$detP"
    conv, Puse = _detect_transport_convention(F[1], Gred, Praw)
    Fred = _transport_packet(F, Puse, conv)
    @assert Fred[1] == Gred "transported F1 does not match reported reduced Gram"
    return (Fred=Fred, Praw=Praw, Puse=Puse, convention=conv,
            source=source, old_bound=old_bound, new_bound=_diagmax(Fred[1]), detP=detP)
end

# ----------------------------
# Exact Hecke isometry helpers
# ----------------------------

function _verify_all_forms(T::Matrix{BigInt}, F::Vector{Matrix{Int}}, G::Vector{Matrix{Int}})
    # Hecke may return either direction/convention depending on setup path.
    Fb = [_big(A) for A in F]
    Gb = [_big(A) for A in G]
    row_G_to_F = all(T * Gb[k] * transpose(T) == Fb[k] for k in 1:4)
    row_F_to_G = all(T * Fb[k] * transpose(T) == Gb[k] for k in 1:4)
    col_G_to_F = all(transpose(T) * Gb[k] * T == Fb[k] for k in 1:4)
    col_F_to_G = all(transpose(T) * Fb[k] * T == Gb[k] for k in 1:4)
    inv_F_to_G = false
    if det(T) == 1 || det(T) == -1
        Ti = _inv_unimodular_bigint(T)
        inv_F_to_G = all(Ti * Fb[k] * transpose(Ti) == Gb[k] for k in 1:4)
    end
    return row_G_to_F || row_F_to_G || col_G_to_F || col_F_to_G || inv_F_to_G
end

function _verify_f1(T::Matrix{BigInt}, F::Matrix{Int}, G::Matrix{Int})
    Fb, Gb = _big(F), _big(G)
    row_G_to_F = T * Gb * transpose(T) == Fb
    row_F_to_G = T * Fb * transpose(T) == Gb
    col_G_to_F = transpose(T) * Gb * T == Fb
    col_F_to_G = transpose(T) * Fb * T == Gb
    inv_F_to_G = false
    if det(T) == 1 || det(T) == -1
        Ti = _inv_unimodular_bigint(T)
        inv_F_to_G = Ti * Fb * transpose(Ti) == Gb
    end
    return row_G_to_F || row_F_to_G || col_G_to_F || col_F_to_G || inv_F_to_G
end

function _choose_input_output(F, G; minbound_direction::Bool=true)
    bF, bG = _diagmax(F[1]), _diagmax(G[1])
    if minbound_direction && bG < bF
        return G, F, :G_input, bG, bF
    else
        return F, G, :F_input, bF, bG
    end
end

function _try_setup(FF, GG; depth::Int=0, bacher_depth::Int=0)
    fl, CF, CG = Hecke._try_iso_setup_small(FF, GG; depth=depth, bacher_depth=bacher_depth)
    if fl
        return CF, CG, true
    end
    CF, CG = Hecke._iso_setup(FF, GG; depth=depth, bacher_depth=bacher_depth)
    return CF, CG, false
end

function hecke_packet_isometric(F::Vector{Matrix{Int}}, G::Vector{Matrix{Int}};
                                depth::Int=0, bacher_depth::Int=0,
                                verify::Bool=true, minbound_direction::Bool=true)
    Fin, Gout, direction, input_bound, output_bound = _choose_input_output(F, G; minbound_direction=minbound_direction)
    tsetup = time()
    CF, CG, setup_small = _try_setup(_zzpacket(Fin), _zzpacket(Gout); depth=depth, bacher_depth=bacher_depth)
    setup_time = time() - tsetup
    tiso = time()
    b, rawT = Hecke.isometry(CF, CG)
    iso_time = time() - tiso
    if Bool(b) && verify
        _verify_all_forms(_mat_big(rawT), Fin, Gout) || error("packet isometry verification failed")
    end
    return (b=Bool(b), direction=direction, input_bound=input_bound, output_bound=output_bound,
            setup_small=setup_small, setup_time=setup_time, iso_time=iso_time,
            total_time=setup_time+iso_time)
end

function hecke_f1_isometric(F1::Matrix{Int}, G1::Matrix{Int};
                            depth::Int=0, bacher_depth::Int=0,
                            verify::Bool=true, minbound_direction::Bool=true)
    F = [F1]
    G = [G1]
    bF, bG = _diagmax(F1), _diagmax(G1)
    if minbound_direction && bG < bF
        Fin, Gout, direction, input_bound, output_bound = G1, F1, :G_input, bG, bF
    else
        Fin, Gout, direction, input_bound, output_bound = F1, G1, :F_input, bF, bG
    end
    tsetup = time()
    CF, CG, setup_small = _try_setup(Hecke.ZZMatrix[_zzmat(Fin)], Hecke.ZZMatrix[_zzmat(Gout)]; depth=depth, bacher_depth=bacher_depth)
    setup_time = time() - tsetup
    tiso = time()
    b, rawT = Hecke.isometry(CF, CG)
    iso_time = time() - tiso
    if Bool(b) && verify
        _verify_f1(_mat_big(rawT), Fin, Gout) || error("F1 isometry verification failed")
    end
    return (b=Bool(b), direction=direction, input_bound=input_bound, output_bound=output_bound,
            setup_small=setup_small, setup_time=setup_time, iso_time=iso_time,
            total_time=setup_time+iso_time)
end

# -----------------------------
# Dynamic classification state
# -----------------------------

mutable struct ClassificationState
    p::Int
    targetH::Int
    Amax::Int
    candidates::Vector{Hermitian2x2}
    raw_packets::Vector{Any}
    red_packets::Vector{Any}
    redinfo::Vector{Any}
    reps::Vector{Hermitian2x2}
    reps_idx::Vector{Int}
    f1_class_reps_idx::Vector{Int}
    f1_bucket_reps_idx::Vector{Vector{Int}}
    old_bounds::Vector{Int}
    new_bounds::Vector{Int}
    packet_checked::Int
    packet_merged::Int
    f1_checked::Int
    f1_merged::Int
    f1_new_classes::Int
    enum_time::Float64
    packet_time::Float64
    reduce_time::Float64
    merge_time::Float64
    packet_setup_time::Float64
    packet_iso_time::Float64
    f1_setup_time::Float64
    f1_iso_time::Float64
    slow_count::Int
end

function _empty_state(p::Int, H::Int)
    return ClassificationState(p, H, 0, Hermitian2x2[], Any[], Any[], Any[], Hermitian2x2[], Int[],
                               Int[], Vector{Int}[], Int[], Int[], 0, 0, 0, 0, 0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0)
end

_candidate_key(H::Hermitian2x2) = (H.u, H.v, _qkey(H.a))

function _old_to_new_index_map(old_cand::Vector{Hermitian2x2}, new_cand::Vector{Hermitian2x2})
    pos = Dict{Tuple{Int,Int,NTuple{4,Int}},Int}()
    for (j,H) in pairs(new_cand)
        pos[_candidate_key(H)] = j
    end
    return [pos[_candidate_key(H)] for H in old_cand]
end

function _next_Amax(A::Int; small_until::Int=12, medium_until::Int=40,
                    small_step::Int=1, medium_step::Int=2, large_step::Int=4)
    if A < small_until
        return A + small_step
    elseif A < medium_until
        return A + medium_step
    else
        return A + large_step
    end
end

function _ordered_new_indices(new_indices::Vector{Int}, old_bounds::Vector{Int}, new_bounds::Vector{Int}, order::Symbol)
    if order == :enumerator
        return new_indices
    elseif order == :small_bound
        return sort(new_indices; by=i -> (new_bounds[i], old_bounds[i], i))
    elseif order == :old_bound
        return sort(new_indices; by=i -> (old_bounds[i], new_bounds[i], i))
    else
        error("unknown process_new_order=$order")
    end
end

function _ordered_reps(reps_idx::Vector{Int}, bounds::Vector{Int}, i::Int, order::Symbol)
    if order == :discovery
        return reps_idx
    elseif order == :small_bound
        return sort(reps_idx; by=r -> (bounds[r], abs(bounds[i]-bounds[r]), r))
    else
        error("unknown rep_order=$order")
    end
end

function _find_f1_class(red_packets, f1_class_reps_idx::Vector{Int}, i::Int;
                        depth::Int=0, bacher_depth::Int=0, verify::Bool=true)
    checked = 0
    merged = 0
    setup = 0.0
    iso = 0.0
    for (cid, rep_i) in enumerate(f1_class_reps_idx)
        checked += 1
        r = hecke_f1_isometric(red_packets[i][1], red_packets[rep_i][1];
                               depth=depth, bacher_depth=bacher_depth, verify=verify)
        setup += r.setup_time
        iso += r.iso_time
        if r.b
            merged += 1
            return cid, false, (checked=checked, merged=merged, setup=setup, iso=iso)
        end
    end
    return length(f1_class_reps_idx) + 1, true, (checked=checked, merged=merged, setup=setup, iso=iso)
end

function _extend_to_Amax!(SP, st::ClassificationState, Amax_new::Int;
                          rep_order::Symbol=:small_bound,
                          process_new_order::Symbol=:small_bound,
                          stop_when_target::Bool=true,
                          depth::Int=0, bacher_depth::Int=0, verify::Bool=true,
                          print_every::Int=100, slow_threshold::Float64=2.0,
                          outfile::Union{Nothing,AbstractString}=nothing)
    p = st.p
    println("\n[$(_now())] p=$p extend Amax $(st.Amax) -> $Amax_new")
    t = time(); new_cand = enumerate_principal_candidates(p; Amax=Amax_new); st.enum_time += time() - t
    println("[$(_now())] p=$p candidates=$(length(new_cand))")
    flush(stdout)

    old_to_new = isempty(st.candidates) ? Int[] : _old_to_new_index_map(st.candidates, new_cand)
    old_set = Set(old_to_new)
    raw = Vector{Any}(undef, length(new_cand))
    red = Vector{Any}(undef, length(new_cand))
    info = Vector{Any}(undef, length(new_cand))

    for old_i in eachindex(st.candidates)
        new_i = old_to_new[old_i]
        raw[new_i] = st.raw_packets[old_i]
        red[new_i] = st.red_packets[old_i]
        info[new_i] = st.redinfo[old_i]
    end

    t = time(); built_raw = 0
    for i in eachindex(new_cand)
        if !isassigned(raw, i)
            raw[i] = canonical_packet_forms(new_cand[i], p)
            built_raw += 1
        end
    end
    st.packet_time += time() - t

    t = time(); built_red = 0
    for i in eachindex(new_cand)
        if !isassigned(red, i)
            R = reduce_packet_by_lll(raw[i])
            info[i] = R
            red[i] = R.Fred
            built_red += 1
            if built_red == 1 || built_red % print_every == 0
                println("[$(_now())] p=$p reduced new packet $built_red old_bound=$(R.old_bound) new_bound=$(R.new_bound)")
                flush(stdout)
            end
        end
    end
    st.reduce_time += time() - t

    old_bounds = [info[i].old_bound for i in eachindex(info)]
    new_bounds = [info[i].new_bound for i in eachindex(info)]

    reps_idx = isempty(st.reps_idx) ? Int[] : [old_to_new[i] for i in st.reps_idx]
    reps = isempty(st.reps) ? Hermitian2x2[] : copy(st.reps)
    f1_class_reps_idx = isempty(st.f1_class_reps_idx) ? Int[] : [old_to_new[i] for i in st.f1_class_reps_idx]
    f1_bucket_reps_idx = Vector{Int}[]
    for bucket in st.f1_bucket_reps_idx
        push!(f1_bucket_reps_idx, [old_to_new[i] for i in bucket])
    end

    new_indices = [i for i in eachindex(new_cand) if !(i in old_set)]
    new_indices = _ordered_new_indices(new_indices, old_bounds, new_bounds, process_new_order)
    println("[$(_now())] p=$p reused=$(length(st.candidates)) built_raw=$built_raw built_red=$built_red process_new=$(length(new_indices))")
    flush(stdout)

    packet_checked = packet_merged = f1_checked = f1_merged = f1_new = slow = 0
    packet_setup = packet_iso = f1_setup = f1_iso = 0.0
    tmerge = time()

    for (count, i) in enumerate(new_indices)
        if stop_when_target && length(reps) >= st.targetH
            println("[$(_now())] p=$p target reached inside Amax=$Amax_new; stopping stage early")
            break
        end

        cid, isnewF1, fs = _find_f1_class(red, f1_class_reps_idx, i;
                                          depth=depth, bacher_depth=bacher_depth, verify=verify)
        f1_checked += fs.checked; f1_merged += fs.merged; f1_setup += fs.setup; f1_iso += fs.iso
        keep = true

        if isnewF1
            push!(f1_class_reps_idx, i)
            push!(f1_bucket_reps_idx, Int[])
            f1_new += 1
        else
            for rep_i in _ordered_reps(f1_bucket_reps_idx[cid], new_bounds, i, rep_order)
                packet_checked += 1
                r = hecke_packet_isometric(red[i], red[rep_i]; depth=depth, bacher_depth=bacher_depth, verify=verify)
                packet_setup += r.setup_time; packet_iso += r.iso_time
                if r.b
                    packet_merged += 1
                    keep = false
                    break
                elseif r.total_time >= slow_threshold
                    slow += 1
                    println("  SLOW p=$p i=$i rep=$rep_i total=$(round(r.total_time; digits=3))s")
                    flush(stdout)
                end
            end
        end

        if keep
            push!(reps, new_cand[i])
            push!(reps_idx, i)
            push!(f1_bucket_reps_idx[cid], i)
            if outfile !== nothing
                open(outfile, "a") do io
                    _write_polarization_row(io, length(reps), new_cand[i], p)
                end
            end
            println("  NEW p=$p rep #$(length(reps)) at i=$i Amax=$Amax_new F1class=$cid newF1=$isnewF1 H=$(new_cand[i])")
            flush(stdout)
            if stop_when_target && length(reps) >= st.targetH
                println("[$(_now())] p=$p target reached by new rep at Amax=$Amax_new")
                break
            end
        end

        if count == 1 || count % print_every == 0 || count == length(new_indices)
            println("[$(_now())] p=$p Amax=$Amax_new processed_new=$count/$(length(new_indices)) reps=$(length(reps))/$(st.targetH) F1classes=$(length(f1_class_reps_idx)) packet_checked_stage=$packet_checked elapsed=$(round(time()-tmerge; digits=1))s")
            flush(stdout)
        end
    end

    st.Amax = Amax_new
    st.candidates = new_cand
    st.raw_packets = raw
    st.red_packets = red
    st.redinfo = info
    st.reps = reps
    st.reps_idx = reps_idx
    st.f1_class_reps_idx = f1_class_reps_idx
    st.f1_bucket_reps_idx = f1_bucket_reps_idx
    st.old_bounds = old_bounds
    st.new_bounds = new_bounds

    st.packet_checked += packet_checked
    st.packet_merged += packet_merged
    st.f1_checked += f1_checked
    st.f1_merged += f1_merged
    st.f1_new_classes += f1_new
    st.merge_time += time() - tmerge
    st.packet_setup_time += packet_setup
    st.packet_iso_time += packet_iso
    st.f1_setup_time += f1_setup
    st.f1_iso_time += f1_iso
    st.slow_count += slow

    println("[$(_now())] p=$p DONE Amax=$Amax_new reps=$(length(st.reps))/$(st.targetH) F1classes=$(length(st.f1_class_reps_idx)) packet_checked_stage=$packet_checked")
    flush(stdout)
    return st
end

function _state_stats(st::ClassificationState, total_time::Float64; Astart, process_new_order, rep_order, depth, bacher_depth)
    sb_old = sort(st.old_bounds); sb_new = sort(st.new_bounds)
    return Dict{Symbol,Any}(
        :Amax => st.Amax, :Astart => Astart, :Astep => "dynamic",
        :depth => depth, :bacher_depth => bacher_depth,
        :process_new_order => process_new_order, :rep_order => rep_order,
        :targetH => st.targetH, :classes => length(st.reps), :candidates => length(st.candidates),
        :total_time => total_time, :enum_time => st.enum_time, :packet_time => st.packet_time,
        :reduce_time => st.reduce_time, :merge_time => st.merge_time,
        :packet_checked => st.packet_checked, :packet_merged => st.packet_merged,
        :f1_checked => st.f1_checked, :f1_merged => st.f1_merged,
        :f1_classes => length(st.f1_class_reps_idx), :f1_new_classes => st.f1_new_classes,
        :packet_setup_time => st.packet_setup_time, :packet_iso_time => st.packet_iso_time,
        :f1_setup_time => st.f1_setup_time, :f1_iso_time => st.f1_iso_time,
        :slow_count => st.slow_count,
        :old_bound_min => isempty(sb_old) ? 0 : first(sb_old),
        :old_bound_median => isempty(sb_old) ? 0 : sb_old[cld(length(sb_old),2)],
        :old_bound_max => isempty(sb_old) ? 0 : last(sb_old),
        :new_bound_min => isempty(sb_new) ? 0 : first(sb_new),
        :new_bound_median => isempty(sb_new) ? 0 : sb_new[cld(length(sb_new),2)],
        :new_bound_max => isempty(sb_new) ? 0 : last(sb_new),
    )
end

function classify_prime_lll_until_target(p::Int;
                                         Astart::Int=1, Amax_cap::Union{Nothing,Int}=nothing,
                                         small_until::Int=12, medium_until::Int=40,
                                         small_step::Int=1, medium_step::Int=2, large_step::Int=4,
                                         rep_order::Symbol=:small_bound,
                                         process_new_order::Symbol=:small_bound,
                                         stop_when_target::Bool=true,
                                         depth::Int=0, bacher_depth::Int=0,
                                         verify::Bool=true, print_every::Int=100,
                                         slow_threshold::Float64=2.0,
                                         max_steps::Union{Nothing,Int}=nothing,
                                         max_wall_seconds::Union{Nothing,Float64}=nothing,
                                         outfile::Union{Nothing,AbstractString}=nothing)
    check_p_11mod12(p)
    st = _empty_state(p, targetH(p))
    t0 = time()
    A = Astart
    steps = 0
    while true
        if Amax_cap !== nothing && A > Amax_cap
            println("[$(_now())] p=$p stopping because next Amax=$A exceeds Amax_cap=$Amax_cap")
            break
        end
        if max_steps !== nothing && steps >= max_steps
            println("[$(_now())] p=$p stopping because max_steps=$max_steps was reached")
            break
        end
        if max_wall_seconds !== nothing && time() - t0 >= max_wall_seconds
            println("[$(_now())] p=$p stopping because max_wall_seconds=$max_wall_seconds was reached")
            break
        end
        steps += 1
        _extend_to_Amax!(SuperspecialPolarizations_ThesisOscar, st, A;
                         rep_order=rep_order, process_new_order=process_new_order,
                         stop_when_target=stop_when_target, depth=depth,
                         bacher_depth=bacher_depth, verify=verify,
                         print_every=print_every, slow_threshold=slow_threshold,
                         outfile=outfile)
        length(st.reps) >= st.targetH && break
        A = _next_Amax(A; small_until=small_until, medium_until=medium_until,
                       small_step=small_step, medium_step=medium_step, large_step=large_step)
    end
    return st, _state_stats(st, time() - t0; Astart=Astart, process_new_order=process_new_order,
                            rep_order=rep_order, depth=depth, bacher_depth=bacher_depth)
end

function _write_summary_row(io, p::Int, stats, outfile::AbstractString, status::String)
    println(io, join((
        p, stats[:targetH], stats[:classes], stats[:Amax], stats[:candidates],
        stats[:packet_checked], stats[:packet_merged], stats[:f1_checked], stats[:f1_classes],
        round(stats[:f1_setup_time]; digits=3), round(stats[:f1_iso_time]; digits=3),
        round(stats[:total_time]; digits=3), round(stats[:enum_time]; digits=3),
        round(stats[:packet_time]; digits=3), round(stats[:reduce_time]; digits=3),
        round(stats[:merge_time]; digits=3), round(stats[:packet_setup_time]; digits=3),
        round(stats[:packet_iso_time]; digits=3), stats[:slow_count],
        stats[:old_bound_max], stats[:new_bound_max], outfile, status,
        stats[:process_new_order], stats[:rep_order]
    ), '\t'))
end

function _write_enhanced_polarizations_file(path::AbstractString, st::ClassificationState, stats)
    open(path, "a") do io
        println(io, "# targetH=$(st.targetH) classes=$(length(st.reps)) Amax=$(st.Amax)")
        println(io, "# timing total=$(round(stats[:total_time]; digits=2)) enum=$(round(stats[:enum_time]; digits=2)) packet=$(round(stats[:packet_time]; digits=2)) reduce=$(round(stats[:reduce_time]; digits=2)) merge=$(round(stats[:merge_time]; digits=2))")
        println(io, "# checks packet_checked=$(st.packet_checked) packet_merged=$(st.packet_merged) f1_checked=$(st.f1_checked) f1_classes=$(length(st.f1_class_reps_idx))")
    end
    return path
end

function run_batch_lll_until_target(; primes::Vector{Int}, outdir::AbstractString="batch_outputs_lll_F1gate",
                                    Astart::Int=1, Amax_cap::Union{Nothing,Int}=nothing,
                                    small_until::Int=12, medium_until::Int=40,
                                    small_step::Int=1, medium_step::Int=2, large_step::Int=4,
                                    rep_order::Symbol=:small_bound,
                                    process_new_order::Symbol=:small_bound,
                                    stop_when_target::Bool=true,
                                    depth::Int=0, bacher_depth::Int=0,
                                    verify::Bool=true, print_every::Int=100,
                                    slow_threshold::Float64=2.0,
                                    max_steps::Union{Nothing,Int}=nothing,
                                    max_wall_seconds::Union{Nothing,Float64}=nothing)
    isdir(outdir) || mkpath(outdir)
    summary = joinpath(outdir, "summary.tsv")
    open(summary, "w") do io
        println(io, "p\ttargetH\tclasses\tAmax\tcandidates\tpacket_checked\tpacket_merged\tf1_checked\tf1_classes\tf1_setup_s\tf1_iso_s\ttotal_s\tenum_s\tpacket_s\treduce_s\tmerge_s\tpacket_setup_s\tpacket_iso_s\tslow_count\told_bound_max\tnew_bound_max\toutfile\tstatus\tprocess_new_order\trep_order")
    end
    results = Dict{Int,Any}()
    for p in primes
        println("\n============================================================")
        println("LLL/F1-gate batch p=$p targetH=$(targetH(p))")
        println("============================================================")
        outfile = joinpath(outdir, "polarizations_p$(p).txt")
        write_polarizations_file(outfile, Hermitian2x2[], p)
        st, stats = classify_prime_lll_until_target(p; Astart=Astart, Amax_cap=Amax_cap,
            small_until=small_until, medium_until=medium_until, small_step=small_step,
            medium_step=medium_step, large_step=large_step, rep_order=rep_order,
            process_new_order=process_new_order, stop_when_target=stop_when_target,
            depth=depth, bacher_depth=bacher_depth, verify=verify, print_every=print_every,
            slow_threshold=slow_threshold, max_steps=max_steps, max_wall_seconds=max_wall_seconds,
            outfile=outfile)
        _write_enhanced_polarizations_file(outfile, st, stats)
        status = length(st.reps) == st.targetH ? "TARGET_PASSED" : (length(st.reps) > st.targetH ? "OVER_TARGET" : "INCOMPLETE")
        open(summary, "a") do io
            _write_summary_row(io, p, stats, outfile, status)
        end
        println("$status p=$p classes=$(length(st.reps))/$(st.targetH), wrote $outfile")
        results[p] = (state=st, stats=stats, outfile=outfile, status=status)
    end
    return (results=results, outdir=outdir, summary=summary)
end


end # module

# If this file is invoked directly as a script, accept a prime p from ARGS
# and run the batch runner for that single prime. This keeps the change
# minimal so the SLURM array script can call `julia RHI.jl <p>`.
if abspath(PROGRAM_FILE) == @__FILE__
    # Module is defined above in this file; refer to it by name to avoid
    # relative `using` scope issues when the file is executed as a script.

    if length(ARGS) < 1
        println("Usage: julia RHI.jl <prime>")
        exit(1)
    end

    p = try
        parse(Int, ARGS[1])
    catch e
        println("Failed to parse prime from ARGS[1]=", ARGS[1])
        rethrow(e)
    end

    outdir = joinpath("batch_outputs_lll_F1gate", string(p))
    mkpath(outdir)
    logfile = joinpath(outdir, "log_p$(p).txt")
    open(logfile, "w") do io
        redirect_stdout(io) do
            redirect_stderr(io) do
                println("Running RHI for p=", p)
                SuperspecialPolarizations_ThesisOscar.run_batch_lll_until_target(
                    primes=[p],
                    outdir=outdir,
                )
            end
        end
    end
end
