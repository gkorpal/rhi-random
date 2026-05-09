module SuperspecialPolarizations_ThesisOscar

using LinearAlgebra
import Oscar

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
       certify_basis_conventions
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
targetH(p::Integer) = classNumbers(p).H

# ----------------------
# Output
# ----------------------

function write_polarizations_file(path::AbstractString, reps::Vector{Hermitian2x2}, p::Int)
    open(path, "w") do io
        println(io, "# principal polarization representatives for p=$p")
        println(io, "# format: u v w x y z where a = w + x*i + y*(i+j)/2 + z*(1+k)/2")
        for (k,H) in enumerate(reps)
            w,x,y,z = order_coordinates(H.a)
            println(io, "$k\t$(H.u)\t$(H.v)\t$w\t$x\t$y\t$z\tdet1=$(is_principal(H,p))")
        end
    end
    return path
end



end # module
