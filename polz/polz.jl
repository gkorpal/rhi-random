using Oscar
using Random
using Printf

# ============================================================
# Z conventions
# ============================================================

const Z    = ZZ
const ZElem = ZZRingElem

# ============================================================
# 0) Primes (p ≡ 11 mod 12)
# ============================================================

"""
    check_p_11mod12(p::Int)

Error unless p is prime and p ≡ 11 (mod 12).
Uses Oscar/Nemo primality via `is_prime(ZZ(p))`.
"""
function check_p_11mod12(p::Int)
    is_prime(ZZ(p)) || error("p must be prime; got p=$p")
    (p % 12 == 11) || error("p must be 11 mod 12; got p=$p")
    return true
end

# ============================================================
# 2) Quaternion algebra B(-1,-p) conventions
# ============================================================

"""
    QElt(A0,A1,A2,A3)

Half-coordinates element (Int):
  (A0/2) + (A1/2)i + (A2/2)j + (A3/2)k,  k=i*j
"""
struct QElt
    A0::Int
    A1::Int
    A2::Int
    A3::Int
end

zeroq() = QElt(0,0,0,0)

"""
    from_wxyz(w,x,y,z) -> QElt

Convert order coordinates (w,x,y,z) in basis
[1, i, (i+j)/2, (1+ij)/2] into half-coordinates.
"""
from_wxyz(w::Int, x::Int, y::Int, z::Int) = QElt(2*w + z, 2*x + y, y, z)

"""
    to_wxyz(a::QElt) -> (w,x,y,z)

Convert half-coordinates back to  order coordinates.
Checks parity conditions for membership in O:
  A0 - A3 even,  A1 - A2 even.
"""
function to_wxyz(a::QElt)
    A0,A1,A2,A3 = a.A0,a.A1,a.A2,a.A3
    iseven(A0 - A3) || error("to_wxyz(QElt): parity fail (A0-A3 odd); not in O")
    iseven(A1 - A2) || error("to_wxyz(QElt): parity fail (A1-A2 odd); not in O")
    z = A3
    y = A2
    w = (A0 - z) ÷ 2
    x = (A1 - y) ÷ 2
    return (w,x,y,z)
end

Base.:+(a::QElt,b::QElt) = QElt(a.A0+b.A0, a.A1+b.A1, a.A2+b.A2, a.A3+b.A3)
Base.:-(a::QElt) = QElt(-a.A0,-a.A1,-a.A2,-a.A3)
Base.:-(a::QElt,b::QElt) = a + (-b)
Base.:(==)(a::QElt,b::QElt) = (a.A0==b.A0 && a.A1==b.A1 && a.A2==b.A2 && a.A3==b.A3)

"""
    conj_q(a::QElt) -> QElt

Quaternion conjugation in half coords: (A0,A1,A2,A3) -> (A0,-A1,-A2,-A3).
"""
conj_q(a::QElt) = QElt(a.A0, -a.A1, -a.A2, -a.A3)

"""
    trd_q(a::QElt) -> Int

Reduced trace. For half coords, trd(a)=A0.
"""
trd_q(a::QElt) = a.A0

"""
    nrd_q(a::QElt, p::Int) -> Int

Reduced norm:
  nrd(a) = (A0^2 + A1^2 + p(A2^2 + A3^2))/4
"""
function nrd_q(a::QElt, p::Int)
    num = a.A0*a.A0 + a.A1*a.A1 + p*(a.A2*a.A2 + a.A3*a.A3)
    (num % 4 == 0) || error("nrd_q(QElt): numerator not divisible by 4")
    return num ÷ 4
end

"""
    mul_q(a::QElt, b::QElt, p::Int) -> QElt

Multiply in B(-1,-p) using half-coordinate formulas.
Ensures result stays in half-lattice (all N* even).
"""
function mul_q(a::QElt, b::QElt, p::Int)
    A0,A1,A2,A3 = a.A0,a.A1,a.A2,a.A3
    B0,B1,B2,B3 = b.A0,b.A1,b.A2,b.A3

    N0 = A0*B0 - A1*B1 - p*A2*B2 - p*A3*B3
    N1 = A0*B1 + A1*B0 + p*A2*B3 - p*A3*B2
    N2 = A0*B2 - A1*B3 + A2*B0 + A3*B1
    N3 = A0*B3 + A1*B2 - A2*B1 + A3*B0

    (N0 % 2 == 0 && N1 % 2 == 0 && N2 % 2 == 0 && N3 % 2 == 0) ||
        error("mul_q: parity fail (product not in half-lattice)")

    return QElt(N0 ÷ 2, N1 ÷ 2, N2 ÷ 2, N3 ÷ 2)
end

# ============================================================
# 3) Hermitian 2x2 forms and unimodularity
# ============================================================

"""
    Hermitian2x2(u,a,v)

Represents H = [ u   a ; conj(a)  v ] with u,v ∈ Z, a ∈ O.
"""
struct Hermitian2x2
    u::Int
    a::QElt
    v::Int
end

"""
    normalize_uv(H) -> Hermitian2x2

Force positive diagonals and u ≤ v by swapping if needed.
If swap happens: a becomes conj(a).
"""
function normalize_uv(H::Hermitian2x2)
    u,a,v = H.u, H.a, H.v
    if u == 0 || v == 0
        error("normalize_uv: diagonal entries must be > 0; got u=$u, v=$v")
    end
    if u < 0 && v < 0
        u = -u; v = -v; a = -a
    end
    (u > 0 && v > 0) || error("normalize_uv: diagonal entries must be > 0; got u=$u, v=$v")
    if u <= v
        return Hermitian2x2(u,a,v)
    else
        return Hermitian2x2(v, conj_q(a), u)
    end
end

"""
    det_unimodular_ok(H,p) -> Bool

For 2×2 Hermitian over quaternions:
  det(H) = u*v - nrd(a)
Unimodular means det(H)=1.
"""
det_unimodular_ok(H::Hermitian2x2, p::Int) = (H.u*H.v - nrd_q(H.a,p) == 1)

"""
    polarization_params(H) -> (u,v,w,x,y,z)

Return parameters in our convention:
- u,v are diagonal integers
- (w,x,y,z) are order coordinates of a in basis [1, i, (i+j)/2, (1+ij)/2]
"""
function polarization_params(H::Hermitian2x2)
    w,x,y,z = to_wxyz(H.a)
    return (H.u, H.v, w, x, y, z)
end

"""
    show_matrix(H) -> String

Pretty-print in the compact form:
  [ u   (w,x,y,z) ; (conj)   v ]
"""
function show_matrix(H::Hermitian2x2)
    u,v,w,x,y,z = polarization_params(H)
    return @sprintf("[ %d   (w,x,y,z)=(%d,%d,%d,%d) ; (conj)   %d ]", u, w,x,y,z, v)
end

# ============================================================
# 4) Norm-1 units in O by Diophantine equation
# ============================================================

"""
    units_norm1_exact(p) -> Vector{QElt}

Enumerate all norm-1 units in O by solving:
  A0^2 + A1^2 + p(A2^2 + A3^2) = 4
with order parity constraints:
  A0-A3 even,  A1-A2 even.
"""
function units_norm1_exact(p::Int)
    check_p_11mod12(p)
    units = QElt[]
    seen = Set{NTuple{4,Int}}()

    maxA2 = isqrt(4 ÷ p)
    for A2 in -maxA2:maxA2, A3 in -maxA2:maxA2
        rhs = 4 - p*(A2*A2 + A3*A3)
        rhs < 0 && continue
        maxA0 = isqrt(rhs)
        for A0 in -maxA0:maxA0
            A1sq = rhs - A0*A0
            A1 = isqrt(A1sq)
            A1*A1 == A1sq || continue
            for A1v in (A1, -A1)
                iseven(A0 - A3) || continue
                iseven(A1v - A2) || continue
                key = (A0,A1v,A2,A3)
                if !(key in seen)
                    push!(seen, key)
                    push!(units, QElt(A0,A1v,A2,A3))
                end
            end
        end
    end
    return units
end

# ============================================================
# 5) Lemma-3 canonical elimination in our coordinates (w,x,y,z mod v)
# ============================================================

modn(a::Int, m::Int) = mod(a, m)

"""
    lemma3_key_cached!(cache, alpha, v, units, p) -> key
"""
function lemma3_key_cached!(cache::Dict{NTuple{5,Int}, Any},
                            alpha::QElt, v::Int,
                            units::Vector{QElt}, p::Int)
    m = abs(v)
    w,x,y,z = to_wxyz(alpha)
    r = (modn(w,m), modn(x,m), modn(y,m), modn(z,m))
    key0 = (m, r[1], r[2], r[3], r[4])
    if haskey(cache, key0)
        return cache[key0]
    end

    best = r
    best == (0,0,0,0) && return (cache[key0] = (m, best))

    wc,xc,yc,zc = to_wxyz(conj_q(alpha))
    rc = (modn(wc,m), modn(xc,m), modn(yc,m), modn(zc,m))
    if rc < best
        best = rc
        best == (0,0,0,0) && return (cache[key0] = (m, best))
    end

    right = Vector{QElt}(undef, length(units))
    @inbounds for j in 1:length(units)
        right[j] = mul_q(alpha, units[j], p)
    end

    for u1 in units
        for t0 in right
            t = mul_q(u1, t0, p)

            wt,xt,yt,zt = to_wxyz(t)
            rr = (modn(wt,m), modn(xt,m), modn(yt,m), modn(zt,m))
            if rr < best
                best = rr
                best == (0,0,0,0) && return (cache[key0] = (m, best))
            end

            tb = conj_q(t)
            wb,xb,yb,zb = to_wxyz(tb)
            rb = (modn(wb,m), modn(xb,m), modn(yb,m), modn(zb,m))
            if rb < best
                best = rb
                best == (0,0,0,0) && return (cache[key0] = (m, best))
            end
        end
    end

    cache[key0] = (m, best)
    return cache[key0]
end

# ----------------------------------------------------------------
# Lift reduction inside a fixed residue class a0 (mod vO)
# ----------------------------------------------------------------

scale_int_q(a::QElt, m::Int) = QElt(m*a.A0, m*a.A1, m*a.A2, m*a.A3)

"""
    reduce_lift_in_coset(a0, v, p; lift_radius=1, Umax=nothing) -> QElt
"""
function reduce_lift_in_coset(a0::QElt, v::Int, p::Int; lift_radius::Int=1, Umax::Union{Nothing,Int}=nothing)
    w0,x0,y0,z0 = to_wxyz(a0)

    center(c) = (2*c <= v ? 0 : -1)
    bw0 = center(w0); bx0 = center(x0); by0 = center(y0); bz0 = center(z0)
    R = lift_radius

    n0 = nrd_q(a0,p)
    best_div = a0
    best_div_n = n0

    found_u = false
    best_u = a0
    best_u_n = n0
    best_u_u = typemax(Int)
    if Umax !== nothing
        ((n0 + 1) % v == 0) || error("reduce_lift_in_coset: base lift not divisible; internal error")
        u0 = (n0 + 1) ÷ v
        if u0 > 0 && u0 <= Umax
            found_u = true
            best_u_u = u0
        end
    end

    for bw in (bw0-R):(bw0+R), bx in (bx0-R):(bx0+R), by in (by0-R):(by0+R), bz in (bz0-R):(bz0+R)
        β = from_wxyz(bw,bx,by,bz)
        a = a0 + scale_int_q(β, v)
        n = nrd_q(a,p)

        ((n + 1) % v == 0) || continue

        if n < best_div_n
            best_div_n = n
            best_div = a
        end

        if Umax !== nothing
            u = (n + 1) ÷ v
            if u > 0 && u <= Umax
                if !found_u || u < best_u_u || (u == best_u_u && n < best_u_n)
                    found_u = true
                    best_u_u = u
                    best_u_n = n
                    best_u = a
                end
            end
        end
    end

    return (Umax === nothing || !found_u) ? best_div : best_u
end

# ============================================================
# 6) Lemma-3 enumeration
# ============================================================

"""
    lemma3_representatives(p; Vmax, Umax, max_reps, dedup_mode, lift_radius)

Returns (reps, units).
"""
function lemma3_representatives(p::Int;
                               Vmax::Int=60,
                               Vmin::Int=1,
                               Umax::Int=60,
                               max_reps::Int=100_000,
                               dedup_mode::Symbol=:units,
                               lift_radius::Int=1)
    check_p_11mod12(p)
    if Vmin < 1 || Vmin > Vmax
        error("lemma3_representatives: require 1 <= Vmin <= Vmax, got Vmin=$Vmin, Vmax=$Vmax")
    end
    units = units_norm1_exact(p)

    if !(dedup_mode in (:none, :units))
        error("Unknown dedup_mode=$dedup_mode. Use :none or :units.")
    end

    reps = Dict{Any, Hermitian2x2}()
    cache = Dict{NTuple{5,Int}, Any}()

    for v in Vmin:Vmax
        M = 4*v

        mapZ = Vector{Dict{Int, Vector{Int}}}(undef, v)
        for z in 0:v-1
            d = Dict{Int, Vector{Int}}()
            for w in 0:v-1
                A0 = 2*w + z
                r0 = mod(A0*A0, M)
                push!(get!(d, r0, Int[]), w)
            end
            mapZ[z+1] = d
        end

        mapY = Vector{Dict{Int, Vector{Int}}}(undef, v)
        for y in 0:v-1
            d = Dict{Int, Vector{Int}}()
            for x in 0:v-1
                A1 = 2*x + y
                r1 = mod(A1*A1, M)
                push!(get!(d, r1, Int[]), x)
            end
            mapY[y+1] = d
        end

        if dedup_mode == :units
            orbit_rep = Dict{Any, QElt}()

            for y in 0:v-1, z in 0:v-1
                s = mod(p*(y*y + z*z), M)
                target = mod(-4 - s, M)

                d0 = mapZ[z+1]
                d1 = mapY[y+1]

                for (r0, wlist) in d0
                    r1 = mod(target - r0, M)
                    xlist = get(d1, r1, nothing)
                    xlist === nothing && continue

                    for w in wlist, x in xlist
                        a0 = from_wxyz(w,x,y,z)
                        n0 = nrd_q(a0,p)
                        (n0 + 1) % v == 0 || continue

                        key0 = lemma3_key_cached!(cache, a0, v, units, p)
                        if !haskey(orbit_rep, key0)
                            orbit_rep[key0] = a0
                        end
                    end
                end
            end

            for (_k0, a0) in orbit_rep
                a = reduce_lift_in_coset(a0, v, p; lift_radius=lift_radius, Umax=Umax)
                n = nrd_q(a,p)
                ((n + 1) % v == 0) || error("internal: reduce_lift_in_coset returned non-divisible lift")
                u = (n + 1) ÷ v
                (u > 0 && u <= Umax) || continue

                H = normalize_uv(Hermitian2x2(u, a, v))
                det_unimodular_ok(H,p) || continue

                key = lemma3_key_cached!(cache, H.a, H.v, units, p)
                if !haskey(reps, key)
                    reps[key] = H
                    length(reps) >= max_reps && return (collect(values(reps)), units)
                end
            end

        else
            for y in 0:v-1, z in 0:v-1
                s = mod(p*(y*y + z*z), M)
                target = mod(-4 - s, M)

                d0 = mapZ[z+1]
                d1 = mapY[y+1]

                for (r0, wlist) in d0
                    r1 = mod(target - r0, M)
                    xlist = get(d1, r1, nothing)
                    xlist === nothing && continue

                    for w in wlist, x in xlist
                        a0 = from_wxyz(w,x,y,z)
                        n0 = nrd_q(a0,p)
                        (n0 + 1) % v == 0 || continue

                        a = reduce_lift_in_coset(a0, v, p; lift_radius=lift_radius, Umax=Umax)
                        n = nrd_q(a,p)
                        ((n + 1) % v == 0) || error("internal: reduce_lift_in_coset returned non-divisible lift")
                        u = (n + 1) ÷ v
                        (u > 0 && u <= Umax) || continue

                        H = normalize_uv(Hermitian2x2(u, a, v))
                        det_unimodular_ok(H,p) || continue

                        key = (H.u, H.v, H.a.A0, H.a.A1, H.a.A2, H.a.A3)
                        if !haskey(reps, key)
                            reps[key] = H
                            length(reps) >= max_reps && return (collect(values(reps)), units)
                        end
                    end
                end
            end
        end
    end

    return (collect(values(reps)), units)
end

# ============================================================
# 7) Big-integer RandomPolarisation with ZZRingElem
#    + probable-prime acceleration + suggest_params + pretty printer
# ============================================================

# uniform random in [0, N] for N::ZZRingElem
function rand_below(rng::AbstractRNG, N::ZElem)::ZElem
    N < 0 && error("rand_below: require N ≥ 0, got $N")
    iszero(N) && return ZZ(0)
    bitlen = Int(flog(N, 2)) + 1
    while true
        x = ZZ(0)
        remaining = bitlen
        while remaining > 0
            take = min(remaining, 64)
            block = rand(rng, UInt64)
            if take < 64
                block &= (UInt64(1) << take) - 1
            end
            x = (x << take) + ZZ(block)
            remaining -= take
        end
        x <= N && return x
    end
end

# Gaussian integer helpers for sum-of-two-squares
gauss_mul(a::ZElem, b::ZElem, c::ZElem, d::ZElem) = (a*c - b*d, a*d + b*c)

function gauss_pow(a::ZElem, b::ZElem, e::Int)
    ra, rb = ZZ(1), ZZ(0)
    ba, bb = a, b
    k = e
    while k > 0
        if isodd(k)
            ra, rb = gauss_mul(ra, rb, ba, bb)
        end
        ba, bb = gauss_mul(ba, bb, ba, bb)
        k >>= 1
    end
    return (ra, rb)
end

# Prime-only Cornacchia (used by factor-based sum_two_squares)
function cornacchia_two_squares_prime(q::ZElem)
    mod(q, 4) == 1 || error("Cornacchia: need q ≡ 1 (mod 4), got q=$q")
    is_prime(q) || error("Cornacchia: modulus must be prime, got q=$q")
    t = sqrtmod(ZZ(-1), q)  # assumes q prime
    r0, r1 = q, t
    while r1*r1 > q
        r0, r1 = r1, mod(r0, r1)
    end
    x = r1
    y2 = q - x*x
    y = isqrt(y2)
    y*y == y2 || error("Cornacchia failed unexpectedly for q=$q")
    return (x, y)
end

# Cornacchia wrapper
function cornacchia_two_squares(q::ZElem; check_prime::Bool=true)
    mod(q, 4) == 1 || error("Cornacchia: need q ≡ 1 (mod 4), got q=$q")
    if check_prime
        is_prime(q) || error("Cornacchia: modulus must be prime, got q=$q")
    end
    t = sqrtmod(ZZ(-1), q)  # if q composite, may throw
    r0, r1 = q, t
    while r1*r1 > q
        r0, r1 = r1, mod(r0, r1)
    end
    x = r1
    y2 = q - x*x
    y = isqrt(y2)
    y*y == y2 || error("Cornacchia failed for q=$q")
    return (x, y)
end

function sum_two_squares(n::ZElem)
    n < 0 && return nothing
    iszero(n) && return (ZZ(0), ZZ(0))
    fac = factor(n)
    a, b = ZZ(1), ZZ(0)
    for (p, e0) in fac
        e = Int(e0)
        if p == 2
            x, y = gauss_pow(ZZ(1), ZZ(1), e)
            a, b = gauss_mul(a, b, x, y)
        else
            r = mod(p, 4)
            if r == 3
                isodd(e) && return nothing
                a, b = gauss_mul(a, b, p^(e ÷ 2), ZZ(0))
            elseif r == 1
                x, y = cornacchia_two_squares_prime(p)
                x, y = gauss_pow(x, y, e)
                a, b = gauss_mul(a, b, x, y)
            else
                error("unexpected prime residue mod 4 for p=$p")
            end
        end
    end
    return (abs(a), abs(b))
end

is_even(x::ZElem) = iszero(mod(x, 2))

struct QEltZZ
    A0::ZElem
    A1::ZElem
    A2::ZElem
    A3::ZElem
end

conj_q(a::QEltZZ) = QEltZZ(a.A0, -a.A1, -a.A2, -a.A3)
in_order(a::QEltZZ) = is_even(a.A0 - a.A3) && is_even(a.A1 - a.A2)

function nrd_q(a::QEltZZ, p::Int)::ZElem
    pZ = ZZ(p)
    num = a.A0*a.A0 + a.A1*a.A1 + pZ*(a.A2*a.A2 + a.A3*a.A3)
    iszero(mod(num, 4)) || error("nrd_q(QEltZZ): numerator not divisible by 4")
    return div(num, 4)
end

"""
    to_wxyz(a::QEltZZ) -> (w,x,y,z)

Convert half-coordinates (ZZRingElem) back to our order coordinates.
Checks parity constraints.
"""
function to_wxyz(a::QEltZZ)
    A0,A1,A2,A3 = a.A0,a.A1,a.A2,a.A3
    is_even(A0 - A3) || error("to_wxyz(QEltZZ): parity fail (A0-A3 odd); not in O")
    is_even(A1 - A2) || error("to_wxyz(QEltZZ): parity fail (A1-A2 odd); not in O")
    z = A3
    y = A2
    w = div(A0 - z, 2)
    x = div(A1 - y, 2)
    return (w,x,y,z)
end

"""
    represent_integer_in_O(N, p; use_probable_prime=true, verify_with_is_prime=false, ...)

Probabilistic solver for r ∈ O with nrd(r)=N.

Fast path:
- random A2,A3
- rhs = 4N - p(A2^2 + A3^2)
- if rhs ≡ 1 (mod 4) and (probable-)prime, try Cornacchia to get A0^2 + A1^2 = rhs

Options:
- `use_probable_prime=true` uses is_probable_prime(rhs) filter (fast)
- `verify_with_is_prime=true` additionally proves rhs prime using is_prime(rhs), but only after passing probable prime
"""
function represent_integer_in_O(N::ZElem, p::Int;
                                rng::AbstractRNG=Random.default_rng(),
                                max_tries::Int=1200,
                                allow_factor::Bool=false,
                                factor_digit_cutoff::Int=50,
                                use_probable_prime::Bool=true,
                                verify_with_is_prime::Bool=false)::QEltZZ
    N > 0 || error("RepresentInteger: need N>0, got N=$N")
    pZ = ZZ(p)
    fourN = 4*N
    B = isqrt(div(fourN, pZ))  # floor sqrt(4N/p)

    for _ in 1:max_tries
        A2 = iszero(B) ? ZZ(0) : (rand_below(rng, 2*B) - B)
        A3 = iszero(B) ? ZZ(0) : (rand_below(rng, 2*B) - B)
        rhs = fourN - pZ*(A2*A2 + A3*A3)
        rhs <= 0 && continue

        # ---- FAST PATH ----
        if mod(rhs, 4) == 1
            passes = use_probable_prime ? is_probable_prime(rhs) : is_prime(rhs)
            if passes
                if verify_with_is_prime && !is_prime(rhs)
                    continue
                end
                try
                    A0, A1 = cornacchia_two_squares(rhs;
                        check_prime = (!use_probable_prime) || verify_with_is_prime)

                    if is_even(A0 - A3) && is_even(A1 - A2)
                        r = QEltZZ(A0, A1, A2, A3)
                        (in_order(r) && nrd_q(r, p) == N) && return r
                    elseif is_even(A1 - A3) && is_even(A0 - A2)
                        r = QEltZZ(A1, A0, A2, A3)
                        (in_order(r) && nrd_q(r, p) == N) && return r
                    end
                catch
                    # rhs composite / sqrtmod failed / Cornacchia failed
                end
            end
        end

        # ---- OPTIONAL FACTOR FALLBACK ----
        if allow_factor || (Int(flog(rhs, 10)) + 1 <= factor_digit_cutoff)
            xy = sum_two_squares(rhs)
            xy === nothing && continue
            A0, A1 = xy

            if is_even(A0 - A3) && is_even(A1 - A2)
                r = QEltZZ(A0, A1, A2, A3)
                nrd_q(r, p) == N && return r
            elseif is_even(A1 - A3) && is_even(A0 - A2)
                r = QEltZZ(A1, A0, A2, A3)
                nrd_q(r, p) == N && return r
            end
        end
    end

    error("RepresentInteger: failed after $max_tries tries. Increase max_tries or restarts (or allow_factor=true).")
end

struct Hermitian2x2ZZ
    u::ZElem
    a::QEltZZ
    v::ZElem
end

det_unimodular_ok(H::Hermitian2x2ZZ, p::Int) = (H.u*H.v - nrd_q(H.a, p) == 1)

"""
    RandomPolarisation(p; sbound=20, ... ) -> Hermitian2x2ZZ
"""
function RandomPolarisation(p::Int; sbound::Int=20,
                            rng::AbstractRNG=Random.default_rng(),
                            max_tries::Int=1200,
                            allow_factor::Bool=false,
                            use_probable_prime::Bool=true,
                            verify_with_is_prime::Bool=false)::Hermitian2x2ZZ
    check_p_11mod12(p)
    pZ = ZZ(p)

    smax = pZ^sbound
    tmax = pZ^(3*sbound)

    s = isone(smax) ? ZZ(1) : (ZZ(1) + rand_below(rng, smax - 1))
    t = isone(tmax) ? ZZ(1) : (ZZ(1) + rand_below(rng, tmax - 1))

    N = s*t - 1
    r = represent_integer_in_O(N, p; rng=rng, max_tries=max_tries, allow_factor=allow_factor,
                               use_probable_prime=use_probable_prime,
                               verify_with_is_prime=verify_with_is_prime)

    H = Hermitian2x2ZZ(s, r, t)
    det_unimodular_ok(H, p) || error("internal: det != 1")
    return H
end

"""
    as_matrix(H::Hermitian2x2ZZ)

Return the actual 2×2 Julia matrix:
[ s  r
  r̄  t ]
"""
as_matrix(H::Hermitian2x2ZZ) = Any[H.u H.a; conj_q(H.a) H.v]

"""
    suggest_params(p; sbound=20, trials=20, mode=:balanced)

Heuristic tuning for the two-level search.
modes:
- :fast        -> smaller max_tries, more restarts
- :balanced    -> moderate max_tries and restarts
- :few_restarts-> larger max_tries, fewer restarts
"""
function suggest_params(p::Int; sbound::Int=20, trials::Int=20, mode::Symbol=:balanced)
    bits = ndigits(p, base=2)
    etries = ceil(Int, 6 * sbound * bits)

    if mode == :fast
        max_tries = max(800, etries)
        restarts_per_trial = 1000
    elseif mode == :few_restarts
        max_tries = max(2000, 5 * etries)
        restarts_per_trial = 80
    else
        max_tries = max(1200, 3 * etries)
        restarts_per_trial = 300
    end

    return (; bits, sbound, trials, etries, max_tries, restarts_per_trial, mode)
end

"""
    stress_test_retry(p; ...)

Retries when RepresentInteger fails for the sampled (s,t).
"""
function stress_test_retry(p::Int; trials::Int=20, sbound::Int=20,
                          restarts_per_trial::Int=500,
                          max_tries::Int=1200, allow_factor::Bool=false,
                          rng::AbstractRNG=Random.default_rng(),
                          use_probable_prime::Bool=true,
                          verify_with_is_prime::Bool=false,
                          use_suggest::Bool=false,
                          suggest_mode::Symbol=:balanced)

    if use_suggest
        pars = suggest_params(p; sbound=sbound, trials=trials, mode=suggest_mode)
        max_tries = pars.max_tries
        restarts_per_trial = pars.restarts_per_trial
    end

    total_restarts = 0
    for k in 1:trials
        H = nothing
        for _ in 1:restarts_per_trial
            try
                H = RandomPolarisation(p; sbound=sbound, rng=rng,
                                       max_tries=max_tries, allow_factor=allow_factor,
                                       use_probable_prime=use_probable_prime,
                                       verify_with_is_prime=verify_with_is_prime)
                break
            catch err
                msg = sprint(showerror, err)
                if occursin("RepresentInteger", msg)
                    total_restarts += 1
                    continue
                else
                    rethrow()
                end
            end
        end
        H === nothing && error("trial $k: failed after $restarts_per_trial restarts")
        det_unimodular_ok(H, p) || error("trial $k: det check failed")
        in_order(H.a)           || error("trial $k: parity fail (r not in O)")
    end

    println("OK: passed $trials trials (RepresentInteger restarts: $total_restarts)")
    return true
end

# ============================================================
# 8) Batch random polarizations to file
# ============================================================

"""
    write_random_polarizations(p, N; sbound=20, ...)

Generate N distinct random polarizations for prime p and write to a text file.
Output format:
  p = 83
  order basis: [1, i, (i+j)/2, (1+ij)/2]
  # Each line: u v w x y z
  [1] u v w x y z
  [2] ...
"""
function write_random_polarizations(p::Int, N::Int;
                                    sbound::Int=20,
                                    rng::AbstractRNG=Random.default_rng(),
                                    max_tries::Int=1200,
                                    allow_factor::Bool=false,
                                    use_probable_prime::Bool=true,
                                    verify_with_is_prime::Bool=false,
                                    max_attempts::Int=100000000)
    check_p_11mod12(p)
    N > 0 || error("write_random_polarizations: N must be positive; got N=$N")
    max_attempts > 0 || error("write_random_polarizations: max_attempts must be positive")

    seen = Set{NTuple{6,String}}()
    attempts = 0

    filepath = "polz_$(p)_$(N).txt"

    open(filepath, "w") do io
        println(io, "p = $p")
        println(io, "order basis: [1, i, (i+j)/2, (1+ij)/2]")
        println(io, "# Each line: u v w x y z")

        while length(seen) < N
            attempts += 1
            attempts <= max_attempts || error("write_random_polarizations: exceeded max_attempts=$max_attempts before collecting N=$N distinct samples")

            try
                H = RandomPolarisation(p; sbound=sbound, rng=rng,
                                       max_tries=max_tries, allow_factor=allow_factor,
                                       use_probable_prime=use_probable_prime,
                                       verify_with_is_prime=verify_with_is_prime)

                w, x, y, z = to_wxyz(H.a)
                u = H.u
                v = H.v

                key = (string(u), string(v), string(w), string(x), string(y), string(z))
                if !(key in seen)
                    push!(seen, key)
                    idx = length(seen)
                    println(io, "[$idx] $(key[1]) $(key[2]) $(key[3]) $(key[4]) $(key[5]) $(key[6])")
                end
            catch
                continue
            end
        end
    end

    return filepath
end

# ---------- Pretty printing for the big-int matrix ----------

# shorten huge ZZRingElem nicely
function shortZ(x::ZZRingElem; maxchars::Int=60)
    s = string(x)
    (maxchars <= 0 || lastindex(s) <= maxchars) && return s
    head = max(10, maxchars ÷ 2 - 3)
    tail = max(10, maxchars - head - 3)
    return s[1:head] * "..." * s[end-tail+1:end] * " (len=$(length(s)))"
end

"""
    pretty_pol_matrix(H::Hermitian2x2ZZ; maxchars=60, show_conj=false)

Prints:
[ s    (w,x,y,z)=(...)
  (conj)    t ]

If show_conj=true, prints explicit conjugate coordinates in the lower-left.
"""
function pretty_pol_matrix(H::Hermitian2x2ZZ; maxchars::Int=60, show_conj::Bool=false)
    s = shortZ(H.u; maxchars=maxchars)
    t = shortZ(H.v; maxchars=maxchars)

    wxyz = to_wxyz(H.a)
    rstr = "(w,x,y,z)=(" * join(string.(wxyz), ",") * ")"

    if show_conj
        wxyzb = to_wxyz(conj_q(H.a))
        rbstr = "(w,x,y,z)=(" * join(string.(wxyzb), ",") * ")"
        return "[ $s    $rstr ;\n  $rbstr    $t ]"
    else
        return "[ $s    $rstr ;\n  (conj)    $t ]"
    end
end


####------ How to use----------
# include("RandomPol2.jl")
# p = 23
# write_random_polarizations(p, p^3)

