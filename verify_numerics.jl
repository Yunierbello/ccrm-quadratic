# verify_numerics.jl
#
# Julia companion to verify_numerics.py for
# "Q-quadratic convergence of the centralized circumcentered-reflection
#  method under a relative interior condition" (Y. Bello-Cruz).
#
# Reproduces, in Float64 (and BigFloat where high precision matters):
#   (1) Table 1: two overlapping discs (isometry to R^2), quad ratio -> 0.5551.
#   (2) Error-bound modulus omega ~ 0.252 and the rate accounting.
#   (3) Lemma 2.7(i): the 1/8 descent constant.
#   (4) Lemma 2.5 application: dist(T,X) / [(1/2) kappa |T-pX|^2] -> 1.
#   (5) Example 4.4 (beta=0 epigraph): asymptotic coefficients and 1 - 1/alpha.
#
# Run:  julia verify_numerics.jl
#
# No external packages required (uses only LinearAlgebra from the stdlib).

using LinearAlgebra
using Printf

# ---------------------------------------------------------------- helpers
vnorm(v) = sqrt(v[1]^2 + v[2]^2)

function circumcenter(a, b, c)
    u = b - a; v = c - a
    uu = dot(u, u); uv = dot(u, v); vv = dot(v, v)
    det2 = uu * vv - uv * uv
    rhs = 0.5 .* [uu, vv]
    if abs(det2) > 0 && det2 / (uu * vv + eps()) > 1e-12
        α = (rhs[1] * vv - rhs[2] * uv) / det2
        β = (uu * rhs[2] - uv * rhs[1]) / det2
        return a .+ α .* u .+ β .* v
    end
    return (a .+ b .+ c) ./ 3
end

# ------------------------------------------------------------------- discs
function disc_example(T = BigFloat)
    r  = T(2)
    cX = T[0, 0]
    cY = T[sqrt(T(15)), 0]

    projball(z, c) = begin
        d = z .- c; n = vnorm(d)
        n <= r ? copy(z) : c .+ (r / n) .* d
    end
    PX(z) = projball(z, cX)
    PY(z) = projball(z, cY)
    RX(z) = 2 .* PX(z) .- z
    RY(z) = 2 .* PY(z) .- z
    zC(z) = begin
        w = PY(PX(z)); pw = PX(w); 0.5 .* (w .+ pw)
    end
    Tmap(z) = (c = zC(z); circumcenter(c, RX(c), RY(c)))
    dball(z, c) = max(T(0), vnorm(z .- c) - r)

    zbar = T[sqrt(T(15)) / 2, T(1) / 2]
    z    = T[sqrt(T(15)) / 2, T(4)]

    println("(1) Table 1 -- disc example")
    @printf("   %2s %14s %12s %12s\n", "k", "|z-zb|", "ratio", "quad ratio")
    prev = vnorm(z .- zbar)
    for k in 0:4
        c = zC(z); Tz = circumcenter(c, RX(c), RY(c))
        e = vnorm(Tz .- zbar)
        @printf("   %2d %14.5e %12.5e %12.5f\n", k, Float64(prev),
                Float64(e / prev), Float64(e / prev^2))
        z = Tz; prev = e
    end

    println("\n(2) omega accounting at the stabilized iterate")
    z = T[sqrt(T(15)) / 2, T(4)]
    for _ in 1:3
        c = zC(z); z = circumcenter(c, RX(c), RY(c))
    end
    c = zC(z); Tz = circumcenter(c, RX(c), RY(c)); pX = PX(c)
    dist = vnorm(z .- zbar); nxt = vnorm(Tz .- zbar)
    maxd = max(dball(Tz, cX), dball(Tz, cY))
    @printf("   max(dist(T,X),dist(T,Y))/dist^2 = %.5f\n", Float64(maxd / dist^2))
    @printf("   dist(next)/dist^2               = %.5f\n", Float64(nxt / dist^2))
    @printf("   implied omega                   = %.5f\n", Float64(maxd / nxt))
    @printf("   contraction |T-pX|/dist         = %.5f\n", Float64(vnorm(Tz .- pX) / dist))

    println("\n(3) Lemma 2.7(i): |T-s|^2 <= |z-s|^2 - (1/8)|z-T|^2")
    s = T[sqrt(T(15)) / 2, 0]; z = T[3, 2]
    ok = true
    for _ in 0:5
        c = zC(z); Tz = circumcenter(c, RX(c), RY(c))
        lhs = vnorm(Tz .- s)^2
        rhs = vnorm(z .- s)^2 - (T(1) / 8) * vnorm(z .- Tz)^2
        ok &= (lhs <= rhs + T(10)^(-40)); z = Tz
    end
    println("   holds at every step: ", ok)

    println("\n(4) Lemma 2.5 application: dist(T,X)/[(1/2)kappa|T-pX|^2] -> 1")
    κ = T(1) / 2; z = T[sqrt(T(15)) / 2, T(4)]
    for k in 0:4
        c = zC(z); Tz = circumcenter(c, RX(c), RY(c)); pX = PX(c)
        rhs = T(1) / 2 * κ * vnorm(Tz .- pX)^2
        @printf("   k=%d: ratio = %.8f\n", k, rhs > 0 ? Float64(dball(Tz, cX) / rhs) : NaN)
        z = Tz
    end
end

# --------------------------------------------------------------- epigraph
function proj_axis(x, α)              # solve x = u + α u^{2α-1}
    u = x
    for _ in 1:200
        f  = u + α * u^(2α - 1) - x
        df = 1 + α * (2α - 1) * u^(2α - 2)
        un = u - f / df
        if abs(un - u) < eps(typeof(x)) * 100
            u = un; break
        end
        u = un
    end
    return u
end

function proj_pt(px, py, α)
    if py >= abs(px)^α
        return (px, py)
    end
    s = px == 0 ? one(px) : sign(px); ax = abs(px)
    t = ax > 0 ? ax : eps(typeof(px))
    for _ in 1:200
        f  = (t - ax) + α * t^(α - 1) * (t^α - py)
        df = 1 + α * (α - 1) * t^(α - 2) * (t^α - py) + α^2 * t^(2α - 2)
        tn = t - f / df
        tn <= 0 && (tn = t / 2)
        if abs(tn - t) < eps(typeof(px)) * 100
            t = tn; break
        end
        t = tn
    end
    return (s * t, t^α)
end

function epigraph_example(T = BigFloat)
    println("\n(5) Example 4.4 (beta=0): coefficients and rate 1 - 1/alpha")
    for α in T[2, 3, 5//2]
        x = T(1) / 1000
        u = proj_axis(x, α); v = proj_axis(u, α)
        a = (u + v) / 2; h = v^α / 2
        (px, _) = proj_pt(a, h, α); p = px
        x1 = p + p^α * (p^α - h) / (p - a)
        # circumcenter equidistance residual for (x1,0) vs z_C, R_X(z_C)
        resid = abs(((x1 - a)^2 + h^2) - ((x1 - (2p - a))^2 + (2 * p^α - h)^2))
        @printf("   alpha=%.2f: circ resid=%.3e | rate %.8f (target %.8f)\n",
                Float64(α), Float64(resid), Float64(x1 / x), Float64(1 - 1 / α))
    end
end

# ------------------------------------------------------------------- main
function main()
    setprecision(BigFloat, 200)
    disc_example(BigFloat)
    epigraph_example(BigFloat)
end

main()
