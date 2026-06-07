# run_experiments.jl
#
# Julia companion to run_experiments.py: runs cCRM, MAP, and DR on
#   (A) equality-constrained ellipsoidal feasibility (reduced to R^2), and
#   (B) fixed-trace spectral feasibility in S^n,
# and prints the distance histories and quadratic ratios.
#
# Uses only the standard library (LinearAlgebra, Printf).  No plotting here;
# Run:  julia run_experiments.jl

using LinearAlgebra
using Printf
using Random

# ---------------------------------------------------------------- projections
function proj_ellipsoid(z, M, center, radius; tol = 1e-15, itmax = 200)
    y = M * (z .- center)
    if norm(y) <= radius
        return copy(z)
    end
    G = Symmetric(M' * M)
    w = eigvals(G); V = eigvecs(G)
    t = V' * (z .- center)
    λ = 0.0
    for _ in 1:itmax
        denom = 1 .+ λ .* w
        f  = sum(w .* t .^ 2 ./ denom .^ 2) - radius^2
        df = sum(-2 .* w .^ 2 .* t .^ 2 ./ denom .^ 3)
        λn = λ - f / df
        λn < 0 && (λn = λ / 2)
        if abs(λn - λ) <= tol * (1 + abs(λ))
            λ = λn; break
        end
        λ = λn
    end
    s = t ./ (1 .+ λ .* w)
    return center .+ V * s
end

proj_ball(z, c, r) = (d = z .- c; n = norm(d); n <= r ? copy(z) : c .+ (r / n) .* d)

function proj_lambda_max(S, a)
    F = eigen(Symmetric(S))
    return F.vectors * Diagonal(min.(F.values, a)) * F.vectors'
end

# ---------------------------------------------------------------- circumcenter
function circumcenter(a, b, c)
    u = b .- a; v = c .- a
    uu = dot(u, u); uv = dot(u, v); vv = dot(v, v)
    det2 = uu * vv - uv * uv
    rhs = 0.5 .* (uu, vv)
    if abs(det2) > 0 && det2 / (uu * vv + eps()) > 1e-12
        α = (rhs[1] * vv - rhs[2] * uv) / det2
        β = (uu * rhs[2] - uv * rhs[1]) / det2
        return a .+ α .* u .+ β .* v
    end
    return (a .+ b .+ c) ./ 3
end

# ------------------------------------------------------------- cCRM/MAP/DR
function operators(PX, PY)
    RX(z) = 2 .* PX(z) .- z
    RY(z) = 2 .* PY(z) .- z
    centralized(z) = (w = PY(PX(z)); 0.5 .* (w .+ PX(w)))
    Tc(z) = (zc = centralized(z); circumcenter(zc, RX(zc), RY(zc)))
    Tm(z) = PY(PX(z))
    Td(z) = z .+ (PY(RX(z)) .- PX(z))
    return Tc, Tm, Td
end

# distance to X ∩ Y via Dykstra
function dist_to_inter(z, PX, PY; itmax = 4000, tol = 1e-14)
    x = copy(z); p = zero(z); q = zero(z)
    for _ in 1:itmax
        y  = PX(x .+ p); p = x .+ p .- y
        xn = PY(y .+ q); q = y .+ q .- xn
        if norm(xn .- x) <= tol * (1 + norm(xn))
            x = xn; break
        end
        x = xn
    end
    return norm(z .- x)
end

function run_iter(T, z0, dist; steps = 30, tol = 1e-14)
    z = copy(z0); ds = [dist(z)]
    for _ in 1:steps
        z = T(z); d = dist(z); push!(ds, d)
        d <= tol && break
    end
    return ds
end

function show_run(name, Tc, z0, dist; steps = 30)
    println("="^60); println(name); println("-"^60)
    ds = run_iter(Tc, z0, dist; steps = steps)
    @printf("%3s %14s %12s\n", "k", "dist", "quad ratio")
    for k in 1:min(length(ds) - 1, 10)
        if ds[k] > 1e-11
            @printf("%3d %14.4e %12.4f\n", k - 1, ds[k], ds[k + 1] / ds[k]^2)
        end
    end
end

# ------------------------------------------------------------------ Problem A
function ellipsoid_problem()
    aX = [2.2, 1.0]; cX = [0.0, 0.0]
    aY = [1.0, 2.2]; cY = [3.0, 0.0]
    MX = Diagonal(1 ./ aX); MY = Diagonal(1 ./ aY)
    PX(z) = proj_ellipsoid(z, MX, cX, 1.0)
    PY(z) = proj_ellipsoid(z, MY, cY, 1.0)
    z0 = [1.5, 3.0]
    return PX, PY, z0, "Equality-constrained ellipsoids (d=2)"
end

# ------------------------------------------------------------------ Problem B
function spectral_problem(n = 8)
    a = 0.35; r = 0.4
    A = randn(n, n); Shat = 0.1 * (A + A') / 2
    I_n = Matrix{Float64}(I, n, n)
    proj_trace(S) = S .- (tr(S) - 1) / n .* I_n
    Shat = proj_trace(Shat)
    ballF(S, C, rad) = (d = S - C; nd = norm(d); nd <= rad ? S : C .+ (rad / nd) .* d)
    function dykstra(S, pA, pB; itmax = 1500, tol = 1e-14)
        x = copy(S); p = zero(S); q = zero(S)
        for _ in 1:itmax
            y  = pA(x + p); p = x + p - y
            xn = pB(y + q); q = y + q - xn
            norm(xn - x) <= tol && (x = xn; break)
            x = xn
        end
        return x
    end
    PXm(S) = dykstra(S, M -> proj_lambda_max(M, a), proj_trace)
    PYm(S) = dykstra(S, M -> ballF(M, Shat, r), proj_trace)
    vecf(S) = vec(S); matf(x) = reshape(x, n, n)
    PX(x) = vecf(PXm(matf(x)))
    PY(x) = vecf(PYm(matf(x)))
    z0 = vecf(proj_trace(0.2 * (randn(n, n) + randn(n, n)') ))
    return PX, PY, z0, "Fixed-trace spectral feasibility (n=$n)"
end

# ------------------------------------------------------------------- main
function main()
    Random.seed!(20260606)
    PX, PY, z0, name = ellipsoid_problem()
    Tc, _, _ = operators(PX, PY)
    show_run(name, Tc, z0, z -> dist_to_inter(z, PX, PY); steps = 25)

    PX, PY, z0, name = spectral_problem(8)
    Tc, _, _ = operators(PX, PY)
    show_run(name, Tc, z0, z -> dist_to_inter(z, PX, PY); steps = 40)
end

main()
