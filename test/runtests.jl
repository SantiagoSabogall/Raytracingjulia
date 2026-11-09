using Test
using Raytracingjulia
using StaticArrays
using LinearAlgebra
using DifferentialEquations
using DelimitedFiles
using Dierckx

# ============================================================================
# Helpers
# ============================================================================

function cov_matrix(bh, x)
    g_tt, g_rr, g_θθ, g_φφ, g_tφ = metric(bh, x)
    return [g_tt    0.0    0.0    g_tφ ;
            0.0     g_rr   0.0    0.0  ;
            0.0     0.0    g_θθ   0.0  ;
            g_tφ    0.0    0.0    g_φφ]
end

function con_matrix(bh, x)
    gtt, grr, gθθ, gφφ, gtφ = inverse_metric(bh, x)
    return [gtt    0.0    0.0    gtφ ;
            0.0    grr    0.0    0.0  ;
            0.0    0.0    gθθ    0.0  ;
            gtφ    0.0    0.0    gφφ]
end

# g^μν g_νλ = δ^μ_λ  (duality of covariant/contravariant metric)
function test_duality(bh, x)
    G = cov_matrix(bh, x)
    Ginv = con_matrix(bh, x)
    return Ginv * G ≈ I(4)
end

# Hamiltonian H = ½ g^μν p_μ p_ν along a solution state
function hamiltonian(bh, u)
    gtt, grr, gθθ, gφφ, gtφ = inverse_metric(bh, u)
    pt = u[5]; pr = u[6]; pθ = u[7]; pφ = u[8]
    return 0.5 * (gtt*pt^2 + grr*pr^2 + gθθ*pθ^2 + gφφ*pφ^2 + 2*gtφ*pt*pφ)
end

# ---------------------------------------------------------------------------
# 1. Metric duality + Hamiltonian conservation on a real photon
# ---------------------------------------------------------------------------
@testset "Metric identities (g g⁻¹ = I)" begin
    sch = SchwarzschildBH()
    kerr = KerrBH(0.7)
    kerrpf = KerrPFDMBH(0.7, 0.05)
    kersc = KerrScalaronBH(0.7, 5.0)

    for bh in (sch, kerr, kerrpf, kersc)
        for (r, θ) in ((5.0, 0.3), (3.2, 1.2), (8.0, 1.5707))
            x = @SVector [0.0, r, θ, 0.0]
            @test test_duality(bh, x)
        end
    end
end

@testset "Hamiltonian conservation along geodesic" begin
    sch = SchwarzschildBH()
    kerr = KerrBH(0.7)
    kerrpf = KerrPFDMBH(0.7, 0.05)
    kersc = KerrScalaronBH(0.7, 5.0)

    for bh in (sch, kerr, kerrpf, kersc)
        # Build a detector and grab one real null photon
        det = Detector(100.0, 85.0*pi/180, 12.0; x_pixels=15, ratio="1:1")
        img = Image(bh, SimpleDisk(bh), det)
        create_photons!(img)
        p = img.photon_list[4]   # off-axis photon (more interesting)

        @test isfinite(hamiltonian(bh, p.iC))  # initial state must be null & finite

        tspan = (0.0, -30.0)
        cb = ContinuousCallback((u,t,i)->u[2]-(bh.EH+1e-5), (i)->terminate!(i))
        sol = solve(ODEProblem(geodesics, p.iC, tspan, bh), Tsit5();
                    reltol=1e-9, abstol=1e-9, saveat=0.5, callback=cb, verbose=false)

        Hmax = 0.0
        for u in sol.u
            h = hamiltonian(bh, u)
            isfinite(h) && (Hmax = max(Hmax, abs(h)))
        end
        # Null geodesics must keep H ≈ 0 along the whole trajectory
        @test Hmax < 1e-5
    end
end

# ---------------------------------------------------------------------------
# 2. Analytic limits
# ---------------------------------------------------------------------------
@testset "Limits: PFDM k→0 recovers Kerr" begin
    a = 0.5
    kerr = KerrBH(a)
    kerrpf = KerrPFDMBH(a, 1e-6)   # tiny k → Kerr

    @test kerrpf.EH ≈ kerr.EH atol=1e-3
    @test kerrpf.ISCOco ≈ kerr.ISCOco atol=1e-3

    for (r, θ) in ((5.0, 0.3), (3.5, 1.0))
        x = @SVector [0.0, r, θ, 0.0]
        @test SVector(metric(kerrpf, x)) ≈ SVector(metric(kerr, x)) rtol=1e-4
        @test SVector(inverse_metric(kerrpf, x)) ≈ SVector(inverse_metric(kerr, x)) rtol=1e-4
    end
end

@testset "Limits: Scalaron Mψ→∞ recovers Kerr" begin
    a = 0.5
    kerr = KerrBH(a)
    kersc = KerrScalaronBH(a, 1000.0)   # large scalaron mass → Kerr

    @test kersc.EH ≈ kerr.EH atol=1e-3
    @test kersc.ISCOco ≈ kerr.ISCOco atol=1e-2

    for (r, θ) in ((5.0, 0.3), (3.5, 1.0))
        x = @SVector [0.0, r, θ, 0.0]
        @test SVector(metric(kersc, x)) ≈ SVector(metric(kerr, x)) rtol=1e-4
        @test SVector(inverse_metric(kersc, x)) ≈ SVector(inverse_metric(kerr, x)) rtol=1e-4
    end
end

@testset "Limits: Schwarzschild ≡ Kerr(a=0)" begin
    sch = SchwarzschildBH()
    kerr = KerrBH(0.0)
    x = @SVector [0.0, 4.0, 0.7, 0.0]
    @test SVector(metric(sch, x)) ≈ SVector(metric(kerr, x))
    @test SVector(inverse_metric(sch, x)) ≈ SVector(inverse_metric(kerr, x))
    @test sch.EH == kerr.EH
    @test sch.ISCOco == kerr.ISCOco
end

# ---------------------------------------------------------------------------
# 3. Kerr sanity: EH and ISCO analytic values
# ---------------------------------------------------------------------------
@testset "Kerr analytic EH / ISCO" begin
    for a in (0.0, 0.3, 0.9)
        bh = KerrBH(a)
        @test bh.EH ≈ 1 + sqrt(1 - a^2)
        # co-rotating ISCO must lie inside the Schwarzschild value 6M
        @test 1 < bh.ISCOco ≤ 6.0
        if a == 0.0
            @test bh.ISCOco == bh.ISCOcounter == 6.0
        else
            @test bh.ISCOco < bh.ISCOcounter
        end
    end
end

# ---------------------------------------------------------------------------
# 4. Schwarzschild shadow: critical impact parameter b = 3√3 M
# ---------------------------------------------------------------------------
@testset "Schwarzschild shadow boundary (b_crit = 3√3)" begin
    sch = SchwarzschildBH()
    D = 100.0

    # Face-on detector → shadow is a circle of radius 3√3 M in screen space.
    det = Detector(D, 0.0, 12.0; x_pixels=15, ratio="1:1")
    img = Image(sch, SimpleDisk(sch), det)
    create_photons!(img)

    function fate(bh, p)
        r_max = 1.2 * D
        cb = CallbackSet(
            ContinuousCallback((u,t,i)->u[2]-(bh.EH+1e-5), (i)->terminate!(i)),
            ContinuousCallback((u,t,i)->u[2]-r_max, (i)->terminate!(i))
        )
        sol = solve(ODEProblem(geodesics, p.iC, (0.0, -1.5*D), bh), Tsit5();
                    reltol=1e-8, abstol=1e-8, save_everystep=false, callback=cb, verbose=false)
        return sol.u[end][2] < bh.EH + 1e-3 ? :captured : :escaped
    end

    # Circular shadow ⇒ boundary in impact-parameter radius ρ = √(α²+β²).
    ρ_captured   = 0.0
    ρ_uncaptured = Inf
    for p in img.photon_list
        ρ = sqrt(p.alpha^2 + p.beta^2)
        if fate(sch, p) == :captured
            ρ_captured = max(ρ_captured, ρ)
        else
            ρ_uncaptured = min(ρ_uncaptured, ρ)
        end
    end

    b_crit = 3.0 * sqrt(3.0)   # ≈ 5.1962 (Schwarzschild photon sphere)
    @test ρ_captured < b_crit < ρ_uncaptured
end

# ---------------------------------------------------------------------------
# 5. Numerical Schwarzschild: interpolated N(r) must reproduce the metric
# ---------------------------------------------------------------------------
@testset "Numerical Schwarzschild reproduces analytic metric" begin
    # Build a temp N(r) = 1 - 2/r, dN/dr = 2/r² table for the standard BH.
    mktempdir() do dir
        rs = range(2.1, 40.0, length=400)
        Np = [1.0 - 2.0/r for r in rs]
        dN = [2.0/r^2 for r in rs]

        fN  = joinpath(dir, "N.txt")
        fdN = joinpath(dir, "dN.txt")
        writedlm(fN,  hcat(collect(rs), Np))
        writedlm(fdN, hcat(collect(rs), dN))

        bh = SchwarzschildNumBH(fN, fdN)
        sch = SchwarzschildBH()

        for (r, θ) in ((4.0, 0.5), (10.0, 1.0), (3.0, 0.9))
            x = @SVector [0.0, r, θ, 0.0]
            @test test_duality(bh, x)
            @test SVector(metric(bh, x)) ≈ SVector(metric(sch, x)) rtol=1e-5
            @test SVector(inverse_metric(bh, x)) ≈ SVector(inverse_metric(sch, x)) rtol=1e-5
        end
    end
end

# ---------------------------------------------------------------------------
# 6. Doppler shift sanity: far, slow emitter → g ≈ 1
# ---------------------------------------------------------------------------
@testset "Doppler shift far from BH ≈ 1" begin
    sch = SchwarzschildBH()
    # A photon state far out on the equatorial plane with no angular momentum.
    u = @SVector [0.0, 80.0, pi/2, 0.0, -1.0, 0.1, 0.0, 0.0]
    I0 = 1.0
    @test doppler_shift(u, I0, sch) ≈ I0 rtol=1e-2
end
