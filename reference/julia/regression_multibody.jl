include("StokesLab.jl")
using .StokesLab, LinearAlgebra, Random, Statistics, Test, Printf

# Julia 1.10+; standard libraries only. Run: julia regression_multibody.jl
# Files are generated in this script's output directory.
out = joinpath(@__DIR__,"output")
mkpath(out)
rng = MersenneTwister(20260913)

@testset "01 tensor and projection" begin
    q = [1.0,2.0,3.0]
    P = projector(q)
    @test P*P ≈ P
    @test norm(P*q) < 1e-14
    @test sort(eigvals(Symmetric(P))) ≈ [0,1,1] atol=1e-14
    A = [2.0 1 3; 1 -1 4; 3 4 -1]
    v = stf5(A)
    @test unstf5(v) ≈ A
    @test dot(v,v) ≈ sum(abs2,A)
end

@testset "02 Oseen derivatives and covariance" begin
    r = [1.3,-0.7,2.1]; f = [0.4,0.8,-0.3]
    G = oseen(r); H = grad_oseen(r)
    h = 1e-5
    for k in 1:3
        e = eye3()[:,k]
        FD = (oseen(r+h*e)-oseen(r-h*e))/(2*h)
        @test isapprox(FD,H[:,:,k];rtol=1e-8,atol=1e-11)
    end
    @test norm([sum(H[i,j,i] for i in 1:3) for j in 1:3]) < 1e-14
    gradp = [(pressure(r+h*eye3()[:,k],f)-pressure(r-h*eye3()[:,k],f))/(2*h) for k in 1:3]
    @test norm(lap_oseen(r)*f-gradp) < 1e-10
    Q = Matrix(qr(randn(rng,3,3)).Q)
    @test oseen(Q*r) ≈ Q*G*transpose(Q)
end

@testset "03 dipole and moving sphere" begin
    r = [1.3,-0.7,2.1]; f = [0.4,0.8,-0.3]; d = [0.2,0.4,0.1]
    D = f*transpose(d)
    e = 1e-3
    finite = (oseen(r-e*d/2)-oseen(r+e*d/2))*f/e
    @test isapprox(finite,dipole(r,D);rtol=1e-7)
    S = (D+transpose(D))/2-tr(D)/3*eye3()
    C = (D-transpose(D))/2
    @test dipole(r,D) ≈ dipole(r,S)+dipole(r,C)
    torque = cross(d,f)
    @test dipole(r,C) ≈ cross(torque,r)/(8*pi*norm(r)^3)
    @test norm(dipole(r,eye3())) < 1e-14
    U = [0.2,0.4,0.6]
    @test sphere_flow(r/norm(r),U) ≈ U
end

@testset "04 two spheres and finite size" begin
    X = [0.0 4; 0 0; 0 0]
    M = mobility(X)
    m0 = 1/(6*pi)
    pair = rpy_pair([4.0,0,0])
    cp = (1-2/(3*4^2))/(4*pi*4)
    ct = (1+2/(3*4^2))/(8*pi*4)
    @test diag(pair) ≈ [cp,ct,ct]
    expected = sort([m0+cp,m0-cp,m0+ct,m0-ct,m0+ct,m0-ct])
    @test eigvals(Symmetric(M)) ≈ expected
    @test isposdef(Symmetric(M))
    @test rpy_pair([2-1e-8,0,0]) ≈ rpy_pair([2+1e-8,0,0])
    F = [0.0 0; 0 0; -1 -1]
    U = velocity(X,F)
    @test U[:,1] ≈ U[:,2]
    @test abs(U[3,1]) > m0
    @test dot(vec(F),M*vec(F)) > 0
end

@testset "05 three-body time step convergence" begin
    X0 = [-3.0 3 0; 0 0 0; 0 0 4]
    F = [0.0 0 0; 0 0 0; -1 -1 -1]
    vfun(X) = velocity(X,F)
    function integrate(n)
        X = copy(X0); dt = 1.0/n
        for k in 1:n
            X = heun_step(X,vfun,dt)
        end
        X
    end
    x1=integrate(10); x2=integrate(20); x3=integrate(40)
    ratio=norm(x1-x2)/norm(x2-x3)
    @printf("Heun convergence ratio (expected about 4): %.5f\n",ratio)
    @test 3.5 < ratio < 4.5
    X=copy(X0); dt=0.05
    open(joinpath(out,"three_spheres.csv"),"w") do io
        println(io,"t,x1,y1,z1,x2,y2,z2,x3,y3,z3")
        for k in 0:100
            println(io,join(vcat(k*dt,vec(X)),","))
            k < 100 && (X=heun_step(X,vfun,dt))
        end
    end
end

@testset "06 normal lubrication limit" begin
    F = [1.0 -1; 0 0; 0 0]
    open(joinpath(out,"lubrication.csv"),"w") do io
        println(io,"gap,speed_rpy,speed_normal_model,asymptotic")
        for h in 10.0 .^ range(-5,-0.5,length=60)
            X = [0.0 2+h; 0 0; 0 0]
            Ur=velocity(X,F)
            Ul=velocity(X,F;lubrication=true)
            vr=Ur[1,1]-Ur[1,2]; vl=Ul[1,1]-Ul[1,2]
            vasym=2*h/(3*pi)
            @test vl <= vr*(1+1e-10)
            @test vl > 0
            if h <= 1e-4
                @test isapprox(vl,vasym;rtol=0.003)
            end
            println(io,join([h,vr,vl,vasym],","))
        end
    end
end

@testset "07 correlated Brownian increments" begin
    X = [0.0 4; 0 0; 0 0]; M=mobility(X)
    kBT=1.0; dt=0.01; nsamples=100000
    L=cholesky(Symmetric(M)).L
    # At fixed X, verify covariance, not a many-step equilibrium trajectory.
    Y=sqrt(2*kBT*dt)*L*randn(rng,6,nsamples)
    Yc=Y .- mean(Y,dims=2)
    C=Yc*transpose(Yc)/(nsamples-1)
    target=2*kBT*dt*M
    err=norm(C-target)/norm(target)
    @printf("Brownian covariance relative error: %.5f\n",err)
    @test err < 0.025
    @test norm(mean(Y,dims=2)) < 0.001
    open(joinpath(out,"brownian_covariance.csv"),"w") do io
        println(io,"i,j,empirical,theory")
        for i in 1:6, j in 1:6
            println(io,join([i,j,C[i,j],target[i,j]],","))
        end
    end
end

@testset "08 sedimentation equilibrium and thermal drift" begin
    # Analytic local zero-flux check: M(z)=1+z; kBT=1; mg=2.
    z=0.7; h=1e-6; c(z)=exp(-2*z); mob(z)=1+z
    cp=(c(z+h)-c(z-h))/(2*h)
    J=mob(z)*(-2*c(z)-cp)
    @test abs(J) < 1e-9
    # If drift M' is omitted, stationary c is exp(-2z)/M, not Boltzmann.
    wrong(z)=exp(-2*z)/mob(z)
    dc=(mob(z+h)*wrong(z+h)-mob(z-h)*wrong(z-h))/(2*h)
    @test abs(-2*mob(z)*wrong(z)-dc) < 1e-9
end

println("All lab checks completed. CSV outputs: ",out)
