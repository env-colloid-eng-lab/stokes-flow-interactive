using LinearAlgebra, Test, DelimitedFiles, Printf
include("FoundationLab.jl");using .FoundationLab
include("TwoSphereLab.jl");using .TwoSphereLab
BLAS.set_num_threads(1)
mkpath("output")
results=String[]
record(s)=(println(s);push!(results,s))
@testset "Surface traction and moments" begin
    a=1.3;mu=.7;U=[.4,-.2,.7];Om=[.2,.8,-.3]
    normals,weights=sphere_rule()
    @test isapprox(sum(weights),4pi;rtol=1e-14)
    mt=surface_moments(r->translation_field(r,U;a=a,mu=mu);a=a)
    mr=surface_moments(r->rotation_field(r,Om;a=a,mu=mu);a=a)
    E=[.2 .4 -.1;.4 -.3 .2;-.1 .2 .1]
    ms=surface_moments(r->strain_field(r,E;a=a,mu=mu);a=a)
    @test mt.F ≈ -6pi*mu*a*U
    @test mr.T ≈ -8pi*mu*a^3*Om
    @test ms.S ≈ 20pi/3*mu*a^3*E
    @test norm(ms.F)<1e-12 && norm(ms.T)<1e-12
    for n in normals
        @test translation_field(a*n,U;a=a,mu=mu).u ≈ U
        @test norm(strain_field(a*n,E;a=a,mu=mu).u)<1e-12
        @test strain_field(a*n,E;a=a,mu=mu).sigma*n ≈ 5mu*E*n
    end
    Fp=zeros(3);Fv=zeros(3)
    for (n,w) in zip(normals,weights)
        fld=translation_field(a*n,U;a=a,mu=mu)
        Fp+=-fld.p*n*a^2*w
        Fv+=mu*(fld.L+transpose(fld.L))*n*a^2*w
    end
    @test Fp≈-2pi*mu*a*U
    @test Fv≈-4pi*mu*a*U
    record(@sprintf("Translation traction integral relative error: %.3e",norm(mt.F+6pi*mu*a*U)/norm(mt.F)))
    record(@sprintf("Strain stresslet relative error: %.3e",norm(ms.S-20pi/3*mu*a^3*E)/norm(ms.S)))
end
@testset "Faxen curvature and Einstein viscosity" begin
    normals,weights=sphere_rule(10,20);a=1.2;mu=.8
    mean_u=sum(w*[(a*n[2])^2,0,0] for (n,w) in zip(normals,weights))/(4pi)
    @test mean_u≈[a^2/3,0,0]
    cubic(r)=(5dot(r,r)*[2r[1],-2r[2],0]-4r*(r[1]^2-r[2]^2))/(42mu)
    S=5mu*a^2*sum(w*stf(cubic(a*n)*transpose(n)) for (n,w) in zip(normals,weights))
    lapE=diagm([2.0,-2.0,0.0])/mu
    @test S≈faxen_stresslet(zeros(3,3),lapE;a=a,mu=mu)
    E=[0.0 .5 0;.5 0 0;0 0 0]
    ms=surface_moments(r->strain_field(r,E;a=a,mu=mu);a=a)
    phi=.02;num=phi/(4pi*a^3/3)
    relative=(mu+num*ms.S[1,2])/mu
    @test relative≈1+2.5phi
    record(@sprintf("Einstein coefficient from surface integral: %.12f",(relative-1)/phi))
end
polys_exact=jo_polynomials(10)
@testset "JO exact coefficients" begin
    targets=[[1],[0,3],[0,9,0],[0,-4,27,-4],[0,-24,81,36,0],[0,0,72,243,72,0]]
    for k in 0:5;@test polys_exact[k+1]==targets[k+1];end
    @test polys_exact[11]==[0,0,2304,20736,42804,115849,76176,39264,20736,2304,0]
end
record("JO coefficients f0-f5 and f10: exact rational agreement.")
polys=jo_polynomials(60;T=Float64)
rows=Vector{Vector{Float64}}()
@testset "JO versus independent boundary collocation" begin
    for (a1,a2,r) in [(1.,1.,4.),(1.,1.,3.),(1.,.5,3.)]
        jo=jo_resistance(a1,a2,r,polys)
        c=axial_collocation([0.,r],[a1,a2];L=24,ncheck=151)
        relative=norm(jo-c.R)/norm(c.R)
        @test relative<2e-9
        @test c.boundary_error<2e-9
        @test norm(jo-transpose(jo))/norm(jo)<2e-13
        @test minimum(eigvals(Symmetric(jo)))>0
        record(@sprintf("a1=%.1f a2=%.1f r=%.1f: JO/collocation error %.3e; boundary %.3e",a1,a2,r,relative,c.boundary_error))
    end
    for gap in [2.,1.,.2],L in [4,8,16,24,32]
        c=axial_collocation([0.,2+gap],[1.,1.];L=L,ncheck=151)
        push!(rows,[gap,L,c.R[1,1]/(6pi),c.R[1,2]/(6pi),c.boundary_error,c.condition])
    end
    one=axial_collocation([0.],[1.3];L=4)
    @test one.R[1,1]≈6pi*1.3
end
writedlm("output/axial_convergence.csv",vcat(reshape(["gap","L","X11","X12","boundary_residual","condition"],1,:),reduce(vcat,transpose.(rows))),',')
seriesrows=[[K,jo_XA(3.,1.,jo_polynomials(K;T=Float64))...] for K in [4,10,20,30,40,60]]
writedlm("output/jo_convergence.csv",vcat(["K" "X11" "X12"],reduce(vcat,transpose.(seriesrows))),',')
tr=pair_motion([0.,4.],[1.,-1.],[1.,1.],polys;dt=.1,steps=100,min_gap=1.)
writedlm("output/jo_pair_motion.csv",vcat(["time" "z1" "z2" "U1" "U2" "gap"],tr),',')
open("output/revision_validation.txt","w") do io
    println(io,"Julia ",VERSION)
    foreach(s->println(io,s),results)
end
println("Revision labs passed.")
