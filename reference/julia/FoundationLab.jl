module FoundationLab
using LinearAlgebra
export stf, sphere_rule, translation_field, strain_field, rotation_field,
       surface_moments, faxen_force, faxen_torque, faxen_stresslet
const I3=Matrix{Float64}(I,3,3)
stf(A)=(A+transpose(A))/2-tr(A)/3*I3

# Gauss-Legendre in t=cos(theta), trapezoid in azimuth; weights integrate dOmega.
function sphere_rule(nt=8,np=16)
    nt>=2 && np>=3 || throw(ArgumentError("increase quadrature order"))
    b=[k/sqrt(4k^2-1) for k in 1:nt-1]
    ev=eigen(SymTridiagonal(zeros(nt),b))
    t=ev.values; w=2*ev.vectors[1,:].^2
    normals=Vector{Vector{Float64}}(); weights=Float64[]
    for j in 1:nt, k in 0:np-1
        phi=2pi*k/np; s=sqrt(1-t[j]^2)
        push!(normals,[s*cos(phi),s*sin(phi),t[j]])
        push!(weights,w[j]*2pi/np)
    end
    normals,weights
end

function translation_field(r,U; a=1.0,mu=1.0)
    R=norm(r); R>=a*(1-1e-12) && a>0 && mu>0 || throw(ArgumentError("exterior points required"))
    n=r/R; q=dot(r,U)
    f=3a/(4R)+a^3/(4R^3); g=3a/(4R^3)-3a^3/(4R^5)
    fp=-3a/(4R^2)-3a^3/(4R^4); gp=-9a/(4R^4)+15a^3/(4R^6)
    u=f*U+g*q*r
    L=fp*U*transpose(n)+gp*q*r*transpose(n)+g*q*I3+g*r*transpose(U)
    p=3mu*a*q/(2R^3)
    sigma=-p*I3+mu*(L+transpose(L))
    (;u,p,L,sigma)
end

function strain_field(r,E; a=1.0,mu=1.0)
    isapprox(E,transpose(E);atol=1e-12) && abs(tr(E))<1e-12 || throw(ArgumentError("E must be symmetric traceless"))
    R=norm(r); R>=a*(1-1e-12) && a>0 && mu>0 || throw(ArgumentError("exterior points required"))
    n=r/R; v=E*r; q=dot(r,v)
    f=1-a^5/R^5; g=-5a^3/(2R^5)+5a^5/(2R^7)
    fp=5a^5/R^6; gp=25a^3/(2R^6)-35a^5/(2R^8)
    u=f*v+g*q*r
    L=fp*v*transpose(n)+f*E+gp*q*r*transpose(n)+g*(2r*transpose(v)+q*I3)
    p=-5mu*a^3*q/R^5
    sigma=-p*I3+mu*(L+transpose(L))
    (;u,p,L,sigma)
end

function rotation_field(r,Omega; a=1.0,mu=1.0)
    R=norm(r); R>=a*(1-1e-12) && a>0 && mu>0 || throw(ArgumentError("exterior points required"))
    v=cross(Omega,r)
    C=[0.0 -Omega[3] Omega[2];Omega[3] 0.0 -Omega[1];-Omega[2] Omega[1] 0.0]
    u=a^3*v/R^3
    L=a^3*(C/R^3-3v*transpose(r)/R^5)
    p=0.0; sigma=mu*(L+transpose(L))
    (;u,p,L,sigma)
end

# t=sigma*n is fluid-on-particle traction. S is the suspension-stress convention.
function surface_moments(field; a=1.0,nt=8,np=16)
    normals,weights=sphere_rule(nt,np)
    F=zeros(3); T=zeros(3); S=zeros(3,3)
    for (n,w) in zip(normals,weights)
        r=a*n; t=field(r).sigma*n; dA=a^2*w
        F+=dA*t; T+=dA*cross(r,t); S+=dA*stf(r*transpose(t))
    end
    (;F,T,S)
end
faxen_force(U,uamb,lapu; a=1.0,mu=1.0)=-6pi*mu*a*(U-uamb-a^2/6*lapu)
faxen_torque(Omega,curlu; a=1.0,mu=1.0)=-8pi*mu*a^3*(Omega-curlu/2)
faxen_stresslet(E,lapE; a=1.0,mu=1.0)=20pi/3*mu*a^3*(E+a^2/10*lapE)
end
