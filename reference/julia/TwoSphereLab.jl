module TwoSphereLab
using LinearAlgebra
export jo_polynomials, jo_XA, jo_resistance, axial_mode, axial_collocation, pair_motion

# Jeffrey-Onishi (1984), eqs. (3.6)-(3.9), (3.15).
# No coefficients are read from a table. T=Rational{BigInt} produces exact polynomials.
function jo_polynomials(K::Int; T=Rational{BigInt})
    K>=0 || throw(ArgumentError("K must be nonnegative"))
    pc=Dict{NTuple{3,Int},T}(); vc=Dict{NTuple{3,Int},T}()
    choose(n,s)=T(binomial(big(n+s),n))
    function P(n,p,q)
        (n<1 || p<0 || q<0) && return zero(T)
        p==0 && q==0 && return T(n==1)
        get!(pc,(n,p,q)) do
            v=zero(T)
            for s in 1:q
                c1=T(n*(2n+1)*(2n*s-n-s+2))/T(2*(n+1)*(2s-1)*(n+s))
                c2=T(n*(2n-1))/T(2*(n+1))
                c3=T(n*(4n^2-1))/T(2*(n+1)*(2s+1))
                v+=choose(n,s)*(c1*P(s,q-s,p-n+1)-c2*P(s,q-s,p-n-1)-c3*V(s,q-s-2,p-n+1))
            end
            v
        end
    end
    function V(n,p,q)
        (n<1 || p<0 || q<0) && return zero(T)
        p==0 && q==0 && return T(n==1)
        get!(vc,(n,p,q)) do
            v=P(n,p,q)
            fac=T(2n)/T((n+1)*(2n+3))
            for s in 1:q
                v-=fac*choose(n,s)*P(s,q-s,p-n-1)
            end
            v
        end
    end
    [[T(big(2)^k)*P(1,k-q,q) for q in 0:k] for k in 0:K]
end

# s=2r/(a1+a2)>2, lambda=a2/a1. This is the truncated FAR-FIELD series.
# Near-contact convergence must be checked; this is not an asymptotically matched evaluator.
function jo_XA(s,lambda,polys)
    s>2 && lambda>0 || throw(ArgumentError("s>2, lambda>0 required"))
    x11=0.0; x12=0.0
    for k in 0:length(polys)-1
        fk=evalpoly(float(lambda),Float64.(polys[k+1]))
        term=fk/((1+lambda)*s)^k
        if iseven(k); x11+=term; else; x12-=2/(1+lambda)*term; end
    end
    x11,x12
end
function jo_resistance(a1,a2,r,polys; mu=1.0)
    a1>0 && a2>0 && mu>0 || throw(ArgumentError("positive radii and viscosity required"))
    s=2r/(a1+a2); lambda=a2/a1
    x11,x12=jo_XA(s,lambda,polys); x22,x21=jo_XA(s,1/lambda,polys)
    [6pi*mu*a1*x11 3pi*mu*(a1+a2)*x12;
     3pi*mu*(a1+a2)*x21 6pi*mu*a2*x22]
end

function legendre_and_derivative(n,t)
    n==0 && return 1.0,0.0
    p0,p1=1.0,t; d0,d1=0.0,1.0
    for k in 1:n-1
        p2=((2k+1)*t*p1-k*p0)/(k+1)
        d2=((2k+1)*(p1+t*d1)-k*d0)/(k+1)
        p0,p1=p1,p2;d0,d1=d1,d2
    end
    p1,d1
end

# Axisymmetric exterior streamfunction mode; k=0 pressure mode, k=1 potential mode.
# Coefficients have units of velocity: psi=c*a^2*(r/a)^q*(1-t^2)*P_n'(t).
function axial_mode(rho,z,a,n,k)
    R=hypot(rho,z);t=z/R;s=rho/R
    P,D=legendre_and_derivative(n,t)
    q=k==0 ? 2-n : -n
    f=a^2*(R/a)^q
    ur=n*(n+1)*f/R^2*P;ut=-q*f/R^2*s*D
    [ur*s+ut*t,ur*t-ut*s]
end
function axial_block(rho,z,centers,radii,L)
    hcat([axial_mode(rho,z-centers[b],radii[b],n,k)
          for b in 1:length(radii) for n in 1:L for k in 0:1]...)
end

# N coaxial spheres, axial translation only. Non-axisymmetric motion is excluded.
# Solve all unit-velocity boundary problems together, without imposing R symmetry.
function axial_collocation(centers,radii; L=16,nc=3L+4,mu=1.0,ncheck=121)
    N=length(radii)
    length(centers)==N && N>=1 && L>=1 && nc>=L+1 && mu>0 || throw(ArgumentError("invalid dimensions/orders"))
    all(radii .> 0) || throw(ArgumentError("positive radii required"))
    for i in 1:N, j in i+1:N
        abs(centers[i]-centers[j])>radii[i]+radii[j] || throw(ArgumentError("spheres overlap or touch"))
    end
    A=zeros(2N*nc,2N*L);rhs=zeros(2N*nc,N)
    for alpha in 1:N, j in 1:nc
        t=cos(pi*(j-.5)/nc);rho=radii[alpha]*sqrt(1-t^2);z=centers[alpha]+radii[alpha]*t
        row=2*((alpha-1)*nc+j)-1
        A[row:row+1,:]=axial_block(rho,z,centers,radii,L)
        rhs[row+1,alpha]=1
    end
    sc=sqrt.(sum(abs2,A;dims=1))
    As=A./sc
    C=(As\rhs)./transpose(sc)
    R=zeros(N,N)
    for alpha in 1:N
        R[alpha,:]=8pi*mu*radii[alpha]*C[2L*(alpha-1)+1,:]
    end
    err=0.0
    for alpha in 1:N,t in range(-.999999,.999999;length=ncheck)
        rho=radii[alpha]*sqrt(1-t^2);z=centers[alpha]+radii[alpha]*t
        target=zeros(2,N);target[2,alpha]=1
        err=max(err,maximum(abs.(axial_block(rho,z,centers,radii,L)*C-target)))
    end
    (;R,C,boundary_error=err,condition=cond(As),L,nc,centers,radii)
end

# Deterministic midpoint stepping. Stop before leaving the series' checked gap range.
function pair_motion(z0,F,radii,polys; dt=.1,steps=100,min_gap=1.0,mu=1.0)
    z=copy(float.(z0));a1,a2=radii
    z[2]>z[1] || throw(ArgumentError("centers must be ordered"))
    rows=Vector{Vector{Float64}}()
    function speed(x)
        R=jo_resistance(a1,a2,x[2]-x[1],polys;mu=mu)
        R\F
    end
    for j in 0:steps
        gap=z[2]-z[1]-a1-a2
        gap>=min_gap || break
        u=speed(z);push!(rows,[j*dt;z;u;gap])
        j==steps && break
        mid=z+dt/2*u
        mid[2]-mid[1]-a1-a2>=min_gap || break
        z=z+dt*speed(mid)
    end
    reduce(vcat,transpose.(rows))
end
end
