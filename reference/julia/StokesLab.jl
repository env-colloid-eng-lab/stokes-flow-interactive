module StokesLab
using LinearAlgebra, Random
export eye3, projector, oseen, grad_oseen, lap_oseen, pressure,
       dipole, sphere_flow, rpy_pair, mobility, normal_resistance,
       velocity, heun_step, brownian_increment, stf_basis, stf5, unstf5

eye3() = Matrix{Float64}(I, 3, 3)

function projector(q)
    dot(q,q) > 0 || throw(ArgumentError("q must be nonzero"))
    eye3() - q * transpose(q) / dot(q,q)
end

function oseen(r; mu=1.0)
    R = norm(r)
    R > 0 && mu > 0 || throw(ArgumentError("r != 0, mu > 0 required"))
    (eye3()/R + r*transpose(r)/R^3)/(8*pi*mu)
end

# H[i,j,k] = partial_k G_ij, not partial_i G_jk.
function grad_oseen(r; mu=1.0)
    R = norm(r)
    R > 0 && mu > 0 || throw(ArgumentError("r != 0, mu > 0 required"))
    H = zeros(3,3,3)
    for i in 1:3, j in 1:3, k in 1:3
        H[i,j,k] = (-(i==j)*r[k]/R^3 + (i==k)*r[j]/R^3 +
                    (j==k)*r[i]/R^3 - 3*r[i]*r[j]*r[k]/R^5)/(8*pi*mu)
    end
    H
end

function lap_oseen(r; mu=1.0)
    R = norm(r)
    R > 0 && mu > 0 || throw(ArgumentError("r != 0, mu > 0 required"))
    (2*eye3()/R^3 - 6*r*transpose(r)/R^5)/(8*pi*mu)
end

pressure(r,f) = dot(f,r)/(4*pi*norm(r)^3) # r != 0

# +f at +d/2 and -f at -d/2: D = f*d', u_i = -H_ijk D_jk.
function dipole(r,D; mu=1.0)
    H = grad_oseen(r; mu=mu)
    [-sum(H[i,j,k]*D[j,k] for j in 1:3, k in 1:3) for i in 1:3]
end

function sphere_flow(r,U; a=1.0, mu=1.0)
    a > 0 && norm(r) >= a || throw(ArgumentError("r must be outside sphere"))
    (oseen(r;mu=mu) + a^2/6*lap_oseen(r;mu=mu))*(6*pi*mu*a*U)
end

# Equal-radius, free-space, translational RPY; overlap extension is not contact physics.
function rpy_pair(r; a=1.0, mu=1.0)
    a > 0 && mu > 0 || throw(ArgumentError("a, mu > 0 required"))
    R = norm(r)
    R == 0 && return eye3()/(6*pi*mu*a)
    nn = (r*transpose(r))/R^2
    if R >= 2*a
        ((1+2*a^2/(3*R^2))*eye3() + (1-2*a^2/R^2)*nn)/(8*pi*mu*R)
    else
        ((1-9*R/(32*a))*eye3() + 3*R/(32*a)*nn)/(6*pi*mu*a)
    end
end

function mobility(X; a=1.0, mu=1.0, model=:rpy)
    size(X,1) == 3 || throw(ArgumentError("X must be 3 by N"))
    a > 0 && mu > 0 || throw(ArgumentError("a, mu > 0 required"))
    model in (:rpy,:oseen,:self) || throw(ArgumentError("unknown model"))
    N = size(X,2)
    M = zeros(3*N,3*N)
    for p in 1:N
        ip = (3*p-2):(3*p)
        M[ip,ip] = eye3()/(6*pi*mu*a)
        for q in (p+1):N
            iq = (3*q-2):(3*q)
            r = X[:,p]-X[:,q]
            block = model == :rpy ? rpy_pair(r;a=a,mu=mu) :
                    model == :oseen ? oseen(r;mu=mu) : zeros(3,3)
            M[ip,iq] = block
            M[iq,ip] = transpose(block)
        end
    end
    M
end

# Teaching model: ONLY a normal leading-order lubrication correction.
# This is not the matched two-body correction of full Stokesian Dynamics.
function normal_resistance(X; a=1.0, mu=1.0, hc=0.2*a)
    hc > 0 || throw(ArgumentError("hc > 0 required"))
    M = mobility(X;a=a,mu=mu)
    n = size(M,1)
    Rmat = cholesky(Symmetric(M)) \ Matrix{Float64}(I,n,n)
    for p in 1:size(X,2), q in (p+1):size(X,2)
        r = X[:,p]-X[:,q]
        gap = norm(r)-2*a
        gap > 0 || throw(ArgumentError("lubrication model requires gap > 0"))
        if gap < hc
            normal = r/norm(r)
            b = zeros(n)
            b[(3*p-2):(3*p)] = normal
            b[(3*q-2):(3*q)] = -normal
            zeta = 3*pi*mu*a^2/2*(1/gap - 1/hc)
            Rmat += zeta*(b*transpose(b))
        end
    end
    Rmat
end

function velocity(X,F; a=1.0, mu=1.0, lubrication=false)
    size(F) == size(X) || throw(DimensionMismatch("F and X must have same shape"))
    u = lubrication ? normal_resistance(X;a=a,mu=mu) \ vec(F) :
                      mobility(X;a=a,mu=mu)*vec(F)
    reshape(u,3,size(X,2))
end

# vfun is a deterministic X -> U function; no stochastic Heun interpretation implied.
function heun_step(X,vfun,dt)
    dt > 0 || throw(ArgumentError("dt > 0 required"))
    k1 = vfun(X)
    k2 = vfun(X + dt*k1)
    X + dt/2*(k1+k2)
end

# Frozen-configuration noise sample. Does not include force or thermal drift.
function brownian_increment(M,kBT,dt; rng=Random.default_rng())
    kBT >= 0 && dt > 0 || throw(ArgumentError("kBT >= 0, dt > 0 required"))
    L = cholesky(Symmetric(M)).L
    sqrt(2*kBT*dt)*(L*randn(rng,size(M,1)))
end

# Orthonormal basis of symmetric traceless 3 by 3 tensors, under A:B=sum(A.*B).
function stf_basis()
    B = zeros(3,3,5)
    B[:,:,1] = diagm([1.0,-1.0,0.0])/sqrt(2)
    B[:,:,2] = diagm([1.0,1.0,-2.0])/sqrt(6)
    B[1,2,3] = B[2,1,3] = 1/sqrt(2)
    B[1,3,4] = B[3,1,4] = 1/sqrt(2)
    B[2,3,5] = B[3,2,5] = 1/sqrt(2)
    B
end
stf5(A) = [sum(A .* stf_basis()[:,:,k]) for k in 1:5]
unstf5(v) = sum(v[k]*stf_basis()[:,:,k] for k in 1:5)
end
