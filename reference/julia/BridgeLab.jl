module BridgeLab
using LinearAlgebra
export reflection_coefficients, reflected_mobility, single_sphere_scatter

"""
Chapter 13.3: generate f0,...,f4 by inverting the mobility POWER SERIES.
lambda is a radius ratio; a Rational input keeps all arithmetic exact.
Removing the stresslet changes f4. This is not an all-order JO generator.
"""
function reflection_coefficients(lambda=1//1; stresslet=true)
    lambda>0 || throw(ArgumentError("positive radius ratio required"))
    T=typeof(lambda/one(lambda))
    h=T(1)/T(2)
    M=[zeros(T,2,2) for k in 0:4]
    M[1]=diagm(T[1,1/lambda])
    M[2]=T[0 3h;3h 0]
    M[4]=T[0 -(1+lambda^2)*h;-(1+lambda^2)*h 0]
    if stresslet
        M[5]=-T(15)/T(4)*diagm(T[lambda^3,1])
    end
    R=[zeros(T,2,2) for k in 0:4]
    R[1]=diagm(T[1,lambda])
    for k in 1:4
        product=zeros(T,2,2)
        for j in 1:k
            product+=M[j+1]*R[k-j+1]
        end
        R[k+1]=-R[1]*product
    end
    [iseven(k) ? T(2)^k*R[k+1][1,1] :
                 -T(2)^k*R[k+1][1,2] for k in 0:4]
end

"""
Far-field axial mobility retaining the displayed r^-4 self correction.
The positive-definiteness of this truncation is NOT guaranteed near contact.
"""
function reflected_mobility(a1,a2,r;mu=1.0,stresslet=true)
    min(a1,a2,mu)>0 && r>a1+a2 || throw(ArgumentError("separated positive spheres required"))
    c=(1-(a1^2+a2^2)/(3r^2))/(4pi*mu*r)
    m1=1/(6pi*mu*a1);m2=1/(6pi*mu*a2)
    if stresslet
        m1-=5a2^3/(8pi*mu*r^4)
        m2-=5a1^3/(8pi*mu*r^4)
    end
    [m1 c;c m2]
end

"""Chapter 12.2: one incoming streamfunction mode and specified surface data."""
function single_sphere_scatter(l,C,D,b,w)
    l>=1 || throw(ArgumentError("l must be positive"))
    [1.0 1.0;2-l -l]\[b-C-D,w-(l+1)*C-(l+3)*D]
end
end
