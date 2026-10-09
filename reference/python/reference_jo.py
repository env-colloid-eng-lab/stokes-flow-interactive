from functools import lru_cache
from fractions import Fraction as Q
from math import comb
@lru_cache(None)
def P(n,p,q):
 if min(p,q)<0 or n<1:return Q(0)
 if p==0 and q==0:return Q(int(n==1))
 return sum((comb(n+s,n)*(Q(n*(2*n+1)*(2*n*s-n-s+2),2*(n+1)*(2*s-1)*(n+s))*P(s,q-s,p-n+1)-Q(n*(2*n-1),2*(n+1))*P(s,q-s,p-n-1)-Q(n*(4*n*n-1),2*(n+1)*(2*s+1))*V(s,q-s-2,p-n+1)) for s in range(1,q+1)),Q(0))
@lru_cache(None)
def V(n,p,q):
 if min(p,q)<0 or n<1:return Q(0)
 if p==0 and q==0:return Q(int(n==1))
 return P(n,p,q)-Q(2*n,(n+1)*(2*n+3))*sum((comb(n+s,n)*P(s,q-s,p-n-1) for s in range(1,q+1)),Q(0))
def polynomials(K):return [[2**k*P(1,k-q,q) for q in range(k+1)] for k in range(K+1)]
def series(s,lam,polys):
 terms=[sum(float(c)*lam**q for q,c in enumerate(cs))/((1+lam)*s)**k for k,cs in enumerate(polys)]
 return sum(terms[::2]),-2/(1+lam)*sum(terms[1::2])
if __name__=='__main__':
 f=polynomials(30)
 for k in range(6):print(k,[str(x) for x in f[k]])
 for s in [4,3,2.2]:print(s,series(s,1,f))
