import numpy as np
from numpy.polynomial.legendre import legval,legder

def mode(rho,z,a,n,k):
 r=np.hypot(rho,z);t=z/r;s=rho/r
 c=np.zeros(n+1);c[-1]=1
 P=legval(t,c);D=legval(t,legder(c));q=2-n if k==0 else -n
 # psi = coefficient (velocity) * a^2 * (r/a)^q * angular
 f=a*a*(r/a)**q
 ur=n*(n+1)*f/r**2*P; ut=-q*f/r**2*s*D
 return np.array([ur*s+ut*t,ur*t-ut*s])
def fit(a=(1,1),gap=1,L=12):
 a=np.array(a);z=np.array([0,sum(a)+gap]);nc=3*L+4
 ts=np.cos(np.pi*(np.arange(nc)+.5)/nc)
 rows=[];rhs=[]
 for al in range(2):
  for t in ts:
   rho=a[al]*np.sqrt(1-t*t);zz=z[al]+a[al]*t
   block=np.column_stack([mode(rho,zz-z[b],a[b],n,k) for b in range(2) for n in range(1,L+1) for k in range(2)])
   rows.extend(block);rhs.extend([[0,0],np.eye(2)[al]])
 A=np.array(rows); B=np.array(rhs);sc=np.linalg.norm(A,axis=0)
 C=np.linalg.lstsq(A/sc,B,rcond=None)[0]/sc[:,None]
 R=8*np.pi*np.array([a[0]*C[0],a[1]*C[2*L]])
 err=0
 for al in range(2):
  for t in np.linspace(-.99999,.99999,200):
   block=np.column_stack([mode(a[al]*np.sqrt(1-t*t),z[al]+a[al]*t-z[b],a[b],n,k) for b in range(2) for n in range(1,L+1) for k in range(2)])
   err=max(err,np.max(np.abs(block@C-np.vstack((np.zeros(2),np.eye(2)[al])))))
 return R,err,np.linalg.cond(A/sc)
if __name__=='__main__':
 for gap in [2,1,.2,.05]:
  for L in [4,8,16,24]:
   R,e,c=fit(gap=gap,L=L);print(gap,L,R[0]/(6*np.pi),e,c)
