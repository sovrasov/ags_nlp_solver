#include <cmath>
#include "Map.hpp"

int n1,nexp,l,iq,iu[10],iv[10];

void xyd ( double *xx, int m, double y[], int n)
{
 /* calculate preimage x  for nearest level  m center to y */
 /* (x - left boundary point of level m interval)          */

 double x,r1;
 double r;
 int iw[11];
 int i,j,it,is;
 void numbr ( int * );

 n1=n-1;
 for ( nexp=1,i=0; i<n; i++ ) {
	nexp*=2; iw[i]=1;
 }
 r=0.5;
 r1=1.0;
 x=0.0;
 it=0;
 for ( j=0;  j<m; j++ ) {
	r*=0.5;
	for ( i=0; i<n; i++ ) {
	  iu[i] = ( y[i]<0 ) ? -1 : 1;
	  y[i]-=r*iu[i];
	  iu[i]*=iw[i];
	}
	i=iu[0];
	iu[0]=iu[it];
	iu[it]=i;
	numbr ( &is );
	i=iv[0];
	iv[0]=iv[it];
	iv[it]=i;
	for ( i=0; i<n; i++ )
	  iw[i]=-iw[i]*iv[i];
	if ( l == 0 ) l=it;
	  else if ( l == it ) l=0;
	it=l;
	r1=r1/nexp;
	x+=r1*is;
 }
 *xx=x;
}
void numbr ( int *iss)
{
 /* calculate s(u)=is,l(u)=l,v(u)=iv by u=iu */

  int i,n,is,iff,k1,k2,l1;

  n=n1+1;
  iff=nexp;
  is=0;
  k1=-1;
  for ( i=0; i<n; i++ ) {
	 iff=iff/2;
	 k2=-k1*iu[i];
	 iv[i]=iu[i];
	 k1=k2;
	 if ( k2<0 ) l1=i;
		else { is+=iff; l=i; }
  }
  if ( is == 0 ) l=n1;
	 else {
		iv[n1]=-iv[n1];
		if ( is == (nexp-1) ) l=n1;
		  else if ( l1 == n1 ) iv[l]=-iv[l];
			 else l=l1;
  }
  *iss=is;
}

__declspec(dllexport) void invmad(int m, double xp[], int kp,
	int *kxx, double p[], int n, int incr)
{
	/* calculate kx preimage p node */
	/*   node type mapping m level  */

	double mne, d1, dd, x, dr, del;// , convers;
	double r, d, u[10], y[10];
	int i, k, kx, nexp;
	void xyd(double *, int, double *, int);

	kx = 0;
	kp--;
	for (nexp = 1, i = 0; i<n; i++) {
		nexp *= 2; u[i] = -1.0;
	}
	dr = nexp;
	for (mne = 1, r = 0.5, i = 0; i<m; i++) {
		mne *= dr; r *= 0.5;
	}
	dr = mne / nexp;

	dr=dr-fmod(dr,1.0);  
	//dr = (dr>0) ? floor(dr) : ceil(dr);

	del = 1. / (mne - dr);
	d1 = del*(incr + 0.5);
	for (kx = -1; kx<kp;) {
		for (i = 0; i<n; i++) {       /* label 2 */
			d = p[i];
			y[i] = d - r*u[i];
		}
		for (i = 0; (i<n) && (fabs(y[i]) < 0.5); i++);
		if (i >= n) {
			xyd(&x, m, y, n);
			dr = x*mne;
			dd=dr-fmod(dr,1.0); 
			//dd = (dr>0) ? floor(dr) : ceil(dr);
			dr = dd / nexp;
			dd=dd-dr+fmod(dr,1.0); 
			//dd = dd - ((dr>0) ? floor(dr) : ceil(dr));
			x = dd*del;
			if (kx>kp) break;
			k = kx++;                     /* label 9 */
			if (kx == 0) xp[0] = x;
			else {
				while (k >= 0) {
					dr = fabs(x - xp[k]);     /* label 11 */
					if (dr <= d1) {
						for (kx--; k<kx; k++, xp[k] = xp[k + 1]);
						goto m6;
					}
					else
						if (x <= xp[k]) {
							xp[k + 1] = xp[k]; k--;
						}
						else break;
				}
				xp[k + 1] = x;
			}
		}
	m6: for (i = n - 1; (i >= 0) && (u[i] = (u[i] <= 0.0) ? 1 : -1)<0; i--);
		if (i<0) break;
	}
	*kxx = ++kx;

}