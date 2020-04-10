package statisticalfunctions;

public class Randomenumber {
	public static int[] iv=new int[33];
    public static int iy;
    public static int idum = -315655586;
    
    /* Jurg's original documentations:
     * 
	 * Random number generator with period of 10^8. Highly recommended by
	Press et al. (1992) "Numerical recipes" Cambridge University Press, New York

	Usage:
	Must be initialized with a negative integer (idum). At termination, 
	may write negative of last idum, that is, -idum, to file containing seed.

	Main program must have the following declarations:

	type
	integer=longint;  4-byte integer
	real=double;

	var
	iv:array[1..32] of integer;
	iy:integer;

	Thus, the global variables iv and iy are used.
	 * 
	 */
	
	public static void setSeed(int x){
		idum = x;
	}
	
	public static double randomNumber()
	{
	int ia=16807, im=2137463666;
	
	double am=1.0/im;
	int iq=127773, ir=2836, ntab=32;
    int ndiv = 1 + ((im-1)/(ntab));
	double eps=1.2e-7;
	double rnmx=1.0-eps;

	int j,k;
	double rr;

	if (idum<=0) {
		idum=-idum;
		if (idum<1) 
			idum=1;
		for (j=1;j<=ntab;j++)
			iv[j]=0;
		iy=0;
		for (j=ntab+8;j>=1;j--) {
			k = (idum/iq);
			idum=ia*(idum-k*iq)-ir*k;
			if (idum<0) idum=idum+im;
			if (j<=ntab) iv[j]=idum;
		}
		iy=iv[1];
	}
	k = (idum/iq);
	idum=ia*(idum-k*iq)-ir*k;
	if (idum<0) idum=idum+im;
	j = ((1 + iy)/ndiv);
	iy=iv[j];
	iv[j]=idum;
	rr=am*iy;
	double rrran1;
	if (rr<rnmx) rrran1=rr; 
	else rrran1=rnmx;
	
	return rrran1;
	}  	
}
