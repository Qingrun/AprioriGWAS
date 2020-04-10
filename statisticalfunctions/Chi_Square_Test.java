package statisticalfunctions;

public class Chi_Square_Test {
	public static double chi2pr(double x2, int ndf){

		/*	From Chi-squire value to probability
		 * 
		 * {INCLUDE file for Linkage Utility programs.
		 *	 Computes the error probability, chipr = 1-f, where f is the
		 *	 chi-square distribution function with ndf degrees of freedom
		 *	 evaluated at x2.       Function required: pnorm.
		 *	 Follows formulas 26.4.4 and 26.4.21 in Abramowitz and Stegun,
		 *	 "Handbook of mathematical functions", Dover, New York 1968}
		 */
		double z,p,x,sum,re,ch,chp, chipr;

		int  n1, n2;	
		if ((x2>0.0) && (ndf>0)){	 
			if (ndf==1){	  
				x=Math.sqrt(x2);
				chipr=2.0*pnorm(x);
			}
			else{  
				if (ndf==2)
				{ chipr=Math.exp(-0.5*x2);} //{ndf=2}
				else{               //{ndf>2}
					n1=(ndf-1)/ 2;
					n2=(ndf-2)/ 2;
					if (n1==n2){
						// {ndf is even and >2}
						sum=1.0;
						re=0.5*x2;
						ch=1.0;
						for (int i=1; i<= n2;i++){
							ch=ch*re/i;
							sum=sum+ch;
						}
						chipr=Math.exp(-re)*sum;
					}
					else{             //{ndf is odd and >1}
						ch=Math.sqrt(x2);
						z=0.39894228*Math.exp(-0.5*x2);
						p=pnorm(ch);
						if (ndf==3)
						{chipr=2.0*(p+z*ch);} //{ndf=3}
						else{             // {ndf odd and >3}
							chp=ch;
							re=1.0;
							for (int i=2; i<= n1;i++){
								re=re+2.0;
								chp=chp*x2/re;
								ch=ch+chp;
							}
							chipr=2.0*(p+z*ch);
						}
					}
				}
			}
		} 
		else	chipr =1.0;
		return chipr;
	}

	public static double chiSquareValueLR(double[][] table){

		//	double table[][] = preCheck4ZeroRowAndColumn(table2);

		int length = table.length;
		//System.out.println("length:"+length);
		int width =  table[0].length;
		double sumtablerow[] = new double[length];
		double sumtablecolum[] = new double[width];
		double sum=0;

		for (int i=0;i<length;i++)
			for (int j=0;j<width;j++){
				if (table[i][j]==0)
					table[i][j]= 0.5;
				sumtablerow[i]=sumtablerow[i]+ table[i][j];
				sumtablecolum[j]=sumtablecolum[j]+ table[i][j];
			}	
		for(int i=0;i<length;i++)
			sum = sum+sumtablerow[i];
		//		  		outputArray(sumtablerow);
		//		  		outputArray(sumfrequencyrow);

		//		  		outputArray(sumtablecolum);
		//		  	 	outputArray(sumfrequencycolum); 			

		double chi2value=0;	
		for (int i=0;i<length;i++){
			for (int j=0;j<width;j++){
				chi2value += 
						table[i][j]*(Math.log(table[i][j])-Math.log(sumtablerow[i])-
								Math.log(sumtablecolum[j]));
			}			
		}
		chi2value+=sum*Math.log(sum);
		return 2*chi2value;
	}

	public static double pnorm(double xori){
		//  {
		//   Copyright (C) Jurg Ott 1988-2006.  Partly programmed by Joe
		//   Terwilliger using continued fractions.  Based on several formulas in
		//   Abramowitz and Stegun, "Handbook of mathematical functions", Dover, May 1968.
		//
		//   Calculates the upper tail probability of the normal
		//   distribution for a given normal deviate, x.
		//
		//		      4 Jan 2001  Last program version
		//   16 Aug 2002  Turned into function to replace old pnorm
		//  }


		double  rp,x,q,logQ,ln10;

		ln10=Math.log(10);
		rp=0.9189385332046727417803296;
		x=Math.abs(xori);
		if (Math.abs(x-0.5)<0.00001){
			q=0.5;
		}
		else if (x>1){ 
			// { Formula 26.2.14 in Abramowitz and Stegun, 1968, page 932}
			int numit=400;  		   
			double inter;  		   
			inter=x+numit+1;
			for (int i=1; i<= numit; i++) 
				inter=(numit-i+1)/inter+x;
			logQ=-Math.log(inter)-0.5*Math.pow(x,2)-rp;
			if (x<20) q=Math.exp(logQ); else q=0;
			logQ =logQ/ln10;
		}
		else {
			//{Implements formula 26.2.12, page 932, in Abramowitz and Stegun}
			double ff,sum,count,x2,two;  		  
			two=2;
			ff=x;
			x2=Math.pow(x,2);
			count=1;
			sum = ff;
			do{
				count = count+two;
				ff = ff*(x2/count);
				sum = sum+ff;
			}while(!(ff<=1.0E-300)); //{until it underflows}
			q=0.5-sum*Math.exp(-0.5*x2-rp);
			logQ=Math.log(q)/ln10;
		}
		if (xori<0) 
			q=1-q;
		double pnorm=q;
		return pnorm;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
//		double[][] table =new double[2][2];
//		table[0][0]=10113; table[0][1]=9271;
//		table[1][0] =32996-table[0][0]; table[1][1]=32996-table[0][1];
//		
//		double che =chiSquareValueLR(table);
//		System.out.println(che);
//		double chi =36.5;
		double p =chi2pr(130.082, 1);
		System.out.println(p);
		

	}


}
