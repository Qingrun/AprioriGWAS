package statisticalfunctions;
import java.util.ArrayList;
public class Permutation {
	public static void randomPermute(double[] indicators){
		java.util.Random r= new java.util.Random();
		int mm = indicators.length;
		for (int ii=0; ii< mm; ii++){
		  // {i1 will be a random element from ii .. mm}
			double randomNumber = r.nextDouble();
			double yy = randomNumber*(mm-1-ii);
			int zz = (int)yy;
			int xx = (yy-zz>=0.5)?(zz+1):(zz);
		    int i1 = ii + (xx);
		   //{Now exchange element ii with element i1 in the array}
		    double temp = indicators[i1];
		    indicators[i1] = indicators[ii];
		    indicators[ii] = temp;
		}
	}
	
	public static void randPermute(char[] indicators){
		int m = indicators.length;
		for (int i=0; i< m; i++){
		  // {i1 will be a random element from ii .. mm}
			double randomNumber = Randomenumber.randomNumber();
			double yy = randomNumber*(m-1-i);
			int zz = (int)yy;
			int xx = (yy-zz>=0.5)?(zz+1):(zz);
		    int i1 = i + (xx);
		   //{Now exchange element ii with element i1 in the array}
		    char temp = indicators[i1];
		    indicators[i1] = indicators[i];
		    indicators[i] = temp;
		}
	}
	public static void randPermute(String[] indicators){
		int m = indicators.length;
		for (int i=0; i< m; i++){
		  // {i1 will be a random element from ii .. mm}
			double randomNumber = Randomenumber.randomNumber();
			double yy = randomNumber*(m-1-i);
			int zz = (int)yy;
			int xx = (yy-zz>=0.5)?(zz+1):(zz);
		    int i1 = i + (xx);
		   //{Now exchange element ii with element i1 in the array}
		    String temp = indicators[i1];
		    indicators[i1] = indicators[i];
		    indicators[i] = temp;
		}
	}
	public static void randPermute(char[] indicators, char[] snp){
		if(indicators.length!=snp.length){
			System.out.println("length for permutation is wrong ");
		}
		int m = indicators.length;
		for (int i=0; i< m; i++){
		  // {i1 will be a random element from ii .. mm}
			double randomNumber = Randomenumber.randomNumber();
			double yy = randomNumber*(m-1-i);
			int zz = (int)yy;
			int xx = (yy-zz>=0.5)?(zz+1):(zz);
		    int i1 = i + (xx);
		   //{Now exchange element ii with element i1 in the array}
		    char temp = indicators[i1];		    
		    indicators[i1] = indicators[i];
		    indicators[i] = temp;		    
		    char tempsnp =snp[i1];
		    snp[i1]=snp[i];
		    snp[i] =tempsnp;
		}
	}

	public static void randPermute(String[] indicators, char[] snp){
		int m = indicators.length;
		for (int i=0; i< m; i++){
		  // {i1 will be a random element from ii .. mm}
			double randomNumber = Randomenumber.randomNumber();
			double yy = randomNumber*(m-1-i);
			int zz = (int)yy;
			int xx = (yy-zz>=0.5)?(zz+1):(zz);
		    int i1 = i + (xx);
		   //{Now exchange element ii with element i1 in the array}
		    String temp = indicators[i1];		    
		    indicators[i1] = indicators[i];
		    indicators[i] = temp;		    
		    char tempsnp =snp[i1];
		    snp[i1]=snp[i];
		    snp[i] =tempsnp;
		}
	}
	
	public static void randPermute(String[] indicators, int[] snp){
		int m = indicators.length;
		for (int i=0; i< m; i++){
		  // {i1 will be a random element from ii .. mm}
			double randomNumber = Randomenumber.randomNumber();
			double yy = randomNumber*(m-1-i);
			int zz = (int)yy;
			int xx = (yy-zz>=0.5)?(zz+1):(zz);
		    int i1 = i + (xx);
		   //{Now exchange element ii with element i1 in the array}
		    String temp = indicators[i1];		    
		    indicators[i1] = indicators[i];
		    indicators[i] = temp;		    
		    int tempsnp =snp[i1];
		    snp[i1]=snp[i];
		    snp[i] =tempsnp;
		}
	}
	public static void condiPermuteSite(String[] indicators, String[]siteinfo, int[] genotype){
//		double pval1 =1; 	
//		System.out.println(siteinfo.length);
		int m =indicators.length;
		ArrayList<String> sites =new ArrayList<String>();
		ArrayList<Integer> genos =new ArrayList<Integer>();
		for(int i=0; i<m; i++){
			if(!sites.contains(siteinfo[i])){
				sites.add(siteinfo[i]);
//				System.out.println(siteinfo[i]);
			}
			if(!genos.contains(genotype[i])){
				genos.add(genotype[i]);
			}
		}	
//		System.out.println(sites.size()+"\t"+genos.size());
		for(int i=0; i<sites.size(); i++){
			String siteid=sites.get(i);
			for(int j=0; j<genos.size(); j++){
				int geno =genos.get(j);
				ArrayList<Integer> toswap =new ArrayList<Integer>();
				for(int k=0; k<m; k++){
					if(siteinfo[k].equals(siteid) && genotype[k]==geno){
						toswap.add(k);
					}
				}
				int num =toswap.size();
				for(int k=0; k<toswap.size(); k++){
					double randomNumber = Randomenumber.randomNumber();	
					double yy = randomNumber*(num-1-k);
					int zz = (int)yy;
					int xx = (yy-zz>=0.5)?(zz+1):(zz);
				    int i1 = k + (xx);
				   //{Now exchange element ii with element i1 in the array}
				    int num1=toswap.get(i1);
				    int num2=toswap.get(k);
				    String temp = indicators[num1];		    
				    indicators[num1] = indicators[num2];
				    indicators[num2] = temp;		    
				    int tempsnp =genotype[num1];
				    genotype[num1]=genotype[num2];
				    genotype[num2] =tempsnp;
				}
				
			}
			
		}
		
	}
}
