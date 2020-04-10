package apriori;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
public class Support {
	public static double array_sum(double[] input){
		double sum =0;
		for(int i=0; i<input.length; i++){
			sum=sum+input[i];
		}
		return sum;
	}
	public static HashSet<String> snplist(String input){
		HashSet<String> extrasigsnp =new HashSet<String>();
		try{
			BufferedReader br =new BufferedReader(new FileReader(input));
			String line =br.readLine();
			while(line!=null){
				extrasigsnp.add(line);				
			}line =br.readLine();
		}catch(Exception e){e.printStackTrace();}
		return extrasigsnp;		
	}

	public static double[][] variant2pheno(String[] snpgeno, String[]pheno){	
		HashMap<String, double[]> newtype =new HashMap<String, double[]>();
		for(int i=0; i<snpgeno.length;i++){
			if(newtype.keySet().contains(snpgeno[i])){
				double[] cccount = newtype.get(snpgeno[i]);
				if(pheno[i].equals("0")){
					cccount[0]++;					
				}else{
					cccount[1]++;
				}
				newtype.put(snpgeno[i], cccount);
			}else if(!snpgeno[i].equals("0")){//don't count missing data, 0 stands for missing, 1,2,3 stands for geno
//				System.out.println("type "+snpgeno[i]);
				double[] cccount =new double[2];
				if(pheno[i].equals("0")){
					cccount[0]++;
				}else{
					cccount[1]++;
				}
				newtype.put(snpgeno[i], cccount);
			}
		}
		double[][] result = new double[newtype.size()][2];
		int rowindex=0;
		for(String key:newtype.keySet()){
			double[] count = newtype.get(key);
			result[rowindex] =count; rowindex++;
		}
		return result;
	}
	public static double[][] variant2pheno(int[] snpgeno, String[]pheno){	
//		System.out.println(snpgeno.length+"\t"+pheno.length);
		HashMap<String, double[]> newtype =new HashMap<String, double[]>();
		for(int i=0; i<snpgeno.length;i++){
			if(newtype.keySet().contains(snpgeno[i]+"")){
				double[] cccount = newtype.get(snpgeno[i]+"");
				if(pheno[i].equals("0")){
					cccount[0]++;					
				}else{
					cccount[1]++;
				}
				newtype.put(snpgeno[i]+"", cccount);
			}else if(snpgeno[i]!=0){
//				System.out.println("type "+snpgeno[i]);
				double[] cccount =new double[2];
				if(pheno[i].equals("0")){
					cccount[0]++;
				}else{
					cccount[1]++;
				}
				newtype.put(snpgeno[i]+"", cccount);
			}
		}
		double[][] result = new double[newtype.size()][2];
		int rowindex=0;
		for(String key:newtype.keySet()){
			double[] count = newtype.get(key);
			result[rowindex] =count; rowindex++;
		}
		return result;
	}
	public static double[][] file2table(String file){
		double[][] result =new double[1][1];
		try{
			BufferedReader br = new BufferedReader(new FileReader(file));
			String lbr =br.readLine(); int col=0; int row=0;
			while(lbr!=null){
				String[]a =lbr.split("\t");
				col=a.length;
				row++;
				lbr =br.readLine();
			}
			result =new double[row][col];
			br = new BufferedReader(new FileReader(file));
			lbr =br.readLine();int rowindex=0;
			while(lbr!=null){
				String[]a =lbr.split("\t");
				for(int i=0; i<a.length; i++){
					result[rowindex][i] =Double.parseDouble(a[i]);
				}
				rowindex++;
				lbr =br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		return result;
	}

	public static String[]phenodata(String input, int permuteindex){
		ArrayList<String> temp =new ArrayList<String>(); 
		try{			
			BufferedReader br =new BufferedReader(new FileReader(input));
			String line =br.readLine(); int rowindex =0;			
			while(line!=null){
				
				String[] array =line.split(",");
				if(permuteindex==rowindex){
					for(int i=0; i<array.length; i++){
						temp.add(array[i]);
					}
				}	
				rowindex++;
				line =br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		String[] result =new String[temp.size()];
		for(int i=0; i<temp.size(); i++){
			result[i] =temp.get(i);
		}
//		System.out.println("geno_"+result.length);
		return result;
	}
	public static int[]genodata(String input, int permuteindex){
		ArrayList<String> temp =new ArrayList<String>(); 
		try{			
			BufferedReader br =new BufferedReader(new FileReader(input));
			String line =br.readLine(); int rowindex =0;			
			while(line!=null){
				
				String[] array =line.split(",");
				if(permuteindex==rowindex){
					for(int i=0; i<array.length; i++){
						temp.add(array[i]);
					}
				}
				rowindex++;
				line =br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		int[] result =new int[temp.size()];
		for(int i=0; i<temp.size(); i++){
			result[i] =Integer.parseInt(temp.get(i));
		}
//		System.out.println("geno_"+result.length);
		return result;
	}
	public static double criticalValue(String input, int pvalcol, double threshold){
		double result =1;
		ArrayList<Double> pvals =new ArrayList<Double>();
		try{
			BufferedReader br =new BufferedReader(new FileReader(input));//read itemset from search_length-1
			String line =br.readLine(); 
			while(line!=null){
				String[] array =line.split("\t");
				pvals.add(Double.parseDouble(array[pvalcol]));
				line =br.readLine(); 			
			}
//			System.out.println(pvals.get(0));
			Collections.sort(pvals);
			System.out.println("smallest "+pvals.get(0));
			System.out.println("Largest "+pvals.get(pvals.size()-1));
			System.out.println("total "+pvals.size());
			int index =(int)(threshold*pvals.size());
			System.out.println("index "+index);
			result =pvals.get(index);
			System.out.println("threshold "+result);
		}catch(Exception e){e.printStackTrace();}
		return result;
	}
	public static ArrayList<String> candidate (String item, String pval,double markernum, double ithresh, double pthresh){
		ArrayList<String> result =new ArrayList<String>();
		try{
			BufferedReader bri =new BufferedReader(new FileReader(item));
			String iline =bri.readLine(); 
			while(iline!=null){
				String[] array =iline.split("\t");	
				if(!result.contains(array[1])){
					result.add(array[1]);
				}
				iline =bri.readLine();
			}
			BufferedReader brp =new BufferedReader(new FileReader(pval));
			String pline =brp.readLine(); int i=0;
			while(pline!=null){
				String[] array =pline.split("\t");
				if(!result.contains(array[1])){
					result.add(array[1]);
				}
				pline =brp.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		System.out.println(result.size());
		return result;
	}
	public static ArrayList<String> candidate (String item, String pval){
		ArrayList<String> result =new ArrayList<String>();
		try{
			BufferedReader bri =new BufferedReader(new FileReader(item));
			String iline =bri.readLine(); 
			while(iline!=null){
				String[] array =iline.split("\t");	
				if(!result.contains(array[1])){
					result.add(array[1]);
				}
				iline =bri.readLine();
			}
			BufferedReader brp =new BufferedReader(new FileReader(pval));
			String pline =brp.readLine(); int i=0;
			while(pline!=null){
				String[] array =pline.split("\t");
				if(!result.contains(array[1])){
					result.add(array[1]);
				}
				pline =brp.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		System.out.println(result.size());
		return result;
	}
	public static ArrayList<String> candidate (String item, String pval, double ithresh, double pthresh){
		ArrayList<String> result =new ArrayList<String>();
		try{
			BufferedReader bri =new BufferedReader(new FileReader(item));
			String iline =bri.readLine(); 
			while(iline!=null){
				String[] array =iline.split("\t");				
				if(Double.compare(ithresh, Double.parseDouble(array[3]))<0){
					if(!result.contains(array[1])){
						result.add(array[1]);
					}
				}
				iline =bri.readLine();
			}
			BufferedReader brp =new BufferedReader(new FileReader(pval));
			String pline =brp.readLine(); int i=0;
			while(pline!=null){
				String[] array =pline.split("\t");
				if(Double.compare(pthresh, Double.parseDouble(array[2]))>0){					
					if(!result.contains(array[1])){
						result.add(array[1]);
					}
				}
				pline =brp.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		System.out.println(result.size());
		return result;
	}
}
