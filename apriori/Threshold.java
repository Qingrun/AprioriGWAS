package apriori;
import java.io.*;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Arrays;

public class Threshold {
	public static HashMap<String, Double>single_pval(String single_result){
//		only snp with pval<=0.01 we store them in hashmap
			HashMap<String, Double> index2pval =new HashMap<String, Double>();
			try{			
				BufferedReader br = new BufferedReader(new FileReader(single_result));
				String line =br.readLine();
				while(line!=null){
					String[] array =line.split("\t");
					double pval =Double.parseDouble(array[2]);
					if(pval<=0.01){
						index2pval.put(array[0],pval);
					}
					line =br.readLine();
				}
			}catch(Exception e){e.printStackTrace();}
			return index2pval;
		}
		public static HashMap<String, Integer>single_level(String single_result){
//			only snp with pval<=0.1 we store them in hashmap
				HashMap<String, Integer> index2level =new HashMap<String, Integer>();
				try{			
					BufferedReader br = new BufferedReader(new FileReader(single_result));
					String line =br.readLine();
					while(line!=null){
						String[] array =line.split("\t");
						double pval =Double.parseDouble(array[2]);
						int level =(int)(-Math.log10(pval));
						if(level>0){
							index2level.put(array[0], level);
						}
						line =br.readLine();
					}
				}catch(Exception e){e.printStackTrace();}
				return index2level;
		}	
		
		public static void item_threshold(String filefolder, String level_snp, int permutetime, int maxlength){
			System.out.println("the output file actually contains n row, m column, row stands for permute time, " +
					"column stands for search length (start from 2), the item chi-square value sort from big to small.");
			try{		
				double[][] bigchi =new double[maxlength-1][permutetime];
				for(int i=0; i<permutetime; i++){
					for(int l=0; l<maxlength-1; l++){				
						double thebig =0;
						String file =filefolder+"/"+level_snp+"/length"+(l+2)+"_p"+i+".item";
						BufferedReader br =new BufferedReader(new FileReader(file));
						String line =br.readLine();
						while(line!=null){
							String [] array =line.split("\t");
							double prop =Double.parseDouble(array[3]);
							if(Double.compare(prop, thebig)>0){
								thebig =prop;
							}
							line =br.readLine();
						}
						bigchi[l][i] =thebig;	
					}
				}
				for(int l=0; l<maxlength-1; l++){
					Arrays.sort(bigchi[l]);
				}
				
				String[] array =level_snp.split("_");
				String out =filefolder+"/level"+array[0]+"_item.thresh";
				BufferedWriter bw = new BufferedWriter(new FileWriter(out));
				for(int i=permutetime-1; i>=0; i--){
					for(int l=0; l<maxlength-1; l++){			
						bw.write(bigchi[l][i]+"\t");
					}
					bw.write("\n");
				}
				bw.flush();
			}catch(Exception e){e.printStackTrace();}
		}
		public static void pval_threshold(String filefolder, String level_snp, int permutetime, int maxlength){
			System.out.println("the output file actually contains n row, m column, row stands for permute time, " +
					"column stands for search length (start from 2), the pvalue value sort from small to big, the most signicant in first");
			try{		
				double[][] smallp =new double[maxlength-1][permutetime];
				for(int i=0; i<permutetime; i++){
					for(int l=0; l<maxlength-1; l++){
						double thesmall =1;
						String file =filefolder+"/"+level_snp+"/length"+(l+2)+"_p"+i+".pval";
						BufferedReader br =new BufferedReader(new FileReader(file));
						String line =br.readLine();
						while(line!=null){
							String [] array =line.split("\t");
							double pval =Double.parseDouble(array[2]);
							if(Double.compare(pval, thesmall)<0){
								thesmall =pval;
							}
							line =br.readLine();
						}
						smallp[l][i] =thesmall;		
					}	
				}			
				for(int l=0; l<maxlength-1; l++){
					Arrays.sort(smallp[l]);
				}
				String[] array =level_snp.split("_");
				String out =filefolder+"/level"+array[0]+"_pval.thresh";
				BufferedWriter bw = new BufferedWriter(new FileWriter(out));
//				for(int i=0; i<permutetime-1;  i++){
				for(int i=permutetime-1; i>=0; i--){
					for(int l=0; l<maxlength-1; l++){				
						bw.write(smallp[l][i]+"\t");
					}
					bw.write("\n");
				}bw.flush();
			}catch(Exception e){e.printStackTrace();}
		}
		public static double[][]threshtable(String threshold){	
			ArrayList<String> thresholds =new ArrayList<String>();		
			try{
				BufferedReader brt =new BufferedReader(new FileReader(threshold));
				String lt =brt.readLine(); 
				while(lt!=null){
					thresholds.add(lt);
					lt =brt.readLine(); 
				}
			}catch(Exception e){e.printStackTrace();}
//			System.out.println(thresholds.size());
			
			int length =thresholds.get(0).split("\t").length;
			double[][] re =new double[thresholds.size()][length]; 
			for(int i=0; i<thresholds.size(); i++){
				String[] array =thresholds.get(i).split("\t");
				for(int j=0; j<array.length; j++){
					re[i][j] =Double.parseDouble(array[j]);
//					System.out.println(re[i][j]);
				}
			}
			return re;
			
		}
		public static void lookuptable_lengthi(String filefolder,int permutetime, int maxlevel,int searchlength,  String itemOrpval){
	// this will generate a lookuptable for search length i;
			try{
				double[][] result =new double[maxlevel+1][permutetime-1];//due to job array submission issue, we don't have value for permute0			
				String thresholdfile =filefolder+"/"+itemOrpval+"_lookup.txt";		
				BufferedWriter bw =new BufferedWriter(new FileWriter(thresholdfile));
				for(int i=0; i<=maxlevel; i++){
					String file =filefolder+"/level"+i+"_"+itemOrpval+".thresh";
					System.out.println(file);
					File temp =new File(file);
					if(temp.exists()){
						BufferedReader br =new BufferedReader(new FileReader(file));
						int rowindex = 0;
						String line =br.readLine();
						while(line!=null){
							String[] array =line.split("\t");
//							System.out.println(line);
							result[i][rowindex]=Double.parseDouble(array[searchlength-2]);	
							line=br.readLine();
							rowindex++;
						}
					}else{
						System.out.println("file"+temp+"_doesn't exist");
					}
				}
				if(itemOrpval.equals("pval")){
					for(int i=0; i<result.length; i++){
						for(int j=result[i].length-1; j>=0; j--){
							bw.write(result[i][j]+"\t");
						}bw.write("\n");
					}bw.flush();
				}else{
					for(int i=0; i<result.length; i++){
						for(int j=0; j<result[i].length; j++){
							bw.write(result[i][j]+"\t");
						}bw.write("\n");
					}bw.flush();
				}
				
			}catch(Exception e){e.printStackTrace();}
		}
		public static void threshtable(String filefolder,int permutetime, int maxlevel,int maxlength, double percent, String itemOrpval){
			try{
				String thresholdfile =filefolder+"/"+itemOrpval+percent+"_threshold.txt";
				BufferedWriter bw =new BufferedWriter(new FileWriter(thresholdfile));
				double[][] result =new double[maxlevel+1][maxlength-1];
				int index =(int)(permutetime*percent);
				for(int i=0; i<=maxlevel; i++){
					String file =filefolder+"/level"+i+"_"+itemOrpval+".thresh";
					File temp =new File(file);
					if(temp.exists()){
						BufferedReader br =new BufferedReader(new FileReader(file));
						int rowindex = 0;
						String line =br.readLine();
						while(line!=null){
							if(rowindex ==index){	
								String[] array =line.split("\t");
								for(int c =0; c<array.length; c++){								
									result[i][c]=Double.parseDouble(line);
								}
							}
							line=br.readLine();
							rowindex++;
						}
					}else{
						System.out.println("file"+temp+"_doesn't exist");
					}
				}
				for(int i=0; i<result.length; i++){
					for(int j=0; j<result[i].length; j++){
						bw.write(result[i][j]+"\t");
					}bw.write("\n");
				}bw.flush();
				
			}catch(Exception e){e.printStackTrace();}
		}

		public static double threshold(HashMap<String, Integer>single_level, double[][] table, String snp_index, int search_length){
			double result=0; 
			String[]array =snp_index.split("_"); int level =0;
			for(int i=0; i<array.length; i++){
				if(single_level.containsKey(array[i])){
					if(level<single_level.get(array[i])){
						level =single_level.get(array[i]);
					}
				}
			}
			result =table[level][search_length-2];
			return result;
		}

}
