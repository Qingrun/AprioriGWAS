package apriori;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import statisticalfunctions.*;
public class Permutation {
//	prepare permuted phenotype and focosnp genotype file, each line stands for one permutation.
	public static String[] siteinfo (String sitecode){
		String[] result =new String[1];
		try{
			BufferedReader br =new BufferedReader(new FileReader(sitecode));
			String lbr =br.readLine();  int samplesize =0;
			while(lbr!=null){
				samplesize++;
				lbr =br.readLine();
			}
			result =new String[samplesize];
			int index =0;
			br =new BufferedReader(new FileReader(sitecode));
			lbr =br.readLine(); 
			while(lbr!=null){
				String[] a =lbr.split("\t");
				if(!a[0].equals("FID")){
					result[index] =a[1];
//					System.out.println(a[1]);
					index++;
				}
				lbr =br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		return result;
	}
	public static void sitepermute(SimpleHDF5 geno, int focosnp, String output, int permute_times, String siteinfo){
		try{
			String[] sitecode =siteinfo(siteinfo);
			int[] focsnpgeno =geno.reader.readIntMatrixBlockWithOffset(geno.index_fast_paths, 1, geno.sample_size, (focosnp-1), 0)[0];
			String[] phenos =geno.pheno_index.clone();
			String pheno =output+"/permute.pheno";
			String focogeno =output+"/permute.geno";
			BufferedWriter bwphe = new BufferedWriter(new FileWriter(pheno));
			BufferedWriter bwgeno = new BufferedWriter(new FileWriter(focogeno));
			for(int i=0; i<permute_times; i++){
				statisticalfunctions.Permutation.condiPermuteSite(phenos, sitecode,focsnpgeno);
				for(int s =0; s<phenos.length-1; s++){
					bwphe.write(phenos[s]+",");
					bwgeno.write(focsnpgeno[s]+",");
				}
				bwphe.write(phenos[phenos.length-1]);
				bwgeno.write(focsnpgeno[focsnpgeno.length-1]+"");
				bwphe.write("\n"); bwgeno.write("\n");
			}bwphe.flush(); bwgeno.flush();			
		}catch(Exception e){e.printStackTrace();}
	}
	public static void permuteData(SimpleHDF5 geno, int focosnp, String output, int permute_times){
		try{			
			int[] focsnpgeno =geno.reader.readIntMatrixBlockWithOffset(geno.index_fast_paths, 1, geno.sample_size, (focosnp-1), 0)[0];
			String[] phenos =geno.pheno_index.clone();
			String pheno =output+"/permute.pheno";
			String focogeno =output+"/permute.geno";
			BufferedWriter bwphe = new BufferedWriter(new FileWriter(pheno));
			BufferedWriter bwgeno = new BufferedWriter(new FileWriter(focogeno));
			for(int i=0; i<permute_times; i++){
				statisticalfunctions.Permutation.randPermute(phenos, focsnpgeno);
				for(int s =0; s<phenos.length-1; s++){
					bwphe.write(phenos[s]+",");
					bwgeno.write(focsnpgeno[s]+",");
				}
				bwphe.write(phenos[phenos.length-1]);
				bwgeno.write(focsnpgeno[focsnpgeno.length-1]+"");
				bwphe.write("\n"); bwgeno.write("\n");
			}bwphe.flush(); bwgeno.flush();			
		}catch(Exception e){e.printStackTrace();}
	}

	public static void summarypermute(String filefolder, int permutetime,int searchLength, double itemthresh, double pvalthresh){
		try{
			String itemout =filefolder+"/length"+searchLength+"item.txt";
			BufferedWriter bwitem = new BufferedWriter(new FileWriter(itemout));
			String pvalout =filefolder+"/length"+searchLength+"pval.txt";
			BufferedWriter bwpval = new BufferedWriter(new FileWriter(pvalout));
			for(int i=0; i<permutetime; i++){
				String item =filefolder+"/length"+searchLength+"_p"+i+".item";
				File fitem =new File(item);
				if(fitem.exists()){
					BufferedReader britem =new BufferedReader(new FileReader(fitem));
					String litem =britem.readLine();
					while(litem!=null){
						String[] array =litem.split("\t");
						double chi =Double.parseDouble(array[3]);
//						double count1=Double.parseDouble(array[4]);
//						double count2 =Double.parseDouble(array[5]);
						if(Double.compare(chi, itemthresh)>0){
							bwitem.write(litem+"\n");
						}
						litem =britem.readLine();
					}
				}
				String pval =filefolder+"/length"+searchLength+"_p"+i+".pval";
				File fpval =new File(pval);
				if(fpval.exists()){
					BufferedReader brpval =new BufferedReader(new FileReader(fpval));
					String lpval =brpval.readLine();
					while(lpval!=null){
						String[] array =lpval.split("\t");
						double p =Double.parseDouble(array[2]);
						if(Double.compare(p, pvalthresh)>0){
							bwitem.write(lpval+"\n");
						}
						lpval =brpval.readLine();
					}
				}
			}
		}catch(Exception e){e.printStackTrace();}
	}
	public static void permute(SimpleHDF5 geno, int focsnp, int permuteindex, String filefolder, double itemthresh, double pvalthresh, int searchLength){	
		System.out.println("for the first permute data, we set a loose itemthresh and pvalthresh, then for other permute data, we will use the threshold we get+" +
				"from the first permute data");
		try{
			String snpgeno =filefolder+"/permute.geno";
			String pheno =filefolder+"/permute.pheno";
			String singletest =filefolder+"/length1_p"+permuteindex+".pval";
			int[] psnpgeno =Support.genodata(snpgeno, permuteindex);
			String[] ppheno =Support.phenodata(pheno, permuteindex);
			double focopval =SingleMarkerTest.snp1test(ppheno, psnpgeno);
			HashSet<String> exsig =SingleMarkerTest.extrasig(singletest, focopval);
			String p_val =filefolder+"/length"+searchLength+"_p"+permuteindex+".pval";
			String item =filefolder+"/length"+searchLength+"_p"+permuteindex+".item";
			BufferedWriter bwp = new BufferedWriter(new FileWriter(p_val));	
			BufferedWriter bwt = new BufferedWriter(new FileWriter(item));	
			if(searchLength==2){			
				for(int i=0; i<geno.num_sites_total; i++){
					if(focsnp!=(i+1) && !exsig.contains((i+1)+"")){
						String snp_index =(i+1)+"";
						Pattern_support pats =new Pattern_support(geno, snp_index, focsnp, psnpgeno,ppheno);
						HashMap<String,double[]> patterns =pats.pat_support;		
						double[][] pattern_counts =new double[patterns.keySet().size()][2]; int row_index =0;
						for(String key:patterns.keySet()){//check each pattern whether they passed proportion test, if past, then write to item file
							double[] count = patterns.get(key);
							pattern_counts[row_index] =count; row_index++;
							double proptest = statisticalfunctions.Proportion_test.Proportiontest(count[0], pats.controlcount, count[1], pats.casecount);
							if(Double.compare(proptest, itemthresh)>0){
								bwt.write(pats.rs+"\t"+pats.snp_seq+"\t"+key+"\t"+proptest+"\t"+count[0]+"\t"+count[1]+"\n");
							}
						}
						double[][] merg =statisticalfunctions.Proportion_test.merged(pattern_counts);
						double chi2 =statisticalfunctions.Chi_Square_Test.chiSquareValueLR(merg);
						double pval =statisticalfunctions.Chi_Square_Test.chi2pr(chi2, (merg[0].length-1)*(merg.length-1));
						if(Double.compare(pval, pvalthresh)<0){
							bwp.write(pats.rs+"\t"+pats.snp_seq+"\t"+pval+"\n");
						}
					}
					bwp.flush(); bwt.flush();
				}
			}else{
				ArrayList<String> candidates=new ArrayList<String>();
				String itembefore =filefolder+"/length"+(searchLength-1)+"_p"+permuteindex+".item";
				String pvalbefore = filefolder+"/length"+(searchLength-1)+"_p"+permuteindex+".pval";
				if(permuteindex==0){
					double criticitem =Support.criticalValue(itembefore, 3,0.95);		
					double criticpval =Support.criticalValue(pvalbefore, 2, 0.05);
					candidates =Support.candidate(itembefore, pvalbefore, criticitem, criticpval);
				}else{
					candidates =Support.candidate(itembefore, pvalbefore);
				}
				System.out.println("candidates number "+candidates.size());
				for(int k=0; k<candidates.size(); k++){			
					String[] pairs =candidates.get(k).split("_");
					String snp ="";
					for(int m=1; m<pairs.length; m++){
						snp=snp+pairs[m]+"_";
					}					
					int lastsnp =Integer.parseInt(pairs[pairs.length-1]);
					for(int i =lastsnp; i<geno.num_sites_total; i++){
						String snpindex =snp+(i+1);
						Pattern_support pats =new Pattern_support(geno, snpindex, focsnp, psnpgeno,ppheno);
						HashMap<String,double[]> patterns =pats.pat_support;		
						double[][] pattern_counts =new double[patterns.keySet().size()][2]; int row_index =0;
						for(String key:patterns.keySet()){//check each pattern whether they passed proportion test, if past, then write to item file
							double[] count = patterns.get(key);
							pattern_counts[row_index] =count; row_index++;
							double proptest = statisticalfunctions.Proportion_test.Proportiontest(count[0], pats.controlcount, count[1], pats.casecount);
							if(Double.compare(proptest, itemthresh)>0){
								bwt.write(pats.rs+"\t"+pats.snp_seq+"\t"+key+"\t"+proptest+"\t"+count[0]+"\t"+count[1]+"\n");
							}
						}
						double[][] merg =statisticalfunctions.Proportion_test.merged(pattern_counts);
						double chi2 =statisticalfunctions.Chi_Square_Test.chiSquareValueLR(merg);
						double pval =statisticalfunctions.Chi_Square_Test.chi2pr(chi2, (merg[0].length-1)*(merg.length-1));
						if(Double.compare(pval, pvalthresh)<0){
							bwp.write(pats.rs+"\t"+pats.snp_seq+"\t"+pval+"\n");
						}
					}
				}
				
				
			}
		}catch(Exception e){e.printStackTrace();}
		
		
		
	}
}
