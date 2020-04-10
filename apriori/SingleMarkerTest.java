package apriori;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashSet;
import statisticalfunctions.*;
public class SingleMarkerTest {
	public static double snp1test(String[] pheno, int[] genotype){		
		double[][] table =Support.variant2pheno(genotype, pheno);
		
		double chi_value =Chi_Square_Test.chiSquareValueLR(table);
		double snp_pval =Chi_Square_Test.chi2pr(chi_value, (table.length-1)*(table[0].length-1));		
		return snp_pval;
	}
	public static double snp1test(String[] pheno, String[] genotype){		
		double[][] table =Support.variant2pheno(genotype, pheno);
		double chi_value =Chi_Square_Test.chiSquareValueLR(table);
		double snp_pval =Chi_Square_Test.chi2pr(chi_value, (table.length-1)*(table[0].length-1));		
		return snp_pval;
	}
	public static void singletest(SimpleHDF5 geno, String out){
		try{
			BufferedWriter bw =new BufferedWriter(new FileWriter(out));
			String[] pheno =geno.pheno_index;			
			for(int i=0; i<geno.num_sites_total; i++){
				int[] genotype =geno.reader.readIntMatrixBlockWithOffset(geno.index_fast_paths, 1, geno.sample_size, i, 0)[0];
				double pval =snp1test(pheno, genotype);
				bw.write(geno.snp_order[i]+"\t"+geno.rs_ids[i]+"\t"+pval+"\n");
				
			}bw.flush(); bw.close();
			
		}catch(Exception e){e.printStackTrace();}
	}
	public static String snp4impute(String singleResult, int level){
		
		String snp ="";
		try{
			double smallestp =1;
			BufferedReader br =new BufferedReader(new FileReader(singleResult));
			String line =br.readLine();
			while(line!=null){
				String[] array =line.split("\t");
				double p =Double.parseDouble(array[2]);
				if((int)(-Math.log10(p))==level){
					if(p<smallestp){
						smallestp=p; snp =array[0];
					}
				}
				line =br.readLine();				
			}
		}catch(Exception e){e.printStackTrace();}
		return snp;
	}
	public static void permute1snp(SimpleHDF5 geno,int focsnp,String filefolder, int permuteindex){		
		try{
			String pheno =filefolder+"/permute.pheno";
			String snpgeno =filefolder+"/permute.geno";
			String[] ppheno =Support.phenodata(pheno, permuteindex);
			int[] psnpgeno =Support.genodata(snpgeno, permuteindex);
			double focopval =snp1test(ppheno, psnpgeno);
			String oneSNP =filefolder+"/length1_p"+permuteindex+".pval";
			BufferedWriter bwp1 = new BufferedWriter(new FileWriter(oneSNP));	
			for(int i=0; i<geno.num_sites_total; i++){	
				if((i+1)!=focsnp){
					int[] genotypes =geno.reader.readIntMatrixBlockWithOffset(geno.index_fast_paths, 1, geno.sample_size, i, 0)[0];						
					double pval =snp1test(ppheno, genotypes);
					bwp1.write(geno.snp_order[i]+"\t"+geno.rs_ids[i]+"\t"+pval+"\n");
				}else{
					double pval =snp1test(ppheno, psnpgeno);
					bwp1.write(geno.snp_order[i]+"\t"+geno.rs_ids[i]+"\t"+pval+"\n");
				}
			}bwp1.flush();
		}catch(Exception e){e.printStackTrace();}
	}
	public static HashSet<String> extrasig(String input, double focopval){
		HashSet<String> result =new HashSet<String>();
		try{
			BufferedReader br = new BufferedReader(new FileReader(input));
			String line =br.readLine();
			while(line!=null){
				String[] array =line.split("\t");
				if(Double.compare(focopval, Double.parseDouble(array[array.length-1]))>0){
					result.add(array[0]);
				}
				line =br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		return result;
	}
	
	public static void genotype_test(String input, String singleresult){
			SimpleHDF5 genotype = new SimpleHDF5(input);
			try{
				String p_val =singleresult+".pval";				
				BufferedWriter bwp = new BufferedWriter(new FileWriter(p_val));	
				for(int i=0; i<genotype.num_sites_total; i++){//go through each SNP
					String[] pheno =genotype.pheno_index;
					int[] snp_geno =genotype.reader.readIntMatrixBlockWithOffset(genotype.index_fast_paths, 1, genotype.sample_size, i, 0)[0];
					double snp_pval =snp1test(pheno,snp_geno);				
					bwp.write(genotype.snp_order[i]+"\t"+genotype.rs_ids[i]+"\t"+snp_pval+"\n");
				}bwp.flush(); bwp.close();
			}catch(Exception e){e.printStackTrace();}
	}
	
	public static void genotype_test(SimpleHDF5 genotype, String singleresult){
		try{
			String p_val =singleresult;				
			BufferedWriter bwp = new BufferedWriter(new FileWriter(p_val));	
			for(int i=0; i<genotype.num_sites_total; i++){//go through each SNP
				String[] pheno =genotype.pheno_index;
				int[] snp_geno =genotype.reader.readIntMatrixBlockWithOffset(genotype.index_fast_paths, 1, genotype.sample_size, i, 0)[0];
				double snp_pval =snp1test(pheno,snp_geno);
				bwp.write(genotype.snp_order[i]+"\t"+genotype.rs_ids[i]+"\t"+snp_pval+"\n");
			}bwp.flush(); bwp.close();
		}catch(Exception e){e.printStackTrace();}
}
}
