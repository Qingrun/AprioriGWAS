package apriori;
import java.io.*;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Arrays;

public class Pattern_support {
	public HashMap<String, double[]> pat_support =new HashMap<String, double[]>();//pattern and support in case and control. [0] for control, [1]for case
//	public double[][] pattern_counts;
	public String rs;
	public String snp_seq;
	public double casecount=0;
	public double controlcount=0;
	
	public Pattern_support(SimpleHDF5 genotype, String snp_index){
		this.snp_seq =snp_index; //snp index if comes from i, then keep it, if come from input number, then should -1
		String[] snp_indexes =snp_index.split("_");
		int[] snps =new int[snp_index.split("_").length];
		for(int i=0; i<snps.length; i++){
			snps[i] =Integer.parseInt(snp_indexes[i]);
		}
		this.rs =genotype.rs_ids[snps[0]-1]; 
		for(int i=1; i<snps.length; i++){
			this.rs=this.rs+"_"+genotype.rs_ids[snps[i]-1];
		}
		int[][] foco_snps =new int[snp_indexes.length][genotype.sample_size];
		int index4snps =0;
		for(int i=0; i<snps.length; i++){		
			foco_snps[index4snps] =genotype.reader.readIntMatrixBlockWithOffset(genotype.index_fast_paths, 1, genotype.sample_size, snps[i]-1, 0)[0];			
			index4snps++;
		}
		for(int sample=0; sample<genotype.sample_size; sample++){
			String pattern =""; boolean missing=false;
			for(int variant =0; variant<foco_snps.length;variant++ ){
				if(foco_snps[variant][sample]==0){
					missing =true;
				}else{
					pattern=pattern+foco_snps[variant][sample];
				}
			}
			if(missing==false){
				if(this.pat_support.keySet().contains(pattern)){
					double[] case_control_support = this.pat_support.get(pattern);
					if(genotype.pheno_index[sample].equals("0")){
						case_control_support[0]++;
						this.controlcount++;
					}else{
						case_control_support[1]++;
						this.casecount++;
					}
					this.pat_support.put(pattern, case_control_support);
				}else{
					double[] case_control_support =new double[2];
					if(genotype.pheno_index[sample].equals("0")){
						case_control_support[0]++;
						this.controlcount++;
					}else{
						case_control_support[1]++;
						this.casecount++;
					}
					this.pat_support.put(pattern, case_control_support);
				}
			}
		}
	}
	
	public Pattern_support(SimpleHDF5 genotype, String snp_index, int focsnp, int[]focsnpgeno, String[] newpheno){
		this.snp_seq =focsnp+"_"+snp_index;//add focusnp
		String[] snp_indexes =snp_index.split("_");
		int[] snps =new int[snp_index.split("_").length];		
		for(int i=0; i<snps.length; i++){
			snps[i] =Integer.parseInt(snp_indexes[i]);
		}
		this.rs =genotype.rs_ids[(focsnp-1)]; // add focusnp
		for(int i=0; i<snps.length; i++){
			this.rs=this.rs+"_"+genotype.rs_ids[snps[i]-1];
		}
		int[][] foco_snps =new int[snps.length+1][genotype.sample_size];
		foco_snps[0] =focsnpgeno;//add focusnp's genotype
		int index4snps =1;		
		for(int i=0; i<snps.length; i++){
			foco_snps[index4snps] =genotype.reader.readIntMatrixBlockWithOffset(genotype.index_fast_paths, 1, genotype.sample_size, (snps[i]-1), 0)[0];			
			index4snps++;
		}
		for(int sample=0; sample<genotype.sample_size; sample++){
			String pattern =""; boolean missing=false;
			for(int variant =0; variant<foco_snps.length;variant++ ){
				if(foco_snps[variant][sample]==0){
					missing =true;
				}else{
					pattern=pattern+foco_snps[variant][sample];
				}
			}
			if(missing==false){
				if(this.pat_support.keySet().contains(pattern)){
					double[] case_control_support = this.pat_support.get(pattern);
					if(newpheno[sample].equals("0")){
						case_control_support[0]++;
						this.controlcount++;
					}else{
						case_control_support[1]++;
						this.casecount++;
					}
					this.pat_support.put(pattern, case_control_support);
				}else{
					double[] case_control_support =new double[2];
					if(newpheno[sample].equals("0")){
						case_control_support[0]++;
						this.controlcount++;
					}else{
						case_control_support[1]++;
						this.casecount++;
					}
					this.pat_support.put(pattern, case_control_support);
				}
			}
		}
	}
	public Pattern_support(SimpleHDF5 genotype, int[] snps){		
		this.rs =genotype.rs_ids[snps[0]-1]; this.snp_seq =snps[0]+"";
		for(int i=1; i<snps.length; i++){
			this.rs=this.rs+"_"+genotype.rs_ids[snps[i]-1];
			this.snp_seq =this.snp_seq+"_"+snps[i];
		}		
		int[][] foco_snps =new int[snps.length][genotype.sample_size];
		int index4snps =0;
		for(int i=0; i<snps.length; i++){
			foco_snps[index4snps] =genotype.reader.readIntMatrixBlockWithOffset(genotype.index_fast_paths, 1, genotype.sample_size, snps[i], 0)[0];
			index4snps++;
		}
		for(int sample=0; sample<genotype.sample_size; sample++){
			String pattern =""; boolean missing=false;
			for(int variant =0; variant<foco_snps.length;variant++ ){
				if(foco_snps[variant][sample]==0){
					missing =true;
				}else{
					pattern=pattern+foco_snps[variant][sample];
				}
			}
			if(missing==false){
				if(this.pat_support.keySet().contains(pattern)){
					double[] case_control_support = this.pat_support.get(pattern);
					if(genotype.pheno_index[sample].equals("0")){
						case_control_support[0]++;
						this.controlcount++;
					}else{
						case_control_support[1]++;
						this.casecount++;
					}
					this.pat_support.put(pattern, case_control_support);
				}else{
					double[] case_control_support =new double[2];
					if(genotype.pheno_index[sample].equals("0")){
						case_control_support[0]++;
						this.controlcount++;
					}else{
						case_control_support[1]++;
						this.casecount++;
					}
					this.pat_support.put(pattern, case_control_support);
				}
			}
		}
	}

}
