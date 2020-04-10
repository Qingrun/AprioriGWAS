package apriori;
import java.io.*;
import java.util.HashMap;
import ch.systemsx.cisd.hdf5.HDF5MDDataBlock;
import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.base.mdarray.MDIntArray;

import statisticalfunctions.*;


public class AprioriGWAS {
	public static void algorithm(SimpleHDF5 genotype, int search_length,  String result_folder, int step_size, int step_index, String pthresh,String singleresult){
//		this algorithm search patterns step by step, first search length=1, then 2, 3, for each step, we do permutation test first to get pro_threshold for 
//		different searching length and different snps's significant level.so pro_threshold is a table each row stands for the most significant level. 
			System.out.println("Search length\t"+search_length);
			long startTime = System.currentTimeMillis();	
			try{
				HashMap<String,Integer> snp2level =Threshold.single_level(singleresult);
//				double[][] itreshtable =Threshold.threshtable(ithresh);	
				double[][] pthreshtable =Threshold.threshtable(pthresh);
				String p_val =result_folder+"/length_"+search_length+"_step_"+step_index+".pval";
//				System.out.println(p_val);
//				String item =result_folder+"/length_"+search_length+"_step_"+step_index+".item";			
//				bwt.write("rs\tsnpseq\tpattern\tcontrol\tcase\tproptest\tthreshold\n");
				if(search_length==2){	
					BufferedWriter bwp = new BufferedWriter(new FileWriter(p_val));	
//					BufferedWriter bwt = new BufferedWriter(new FileWriter(item));	
					for(int i=step_size*step_index; i<step_size*(step_index+1); i++){
						for(int j =i+1; j<genotype.num_sites_total; j++){//get snp pairs
							String snp_indexs =(i+1)+"_"+(j+1);//System.out.println("snp_index "+snp_indexs);
//							double newithresh =Threshold.threshold(snp2level, itreshtable, snp_indexs, search_length);
							double newpthresh =Threshold.threshold(snp2level, pthreshtable, snp_indexs, search_length);
							Pattern_support pats =new Pattern_support(genotype, snp_indexs);//get pattern support for involved snps
							HashMap<String,double[]> patterns =pats.pat_support;		//System.out.println(patterns.size());
//							boolean meet =false;//to check whether pairs of SNP has genotype patterns which passed proportion test,for further whole pattern test.
							double[][] pattern_counts =new double[patterns.keySet().size()][2]; 
							int row_index =0;
							for(String key:patterns.keySet()){//check each pattern whether they passed proportion test, if past, then write to item file
								double[] count = patterns.get(key);
								pattern_counts[row_index] =count; row_index++;
//								double proptest = statiticalfunctions.Proportion_test.Proportiontest(count[0], pats.controlcount, count[1], pats.casecount);
//								if(Double.compare(proptest, newithresh)>0 ){
//									meet =true;
//									bwt.write(pats.rs+"\t"+pats.snp_seq+"\t"+key+"\t"+count[0]+"\t"+count[0]/pats.controlcount+"\t"+count[1]+"\t"+count[1]/pats.casecount+"\t"+proptest+"\t"+newithresh+"\n");
//								}
							}						
//							if(meet==true){
								double[][] merg =statisticalfunctions.Proportion_test.merged(pattern_counts);
								double chi2 =statisticalfunctions.Chi_Square_Test.chiSquareValueLR(merg);
								double pval =statisticalfunctions.Chi_Square_Test.chi2pr(chi2, (merg[0].length-1)*(merg.length-1));
								if(Double.compare(pval, newpthresh)<0){
									bwp.write(pats.rs+"\t"+pats.snp_seq+"\t"+pval+"\t"+newpthresh+"\n");
								}
//							}
						}
						bwp.flush(); 
//						bwt.flush();					
					}System.out.println("new_apriori");
				}else{
//					BufferedWriter bwp = new BufferedWriter(new FileWriter(p_val));	
//					BufferedWriter bwt = new BufferedWriter(new FileWriter(item));	
					System.out.println("To be updated");
					
				}
			}catch(Exception e){e.printStackTrace();}
			long endTime = System.currentTimeMillis();	
			System.out.println("Time consume\t"+(endTime-startTime)/1000+"seconds");
		}
		


}
