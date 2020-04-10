package apriori;

import java.io.IOException;

public class Start {
	static String mit="Copyright (C) 2015, QINGRUN ZHANG and QUAN LONG" +
			"\n\nPermission is hereby granted, free of charge, to any person obtaining a copy " +
			"\nof this software and associated documentation files (the \"Software\"), " +
			"\nto deal in the Software without restriction, including without limitation the " +
			"\nrights to use, copy, modify, merge, publish, distribute, sublicense, " +
			"\nand/or sell copies of the Software, and to permit persons to whom the " +
			"\nSoftware is furnished to do so, subject to the following conditions:" +
			"\n\nThe above copyright notice and this permission notice shall be included " +
			"\nin all copies or substantial portions of the Software." +
			"\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS " +
			"\nOR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, " +
			"\nFITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE " +
			"\nAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER " +
			"\nLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING " +
			"\nFROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS " +
			"\nIN THE SOFTWARE.";


	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		// TODO Auto-generated method stub
		if(args.length==0){
			System.out.println("==========================================\n"+mit+"\n==========================================\n");
			System.out.println("AprioriGWAS: Java implementation of gene-gene interaction from case control data\n" +
					"Developer: QINGRUN ZHANG & QUAN LONG\n" +
					"Usage: java -Xmx2g -jar AprioriGWAS.jar [function]");
			System.out.println("Supported functions:" +
					"\n\ttped2csv"+
					"\n\timporthdf5" +
					"\n\tSingleMarkerTest" +
					"\n\tSNP4Permutaton" +
					"\n\tPreparePermute"+
					"\n\tPermutation" +
					"\n\tThreshold" +
					"\n\tThresholdtable"+
					"\n\tApriori"+
					"\n\tGetpattern");
			System.exit(0);
		}
		String function=args[0];
		if(function.equals("tped2csv")){
			if(args.length<8){
				System.out.println("transfer tepd file to aprioriGWAS csv file");
				System.out.println("Usage: \n\t[-p\ttped_file]\n\t" +
						"[-f\ttfam file]\n\t" +
						"[-c\tcsv file]\n\t" +
						"[-a\tfrquency file]\n\t");
				System.exit(0);
			}else{		
				//			String prefix ="/Users/qingrunzhang/Documents/Working/AprioriGWAS/iccbd/";
				String tped =null;
				String tfam=null;
				String frq =null;
				String csv =null;

				for(int k=0;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-p")) tped=args[k+1];
						else if(args[k].equals("-f")) tfam=args[k+1];
						else if(args[k].equals("-a")) frq=args[k+1];
						else if(args[k].equals("-c")) csv=args[k+1];
					}
				}
				
				Tped2csv.tped2csv(tped, tfam, frq, csv);
			}
		}
		if(function.equals("importhdf5")){
			if(args.length==1){
				System.out.println("Import genotype from .csv file to .hdf5 file");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-o\tout_put_hdf5_file>\n\t" +
						"[-b\tblock_size (df=2000)]\n\t");
				System.exit(0);
			}else{
				String file_csv=null, file_hdf5=null;
				int block_size=2000;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))file_csv=args[k+1];
						else if(args[k].equals("-o"))file_hdf5=args[k+1];
						else if(args[k].equals("-b"))block_size=Integer.parseInt(args[k+1]);
					}
				}if(file_csv==null||file_hdf5==null){
					System.out.println("file_csv or file_hdf5 can't be null!");
				}else{
					SimpleHDF5.importCSV(file_csv, file_hdf5, block_size);
				}
			}
		}else if(function.equals("SingleMarkerTest")){
			if(args.length==1){
				System.out.println("Do single marker test on original data");
				System.out.println("Usage: \n\t<-g\thdf_file>\n\t" +
						"<-o\tout_put>\n\t");
				System.exit(0);
			}else{
				String  file_hdf5=null; String  result =null;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-g"))file_hdf5=args[k+1];
						else if(args[k].equals("-o"))result=args[k+1];
					}
				}if(result==null||file_hdf5==null){
					System.out.println("result or file_hdf5 can't be null!");
				}else{
					long startTime = System.currentTimeMillis();	
					SimpleHDF5 geno =new SimpleHDF5(file_hdf5);
					SingleMarkerTest.singletest(geno, result);
					long endTime = System.currentTimeMillis();	
					System.out.println("Time consume\t"+(endTime-startTime)/1000+"seconds");
				}
			}
		}else if(function.equals("SNP4Permutaton")){
			if(args.length==1){
				System.out.println("chose snp of each level for permutation");
				System.out.println("Usage: \n\t<-s\tonesnpresult>\n\t" +
						"<-l\tMax_level4impute>\n\t");
				System.exit(0);
			}else{
				String  onesnpresult=null; int level=0;;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-s"))onesnpresult=args[k+1];
						else if(args[k].equals("-l"))level=Integer.parseInt(args[k+1]);
					}
				}if(onesnpresult==null || level==0){
					System.out.println("onesnpresult can't be null! level should bigger than 0");
				}else{
					java.lang.Runtime rt = java.lang.Runtime.getRuntime();
					java.lang.Process p =rt.exec("mkdir permutation");// rt.exec("mkdir Permutation");
					for(int le=0; le<level; le++){
						String snpid =SingleMarkerTest.snp4impute(onesnpresult, le);
						java.lang.Process g = rt.exec("mkdir Permutation/"+le+"_"+snpid);
						System.out.println(le+"\t"+snpid);
					}
				}
			}
		}else if(function.equals("PreparePermute")){
			if(args.length==1){
				System.out.println("Prepare conditional Permutation for selected SNP");
				System.out.println("Usage: \n\t<-g\thdf_file>\n\t" +"<-s\tsiteinfo\n\t>"+
						"<-f\tfocosnp>\n\t"+"<-o\toutputprefix>\n\t"+"<-n\tpermute_times>\n\t");
				System.exit(0);
			}else{
				String  file_hdf5=null; int focosnp=0; String outprefix =null; int permute_times=0;
				String siteinfo =null;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-g"))file_hdf5=args[k+1];
						else if(args[k].equals("-f"))focosnp=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-o"))outprefix=args[k+1];
						else if(args[k].equals("-n"))permute_times=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-s"))siteinfo=args[k+1];
					}
				}if(file_hdf5==null || focosnp==0 ||outprefix==null|| permute_times==0){
					System.out.println("file_hdf5 or pheno or focogeno can not be null or permute_times can't be 0!");
				}else{
					SimpleHDF5 geno =new SimpleHDF5(file_hdf5);
//					Permutation.permuteData(geno, focosnp, outprefix, permute_times);
					Permutation.sitepermute(geno, focosnp, outprefix, permute_times, siteinfo);
				}
			}

		}else if(function.equals("Permutation")){
			if(args.length==1){
				System.out.println("Permutation for selected SNP");
				System.out.println("Usage: \n\t<-g\thdf_file>\n\t" +
						"<-f\tfocosnp>\n\t"+
						"<-p\tpermuteindex>\n\t"+
						"<-i\tfilefolder>\n\t"+
						"<-it\titemthresh>\n\t"+
						"<-pt\tpvalthresh>\n\t"+
						"<-s\tsearchlength>\n\t"						);
				System.exit(0);
			}else{
				String  file_hdf5=null; int focsnp=0; String filefolder=null; int permuteindex =0; int searchlength=0; 
				double itemthresh=0; double pvalthresh=1;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-g"))file_hdf5=args[k+1];
						else if(args[k].equals("-f"))focsnp=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-p"))permuteindex=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-i"))filefolder=args[k+1];
						else if(args[k].equals("-it"))itemthresh=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-pt"))pvalthresh=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-s"))searchlength=Integer.parseInt(args[k+1]);
					}
				}if(file_hdf5==null|| filefolder==null ||searchlength==0||focsnp==0){
					System.out.println("outprefix or file_hdf5 or inprefix can't be null! searchlength focosnp permuteindex can't be 0");
				}else{
					long startTime = System.currentTimeMillis();	
					SimpleHDF5 geno =new SimpleHDF5(file_hdf5);
					
					SingleMarkerTest.permute1snp(geno, focsnp, filefolder, permuteindex);
					Permutation.permute(geno, focsnp, permuteindex, filefolder, itemthresh, pvalthresh, searchlength);
//					Permutation.permute(geno, focsnp, permuteindex, filefolder, searchlength);
					long endTime = System.currentTimeMillis();	
					System.out.println("Time consume\t"+(endTime-startTime)/1000+"seconds");
				}
			}
		}else if(function.equals("Threshold")){
			if(args.length==1){
				System.out.println("Get threshold from permutation data");
				System.out.println("Usage: \n\t<-f\tfilefolder>\n\t" +
						"<-s\tlevel_snp>\n\t"+
						"<-p\tpermute_times>\n\t"+					
						"<-l\tmaxlength>\n\t"
						
										);
				System.exit(0);
			}else{
				String  filefolder=null; String level_snp =null; int permutetime=0; int maxlength =0;
//				String itemOrpval =null; double percent=0; int maxlevel =0; int maxlength=0;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-f"))filefolder=args[k+1];
						else if(args[k].equals("-s"))level_snp=args[k+1];
						else if(args[k].equals("-p"))permutetime=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-l"))maxlength=Integer.parseInt(args[k+1]);
						
					}
				}if(filefolder==null||level_snp==null ||maxlength==0||permutetime==0){
					System.out.println("outprefix or file_hdf5 or inprefix can't be null! searchlength focosnp permuteindex can't be 0");
				}else{
					Threshold.item_threshold(filefolder, level_snp, permutetime, maxlength);
					Threshold.pval_threshold(filefolder, level_snp, permutetime, maxlength);
				}
			}
		}else if(function.equals("Thresholdtable")){
			if(args.length==1){
				System.out.println("Get threshold from permutation data");
				System.out.println("Usage: \n\t<-f\tfilefolder>\n\t" +				
						"<-p\tpermute_times>\n\t"+										
						"<-m\tmaxlevel>\n\t"+
						"<-n\tmaxlength>\n\t"+
						"<-q\tpercent>\n\t"+
						"<-o\titemOrpval>\n\t"
										);
				System.exit(0);
			}else{
				String  filefolder=null; int permutetime=0; 
				String itemOrpval =null; double percent=0; int maxlevel =0; int maxlength=0;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-f"))filefolder=args[k+1];						
						else if(args[k].equals("-p"))permutetime=Integer.parseInt(args[k+1]);						
						else if(args[k].equals("-m"))maxlevel=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-n"))maxlength=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-q"))percent=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-o"))itemOrpval=args[k+1];
					}
				}if(filefolder==null||permutetime==0){
					System.out.println("outprefix or file_hdf5 or inprefix can't be null! ");
				}else{
					Threshold.threshtable(filefolder, permutetime, maxlevel, maxlength, percent, itemOrpval);
				}
			}
		}else if(function.equals("Apriori")){
			if(args.length==1){
				System.out.println("Start AprioriGWAS");
				System.out.println("Usage: \n\t<-g\thdf5_file>\n\t" +
						"<-l\tsearch_length>\n\t"+
						"<-o\tresult_folder>\n\t"+
						"<-s\tstep_size>\n\t"+
						"<-i\tstep_index>\n\t"+
//						"<-it\titemthreshtable>\n\t"+	
						"<-pt\tpvalthreshtable>\n\t"+	
						"<-v\tsingleresult>\n\t"										);
				System.exit(0);
			}else{
				String  file_hdf5=null; int search_length=0; String result_folder=null; int step_size =0;  int step_index=0; 
				String itemthresh =null; String singleresult =null; String pvalthresh =null;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-g"))file_hdf5=args[k+1];
						else if(args[k].equals("-l"))search_length=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-o"))result_folder=args[k+1];
						else if(args[k].equals("-s"))step_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-i"))step_index=Integer.parseInt(args[k+1]);	
//						else if(args[k].equals("-it"))itemthresh=args[k+1];
						else if(args[k].equals("-pt"))pvalthresh=args[k+1];
						else if(args[k].equals("-v"))singleresult=args[k+1];
					}
				}if(file_hdf5==null||result_folder==null ||search_length==0||step_size==0 ||itemthresh==null||singleresult==null){
					System.out.println("file_hdf5 or result_folder or threshtable or singleresult can't be null! searchlength step_size step_index can't be 0");
				}else{
					SimpleHDF5 genotype =new SimpleHDF5(file_hdf5);
//					AprioriGWAS.algorithm(genotype, search_length, result_folder, step_size, step_index, threshtable, singleresult)
					AprioriGWAS.algorithm(genotype, search_length, result_folder, step_size, step_index,  pvalthresh, singleresult);
				}
			}
		}else if(function.equals("Getpattern")){
			if(args.length==1){
				System.out.println("Get patterns for significant pairs");
				System.out.println("Usage: \n\t<-g\thdf_file>\n\t" +"<-p\tpval_file>\n\t"+"<-s\tSinglePval_file>\n\t"+
						"<-o\tout_put>\n\t");
				System.exit(0);
			}else{
				String  file_hdf5=null; String  out =null;String pval=null; String singleResult=null;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-g"))file_hdf5=args[k+1];
						else if(args[k].equals("-o"))out=args[k+1];
						else if(args[k].equals("-p"))pval=args[k+1];
						else if(args[k].equals("-s"))singleResult=args[k+1];
					}
				}if(out==null||file_hdf5==null){
					System.out.println("result or file_hdf5 can't be null!");
				}else{
					
					SimpleHDF5 geno =new SimpleHDF5(file_hdf5);
					GetPattern.getsigpatterns(pval, geno, singleResult,out);					
				}
			}
		}

	}


	}


