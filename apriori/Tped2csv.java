package apriori;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import dataprocess.Miscellaneous;
public class Tped2csv {
	public static void tped2csv(String tped, String tfam, String frq, String csv){
		try{
			BufferedWriter bw =new BufferedWriter(new FileWriter(csv));
			bw.write("rsID,Chr,Location");
			
			BufferedReader brfam =new BufferedReader(new FileReader(tfam));
			String lfam =brfam.readLine(); int rowindex=0;
			while(lfam!=null){
				String[]a =lfam.split(" ");
//				System.out.println(a[a.length-1]);
				if(a[a.length-1].equals("-9")){
					System.out.println(lfam);
					System.out.println("missing phenotype is not allowed");
				}else if(a[a.length-1].equals("1")){
					bw.write(",0");
				}else if(a[a.length-1].equals("2")){
					bw.write(",1");
				}
				rowindex++;
				lfam =brfam.readLine();
			}
			bw.write("\n");
			
			BufferedReader brfrq =new BufferedReader(new FileReader(frq));
			String lfrq =brfrq.readLine();
			lfrq =brfrq.readLine();			
			BufferedReader brped =new BufferedReader(new FileReader(tped));
			String lped =brped.readLine();
			while(lped!=null){
				String[]ageno =lped.split(" ");
				bw.write(ageno[1]+","+ageno[0]+","+ageno[3]);
				String[]bfrq =lfrq.split(" ");
				String[] snpfrq =dataprocess.Miscellaneous.array_no_space(bfrq);
//				System.out.println(snpfrq[1]);
				
				if(!snpfrq[1].equals(ageno[1])){
					System.out.println("snp not match ");
				}else {
					for(int i=4; i<ageno.length; i=i+2){
						if(!ageno[i].equals("0")){
							if(ageno[i].equals(ageno[i+1])){
								if(ageno[i].equals(snpfrq[2])){
									bw.write(","+3);
								}else{
									bw.write(","+1);
								}
							}else{
								bw.write(","+2);
							}
						}else{
							bw.write(","+0);
						}
					}
					bw.write("\n");
				}
				lped =brped.readLine();
				lfrq =brfrq.readLine();		
			}
			bw.flush();
		}catch(Exception e){e.printStackTrace();}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
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
			
			tped2csv(tped, tfam, frq, csv);

		}
	}

}
