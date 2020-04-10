package apriori;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
public class GetPattern {
	public static void getsigpatterns(String pval, SimpleHDF5 genotype, String singleResult,String out){	
		try{
			BufferedReader brs =new BufferedReader(new FileReader(singleResult));
			HashMap<String, Double> singlePval =new HashMap<String, Double>();
			String lbrs =brs.readLine();
			while(lbrs!=null){
				String[] a =lbrs.split("\t");
				double p =Double.parseDouble(a[2]);
				singlePval.put(a[1], p);
				lbrs =brs.readLine();
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));		
			BufferedReader br =new BufferedReader(new FileReader(pval));
			String brl =br.readLine();
			while(brl!=null){
				String[] a =brl.split("\t");
				String[] b =a[0].split("_");
				String[] s =a[1].split("_");
				bw.write(brl); String snpindex =(Integer.parseInt(s[0])-1)+"";
				bw.write("\t"+b[0]+"\t"+singlePval.get(b[0]));
				for(int i=1; i<b.length; i++){
					snpindex=snpindex+"_"+(Integer.parseInt(s[i])-1);
					bw.write("\t"+b[i]+"\t"+singlePval.get(b[i]));
				}bw.write("\n");
				bw.write("pattern\tcontrol\tcase\n");
				System.out.println(a[1]+"\t"+snpindex);
				Pattern_support pats =new Pattern_support(genotype, snpindex);
				HashMap<String,double[]> patterns =pats.pat_support;		
				for(String key:patterns.keySet()){//check each pattern whether they passed proportion test, if past, then write to item file
					double[] count = patterns.get(key);
					bw.write(key+"\t"+count[0]+"\t"+count[1]+"\n");
				}
				bw.write("_______________________________\n");
				brl =br.readLine();
			}	
			bw.flush();
		}catch(Exception e){e.printStackTrace();}
		
	}

}
