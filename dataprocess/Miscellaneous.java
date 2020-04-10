package dataprocess;

public class Miscellaneous {
	public static String[] array_no_space(String[] array){
		int elements =0;
		for(int i=0; i<array.length; i++){
			if(!array[i].equals("")){
				elements++;
			}
		}
		String[] result = new String[elements];
		int index =0;
		for(int i=0; i<array.length; i++){
			if(!array[i].equals("")){
				result[index]=array[i];
				index++;
			}
		}
		return result;
	}

//	public static double[][] case_control_genotype(char[] ccindex, char[] genotype){
//		double[][] genotype_contingency_table = new double[2][3];
//		for(int i=0; i<ccindex.length; i++){
//			if(ccindex[i]=='1' && genotype[i]=='1'){
//				genotype_contingency_table[0][0]++;
//			}else if(ccindex[i]=='1' && genotype[i]=='2'){
//				genotype_contingency_table[0][1]++;
//			}else if(ccindex[i]=='1' && genotype[i]=='3'){
//				genotype_contingency_table[0][2]++;
//			}else if(ccindex[i]=='2' && genotype[i]=='1'){
//				genotype_contingency_table[1][0]++;
//			}else if(ccindex[i]=='2' && genotype[i]=='2'){
//				genotype_contingency_table[1][1]++;
//			}else if(ccindex[i]=='2' && genotype[i]=='3'){
//				genotype_contingency_table[1][2]++;
//			}
//		}
//		return genotype_contingency_table;		
//	}
//	public static double[][] case_control_genotype(String[] ccindex, char[] genotype){
//		double[][] genotype_contingency_table = new double[2][3];
//		for(int i=0; i<ccindex.length; i++){
//			if(ccindex[i].equals("1") && genotype[i]=='1'){
//				genotype_contingency_table[0][0]++;
//			}else if(ccindex[i].equals("1") && genotype[i]=='2'){
//				genotype_contingency_table[0][1]++;
//			}else if(ccindex[i].equals("1")&& genotype[i]=='3'){
//				genotype_contingency_table[0][2]++;
//			}else if(ccindex[i].equals("2") && genotype[i]=='1'){
//				genotype_contingency_table[1][0]++;
//			}else if(ccindex[i].equals("2") && genotype[i]=='2'){
//				genotype_contingency_table[1][1]++;
//			}else if(ccindex[i].equals("2") && genotype[i]=='3'){
//				genotype_contingency_table[1][2]++;
//			}
//		}
//		return genotype_contingency_table;		
//	}
	
	public static double[][] case_control_genotype(String[] ccindex, String[] genotype){
		for(int i=0; i<ccindex.length; i++){
			System.out.print(ccindex[i]+"\t");
		}System.out.println();
		double[][] genotype_contingency_table = new double[2][3];
		for(int i=0; i<ccindex.length; i++){
			if(ccindex[i].equals("1") && genotype[i].equals("1")){
				genotype_contingency_table[0][0]++;
			}else if(ccindex[i].equals("1") && genotype[i].equals("2")){
				genotype_contingency_table[0][1]++;
			}else if(ccindex[i].equals("1") && genotype[i].equals("3")){
				genotype_contingency_table[0][2]++;
			}else if(ccindex[i].equals("2") && genotype[i].equals("1")){
				genotype_contingency_table[1][0]++;
			}else if(ccindex[i].equals("2") && genotype[i].equals("2")){
				genotype_contingency_table[1][1]++;
			}else if(ccindex[i].equals("2") && genotype[i].equals("3")){
				genotype_contingency_table[1][2]++;
			}
		}
		return genotype_contingency_table;		
	}

}
