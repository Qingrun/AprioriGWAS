package statisticalfunctions;

public class Proportion_test {
	public static double Proportiontest(double controlcount, double control_size, double casecount, double case_size){
		//		((p1 - p2)^2) / (p*q(1/n1+1/n2))
		double chival =0;
		double pcontrol = controlcount/control_size;
		double pcase = casecount/case_size;
		double pboth = (controlcount+casecount)/(control_size+case_size);
		chival = (Math.pow((pcontrol - pcase), 2))/((pboth*(1-pboth)*((1/case_size)+(1/control_size))));
		return chival;	
	}
	
	public static double[][] merged(double[][] input){
//		input is nx2 table, 2 stands for case and control, n stands for multiple patterns, merge rare count into one col.		
		
		int small =0; int rare =0;
		for(int row =0; row<input.length; row++){
			if(input[row][0]<5 && input[row][1]<5){
				small++;
				if(input[row][0]>0 || input[row][1]>0){
					rare =1;
				}
				
			}
		}
		double[][] result =new double[input.length-small+rare][input[0].length];
		int rowindex =0;
		for(int row=0; row<input.length; row++){
			if(input[row][0]>=5 || input[row][1]>=5){
				result[rowindex][0] =input[row][0]; 
				result[rowindex][1] =input[row][1]; 
				rowindex++;
			}else if(rare==1){
				result[result.length-1][0] =result[result.length-1][0]+input[row][0]; 
				result[result.length-1][1] =result[result.length-1][1]+input[row][1]; 
			}
		}
		return result;
	}


}
