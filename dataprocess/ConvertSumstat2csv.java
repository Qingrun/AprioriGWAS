package dataprocess;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;

public class ConvertSumstat2csv {
	public static void convert(String sumstat, String csv) {
		try {
			BufferedReader br =new BufferedReader(new FileReader(sumstat));
			String lbr =br.readLine();
			int line =0; int col=0;
			while(lbr!=null) {
				String[] a =lbr.split(" ");
				String[] b =Miscellaneous.array_no_space(a);
				line++; col=b.length;
				lbr =br.readLine();
			}
			String[][] table =new String[line][col];
			br =new BufferedReader(new FileReader(sumstat));
			lbr =br.readLine();
			int lineindex =0;
			while(lbr!=null) {
				String[] a =lbr.split(" ");
				String[] b =Miscellaneous.array_no_space(a);
				for(int i=0; i<b.length; i++) {
					table[lineindex][i] =b[i];
				}
				lineindex++; col=b.length;
				lbr =br.readLine();
			}br.close();
			BufferedWriter bw =new BufferedWriter(new FileWriter(csv));
			bw.write("rs,chr,loc");
			for(int j=0; j<table[0].length-3; j++) {
				bw.write(","+table[table.length-1][j]);
			}bw.write("\n");
			for(int i=0; i<table.length-1; i++) {
				bw.write(table[i][table[i].length-1]+","+table[i][table[i].length-3]+","+table[i][table[i].length-2]);
				for(int j=0; j<table[i].length-3; j++) {
					bw.write(","+table[i][j]);
				}bw.write("\n");
			}
			bw.flush();bw.close();
		}catch(Exception e){e.printStackTrace();}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String prefix="/Users/qingrunzhang/Dropbox (SH Corp)/Qingrun/My_Projects/AprioriGWAS/";
		String sumstat =prefix+"sumstatAMDchromsort.dat";
		String csv =prefix+"AMD.csv";
		convert(sumstat, csv);

	}

}
