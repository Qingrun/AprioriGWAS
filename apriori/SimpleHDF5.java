package apriori;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;


public class SimpleHDF5 {
//	data format csv, first row sample ID, the second row with 
	public int sample_size;
	public int num_sites_total;
	public int block_size=1000;
	public String[] locations; //chr_location
	public String[] rs_ids; //rs_id, chr, location
	public int[] snp_order;
//	String[] sample_ids;
	public String[] pheno_index;//for case control sample 0 stands for control, 1 stands for case, the first row of data for pheno index
	public int block_number;
	public IHDF5Reader reader;
//	Iterable<HDF5MDDataBlock<MDDoubleArray>>[] index_fast_blocks;
	public String index_fast_paths;
	
	public static SimpleHDF5 importCSV(String file_csv,String file_hdf5,int block_size){
		long startTime = System.currentTimeMillis();	
		System.out.println("Start importing data from "+file_csv+" to HDF5 format.");
		try{
			BufferedReader br= new BufferedReader(new FileReader(file_csv));
			String line=br.readLine();//header
			String[] temp=line.split(",");
			int snp_num =0;
			int sample_size=temp.length-3;//rs_id, chr, location, then followed by case control index
			
			String[] pheno_index =new String[sample_size];
			for(int k=0;k<sample_size;k++)pheno_index[k]=temp[k+3]; 
			
			ArrayList<String> var_location=new ArrayList<String>();
			ArrayList<String> var_rs_list = new ArrayList<String>();
			ArrayList<Integer> var_order = new ArrayList<Integer>();
			line=br.readLine(); // the first data line
			temp=line.split(",");
			while(line!=null){
				snp_num++;
				temp=line.split(",");
				var_location.add(temp[1]+"_"+temp[2]);
				var_rs_list.add(temp[0]);	
				var_order.add(snp_num);
				line=br.readLine();
			}
			String[] location_info = new String[var_location.size()];
			String[] rs_info =new String[var_rs_list.size()];
			int[] rs_order_info = new int[var_order.size()];
			for(int i=0; i<var_location.size(); i++){
				location_info[i] =var_location.get(i);
				rs_info[i] =var_rs_list.get(i);
				rs_order_info[i] =var_order.get(i);
			}
			int bloc_num =rs_info.length/block_size+1;
			System.out.println("Finished scanning the .csv file and found "+snp_num+" snps.\n");
	
			File hdf5_old=new File(file_hdf5);
			if(hdf5_old.exists()){hdf5_old.delete();System.out.println(file_hdf5+" alreday exist. Deleted!");}
			IHDF5Writer writer = HDF5Factory.open(file_hdf5);
			// CODE FOR constructing above structure (all groups and tables and importing data)
			writer.createGroup("/genotype");
			writer.setIntAttribute("/genotype", "sample_size", sample_size);
			writer.setIntAttribute("/genotype", "block_size", block_size);
			writer.setIntAttribute("/genotype", "num_sites_total", snp_num);
			writer.setIntAttribute("/genotype", "block_number", bloc_num);			
			writer.writeIntArray("/genotype/snp_order", rs_order_info);
			writer.writeStringArray("/genotype/rs_ids", rs_info); 
			writer.writeStringArray("/genotype/location", location_info);
			writer.writeStringArray("/genotype/pheno_index",pheno_index);
			writer.createIntMatrix("/genotype/genotype_data",  snp_num, sample_size, block_size, sample_size);
							
			// Read data again to fill matrix:
			System.out.println("Reading data again to create HDF5 file for genotypes.");
			br= new BufferedReader(new FileReader(file_csv));
			line=br.readLine();//header
			int[][] data4thisblock=new int[block_size][sample_size];
			line=br.readLine();//data line
			int current_block=0;
			int num_of_variants_in_block=0;
			while(line!=null){
				temp=line.split(",");
				if(num_of_variants_in_block==block_size){//output one block to h5	
					
					writer.writeIntMatrixBlock("/genotype/genotype_data",data4thisblock,current_block,0);
					data4thisblock =new int[block_size][sample_size];
					num_of_variants_in_block=0;
					current_block++;				
				}
				for(int k=0;k<sample_size;k++){
					data4thisblock[num_of_variants_in_block][k]=Integer.parseInt(temp[k+3]); //different from a.thaliana
				}
				num_of_variants_in_block++;
				line=br.readLine();
			}// write the final block:
			writer.writeIntMatrixBlockWithOffset("/genotype/genotype_data",
					data4thisblock, num_of_variants_in_block, sample_size, current_block*block_size, 0L);
			//write the indexes files		
			System.out.println("Finished writing genotypes. Now the HDF5 data is ready to use.");
			System.out.println("sample_size\t"+sample_size);
			writer.close();			
			long endTime = System.currentTimeMillis();	
			System.out.println("Time consume\t"+(endTime-startTime)/1000+"seconds");
		}catch(Exception e){e.printStackTrace();}		
		return new SimpleHDF5(file_hdf5);
		
	}
	
	/*
	 * constructor: load an hdf5 object
	 */
	public SimpleHDF5(String file_hdf5){
		this.reader=HDF5Factory.openForReading(file_hdf5);	
		this.block_number =this.reader.getIntAttribute("/genotype","block_number");
		this.sample_size = this.reader.getIntAttribute("/genotype", "sample_size");
		this.num_sites_total=this.reader.getIntAttribute("/genotype","num_sites_total");
		this.block_size=this.reader.getIntAttribute("/genotype", "block_size");		
		this.locations=this.reader.readStringArray("/genotype/location");
		this.rs_ids =this.reader.readStringArray("/genotype/rs_ids");
		this.snp_order =this.reader.readIntArray("/genotype/snp_order");
		this.pheno_index =this.reader.readStringArray("/genotype/pheno_index");
		this.index_fast_paths="/genotype/genotype_data";
		System.out.println("Finsehd reading "+file_hdf5);
	}
	
	public int[] load_snp_order_by_index(int snp_index){
		int block_index = (snp_index-1)/this.block_size;
		if(block_index!=(this.block_number-1)){
			int[] snpindex =new int[this.block_size]; 
			for(int i =block_index*this.block_size; i<(block_index+1)*this.block_size; i++){
				snpindex[i-block_index*this.block_size] =this.snp_order[i];
			}
			return snpindex;
		}else{
			int[] snpindex =new int[(int)(this.num_sites_total-(this.block_number-1)*this.block_size)];
			for(int i =(this.block_number-1)*this.block_size; i<this.num_sites_total; i++){
				snpindex[i-(this.block_number-1)*this.block_size] =this.snp_order[i];
			}
			return snpindex;
		}
	}
	public String[] load_rs_id_by_index(int snp_index){
		int block_index = (snp_index-1)/this.block_size;
		if(block_index!=(this.block_number-1)){
			String[] rsids =new String[this.block_size]; 
			for(int i =block_index*this.block_size; i<(block_index+1)*this.block_size; i++){
				rsids[i-block_index*this.block_size] =this.rs_ids[i];
			}
			return rsids;
		}else{
			String[] rsids =new String[(int)(this.num_sites_total-(this.block_number-1)*this.block_size)];
			for(int i =(this.block_number-1)*this.block_size; i<this.num_sites_total; i++){
				rsids[i-(this.block_number-1)*this.block_size] =this.rs_ids[i];
			}
			return rsids;
		}
	}
	public int[][] load_variants_by_index(int snp_index){
		int block_index = (snp_index-1)/this.block_size;
		System.out.println(block_index+"\t"+this.block_number);
		if(block_index!=(this.block_number-1)){
			return reader.readIntMatrixBlock(index_fast_paths, block_size, sample_size, block_index, 0);
		}else{
			return reader.readIntMatrixBlockWithOffset(index_fast_paths, (int)(this.num_sites_total-(this.block_number-1)*this.block_size),
					this.sample_size, this.block_size*(this.block_number-1), 0);
		}
	}
	

}
