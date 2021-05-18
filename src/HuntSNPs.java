

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

public class HuntSNPs {

	public static void main(String[] args){
		
			
			
			
			
			String help = 
					"Usage:\n "+
					"java -jar HuntSNPs.jar -i inputFiles\n" +
					"Options:\n "+
					"-i\tString\tInput Files. There files were derived from QuickSNP\n "+
					"-w\tfloat\twindow size. Size of interval to find mutations (default 10000)\n "+
					"-n\tint\tMinimum number of mutants with mutation in an interval (default: 4)\n " + 
					"-e\tint\tMinimum number of mutants with canonical Sodium Azide mutation in an interval (default: 4)\n " +
					"-a\tdouble\tMaxium reference allele frequency to consider a SNP (default: 0.1)\n " +
					"-c\tdouble\tMinimum coverage to consider a SNP (default: 20)\n " +
					"-s\tdouble\tMaximum mutants allowed to share a position (default: 1)\n " 
					;
			
			try {
			
				CLI cli = new CLI();
				cli.parseOptions(args);
				
				if(!cli.hasOption("i")) {
					throw new CLIParseException("parameter -i is mandatory. Please specify input files.");
				}
				
				Vector<File> snpFiles = new Vector<File>();
				for( Iterator<String> iterator = cli.getArgs("i").iterator(); iterator.hasNext();) {
					String filename = iterator.next();
					File file = new File(filename);
					snpFiles.add(file);
				}
				
				int windowSize = 10000;
				int minSNPs = 4;
				int minEMSSNPs = 4;
				double maxReferenceAlleleFrequency = 0.1;
				int maxMutantsSharePosition = 1;
				int coverage = 20;
				
				
				if(cli.hasOption("w")) {
					windowSize = Integer.parseInt(cli.getArg("w"));
				}
				if(cli.hasOption("n")) {
					minSNPs = Integer.parseInt(cli.getArg("n"));
				}
				if(cli.hasOption("e")) {
					minEMSSNPs = Integer.parseInt(cli.getArg("e"));
				}
				if(cli.hasOption("a")) {
					maxReferenceAlleleFrequency = Double.parseDouble(cli.getArg("a"));
				}
				if(cli.hasOption("c")) {
					coverage = Integer.parseInt(cli.getArg("c"));
				}
				if(cli.hasOption("s")) {
					maxMutantsSharePosition = Integer.parseInt(cli.getArg("s"));
				}
				
				
				hunt4(snpFiles, windowSize, minSNPs, minEMSSNPs, maxReferenceAlleleFrequency, maxMutantsSharePosition, coverage);
				
				
			
			}
			catch (CLIParseException e) {
				e.printStackTrace();
				System.out.println(help);
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			
			
			
		
	}
		
	
	
	
	public static void hunt4(Vector<File> snpFiles, int windowSize, int minMutations, int minEMSMutations, double maxRefAlleleFrequency,int maxMutantsSharePosition, int cov)throws IOException{
		Hashtable<String, Hashtable<String,Hashtable<Integer, String[]>>> mainTable = new Hashtable<String, Hashtable<String,Hashtable<Integer, String[]>>>();
		
		
		//read in SNP files
		for(Enumeration<File> enum_inputFiles = snpFiles.elements(); enum_inputFiles.hasMoreElements();){
			File file = enum_inputFiles.nextElement();
			String fileName = file.getName();
			BufferedReader in = new BufferedReader(new FileReader(file));
			String currentContig = "";
			Hashtable<Integer, String[]> h = new Hashtable<Integer, String[]>();
			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				if(inputline.startsWith("#")){
					continue;
				}
				String[] split = inputline.split("\t"); //contig	position	referenceAllele	alternativeAllele	refAlleleFreq largestAltAlleleFreq	Coverage
				
				if(Integer.parseInt(split[6])<cov) {
					continue;
				}
				
				if(!currentContig.equalsIgnoreCase(split[0])){
					if(!mainTable.containsKey(currentContig)){
						mainTable.put(currentContig, new Hashtable<String,Hashtable<Integer, String[]>>());
					}
					
					mainTable.get(currentContig).put(fileName, h);
					h = new Hashtable<Integer, String[]>();
					currentContig = split[0];
				}
				
				if( 	(split[2].equalsIgnoreCase("A") || split[2].equalsIgnoreCase("T") || split[2].equalsIgnoreCase("G") || split[2].equalsIgnoreCase("C") ) &&
						(split[3].equalsIgnoreCase("A") || split[3].equalsIgnoreCase("T") || split[3].equalsIgnoreCase("G") || split[3].equalsIgnoreCase("C") )){
					
					try{
					String[] a ={split[4], split[6], split[2], split[3]};  //allele frequency, coverage, reference allele, alternative allele
					
					h.put(Integer.parseInt(split[1]), a);
					}catch(ArrayIndexOutOfBoundsException e){
						System.out.println(fileName + "\t" +inputline);
					}
				}	
			}

			in.close();
		}
		
		
		
		
		
		for(Enumeration<String> enumContig = mainTable.keys(); enumContig.hasMoreElements();){
			String contig = enumContig.nextElement();
		
			
			
			//mark unusable mutants
			Hashtable<Integer, Integer> wasThere = new Hashtable<Integer, Integer>();
			HashSet<Integer> remove = new HashSet<Integer>();
			for(Enumeration<String> enumFileNames = mainTable.get(contig).keys(); enumFileNames.hasMoreElements();){
				String fileName = enumFileNames.nextElement();
				for(Enumeration<Integer> enumPositions = mainTable.get(contig).get(fileName).keys();enumPositions.hasMoreElements();){
					Integer position = enumPositions.nextElement();
					int num = 0;
					if(wasThere.containsKey(position)){
						num = wasThere.get(position).intValue();
					}
					num++;
					wasThere.put(position, num);
				}
			}
			
			for(Enumeration<Integer> enumPositions = wasThere.keys(); enumPositions.hasMoreElements();){
				Integer position = enumPositions.nextElement();
				if( wasThere.get(position).intValue() > maxMutantsSharePosition){
					remove.add(position);
				}else{
					
				}
			}
			
			
			
			
			//create a data structure that can quickly be searched with sliding window.
			Hashtable<String, Vector<Integer>> positions = new Hashtable<String,Vector<Integer>>();
			int maxPosition = 0; //need that for later
			for(Enumeration<String> enumFileNames = mainTable.get(contig).keys(); enumFileNames.hasMoreElements(); ){
				String fileName = enumFileNames.nextElement();
				Vector<Integer> v = new Vector<Integer>();
				for(Enumeration<Integer> enumPositions = mainTable.get(contig).get(fileName).keys(); enumPositions.hasMoreElements();){
					Integer position = enumPositions.nextElement();
					if( !remove.contains(position)  && Double.parseDouble( mainTable.get(contig).get(fileName).get(position)[0]) <=maxRefAlleleFrequency ){
						v.add(position);
						if(position>maxPosition){
							maxPosition = position;
						}
					}
					
				}
				Collections.sort(v);
				if(v.size()>0){	
					positions.put(fileName, v);
				}
			}
			
			
			
			if(positions.size() < minMutations){
				continue;
			}
			
			
			
			
			
			
			
			
			Hashtable<String, Integer> indices = new Hashtable<String, Integer>();
			for(Enumeration<String> enumFileNames = positions.keys(); enumFileNames.hasMoreElements();){
				String fileName = enumFileNames.nextElement();
				indices.put(fileName, 0);
			}
			
			while(true){
				int minPosition = maxPosition+1; // only to initialize;
				String fileNameWithLowestPosition = "";
				for(Enumeration<String> enumFileNames = positions.keys(); enumFileNames.hasMoreElements();){
					String fileName = enumFileNames.nextElement();
					if(positions.get(fileName).size() > indices.get(fileName)){
						int position = positions.get(fileName).get(indices.get(fileName));
						if(position < minPosition){
							minPosition = position;
							fileNameWithLowestPosition = fileName;
						}
					}
				}
				
				if(fileNameWithLowestPosition.equalsIgnoreCase("")){ //stop condition
					break;
				}
				
				int numSNPs = 0; 
				int numEMSSNPs = 0;
				
				String report = contig + "\t" + contig + ":"+ (minPosition-1000) +"-" + (minPosition+ windowSize + 1000) ;
				for(Enumeration<String> enumFileNames = positions.keys(); enumFileNames.hasMoreElements();){
					String fileName = enumFileNames.nextElement();
					
					int index = indices.get(fileName).intValue();
					if(positions.get(fileName).size()>index){
						int position = positions.get(fileName).get(indices.get(fileName));
						if( position >= minPosition && position <= minPosition + windowSize){
							numSNPs++;
							String refAllele = mainTable.get(contig).get(fileName).get(position)[2];
							String altAllele = mainTable.get(contig).get(fileName).get(position)[3];
							report = report + "\t" + fileName +"(" +  position+"):"+refAllele+"->"+altAllele;
							//System.out.println(contig + "\t" + fileName + "\t" + index + "\t" + minPosition + "\t"+ position + "\t" + numSNPs);
							if( (refAllele.equalsIgnoreCase("G") && altAllele.equalsIgnoreCase("A")) || (refAllele.equalsIgnoreCase("C") && altAllele.equalsIgnoreCase("T") )){
								numEMSSNPs++;
							}
						}
					}
					
				}
				
				
				if(numSNPs >= minMutations && numEMSSNPs >= minEMSMutations ){
					System.out.println(report);
				}
				
				
				int index = indices.get(fileNameWithLowestPosition).intValue();
				index ++;
				indices.put(fileNameWithLowestPosition, index);
				
			}
			
			
			
			
			
			
			
		}
		
		
		
	}
	
	
	
	
	
		
		
		
		
		
	
	

}
