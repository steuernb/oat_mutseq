

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;


public class QuickSNP {

	public static void main(String[] args) {
		
	
		
			
			
			String help = 
					"Usage:\n "+
					"java -jar QuickSNP.jar -i inputFile\n" +
					"Options:\n "+
					"-i\tString\tInput File in mpileup format\n "+
					"-o\tString\tOutput File (default:stdout)\n " +
					"-a\tfloat\tMaximum reference allele frequency (default 0.1)\n "+
					"-c\tint\tMinimum coverage to record (default: 5)\n " + 
					"-e\t\tonly output EMS specific SNPs. These are C->T and G->A";
			
			
		
		try {	
			
			
			
			
			
			CLI cli = new CLI();
			cli.parseOptions(args);
			
			File inputFile=new File("");
			File outputFile=new File("");
			double maxRefAlleleFrequency = 0.1;
			int minCov = 5;
			boolean ems = false;
			boolean writeStdout = true;
			
			if(!cli.hasOption("i") ) {
				throw new CLIParseException("No input file specified.");
			}else {
				inputFile = new File(cli.getArg("i"));
			}
			
			if( cli.hasOption("o")) {
				writeStdout = false;
				outputFile = new File(cli.getArg("o"));
			}
			
			if( cli.hasOption("e")) {
				ems = true;
			}
			
			
			if( cli.hasOption("c")) {
				minCov = Integer.parseInt(cli.getArg("c"));
			}
			
			if( cli.hasOption("a")) {
				maxRefAlleleFrequency = Double.parseDouble(cli.getArg("a"));
			}
			
			
			if(writeStdout) {
				
			}else {
				writeSNPs(inputFile, outputFile, minCov, maxRefAlleleFrequency, ems, writeStdout);
			}
			
		}
		catch (CLIParseException e) {
			e.printStackTrace();
			System.out.println(help);
		}
		catch (IOException e) {
			e.printStackTrace();
		}
			
			
	}
		
	
	
	
	
	
		public static void writeSNPs(File inputFile, File outputFile, int cov, double refAlleleFreq, boolean ems, boolean stdout)throws IOException{
		
			BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
			
			if(stdout) {
				System.out.println("#contig\tposition\trefBase\taltBase\trefAlleleFrequency\tlargestAlternativeAlleleFrequency\tcoverage");
			}else {
				out.write("#contig\tposition\trefBase\taltBase\trefAlleleFrequency\tlargestAlternativeAlleleFrequency\tcoverage");
				out.newLine();
					
			}
			
			BufferedReader in;
			
			FileInputStream fis = new FileInputStream(inputFile);
			byte[] bytes = new byte[2];
			fis.read(bytes);
			int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
			boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
			fis.close();

			
			if(gzip){
				in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile))));
			}else{
				in = new BufferedReader((new FileReader(inputFile)));
			}
			
			
		
			
			String inputline = in.readLine();
			while (inputline != null) {
				MPileupLine line = new MPileupLine(inputline);
				String alt = line.getLargestAlternativeAllele();
				String ref = line.getRefBase()+"";
				if( line.getRealCoverage() >= cov  && line.getReferenceAlleleFrequency()<= refAlleleFreq ){
					if(ems &&  (   (!ref.equalsIgnoreCase("C") || !alt.equalsIgnoreCase("T")   )  &&    (!ref.equalsIgnoreCase("G") || !alt.equalsIgnoreCase("A")          )  )      ){  //hope this is right!
						inputline = in.readLine();
						continue;
					}
					
					if(stdout) {
						System.out.println(line.getChromosome() + "\t" + line.getPosition() +"\t" + line.getRefBase() +"\t" + line.getLargestAlternativeAllele() + "\t" + line.getReferenceAlleleFrequency() +"\t"+line.getLargestAlternativeAlleleFrequency() +"\t" + line.getCoverage());
					}else{
						out.write(line.getChromosome() + "\t" + line.getPosition() +"\t" + line.getRefBase() +"\t" + line.getLargestAlternativeAllele() + "\t" + line.getReferenceAlleleFrequency() +"\t"+line.getLargestAlternativeAlleleFrequency() +"\t" + line.getCoverage());
						out.newLine();
					}
					
					
				}
				
				
				inputline = in.readLine();
			}
			in.close();
			
			out.close();
			
		}	
			
			
			

			
		
	

}
