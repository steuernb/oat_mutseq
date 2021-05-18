

import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MPileupLine {
	//String line;
	String[] split;
	Hashtable<String,Integer> alleles;
	
	public MPileupLine (String line){
		//this.line = line;
		split = line.split("\t");
	}
	
	
	public String getChromosome(){
		//return line.split("\t")[0];
		return split[0];
	}
	public char getRefBase(){
		return split[2].toCharArray()[0];
		
	}
	public int getPosition(){
		return Integer.parseInt(split[1]);
	}
	public int getCoverage(){
		return Integer.parseInt(split[3]);
	}
	
	
	/**
	 * 
	 * calculate coverage from column 4 by counting number of , . A T G C a t g c
	 * 
	 * @return
	 */
	public int getRealCoverage(){
		
		int cov = 0;
		
		if( getAlleles().containsKey("A")) {cov = cov + getAlleles().get("A");}
		if( getAlleles().containsKey("T")) {cov = cov + getAlleles().get("T");}
		if( getAlleles().containsKey("G")) {cov = cov + getAlleles().get("G");}
		if( getAlleles().containsKey("C")) {cov = cov + getAlleles().get("C");}
		
		
		return cov;
	}
	
	
	
	public Hashtable<String,Integer> getAlternativeAlleles(){
		Hashtable<String,Integer> h =new Hashtable<String,Integer> ();
		Hashtable<String,Integer> alleles = this.getAlleles();
		for(Enumeration<String> myenum = alleles.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			if(!key.equalsIgnoreCase(""+this.getRefBase())){
				h.put(key,alleles.get(key) );
			}
		}
		return h;
	}
	
	/**
	 * returns the original String representation of the line in the MPileup file, without the line terminator.
	 * 
	 */
	public String toString(){
		String s = split[0];
		for( int i = 1; i< split.length; i++) {
			s = s + "\t" +split[i];
		}
		
		return s;
	}
	
	
	
	public Hashtable<String,Integer>  getAlleles()throws IllegalStateException{
		if(this.alleles != null){
			return alleles;
		}
		
		this.alleles = new Hashtable<String,Integer>();
		String s = "";
		try{
			 s = split[4];
		}catch(ArrayIndexOutOfBoundsException e){return alleles;}
		char[] c = s.toCharArray();
		
		
		
		int position = 0;
		while(position<c.length){
			try {
				
			
				String allele = "";
				
				
				
				if      ( c[position]== '.' || c[position] == ',' ){ allele = this.getRefBase()+"";}
				else if ( c[position]== 'A' || c[position] == 'a' ){ allele = 'A' + "";}
				else if ( c[position]== 'T' || c[position] == 't' ){ allele = 'T' + "";}
				else if ( c[position]== 'G' || c[position] == 'g' ){ allele = 'G' + "";}
				else if ( c[position]== 'C' || c[position] == 'c' ){ allele = 'C' + "";}
				else if ( c[position]== '*' ){ allele = "*" ;}      //this one is new it is a deletion in the read
				else if ( c[position]== '$'                       ){ }  //end of read symbol; do nothing
				else if ( c[position]== '^'                       ){ position++;}// beginning of read symbol. Skip the next symbol this is read mapping quality
				else if ( c[position]== '+'){
					position ++; 
					
					Pattern p = Pattern.compile("([\\d]+)([\\w]+)");
					Matcher m = p.matcher(s.substring(position));
					m.find();
					String mydigit = m.group(1);
					int digit = Integer.parseInt(mydigit);
					allele = "+" + mydigit + s.substring(position+mydigit.length(), position + mydigit.length() + digit).toUpperCase();
					position = position + mydigit.length() + digit -1;
					
				}else if( c[position]== '-'){
					position ++; 
					Pattern p = Pattern.compile("([\\d]+)([\\w]+)");
					Matcher m = p.matcher(s.substring(position));
					m.find();
					String mydigit = m.group(1);
					int digit = Integer.parseInt(mydigit);
					allele = "-" + mydigit + s.substring(position+mydigit.length(), position + mydigit.length() + digit).toUpperCase();
					position = position + mydigit.length() + digit -1;
				}
					
				int num = 0;
				if( alleles.containsKey(allele)){
					num = alleles.get(allele);
				}
				num++;
				if(!allele.equalsIgnoreCase("")){
					alleles.put(allele, num);
				}
				
				
			
				
				position++;
			} catch (IllegalStateException e) {
				System.out.println("\n\n\n\n"+this.toString());
				throw new IllegalStateException (e.getMessage());
			}
		}
		
		
		return this.alleles;
	}
	
	
	
	
	
	
	public boolean hasInDels(int minNumber){
		for(Enumeration<String> myenum = this.getAlleles().keys(); myenum.hasMoreElements();){
			String allele = myenum.nextElement();
			if( allele.startsWith("+") || allele.startsWith("-")){
				if(alleles.get(allele).intValue()>=minNumber){
					return true;
				}
			}
		}
		return false;
	}
	
	public boolean hasInDels(double alleleFrequency){
		for(Enumeration<String> myenum = this.getAlleles().keys(); myenum.hasMoreElements();){
			String allele = myenum.nextElement();
			if( allele.startsWith("+") || allele.startsWith("-")){
				if(alleles.get(allele).doubleValue()/(double)this.getCoverage()>=alleleFrequency){
					return true;
				}
			}
		}
		return false;
	}
	
	public boolean hasInDels(){
		return hasInDels(0);
	}
	
	public String getLargestAlternativeAllele(){
		Hashtable<String,Integer>  h =  this.getAlternativeAlleles();
		String allele = "";
		int i = 0;
		for(Enumeration<String> myenum = h.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			int j = h.get(key);
			if(j>i){
				i = j;
				allele = key;
			}
		}
		
		return allele;
		
	}
	
	public int getReferenceAlleleCount(){
		if( this.getAlleles().containsKey(this.getRefBase()+"")){	
			return this.alleles.get(this.getRefBase()+"");
		}else{
			return 0;
		}
	}
	
	public int getLargestAlternativeAlleleCount(){
		Hashtable<String,Integer>  h =  this.getAlternativeAlleles();
		//String allele = "";
		int i = 0;
		for(Enumeration<String> myenum = h.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			int j = h.get(key);
			if(j>i){
				i = j;
				//allele = key;
			}
		}
		
		return i;
	}
	
	public double getAlleleFrequency(String allele) {
		this.getAlleles();
		
		if(this.alleles.containsKey(allele) || this.alleles.containsKey(allele.toUpperCase())) {
			
			return (double)this.alleles.get(allele.toUpperCase()) / (double)this.getRealCoverage();
			
		}else {
			return 0;
		}
		
	}
	
	public double getLargestAlternativeAlleleFrequency(){
		double cov = this.getCoverage();
		try{
			double allele = this.getAlternativeAlleles().get(getLargestAlternativeAllele());
			return allele / cov;
		}catch(NullPointerException e){return 0.0;}
	}
	
	public double getReferenceAlleleFrequency(){
		int ref = 0;
		if(this.getAlleles().containsKey(this.getRefBase()+"")){
			ref = this.getAlleles().get(this.getRefBase()+"");
		}
		
		return ((double) ref) / ((double)this.getRealCoverage()) ;
	}
	
	public Vector<Integer> getAlternativeAlleleFrequencies(){
		Vector<Integer>  v = new Vector<Integer>();
		for(Enumeration<Integer> myenum = this.getAlternativeAlleles().elements(); myenum.hasMoreElements();){
			v.add(myenum.nextElement());
		}
		Collections.sort(v);
		return v;
	}
	
	
	
	
	
}
