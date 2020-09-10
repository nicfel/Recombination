package recombination.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import recombination.network.BreakPoints;

public class Partials {
	

	double[][][][] fancyPatterns;
	
//	double[][][] patterns;
//	double[][][] storedPatterns;

	
	List<Integer> ID;
	List<List<BreakPoints>> breaks;
	
	List<Integer> storeID;
	List<List<BreakPoints>> storedBreaks;	
	
	boolean[][] patternIndices;
	boolean[][] storedPatternIndices;
	
	

	
	public Partials(int nrNodes, int nrBreaks, int nrPatterns){
		ID = new ArrayList<>();
		breaks = new ArrayList<>();
	
		patternIndices = new boolean[nrNodes][nrBreaks];
		storedPatternIndices = new boolean[nrNodes][nrBreaks];
		fancyPatterns = new double[nrNodes][nrBreaks][2][nrPatterns];

	}
	
	public void addNode(Integer new_id) {
		int i=ID.indexOf(null);
		if (i!=-1) {
			ID.set(i,new_id);
			breaks.set(i,new ArrayList<>());			
		}else {
			ID.add(new_id);
			breaks.add(new ArrayList<>());
		}
	}
	
	public void addNewNode(Integer new_id) {
		if (ID.contains(new_id))
			return;

		int i=ID.indexOf(null);
		if (i!=-1) {
			ID.set(i,new_id);
			breaks.set(i,new ArrayList<>());			
		}else {
			ID.add(new_id);
			breaks.add(new ArrayList<>());
		}		
	}

	
	public void remove(Integer new_id) {
		if (new_id==null)
			return;

		int i = ID.indexOf(new_id);
		ID.set(i, null);
		breaks.set(i, null);
	}
	
	public void removeBreaks(Integer new_id) {
		int i = ID.indexOf(new_id);
		breaks.set(i, new ArrayList<>());
	}
	

	public void addBreaks(Integer new_id, BreakPoints bp) {
		int i=ID.indexOf(null);
		int j=breaks.get(i).indexOf(new BreakPoints());
		if (j!=-1) {
			breaks.get(i).set(j,bp);			
		}else {
			breaks.get(i).add(bp);			
		}		
	}
	
	public void addBreaksFast(Integer i, BreakPoints bp) {
		int j=breaks.get(i).indexOf(new BreakPoints());
		if (j!=-1) {
			breaks.get(i).set(j,bp);			
		}else {
			breaks.get(i).add(bp);			
		}		
	}

	
	public void replaceBreaks(Integer new_id, BreakPoints old_bp, BreakPoints new_bp) {
		int i = ID.indexOf(new_id);
		int j = breaks.get(i).indexOf(old_bp);
		breaks.get(i).set(j, new_bp);
		
	}
	
	public List<BreakPoints> getBreaks(Integer new_id) {
		return breaks.get(ID.indexOf(new_id));
	}

	
	public double[] getPartials(Integer new_id, BreakPoints bp) {
		int i = ID.indexOf(new_id);
		int j = breaks.get(i).indexOf(bp);
		
		
		if (patternIndices[i][j]) {
			return fancyPatterns[i][j][0];
		}else
			return fancyPatterns[i][j][1];
	}
	
	public double[] getPartialsAdd(Integer new_id, BreakPoints bp) {
		int i = ID.indexOf(new_id);
		int j = breaks.get(i).indexOf(bp);
		if (j==-1) {
			addBreaksFast(i, bp);
			j = breaks.get(i).indexOf(bp);
		}

		if (patternIndices[i][j])
			return fancyPatterns[i][j][0];
		else
			return fancyPatterns[i][j][1];
	}
	
	public double[] getPartialsOperation(Integer new_id, BreakPoints bp) {
		int i = ID.indexOf(new_id);
		int j = breaks.get(i).indexOf(bp);
		if (j==-1) {
			addBreaksFast(i, bp);
			j = breaks.get(i).indexOf(bp);
		}
		
		patternIndices[i][j] = !patternIndices[i][j];

		if (patternIndices[i][j])
			return fancyPatterns[i][j][0];
		else
			return fancyPatterns[i][j][1];
	}

	
	public List<Integer> keySet(){
		return ID;
	}

	public void store() {
		storeID = new ArrayList<>(ID);
		storedBreaks = new ArrayList<>();
		for (int i = 0; i < breaks.size(); i++) {
			if (breaks.get(i)==null) {
				storedBreaks.add(new ArrayList<>());
			}else {
				storedBreaks.add(new ArrayList<>());
				for (BreakPoints bp : breaks.get(i))
					storedBreaks.get(i).add(bp.copy());
			}
		}
		for (int i = 0; i < storeID.size(); i++) {
			if (storeID.get(i)!=null) {
				for (int j = 0; j < storedBreaks.get(i).size();j++) {
					System.arraycopy(patternIndices[i], 0, storedPatternIndices[i], 0, patternIndices[i].length);
				}
			}
		}
	}
	
	
	public void restore() {
//        double[][][] tmp1 = patterns;
//        patterns = storedPatterns;
//        storedPatterns = tmp1;
                
    	List<Integer> tmp2 = ID;
    	ID = storeID;
    	storeID = tmp2;

        List<List<BreakPoints>> tmp3 = breaks;
        breaks = storedBreaks;
        storedBreaks = tmp3;
        
        
		for (int i = 0; i < ID.size(); i++) {
			if (storeID.get(i)!=null) {
				for (int j = 0; j < breaks.get(i).size();j++) {
					System.arraycopy(storedPatternIndices[i], 0, patternIndices[i], 0, patternIndices[i].length);
				}
			}
		}

        
//		System.out.println("restore");
//		int i = ID.indexOf(661898545);
//		if (i!=-1) {
//			int j = breaks.get(i).indexOf(new BreakPoints(0,9999));
//			if (j!=-1)
//				System.out.println(patterns[i][j][0]);
//		}
	}

	public boolean containsKey(Integer new_id, BreakPoints computeFor) {
		return breaks.get(ID.indexOf(new_id)).contains(computeFor);
	}

	
	
}
