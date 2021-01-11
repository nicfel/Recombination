package recombination.util;

import java.util.ArrayList;
import java.util.List;

import recombination.network.BreakPoints;

public class Partials {
	

	double[][][] fancyPatterns;
	
	List<Integer> ID;
	List<List<BreakPoints>> breaks;
	List<List<Integer>> startPoint;
	
	List<Integer> storeID;
	List<List<BreakPoints>> storedBreaks;
	List<List<Integer>> storedStartPoint;
	
	boolean[] patternIndices;
	boolean[] storedPatternIndices;
	
	boolean[] nextNr;
	boolean[] storedNextNr;
	
	int currentSize;
	int nrPatterns;
	
	
	public Partials(int initSize, int nrPatterns){
		ID = new ArrayList<>();
		breaks = new ArrayList<>();
		startPoint = new ArrayList<>();
		
		currentSize = initSize;
		this.nrPatterns = nrPatterns;
	
		patternIndices = new boolean[currentSize];
		storedPatternIndices = new boolean[currentSize];
		
		nextNr = new boolean[currentSize];	
		storedNextNr = new boolean[currentSize];
		
		fancyPatterns = new double[currentSize][2][nrPatterns];		
	}
	
	public void addNode(Integer new_id) {
		int i=ID.indexOf(null);
		if (i!=-1) {
			ID.set(i,new_id);
			breaks.set(i,new ArrayList<>());		
			startPoint.set(i,new ArrayList<>());			
		}else {
			ID.add(new_id);
			breaks.add(new ArrayList<>());
			startPoint.add(new ArrayList<>());
		}
	}
	
	public void addNewNode(Integer new_id) {
		if (ID.contains(new_id))
			return;

		int i=ID.indexOf(null);
		if (i!=-1) {
			ID.set(i, new_id);
			breaks.set(i, new ArrayList<>());
			startPoint.set(i, new ArrayList<>());			
		}else {
			ID.add(new_id);
			breaks.add(new ArrayList<>());
			startPoint.add(new ArrayList<>());
		}		
	}

	
	public void remove(Integer new_id) {
		if (new_id==null)
			return;

		int i = ID.indexOf(new_id);
		ID.set(i, null);		
		breaks.set(i, null);

		for (Integer rem : startPoint.get((i)))
			if (rem!=null)
				nextNr[rem] = false;
				
		startPoint.set(i, null);
	}
	
	public void removeBreaks(Integer new_id) {
		int i = ID.indexOf(new_id);
		
		breaks.set(i, new ArrayList<>());
		
		if (startPoint.get(i)!=null)
			for (Integer rem : startPoint.get((i)))
				if (rem!=null)
					nextNr[rem] = false;

		startPoint.set(i, new ArrayList<>());
	}
	

	public void addBreaks(Integer new_id, BreakPoints bp) {
		int i=ID.indexOf(null);
		int j=breaks.get(i).indexOf(new BreakPoints());
		if (j!=-1) {
			breaks.get(i).set(j,bp);
			startPoint.get(i).set(j, nextEmpty());
		}else {
			breaks.get(i).add(bp);			
			startPoint.get(i).add(nextEmpty());
		}		
	}
	
	public void addBreaksFast(Integer i, BreakPoints bp) {
		int j=breaks.get(i).indexOf(new BreakPoints());
		if (j!=-1) {
			breaks.get(i).set(j,bp);
			startPoint.get(i).set(j,nextEmpty());
		}else {
			breaks.get(i).add(bp);			
			startPoint.get(i).add(nextEmpty());
		}		
	}

	
	public void replaceBreaks(Integer new_id, BreakPoints old_bp, BreakPoints new_bp) {
		int i = ID.indexOf(new_id);
		int j = breaks.get(i).indexOf(old_bp);
		breaks.get(i).set(j, new_bp);
		if (new_bp.isEmpty()) {
			nextNr[startPoint.get(i).get(j)] = false;
			startPoint.get(i).set(j, null);
		}
	
		
	}
	
	public List<BreakPoints> getBreaks(Integer new_id) {
		if (ID.indexOf(new_id)==-1)
			return null;
		return breaks.get(ID.indexOf(new_id));
	}

	
	public double[] getPartials(Integer new_id, BreakPoints bp) {
		int i = ID.indexOf(new_id);
		int j = breaks.get(i).indexOf(bp);
		
		if (patternIndices[startPoint.get(i).get(j)]) {
			return fancyPatterns[startPoint.get(i).get(j)][0];
		}else
			return fancyPatterns[startPoint.get(i).get(j)][1];
	}
	
	public double[] getPartialsAdd(Integer new_id, BreakPoints bp) {
		int i = ID.indexOf(new_id);
		int j = breaks.get(i).indexOf(bp);
		if (j==-1) {
			addBreaksFast(i, bp);
			j = breaks.get(i).indexOf(bp);
		}

		if (patternIndices[startPoint.get(i).get(j)])
			return fancyPatterns[startPoint.get(i).get(j)][0];
		else
			return fancyPatterns[startPoint.get(i).get(j)][1];
	}
	
	public double[] getPartialsOperation(Integer new_id, BreakPoints bp) {
		int i = ID.indexOf(new_id);
		int j = breaks.get(i).indexOf(bp);
		if (j==-1) {
			addBreaksFast(i, bp);
			j = breaks.get(i).indexOf(bp);
		}
		
		patternIndices[startPoint.get(i).get(j)] = !patternIndices[startPoint.get(i).get(j)];

		if (patternIndices[startPoint.get(i).get(j)])
			return fancyPatterns[startPoint.get(i).get(j)][0];
		else
			return fancyPatterns[startPoint.get(i).get(j)][1];
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
		
		
		storedStartPoint = new ArrayList<>();
		for (int i = 0; i < startPoint.size(); i++) {
			if (startPoint.get(i)==null) {
				storedStartPoint.add(new ArrayList<>());
			}else {
				storedStartPoint.add(new ArrayList<>());
				for (Integer bp : startPoint.get(i))
					storedStartPoint.get(i).add(bp);
			}
		}

		
		System.arraycopy(patternIndices, 0, storedPatternIndices, 0, patternIndices.length);
		
		System.arraycopy(nextNr, 0, storedNextNr, 0, nextNr.length);

		
	}
	
	
	public void restore() {

    	List<Integer> tmp2 = ID;
    	ID = storeID;
    	storeID = tmp2;

        List<List<BreakPoints>> tmp3 = breaks;
        breaks = storedBreaks;
        storedBreaks = tmp3;
        
        List<List<Integer>> tmp4 = startPoint;
        startPoint = storedStartPoint;
        storedStartPoint = tmp4;        
        
		System.arraycopy(storedPatternIndices, 0, patternIndices, 0, patternIndices.length);		
		System.arraycopy(storedNextNr, 0, nextNr, 0, storedNextNr.length);        
	}

	public boolean containsKey(Integer new_id, BreakPoints computeFor) {
		return breaks.get(ID.indexOf(new_id)).contains(computeFor);
	}

	
	private int nextEmpty() {
		for (int i = 0; i < nextNr.length; i++) {
			if (!nextNr[i]) {
				nextNr[i] = true;
				return i;
			}
		}
		resize();
		return nextEmpty();		
	}
	
	/**
	 * resizes the partials etc.
	 */
	private void resize() {
//		System.err.println("resize matrices to allow for more nodes curr size = " + currentSize + " " + c);
		int newsize = (int) (currentSize*1.5);
		currentSize = newsize;
		
		boolean[] newPatternIndices = new boolean[newsize];
		System.arraycopy(patternIndices, 0, newPatternIndices, 0, patternIndices.length);
		patternIndices = newPatternIndices;
		
		
		boolean[] newStoredPatternIndices = new boolean[newsize];
		System.arraycopy(storedPatternIndices, 0, newStoredPatternIndices, 0, storedPatternIndices.length);
		storedPatternIndices = newStoredPatternIndices;

		boolean[] newNextNr = new boolean[newsize];
		System.arraycopy(nextNr, 0, newNextNr, 0, nextNr.length);
		nextNr = newNextNr;

		boolean[] newStoredNextNr = new boolean[newsize];
		System.arraycopy(storedNextNr, 0, newStoredNextNr, 0, storedNextNr.length);
		storedNextNr = newStoredNextNr;

		double[][][] newFancyPatterns = new double[newsize][2][nrPatterns];
		for (int i = 0; i < fancyPatterns.length;i++) {
			for (int j= 0; j < fancyPatterns[i].length;j++) {
				System.arraycopy(fancyPatterns[i][j], 0, newFancyPatterns[i][j], 0, fancyPatterns[i][j].length);
			}
		}			
		fancyPatterns = newFancyPatterns;
	}
	
	public void purge() {
		ID = new ArrayList<>();
		breaks = new ArrayList<>();
		startPoint = new ArrayList<>();
		patternIndices = new boolean[currentSize];
		nextNr = new boolean[currentSize];	
	}
	
}
