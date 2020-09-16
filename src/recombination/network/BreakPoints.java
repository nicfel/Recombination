package recombination.network;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import beast.evolution.tree.Node;

/**
 * Keep track of recombination break points in the network
 * @author Nicola Müller
 */
public class BreakPoints {
	
	public List<Range> breakPoints;	
	private BreakPoints leftBreakPoints;
	private BreakPoints rightBreakPoints;   
//	private int hashCode = -1;
	
	final static RangeComparator rc = new RangeComparator();

		
	public BreakPoints(int totalLength) { 
		breakPoints = new ArrayList<>();
		breakPoints.add(new Range(0, totalLength-1));
		
//		updateHashCode();
	}

	
	public BreakPoints(int start, int end) { 
		breakPoints = new ArrayList<>();
		breakPoints.add(new Range(start, end));
		
//		updateHashCode();
	}
	
	public BreakPoints() { }
	
	public void init(List<Integer> breakPoints) { 
		this.breakPoints = new ArrayList<>();
		for (int i = 0; i < breakPoints.size(); i=i+2)
			this.breakPoints.add(new Range(breakPoints.get(i), breakPoints.get(i+1)));
		
//		updateHashCode();
	}
	
	public BreakPoints(List<Range> breakPoints) { 
		if (breakPoints!=null)
			this.breakPoints = new ArrayList<>(breakPoints);
		
//		updateHashCode();
	}

		
	
	/**
	 * computes new break points list based on a new break point introduced
	 * @param breakpoint
	 */
	public void computeLeftAndRight(int breakpoint) {
		List<Range> leftBreakPointsList = new ArrayList<>();
		List<Range> rightBreakPointsList = new ArrayList<>();
			
		
		int i=0;
		
		while (breakPoints.get(i).to <= breakpoint) {
			leftBreakPointsList.add(breakPoints.get(i));
			i++;
		}
		if (breakPoints.get(i).from <= breakpoint) {
			leftBreakPointsList.add(new Range(breakPoints.get(i).from, breakpoint));
			rightBreakPointsList.add(new Range(breakpoint+1, breakPoints.get(i).to));
			i++;
		}
		while (i < breakPoints.size()) {
			rightBreakPointsList.add(breakPoints.get(i));
			i++;
		}
		
		leftBreakPoints= new BreakPoints(leftBreakPointsList);
		rightBreakPoints= new BreakPoints(rightBreakPointsList);
	}

	/**
	 * gets all the info from left of a new breakpoint 
	 * @param position
	 * @return
	 */
	public BreakPoints getLeft() {
		return leftBreakPoints.copy();
	}
	
	/**
	 * gets all the info from the right of a new breakpoint 
	 * @param position
	 * @return
	 */
	public BreakPoints getRight() {
		return rightBreakPoints.copy();
	}

	/**
	 * gets the difference between the first and last event
	 * @return
	 */
	public double getLength() {
		if (!isEmpty())
			return breakPoints.get(breakPoints.size()-1).to-breakPoints.get(0).from+1;
		
		return -1;
	}
	
	/**
	 * gets the total amount of genetic material encompassed
	 * by this break point as defined by 
	 * @return
	 */
	public double getGeneticLength() {
		if (isEmpty())
			return 0;

		int l = 0;
		for (int i = 0; i < breakPoints.size(); i++)
			l += breakPoints.get(i).size();
		return l;
	}
	
	public int getGeneticLengthInt() {
		if (isEmpty())
			return 0;

		int l = 0;
		for (int i = 0; i < breakPoints.size(); i++)
			l += breakPoints.get(i).size();
		return l;
	}
	
	/**
	 * gets the number of different breakpoints	 
	 * */
	public int size() {
		if (isEmpty())
			return 0;
		return breakPoints.size();
	}
	
	/**
	 * gets the i-th range	 
	 * */
	public Range getRange(int i) {
		return breakPoints.get(i);
	}

	
	public boolean isEmpty() {
		if (breakPoints==null)
			return true;
		
		return breakPoints.size()==0;
	}
	

	public boolean withinLimits(int breakpoint) {
		if (this.breakPoints.get(0).from<=breakpoint && this.breakPoints.get(this.breakPoints.size()-1).to>breakpoint)
			return true;
		
		return false;		
	} 
	
	public boolean contains(int point) {
		for (int i=0; i < size();i++) {
			if (point>=getRange(i).from && point <= getRange(i).to)
				return true;			
		}
		return false;
	}
	
	public BreakPoints copy() {
		BreakPoints newBreakPoints = new BreakPoints(breakPoints);
		return newBreakPoints;
	}
		
	public String toString() {
		if (isEmpty())
			return "null";
		
		String val="";
		for (int i = 0; i < this.breakPoints.size(); i++) {
			val = val + "," + this.breakPoints.get(i).toString();
		}
				
		return val.substring(1);
		
	}
	
	
	
	public static class RangeComparator implements Comparator<Range> {
	    final int lessThan = -1;
	    final int greaterThan = 1;

	    @Override
	    public int compare(Range ra, Range rb) {
	    	if (ra.from<rb.from)
	    		return lessThan;
	    	else
	    		return greaterThan;
	    	
	    }

	}

	/**
	 * compute the intersection between this.breakpoonts and breakpoints
	 * @param breakPoints
	 */
	public void and(BreakPoints breakPoints) {
		List<Range> newBreaks = new ArrayList<>();
		if (breakPoints==null || breakPoints.isEmpty() || isEmpty()) {
			this.breakPoints = null;
			return;
		}

		
		int j = 0;
		

		for (int i = 0; i < this.breakPoints.size(); i++) {
			while (breakPoints.breakPoints.get(j).to < this.breakPoints.get(i).from) {
				j++;
				if (j==breakPoints.breakPoints.size()) {
					setBreakPointsAnd(newBreaks);
					return;
				}

			}

			while (breakPoints.breakPoints.get(j).from <= this.breakPoints.get(i).to) {
				Range newR = this.breakPoints.get(i).getOverlap(breakPoints.breakPoints.get(j));
				if (newR!=null)
					newBreaks.add(newR);
				 
				if (breakPoints.breakPoints.get(j).to <= this.breakPoints.get(i).to)
					j++;
				else
					break;
				
				if (j==breakPoints.breakPoints.size()) {
					setBreakPointsAnd(newBreaks);
					return;
				}				
			}	
		}
		setBreakPointsAnd(newBreaks);
	}

	
	/**
	 * remove breakpoints from this.breakpoints
	 * @param breakPoints
	 */
	public void andNot(BreakPoints breakPoints) {		
		if (breakPoints==null || breakPoints.isEmpty() || isEmpty())
			return;
		
		List<Range> newBreaks = new ArrayList<>();
		int j = 0;
		

		
		for (int i = 0; i < this.breakPoints.size(); i++) {
			boolean rangeAdded = false;
			while (breakPoints.breakPoints.get(j).to < this.breakPoints.get(i).from) {
				j++;
				if (j==breakPoints.breakPoints.size()) {
					if (!rangeAdded)
						newBreaks.add(this.breakPoints.get(i));
					i++;
					while (i < this.breakPoints.size()) {
						newBreaks.add(this.breakPoints.get(i));
						i++;
					}
					setBreakPoints(newBreaks);
					return;
				}

			}

			while (breakPoints.breakPoints.get(j).from <= this.breakPoints.get(i).to) {
				List<Range> newR = this.breakPoints.get(i).getRemoved(breakPoints.breakPoints.get(j));
				if (newR!=null) {
					newBreaks.addAll(newR);
					rangeAdded = true;
				}
				 
				if (breakPoints.breakPoints.get(j).to <= this.breakPoints.get(i).to)
					j++;
				else
					break;
				
				if (j==breakPoints.breakPoints.size()) {
					if (!rangeAdded)
						newBreaks.add(this.breakPoints.get(i));
					i++;
					while (i < this.breakPoints.size()) {
						newBreaks.add(this.breakPoints.get(i));
						i++;
					}
					setBreakPoints(newBreaks);
					return;
				}	
			}	
			
			if (!rangeAdded)
				newBreaks.add(this.breakPoints.get(i));
				
		}
		setBreakPoints(newBreaks);
	}
	
	/**
	 * newly sets the BreakPoints after checking that the overlaps are correct
	 * @param newBreaks
	 */
	public void setBreakPoints(List<Range> newBreaks) {
		for (int i = newBreaks.size()-1; i > 0 ; i--) {
			if (newBreaks.get(i-1).from > newBreaks.get(i).from) {
				newBreaks.get(i-1).to = newBreaks.get(i).to;
				newBreaks.remove(i);
			}
		}
		this.breakPoints = new ArrayList<>(newBreaks);
//		updateHashCode();
	}
	
	/**
	 * newly sets the BreakPoints after and operation checking that the overlaps are correct
	 * @param newBreaks
	 */
	public void setBreakPointsAnd(List<Range> newBreaks) {
		for (int i = newBreaks.size()-1; i > 0 ; i--) {
			if (newBreaks.get(i-1).to > newBreaks.get(i).from) {
				newBreaks.get(i-1).to = newBreaks.get(i).to;
				newBreaks.remove(i);
			}
		}
		this.breakPoints = new ArrayList<>(newBreaks);
//		updateHashCode();
	}

	@Override
	public boolean equals(Object obj) {
        if (this == obj)
    		return true;
        
        
        
        BreakPoints bp = (BreakPoints) obj;
        
        if (isEmpty() && bp.isEmpty())
        	return true;
        	
        if (bp.size()!=this.size())
        	return false;
        
        if (bp.getRange(0).from!=getRange(0).from)
        	return false;
        
        if (bp.getRange(0).to!=getRange(0).to)
        	return false;
        
        for (int i=1; i < size(); i++) {
            if (bp.getRange(i).from!=getRange(i).from)
            	return false;
            
            if (bp.getRange(i).to!=getRange(i).to)
            	return false;
        }        
        return true;
	}
	
//	public void updateHashCode() {
//		// TODO more unique hashcode
//		hashCode = 0;
//		for (int i = 0; i < size(); i++) {
//			hashCode += (getRange(i).to*getLength() + getRange(i).from*getLength())*(i+13)+ getLength();
//		}
//	}

	
//	@Override
//	public int hashCode() {
//		if (hashCode==-1)
//			updateHashCode();
//		
//		return hashCode;
//	}



	/**
	 * combine two lists of breakpoints
	 * @param breakPoints
	 */
	public void or(BreakPoints breakPoints) {
		if (this.breakPoints==null) {
			if (breakPoints==null || breakPoints.isEmpty()) {
				return;
			}
			
			this.breakPoints = new ArrayList<>(breakPoints.breakPoints);
			return;
		}
		if (breakPoints==null || breakPoints.isEmpty()) {
			return;
		}
			
		// make a new list containing all breakpoints
		List<Range> newBreaks = new ArrayList<>();
		this.breakPoints.addAll(breakPoints.breakPoints);
		Collections.sort(this.breakPoints, rc);
				
		int nextto = this.breakPoints.get(0).to;
		int lastfrom = this.breakPoints.get(0).from;
		for (int i = 1; i < this.breakPoints.size(); i++) {
			if (this.breakPoints.get(i).from>(nextto+1)) {
				newBreaks.add(new Range(lastfrom, nextto));
				nextto = this.breakPoints.get(i).to;
				lastfrom = this.breakPoints.get(i).from;
			}else {
				nextto = Math.max(nextto, this.breakPoints.get(i).to);
			}
				
		}
		newBreaks.add(new Range(lastfrom, nextto));				
					
		this.breakPoints = new ArrayList<>(newBreaks);	
//		updateHashCode();
	}
	
	public int getMin() {
		if (isEmpty())
			return -1;
		
		return breakPoints.get(0).from;
	}
	
	public int getMax() {
		if (isEmpty())
			return -1;
		
		return breakPoints.get(breakPoints.size()-1).to;
	}

	
	public Range getNewRange(int from, int to) {
		return new Range(from, to);
	}

	public class Range{
		public int from;
		public int to;
	
		public Range(int from, int to) {
			this.from = from;
			this.to = to;
		}
				
		public void combine(Range range) {
			this.from = Math.min(from, range.from);
			this.to = Math.max(to, range.to);
		}
		
		public int size() {
			return to-from+1;
		}
		
		public boolean isSmaller(Range range) {
			if (to < range.from)
				return true;
			return false;
		}
		
		public boolean isLarger(Range range) {
			if (from > range.to)
				return true;
			return false;
		}
		
		public Range getOverlap(Range range) {
			int newfrom = Math.max(from, range.from);
			int newto = Math.min(to, range.to);
			if (newto<newfrom)
				return null;
			
			return new Range(newfrom, newto);			
		}
		
		public List<Range> getRemoved(Range range) {
			// get the overlap between the two
			Range overlap = getOverlap(range);
			if (overlap==null)
				return null;
			
			List<Range> removed = new ArrayList<>();
			
			// remove the overlap from this.range
			if (from < overlap.from) {
				removed.add(new Range(from, overlap.from-1));
			}
			
			if (to > overlap.to) {
				removed.add(new Range(overlap.to+1, to));
			}
				
			return removed;
		}

		
		public String toString() {
			return from +"-"+to;
			
		}
	}

	
}
