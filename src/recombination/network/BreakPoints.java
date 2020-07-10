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
	
	final static RangeComparator rc = new RangeComparator();

		
	public BreakPoints(int totalLength) { 
		breakPoints = new ArrayList<>();
		breakPoints.add(new Range(0, totalLength-1));
	}

	
	public BreakPoints(int start, int end) { 
		breakPoints = new ArrayList<>();
		breakPoints.add(new Range(start, end));
	}
	
	public BreakPoints() { }
	
	public void init(List<Integer> breakPoints) { 
		this.breakPoints = new ArrayList<>();
		for (int i = 0; i < breakPoints.size(); i=i+2)
			this.breakPoints.add(new Range(breakPoints.get(i), breakPoints.get(i+1)));
	}
	
	public BreakPoints(List<Range> breakPoints) { 
		this.breakPoints = new ArrayList<>(breakPoints);
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
		return breakPoints.get(breakPoints.size()-1).to-breakPoints.get(0).from;
	}
	
	/**
	 * gets the total amount of genetic material encompassed
	 * by this break point as defined by 
	 * @return
	 */
	public double getGeneticLength() {
		int l = 0;
		for (int i = 0; i < breakPoints.size(); i++)
			l += breakPoints.get(i).size();
		return l;
	}
	
	public boolean isEmpty() {
		return breakPoints.size()==0;
	}
	

	public boolean withinLimits(int breakpoint) {
		if (this.breakPoints.get(0).from<=breakpoint && this.breakPoints.get(this.breakPoints.size()-1).to>breakpoint)
			return true;
		
		return false;
		
	} 
	
	public BreakPoints copy() {
		BreakPoints newBreakPoints = new BreakPoints(breakPoints);
		return newBreakPoints;
	}
		
	public String toString() {
		String val="";
		for (int i = 0; i < this.breakPoints.size(); i++) {
			val = val + "," + this.breakPoints.get(i).toString();
		}
				
		return val.substring(1);
		
	}
	
	
	
	public class Range{
		int from,to;
	
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
		
		public String toString() {
			return from +":"+to;
			
		}


		
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

	public void and(BreakPoints breakPoints) {
		// TODO Auto-generated method stub
		
	}


	public void andNot(BreakPoints breakPoints) {
		// TODO Auto-generated method stub
		
	}


	public void or(BreakPoints breakPoints) {
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
	}


	
}
