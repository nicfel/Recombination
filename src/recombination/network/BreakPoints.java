package recombination.network;

import java.util.ArrayList;
import java.util.List;

/**
 * Keep track of recombination break points in the network
 * @author Nicola Müller
 */
public class BreakPoints {
	
	public List<Integer> breakPoints;	
	private List<Integer> leftBreakPoints;
	private List<Integer> rightBreakPoints;
	
	private int totalLength;
	
	public BreakPoints(int totalLength) { 
		breakPoints = new ArrayList<>();
		breakPoints.add(0);
		breakPoints.add(totalLength);
		this.totalLength = totalLength;
	}

	
	public BreakPoints(int start, int end, int totalLength) { 
		breakPoints = new ArrayList<>();
		breakPoints.add(start);
		breakPoints.add(end);
		this.totalLength = totalLength;
	}
	
	public BreakPoints(List<Integer> breakPoints, int totalLength) { 
		this.breakPoints = new ArrayList<>(breakPoints);
		this.totalLength = totalLength;
	}

	/**
	 * combines breakpoints with another breakpoint
	 * @param breakPoint
	 */
	public void combine(BreakPoints breakPoints){
		int i=0,j=0;
		while (i<this.breakPoints.size() || j < breakPoints.breakPoints.size()) {
			
		}
	}
	
	/**
	 * computes new break points list based on a new break point introduced
	 * @param breakpoint
	 */
	public void computeLeftAndRight(int breakpoint) {
		List<Integer> leftBreakPoints = new ArrayList<>();
		List<Integer> rightBreakPoints = new ArrayList<>();
		int i=0;
		while (leftBreakPoints.get(i) < breakpoint) {
			leftBreakPoints.add(breakPoints.get(i));
			i++;
		}
		
		// if the new breakpoints falls within an area that is carried by the edge.
		if (leftBreakPoints.size() % 2 == 1) {
			leftBreakPoints.add(breakpoint);
			rightBreakPoints.add(breakpoint);
		}
		
		while (i < breakPoints.size()) {
			rightBreakPoints.add(breakPoints.get(i));
			i++;
		}			
	}

	/**
	 * gets all the info from left of a new breakpoint 
	 * @param position
	 * @return
	 */
	public List<Integer> getLeft() {
		return leftBreakPoints;
	}
	
	/**
	 * gets all the info from the right of a new breakpoint 
	 * @param position
	 * @return
	 */
	public List<Integer> getRight(int breakpoint) {
		return rightBreakPoints;
	}

	/**
	 * gets the difference between the first and last event
	 * @return
	 */
	public double getLength() {
		return (breakPoints.get(breakPoints.size()-1)-breakPoints.get(0))/((double) totalLength);
	}
	
	/**
	 * gets the total amount of genetic material encompassed
	 * by this break point as defined by 
	 * @return
	 */
	public double getGeneticLength() {
		int l = 0;
		for (int i = 0; i < breakPoints.size(); i=i+2)
			l += breakPoints.get(i+1)-breakPoints.get(i);
		return l/((double) totalLength);
	}

	public BreakPoints copy() {
		BreakPoints newBreakPoints = new BreakPoints(breakPoints, totalLength);
		return newBreakPoints;
	}
	
}
