package recombination.network;

import java.util.*;

import recombination.network.BreakPoints.Range;

public class RecombinationNetworkEdge {

    public RecombinationNetworkNode parentNode, childNode;
    public BreakPoints breakPoints;
    public BreakPoints passingRange;
    public BreakPoints carryingRange;
    
    // keeps track of the matrices for the likelihood calculations
    public double[] matrixList;

    public RecombinationNetworkEdge() { }

    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		BreakPoints breakPoints, double[] matrixList) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = breakPoints;
        if (matrixList!=null) {
        	this.matrixList = new double[matrixList.length];
        	System.arraycopy(matrixList, 0, this.matrixList, 0, this.matrixList.length);
        }
    }
    
    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		BreakPoints breakPoints, BreakPoints passingRange, double[] matrixList) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = breakPoints;
        this.passingRange = passingRange;
        if (matrixList!=null) {
        	this.matrixList = new double[matrixList.length];
        	System.arraycopy(matrixList, 0, this.matrixList, 0, this.matrixList.length);
        }
    }
   
    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		int totalLength) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = new BreakPoints(totalLength);
   }
    
    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		List<Integer> breakPointsList) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = new BreakPoints();
        this.breakPoints.init(breakPointsList);
    }

    public double getRecombinationLength() {
        // There are always two reassortment configurations that
        // produce an unobserved reassortment: 1111 and 0000
        // (assuming 4 segs on lineage)
        return breakPoints.getLength();
    }

    public double getLength() {
        return parentNode.getHeight() - childNode.getHeight();
    }

    public boolean isRootEdge() {
        return parentNode == null;
    }
    
    public boolean isLeafEdge() {
    	return childNode.isLeaf();
    }

    public RecombinationNetworkEdge getCopy() {
        return getCopy(new HashMap<>());
    }

    public RecombinationNetworkEdge getCopy(Map<RecombinationNetworkNode,RecombinationNetworkNode> seenNodes) {
        RecombinationNetworkEdge edgeCopy;
        if (passingRange!=null)
        	edgeCopy = new RecombinationNetworkEdge(null, null, breakPoints.copy(), passingRange.copy(), matrixList);
        else
        	edgeCopy = new RecombinationNetworkEdge(null, null, breakPoints.copy(), matrixList);
        
        RecombinationNetworkNode childNodeCopy;
        boolean traverse = true;
        if (seenNodes.containsKey(childNode)) {
            childNodeCopy = seenNodes.get(childNode);
            traverse = false;
        } else {
            childNodeCopy = new RecombinationNetworkNode();
            childNodeCopy.setHeight(childNode.getHeight());
            childNodeCopy.setTaxonLabel(childNode.getTaxonLabel());
            childNodeCopy.setTaxonIndex(childNode.getTaxonIndex());
            childNodeCopy.setTypeIndex(childNode.typeIndex);
            childNodeCopy.setTypeLabel(childNode.typeLabel);
            childNodeCopy.setTypeLabel(childNode.typeLabel);
            if (childNode.states!=null) {
            	childNodeCopy.states = new int[childNode.states.length];
            	System.arraycopy(childNode.states, 0,  childNodeCopy.states, 0, childNode.states.length);
            }
            if (childNode.partials!=null) {
            	childNodeCopy.partials = new double[childNode.partials.length];
            	System.arraycopy(childNode.partials, 0,  childNodeCopy.partials, 0, childNode.partials.length);
            }
            seenNodes.put(childNode, childNodeCopy);
        }

        childNodeCopy.addParentEdge(edgeCopy);

        if (traverse) {
            for (RecombinationNetworkEdge childEdge : childNode.getChildEdges()) {
                RecombinationNetworkEdge childEdgeCopy = childEdge.getCopy(seenNodes);
                childNodeCopy.addChildEdge(childEdgeCopy);
            }
        }

        return edgeCopy;
    }
    
    /**
     * set the range of loci that goes left a this recombination node, only if 
     * @param from
     * @param to
     */
    public void setPassingRange(int from, int to) {
    	passingRange = new BreakPoints(from, to);
    }
    
    public BreakPoints getPassingRange() {
    	return passingRange;
    }

	public void setPassingRange(BreakPoints lociToDivert) {
		passingRange = lociToDivert;
	}

	
	public void setMatrix(double[] matrixList) {
		// TODO, allow for multiple
		System.arraycopy(this.matrixList, 0, matrixList,
                0, matrixList.length);
	}
	
	
	public double[] getMatrix() {
		return matrixList;
	}


	
}
