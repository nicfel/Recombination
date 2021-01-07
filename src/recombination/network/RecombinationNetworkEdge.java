package recombination.network;

import java.util.*;

import beast.evolution.tree.Tree;
import recombination.util.NodeEdgeID;

public class RecombinationNetworkEdge {

    public RecombinationNetworkNode parentNode, childNode;
    public BreakPoints breakPoints;
    public BreakPoints passingRange;
    public BreakPoints carryingRange;
    
    
    public Integer ID;
    
    /**
     * status of this node after an operation is performed on the state *
     */
    int isDirty = Tree.IS_CLEAN;

    public boolean visited;
   
    public RecombinationNetworkEdge(NodeEdgeID nodeEdgeIDs) { 
        ID = nodeEdgeIDs.getNewEdgeID();
        isDirty = Tree.IS_DIRTY;
    }
    
    public RecombinationNetworkEdge(Integer ID) { 
        this.ID = ID;
        isDirty = Tree.IS_DIRTY;
    }


    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		BreakPoints breakPoints, BreakPoints passingRange, NodeEdgeID nodeEdgeIDs) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = breakPoints;
        this.passingRange = passingRange;
        ID = nodeEdgeIDs.getNewEdgeID();
        isDirty = Tree.IS_DIRTY;
    }
    
    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		BreakPoints breakPoints, BreakPoints passingRange, int id) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = breakPoints;
        this.passingRange = passingRange;
        ID = id;
        isDirty = Tree.IS_DIRTY;
    }

    

       
//    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
//    		int totalLength) {
//        this.parentNode = parentNode;
//        this.childNode = childNode;
//        this.breakPoints = new BreakPoints(totalLength);
//        isDirty = Tree.IS_FILTHY;
//   }
    
//    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
//    		List<Integer> breakPointsList) {
//        this.parentNode = parentNode;
//        this.childNode = childNode;
//        this.breakPoints = new BreakPoints();
//        this.breakPoints.init(breakPointsList);
//        isDirty = Tree.IS_FILTHY;
//    }

    public double getRecombinationLength() {
        // There are always two reassortment configurations that
        // produce an unobserved reassortment: 1111 and 0000
        // (assuming 4 segs on lineage)
        return breakPoints.getLength()-1;
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
        	edgeCopy = new RecombinationNetworkEdge(null, null, breakPoints.copy(), passingRange.copy(), ID);
        else
        	edgeCopy = new RecombinationNetworkEdge(null, null, breakPoints.copy(), null, ID);
        

       
        RecombinationNetworkNode childNodeCopy;
        boolean traverse = true;
        if (seenNodes.containsKey(childNode)) {
            childNodeCopy = seenNodes.get(childNode);
            traverse = false;
        } else {
            childNodeCopy = new RecombinationNetworkNode(childNode.ID);
            
            childNodeCopy.setHeight(childNode.getHeight());
            childNodeCopy.setTaxonLabel(childNode.getTaxonLabel());
            childNodeCopy.setTaxonIndex(childNode.getTaxonIndex());
            childNodeCopy.setTypeIndex(childNode.typeIndex);
            childNodeCopy.setTypeLabel(childNode.typeLabel);
            childNodeCopy.setTypeLabel(childNode.typeLabel);

            seenNodes.put(childNode, childNodeCopy);
        }

        childNodeCopy.addParentEdge(edgeCopy);
       	edgeCopy.isDirty = isDirty;

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

		
    public int isDirty() {
        return isDirty;
    }

    public void makeDirty(final int dirty) {
    	if (isDirty!=Tree.IS_FILTHY)
    		isDirty = dirty;
    }
    
    public void forceDirty(final int dirty) {
   		isDirty = dirty;
    }




	
}
