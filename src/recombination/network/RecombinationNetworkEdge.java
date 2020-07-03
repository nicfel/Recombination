package recombination.network;

import java.util.*;

public class RecombinationNetworkEdge {

    public RecombinationNetworkNode parentNode, childNode;
    public BreakPoints breakPoints;

    public RecombinationNetworkEdge() { }

    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		BreakPoints breakPoints) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = breakPoints;
    }
    
    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		int totalLength) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = breakPoints;
    }


    public double getRecombinationObsProb(double p) {
        // There are always two reassortment configurations that
        // produce an unobserved reassortment: 1111 and 0000
        // (assuming 4 segs on lineage)
        return 1.0 - breakPoints.getLength();
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

        RecombinationNetworkEdge edgeCopy = new RecombinationNetworkEdge(null, null, breakPoints.copy());
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

}
