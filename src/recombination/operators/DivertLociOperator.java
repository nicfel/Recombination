package recombination.operators;

import beast.core.Input;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class DivertLociOperator extends EmptyEdgesRecombinationNetworkOperator {

    public Input<Double> scaleFactorInput = new Input<>(
            "scaleFactor",
            "Scale factor tuning parameter.",
            1.0);

    double lambdaDiversion;
    public int totalLength;
    
    boolean stop = false;
    
    public void initAndValidate() {
    	lambdaDiversion = scaleFactorInput.get();    	
    	super.initAndValidate();
    	totalLength = network.totalLength;
    }
	
    @Override
    public double networkProposal() {    	    	
        double logHR = 0.0;
        
        List<RecombinationNetworkNode> sourceNodes = network.getNodes().stream()
                .filter(e -> e.isRecombination())
                .filter(e -> !e.getParentEdges().get(0).breakPoints.isEmpty())
                .filter(e -> !e.getParentEdges().get(1).breakPoints.isEmpty())
                .collect(Collectors.toList());

        if (sourceNodes.isEmpty())
            return Double.NEGATIVE_INFINITY;

        logHR -= Math.log(1.0/sourceNodes.size());

        network.startEditing(this);
        
    	int newBreakPoint = Randomizer.nextInt(totalLength-1)+1;
       
        // get both edges
        RecombinationNetworkNode node = sourceNodes.get(Randomizer.nextInt(sourceNodes.size()));
        RecombinationNetworkEdge edge1,edge2;
        if (node.getParentEdges().get(0).breakPoints.getMax()>node.getParentEdges().get(1).breakPoints.getMax()) {
        	edge1 = node.getParentEdges().get(1);
        	edge2 = node.getParentEdges().get(0);
        }else {
        	edge1 = node.getParentEdges().get(0);
        	edge2 = node.getParentEdges().get(1);       	
        }
              

        edge1.makeDirty(Tree.IS_DIRTY);
        edge2.makeDirty(Tree.IS_DIRTY);

    	edge2.passingRange = new BreakPoints(newBreakPoint, totalLength-1);
    	edge1.passingRange = new BreakPoints(0,newBreakPoint-1);

    	if (edge2.breakPoints.getMin()>newBreakPoint) {
        	BreakPoints rangeToDivert = new BreakPoints(newBreakPoint, edge2.breakPoints.getMin()-1);
        	rangeToDivert.and(edge1.breakPoints);
	        logHR -= addLociToAncestors(edge2, rangeToDivert);
	        logHR += removeLociFromAncestors(edge1, rangeToDivert);
        }else {
        	BreakPoints rangeToDivert = new BreakPoints(edge2.breakPoints.getMin(), newBreakPoint-1);
        	rangeToDivert.and(edge2.breakPoints);
	        logHR -= addLociToAncestors(edge1, rangeToDivert);
	        logHR += removeLociFromAncestors(edge2, rangeToDivert);
        }
        
        int reverseSourceEdgeCount = (int)(network.getNodes().stream()
                .filter(e -> e.isRecombination())
                .filter(e -> !e.getParentEdges().get(0).breakPoints.isEmpty())
                .filter(e -> !e.getParentEdges().get(1).breakPoints.isEmpty())
                .count());

        logHR += Math.log(1.0/reverseSourceEdgeCount);
        return logHR;
    }


    /**
     * Remove segments from this edge and ancestors.
     *
     * @param edge edge at which to start removal
     * @param segsToRemove segments to remove from edge and ancestors
     * @return log probability of reverse operation
     */
    public double removeLociFromAncestors(RecombinationNetworkEdge edge, BreakPoints rangeToRemove) {
        double logP = 0.0;
        

        rangeToRemove = rangeToRemove.copy();
        
        rangeToRemove.and(edge.breakPoints);
        
        if (rangeToRemove.isEmpty())
            return logP;
        
        if (edge.isRootEdge())
            return logP;

        edge.childNode.dirtyBreakPoints = new BreakPoints(0,totalLength-1);

        edge.breakPoints.andNot(rangeToRemove);                      
        
        edge.makeDirty(Tree.IS_FILTHY); 

        if (edge.parentNode.isRecombination()) {        	
            logP += removeLociFromAncestors(edge.parentNode.getParentEdges().get(0), rangeToRemove);
            logP += removeLociFromAncestors(edge.parentNode.getParentEdges().get(1), rangeToRemove);
        } else {
        	rangeToRemove.andNot(getSisterEdge(edge).breakPoints);
            logP += removeLociFromAncestors(edge.parentNode.getParentEdges().get(0), rangeToRemove);
        }

        return logP;
    }

    /**
     * Add segments to this edge and ancestors.
     *
     * @param edge edge at which to start addition
     * @param segsToAdd segments to add to the edge and ancestors
     * @return log probability of operation
     */
    public double addLociToAncestors(RecombinationNetworkEdge edge, BreakPoints rangeToAdd) {
        double logP = 0.0;

        rangeToAdd = rangeToAdd.copy();
               
        if (rangeToAdd.isEmpty())
            return logP;

        rangeToAdd.andNot(edge.breakPoints);

        if (rangeToAdd.isEmpty())
            return logP;        


        edge.breakPoints.or(rangeToAdd);
        
        edge.makeDirty(Tree.IS_FILTHY); 


        if (edge.isRootEdge())
            return logP;        
        
        if (edge.parentNode.isRecombination()) {        	
            BreakPoints rangeToAddLeft = rangeToAdd.copy();
            BreakPoints rangeToAddRight = rangeToAdd.copy();
                        
            rangeToAddLeft.and(edge.parentNode.getParentEdges().get(0).passingRange);
            rangeToAddRight.and(edge.parentNode.getParentEdges().get(1).passingRange);
                        
            logP += addLociToAncestors(edge.parentNode.getParentEdges().get(0), rangeToAddLeft);
            logP += addLociToAncestors(edge.parentNode.getParentEdges().get(1), rangeToAddRight);
        } else {
            logP += addLociToAncestors(edge.parentNode.getParentEdges().get(0), rangeToAdd);
        }
        return logP;
    }       
}
