package recombination.operators;

import beast.core.Input;
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
            1000.0);

    double lambdaDiversion;
    
    public void initAndValidate() {
    	lambdaDiversion = scaleFactorInput.get();
    	super.initAndValidate();
    }
	
    @Override
    public double networkProposal() {
    	
        double logHR = 0.0;
        
        List<RecombinationNetworkEdge> sourceEdges = network.getEdges().stream()
                .filter(e -> e.childNode.isRecombination())
                .filter(e -> !e.breakPoints.isEmpty())
                .collect(Collectors.toList());

        if (sourceEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        logHR -= Math.log(1.0/sourceEdges.size());

        RecombinationNetworkEdge sourceEdge = sourceEdges.get(Randomizer.nextInt(sourceEdges.size()));
        RecombinationNetworkEdge destEdge = getSpouseEdge(sourceEdge);

        network.startEditing(this);
                
        BreakPoints rangeToDivert = getRangeToDivert(sourceEdge, destEdge);
        
        if (rangeToDivert==null)
        	return Double.NEGATIVE_INFINITY;
        
        // TODO account for HR contribution
        BreakPoints lociToDivert = rangeToDivert.copy();
        
    	lociToDivert.and(sourceEdge.breakPoints);

        if (lociToDivert.isEmpty())
        	return Double.NEGATIVE_INFINITY;

        

        logHR -= addLociToAncestors(destEdge, lociToDivert);
        logHR += removeLociFromAncestors(sourceEdge, lociToDivert);
        
    	sourceEdge.passingRange.andNot(rangeToDivert);
    	
//        logHR += getLogUnconditionedSubsetProb(destEdge.breakPoints);

        int reverseSourceEdgeCount = (int)(network.getEdges().stream()
                .filter(e -> e.childNode.isRecombination())
                .filter(e -> !e.breakPoints.isEmpty())
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

        edge.breakPoints.andNot(rangeToRemove.copy());
                
        if (edge.isRootEdge())
            return logP;

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
        
    public BreakPoints getRangeToDivert(RecombinationNetworkEdge sourceEdge, RecombinationNetworkEdge destEdge) {
    	
    	BreakPoints rangeToDivert = new BreakPoints();
    	BreakPoints range = sourceEdge.getPassingRange();  
    	    	    	
    	int step = (int) Randomizer.nextPoisson(lambdaDiversion);
    	
    	if (range.breakPoints.get(0).from!=0) {
    		if ((range.breakPoints.get(0).from+step)>range.breakPoints.get(0).to)
    			return null;
    		
    		rangeToDivert = new BreakPoints(range.breakPoints.get(0).from, range.breakPoints.get(0).from+step);
    	}else {
    		if ((range.breakPoints.get(0).to-step)<0)
    			return null;
    		
			rangeToDivert = new BreakPoints(range.breakPoints.get(0).to-step , range.breakPoints.get(0).to);
    	}
    	    	
    	
    	// change the passing Range that goes through the different edges
    	if (destEdge.getPassingRange()!=null) {
    		destEdge.passingRange.or(rangeToDivert);
    	}else {
    		destEdge.setPassingRange(rangeToDivert.copy());
    	}
    	

    	return rangeToDivert;
    }


    public BreakPoints getNewRangeToDivert(RecombinationNetworkEdge sourceEdge, RecombinationNetworkEdge destEdge) {
    	
    	BreakPoints rangeToDivert = new BreakPoints();
    	
    	int newBreakPoint = Randomizer.nextInt(network.totalLength);

    	
    	if (Randomizer.nextBoolean()) {
    		rangeToDivert = new BreakPoints(0,newBreakPoint);
    	}else {
    		rangeToDivert = new BreakPoints(newBreakPoint, network.totalLength-1);
    	}
    	    	    	    	    	
		destEdge.setPassingRange(rangeToDivert.copy());	

    	return rangeToDivert;
    }



    
//    protected BitSet getRandomUnconditionedSubset(BreakPoints sourceSegments) {
//    	
//    	
//    	
//    	
//        BitSet destSegments = new BitSet();
//
//        destSegments.clear();
//
//        for (int segIdx = sourceSegments.nextSetBit(0); segIdx != -1;
//             segIdx = sourceSegments.nextSetBit(segIdx + 1)) {
//
//            if (Randomizer.nextBoolean())
//                destSegments.set(segIdx);
//        }
//
//        return destSegments;
//    }
//
//    protected double getLogUnconditionedSubsetProb(BreakPoints sourceSegments) {
//        return sourceSegments.cardinality()*Math.log(0.5);
//    }
//

}
