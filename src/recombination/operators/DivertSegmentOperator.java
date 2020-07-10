package recombination.operators;

import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetworkEdge;

import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class DivertSegmentOperator extends EmptyEdgesRecombinationNetworkOperator {

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
        
        BreakPoints breaksToDivert = getRandomUnconditionedSubset(sourceEdge.breakPoints);
        logHR -= getLogUnconditionedSubsetProb(sourceEdge.breakPoints);
        
        if (breaksToDivert.isEmpty())
        	return Double.NEGATIVE_INFINITY;

        network.startEditing(this);
        

        logHR -= addSegmentsToAncestors(destEdge, segsToDivert);
        logHR += removeSegmentsFromAncestors(sourceEdge, segsToDivert);

        logHR += getLogUnconditionedSubsetProb(destEdge.hasSegments);

        int reverseSourceEdgeCount = (int)(network.getEdges().stream()
                .filter(e -> e.childNode.isReassortment())
                .filter(e -> e.hasSegments.cardinality()>0)
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
    double removeSegmentsFromAncestors(RecombinationNetworkEdge edge, BreakPoints rangeToRemove) {
        double logP = 0.0;

        rangeToRemove = rangeToRemove.copy();
        
        rangeToRemove.and(edge.breakPoints);

        if (rangeToRemove.isEmpty())
            return logP;

        edge.breakPoints.andNot(rangeToRemove);

        if (edge.isRootEdge())
            return logP;

        if (edge.parentNode.isRecombination()) {

            logP += Math.log(0.5)*rangeToRemove.cardinality();

            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), rangeToRemove);
            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(1), rangeToRemove);

        } else {

        	rangeToRemove.andNot(getSisterEdge(edge).breakPoints);
            logP += removeSegmentsFromAncestors(edge.parentNode.getParentEdges().get(0), rangeToRemove);

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
    double addSegmentsToAncestors(RecombinationNetworkEdge edge, BreakPoints rangeToAdd) {
        double logP = 0.0;

        rangeToAdd = rangeToAdd.copy();
        rangeToAdd.andNot(edge.breakPoints);

        if (rangeToAdd.isEmpty())
            return logP;

        edge.breakPoints.or(rangeToAdd);

        if (edge.isRootEdge())
            return logP;

        if (edge.parentNode.isRecombination()) {

            BitSet rangeToAddLeft = new BitSet();
            BitSet rangeToAddRight = new BitSet();

            for (int segIdx=segsToAdd.nextSetBit(0); segIdx != -1;
                    segIdx=segsToAdd.nextSetBit(segIdx+1)) {
                if (Randomizer.nextBoolean())
                    segsToAddLeft.set(segIdx);
                else
                    segsToAddRight.set(segIdx);

                logP += Math.log(0.5);
            }

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), segsToAddLeft);
            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(1), segsToAddRight);

        } else {

            logP += addSegmentsToAncestors(edge.parentNode.getParentEdges().get(0), breakPoints);
        }

        return logP;
    }
    
    protected BitSet getRandomUnconditionedSubset(BreakPoints sourceSegments) {
    	
    	
    	
    	
        BitSet destSegments = new BitSet();

        destSegments.clear();

        for (int segIdx = sourceSegments.nextSetBit(0); segIdx != -1;
             segIdx = sourceSegments.nextSetBit(segIdx + 1)) {

            if (Randomizer.nextBoolean())
                destSegments.set(segIdx);
        }

        return destSegments;
    }

    protected double getLogUnconditionedSubsetProb(BreakPoints sourceSegments) {
        return sourceSegments.cardinality()*Math.log(0.5);
    }


}
