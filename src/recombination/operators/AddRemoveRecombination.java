package recombination.operators;

import beast.core.Input;
import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.util.*;
import java.util.stream.Collectors;

public class AddRemoveRecombination extends DivertLociOperator {

    public Input<Double> alphaInput = new Input<>("alpha",
            "Mean of exponential used for choosing root attachment times.",
            Input.Validate.REQUIRED);

    private double alpha;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        alpha = alphaInput.get();
    }

    @Override
    public double networkProposal() {
        double logHR;
        if (Randomizer.nextBoolean()) {
            logHR = addRecombination();
        }else {
            logHR = removeRecombination();
        }

        return logHR;
    }

    double addRecombination() {
        double logHR = 0.0;

        List<RecombinationNetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        List<RecombinationNetworkEdge> possibleSourceEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.breakPoints.getLength()>1)
                .collect(Collectors.toList());
        
        RecombinationNetworkEdge sourceEdge = possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size()));
        double sourceTime = Randomizer.nextDouble() * sourceEdge.getLength() + sourceEdge.childNode.getHeight();

        logHR -= Math.log(1.0/(double)possibleSourceEdges.size())
                + Math.log(1.0/sourceEdge.getLength());

        RecombinationNetworkEdge destEdge = networkEdges.get(Randomizer.nextInt(networkEdges.size()));
        logHR -= Math.log(1.0/networkEdges.size());

        if (!destEdge.isRootEdge() && destEdge.parentNode.getHeight() < sourceTime)
            return Double.NEGATIVE_INFINITY;

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        double destTime;
        if (destEdge.isRootEdge()) {
            destTime = minDestTime + Randomizer.nextExponential(1.0/alpha);
            logHR -= -(1.0/alpha)*(destTime-minDestTime) + Math.log(1.0/alpha);
        } else {
            destTime = Randomizer.nextDouble()*(destEdge.parentNode.getHeight()-minDestTime) + minDestTime;
            logHR -= Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));
        }
        

        // Create new reassortment edge
        logHR += addRecombinationEdge(sourceEdge, sourceTime, destEdge, destTime);
       
        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;        

        // HR contribution for reverse move
        int nRemovableEdges = (int) network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> !e.breakPoints.isEmpty())
                .filter(e -> e.childNode.isRecombination())
                .filter(e -> e.parentNode.isCoalescence())
                .count();
        
        logHR += Math.log(1.0/nRemovableEdges);

        return logHR;
    }
    
    public double addRecombinationEdge(RecombinationNetworkEdge sourceEdge, double sourceTime,
    		RecombinationNetworkEdge destEdge, double destTime) {

		double logHR = 0.0;	

		network.startEditing(this);
		
        BreakPoints rangeToDivert = getNewRangeToDivert(sourceEdge);
        BreakPoints lociToDivert = rangeToDivert.copy();
    	lociToDivert.and(sourceEdge.breakPoints);
    	
    	logHR -= Math.log(0.5) + Math.log(1/(sourceEdge.breakPoints.getLength()-1));
    	    	    	    	    	
		RecombinationNetworkNode sourceNode = new RecombinationNetworkNode(network.nodeEdgeIDs);
		sourceNode.setHeight(sourceTime);		
		
		RecombinationNetworkNode oldSourceEdgeParent = sourceEdge.parentNode;
		oldSourceEdgeParent.removeChildEdge(sourceEdge);
		sourceNode.addChildEdge(sourceEdge);		
		
		RecombinationNetworkEdge newEdge1 = new RecombinationNetworkEdge(network.nodeEdgeIDs);
		sourceNode.addParentEdge(newEdge1);
		oldSourceEdgeParent.addChildEdge(newEdge1);
		
		newEdge1.breakPoints = sourceEdge.breakPoints.copy();
			
		if (destEdge == sourceEdge)
			destEdge = newEdge1;
		
		RecombinationNetworkNode destNode = new RecombinationNetworkNode(network.nodeEdgeIDs);
		destNode.setHeight(destTime);
		
		RecombinationNetworkNode oldDestEdgeParent = destEdge.parentNode;
		if (oldDestEdgeParent != null) {
			oldDestEdgeParent.removeChildEdge(destEdge);
		}
		
		destNode.addChildEdge(destEdge);
		
		RecombinationNetworkEdge newEdge2 = new RecombinationNetworkEdge(network.nodeEdgeIDs);
		destNode.addParentEdge(newEdge2);
		
		if (oldDestEdgeParent == null) {
			network.setRootEdge(newEdge2);
		} else {
			oldDestEdgeParent.addChildEdge(newEdge2);
		}
		
		newEdge2.breakPoints = destEdge.breakPoints.copy();
		
		RecombinationNetworkEdge reassortmentEdge = new RecombinationNetworkEdge(network.nodeEdgeIDs);
		sourceNode.addParentEdge(reassortmentEdge);
		destNode.addChildEdge(reassortmentEdge);
		reassortmentEdge.breakPoints = new BreakPoints();	
        reassortmentEdge.passingRange = rangeToDivert.copy();  
        
        logHR -= addLociToAncestors(reassortmentEdge, lociToDivert);
        logHR += removeLociFromAncestors(newEdge1, lociToDivert);
        
        markDirty(reassortmentEdge);
        markDirty(newEdge1);

        
        newEdge1.passingRange = new BreakPoints(0,network.totalLength-1);        
        newEdge1.passingRange.andNot(rangeToDivert);        

        return logHR;
    }

    double removeRecombination() {
        double logHR = 0.0;
        
        List<RecombinationNetworkEdge> removableEdges = network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> !e.breakPoints.isEmpty())
                .filter(e -> e.childNode.isRecombination())
                .filter(e -> e.parentNode.isCoalescence())
                .collect(Collectors.toList());
        
        
        if (removableEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        RecombinationNetworkEdge edgeToRemove = removableEdges.get(Randomizer.nextInt(removableEdges.size()));
        logHR -= Math.log(1.0/(removableEdges.size()));
        

        double sourceTime = edgeToRemove.childNode.getHeight();
        RecombinationNetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
        RecombinationNetworkEdge destEdge = getSisterEdge(edgeToRemove);
        if (destEdge.childNode == edgeToRemove.childNode)
            destEdge = sourceEdge;
        double destTime = edgeToRemove.parentNode.getHeight();
        
        
        // Remove reassortment edge
        logHR += removeRecombinationEdge(edgeToRemove);
        

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;

        // HR contribution for reverse move

        Set<RecombinationNetworkEdge> finalNetworkEdges = network.getEdges();

        int nPossibleSourceEdges = (int)finalNetworkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.breakPoints.getLength()>1)
                .count();

        logHR += Math.log(1.0/(double)nPossibleSourceEdges)
                + Math.log(1.0/sourceEdge.getLength());

        logHR += Math.log(1.0/finalNetworkEdges.size());

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        if (destEdge.isRootEdge()) {
            logHR += -(1.0/alpha)*(destTime-minDestTime) + Math.log(1.0/alpha);
        } else {
            logHR += Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));
        }

        return logHR;
    }


    public double removeRecombinationEdge(RecombinationNetworkEdge edgeToRemove) {
        double logHR = 0.0;

        network.startEditing(this);

        RecombinationNetworkNode nodeToRemove = edgeToRemove.childNode;
        RecombinationNetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
        RecombinationNetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;

        if (edgeToRemoveSpouse.breakPoints.isEmpty())
        	return Double.NEGATIVE_INFINITY;
        
        // Divert segments away from chosen edge
        BreakPoints lociToDivert = edgeToRemove.breakPoints.copy();
        
        logHR -= addLociToAncestors(edgeToRemoveSpouse, lociToDivert);
        logHR += removeLociFromAncestors(edgeToRemove, lociToDivert);
        logHR += Math.log(0.5) + Math.log(1/(edgeToRemoveSpouse.breakPoints.getLength()-1));
        
        markDirty(edgeToRemoveSpouse);
        markDirty(edgeToRemove);
               
        
        // Remove edge and associated nodes
        RecombinationNetworkEdge edgeToExtend = nodeToRemove.getChildEdges().get(0);
        nodeToRemove.removeChildEdge(edgeToExtend);
        nodeToRemove.removeParentEdge(edgeToRemove);
        nodeToRemove.removeParentEdge(edgeToRemoveSpouse);
        edgeToRemoveSpouseParent.removeChildEdge(edgeToRemoveSpouse);
        edgeToRemoveSpouseParent.addChildEdge(edgeToExtend);

        RecombinationNetworkNode secondNodeToRemove = edgeToRemove.parentNode;
        RecombinationNetworkEdge secondEdgeToExtend = getSisterEdge(edgeToRemove);

        secondNodeToRemove.removeChildEdge(secondEdgeToExtend);
        secondNodeToRemove.removeChildEdge(edgeToRemove);

        if (secondNodeToRemove.getParentEdges().get(0).isRootEdge()) {
            network.setRootEdge(secondEdgeToExtend);
        } else {
        	RecombinationNetworkEdge secondNodeToRemoveParentEdge = secondNodeToRemove.getParentEdges().get(0);
        	RecombinationNetworkNode secondNodeToRemoveParent = secondNodeToRemoveParentEdge.parentNode;
            secondNodeToRemoveParent.removeChildEdge(secondNodeToRemoveParentEdge);
            secondNodeToRemove.removeParentEdge(secondNodeToRemoveParentEdge);
            secondNodeToRemoveParent.addChildEdge(secondEdgeToExtend);
        }
    	
        
        if (!networkTerminatesAtMRCA())
            return Double.NEGATIVE_INFINITY;

        return logHR;
    }
    
    
    public BreakPoints getNewRangeToDivert(RecombinationNetworkEdge sourceEdge) {
    	
    	BreakPoints rangeToDivert = new BreakPoints();
    	
    	int newBreakPoint = Randomizer.nextInt((int) sourceEdge.breakPoints.getLength()-1) + 
    			sourceEdge.breakPoints.getMin();
    	
//    	System.out.println(newBreakPoint);
//    	System.out.println(sourceEdge.breakPoints.getLength());
//    	System.out.println(sourceEdge.breakPoints);
//    	System.out.println(newBreakPoint);
    	
    	
    	if (Randomizer.nextBoolean()) {
    		rangeToDivert = new BreakPoints(0,newBreakPoint);
    	}else {
    		rangeToDivert = new BreakPoints(newBreakPoint+1, totalLength-1);
    	}    
    	return rangeToDivert;
    }
    
//    /**
//     * automatic parameter tuning *
//     */
//    @Override
//    public void optimize(final double logAlpha) {
//        if (optimiseInput.get()) {
//            double delta = calcDelta(logAlpha);
//            delta += Math.log(alpha);
//            final double f = Math.exp(delta);
//            if( alpha > 0 ) {
//                final RecombinationNetwork network = networkInput.get();
//                final double h = network.getRootEdge().childNode.getHeight();
//                final double k = Math.log(network.getLeafNodes().size()) / Math.log(2);
//                final double lim = (h / k) * alpha;
//                if( f <= lim ) {
//                	alpha = f;
//                }
//            } else {
//            	alpha = f;
//            }
//        }
//    }

    
}
