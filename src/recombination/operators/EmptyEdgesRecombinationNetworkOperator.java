package recombination.operators;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import beast.core.Input;
import beast.core.Operator;
import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

public abstract class EmptyEdgesRecombinationNetworkOperator extends RecombinationNetworkOperator {

    public Input<Double> emptyAlphaInput = new Input<>("emptyAlpha",
            "Mean of exponential used for choosing root attachment times.",
            0.1);
    
    public Input<Double> lambdaInput = new Input<>("lambda",
            "lambda of the poisson distribution for how many empty edges to add.",
            0.1);
    
    public Input<Boolean> addRemoveEmptyEdgesInput = new Input<>("addRemoveEmptyEdges",
            "adds empty edges before calling the networkproposal and then removes all empty edges at the end again",
            true);

    private double emptyAlpha;
    private double lambda;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        emptyAlpha = emptyAlphaInput.get();
        lambda = lambdaInput.get();
    }

    @Override
    public double proposal() {
//    	System.out.println();
       
    	double logHR = 0.0;
        // Adds empty network edges
        if (addRemoveEmptyEdgesInput.get()){
        	logHR += addEmptyNetworkEdges();
        	
            if (logHR == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;
        }     
                
        
        // calls the operator
        logHR += networkProposal();
        
    	

        
        
        if (logHR == Double.NEGATIVE_INFINITY)
        	return Double.NEGATIVE_INFINITY;


        // removes all the empty network edges in the network again
        if (addRemoveEmptyEdgesInput.get()){
            if (logHR == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;
            
        	logHR += RemoveAllEmptyNetworkSegments();
        }
        

        // case there are empty edges, which can happen when addRemoveEmptyEdges is false
		if (!allEdgesAncestral()){
//			System.out.println("/////////////////");
//			System.out.println(network);
//			System.out.println(logHR);
//			throw new IllegalArgumentException("ancestral issues");
            return Double.NEGATIVE_INFINITY;
		}
		
//        // case there are empty edges, which can happen when addRemoveEmptyEdges is false
//		if (!allRangesMatch()){
//			System.out.println(network);
////			throw new IllegalArgumentException("Range issues");
//
////            return Double.NEGATIVE_INFINITY;
//		}
		
		
        // case there are empty edges, which can happen when addRemoveEmptyEdges is false
		if (!allBreaksDiffer()){
			System.out.println(network);
			throw new IllegalArgumentException("Break error");
//            return Double.NEGATIVE_INFINITY;
		}
		
//        // case there are empty edges, which can happen when addRemoveEmptyEdges is false
//		if (!splitsMakeSense()){
//			System.out.println(network);
//			throw new IllegalArgumentException("Split error");
////            return Double.NEGATIVE_INFINITY;
//		}

		
		if (!coalBPcheck()){
			System.out.println(network);
			throw new IllegalArgumentException("BreakPoints error");
		}
		
//
//    	System.out.println(network);
//    	System.out.println(logHR);
		
        return logHR;
    }
    
    private double addEmptyNetworkEdges(){
    	double logHR = 0.0;
    	
    	// randomly sample the number of edges to add
    	int nrEmptyEdges = (int) Randomizer.nextPoisson(lambda);
    	    	
    	for (int i = 0; i < nrEmptyEdges; i ++){
    		logHR += addEmptyRecombination();
    	}      	
    	
    	logHR -= Math.log(Math.pow(lambda, nrEmptyEdges)) - lambda - Math.log(factorial(nrEmptyEdges));
    	
    	return logHR;
    }
    

    
    double addEmptyRecombination() {
        double logHR = 0.0;

        List<RecombinationNetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        // add empty reassortment edges to non empty edges
        List<RecombinationNetworkEdge> possibleSourceEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .collect(Collectors.toList());

        RecombinationNetworkEdge sourceEdge = possibleSourceEdges.get(Randomizer.nextInt(possibleSourceEdges.size()));
        double sourceTime = Randomizer.nextDouble()*sourceEdge.getLength() + sourceEdge.childNode.getHeight();

        logHR -= Math.log(1.0/(double)possibleSourceEdges.size())
                + Math.log(1.0/sourceEdge.getLength());
               
        List<RecombinationNetworkEdge> possibleDestEdges = networkEdges.stream()
                .collect(Collectors.toList());

        RecombinationNetworkEdge destEdge = possibleDestEdges.get(Randomizer.nextInt(possibleDestEdges.size()));
    	// works
        logHR -= Math.log(1.0/possibleDestEdges.size());

        
        if (!destEdge.isRootEdge() && destEdge.parentNode.getHeight() < sourceTime)
            return Double.NEGATIVE_INFINITY;

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        double destTime;
        if (destEdge.isRootEdge()) {
        	// works
            destTime = minDestTime + Randomizer.nextExponential(1.0/emptyAlpha);
            logHR -= -(1.0/emptyAlpha)*(destTime-minDestTime) + Math.log(1.0/emptyAlpha);
        } else {
            destTime = Randomizer.nextDouble()*(destEdge.parentNode.getHeight()-minDestTime) + minDestTime;
            logHR -= Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));
        }
        
      

        // Create new reassortment edge
        logHR += addEmptyRecombinationEdge(sourceEdge, sourceTime, destEdge, destTime);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
                
        
        // HR contribution for reverse move
        int nRemovableEdges = (int) network.getEdges().stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.breakPoints.isEmpty())
                .filter(e -> e.childNode.isRecombination())
                .filter(e -> e.parentNode.isCoalescence())
                .count();

        // works
        logHR += Math.log(1.0/nRemovableEdges);       
        
        if(!networkTerminatesAtMRCA())
        	return Double.NaN;


        return logHR;
    }
    
    // Only adds the recombination edge, but does not diverge segments
    double addEmptyRecombinationEdge(RecombinationNetworkEdge sourceEdge, double sourceTime,
    		RecombinationNetworkEdge destEdge, double destTime) {

		double logHR = 0.0;
		
		network.startEditing(this);
		
		RecombinationNetworkNode sourceNode = new RecombinationNetworkNode(network.nodeEdgeIDs);
		sourceNode.setHeight(sourceTime);
		
		RecombinationNetworkNode oldSourceEdgeParent = sourceEdge.parentNode;
		oldSourceEdgeParent.removeChildEdge(sourceEdge);
		sourceNode.addChildEdge(sourceEdge);
		
		RecombinationNetworkEdge newEdge1 = new RecombinationNetworkEdge(network.nodeEdgeIDs);
		sourceNode.addParentEdge(newEdge1);
		oldSourceEdgeParent.addChildEdge(newEdge1);
		
		newEdge1.breakPoints = sourceEdge.breakPoints.copy();
		newEdge1.setPassingRange(0,network.totalLength-1);

			
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
		
	
		

		

		
		return logHR;
	}
    
    public double RemoveRandomEmptyNetworkSegments() {
    	double logHR = 0.0;
    	
        List<RecombinationNetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        List<RecombinationNetworkEdge> removableEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.isRecombination())
                .filter(e -> e.breakPoints.isEmpty())
                .filter(e -> e.parentNode.isCoalescence())
                .collect(Collectors.toList());
        
        if (removableEdges.size()>0){
	        logHR -= Math.log(1.0/(removableEdges.size()));
	        int edgeInd = Randomizer.nextInt(removableEdges.size());            
	    	logHR += removeEmptyRecombinationEdge(removableEdges.get(edgeInd));
	    	networkEdges = new ArrayList<>(network.getEdges());
    	}
	    	
        if (logHR == Double.NEGATIVE_INFINITY)
        	return Double.NEGATIVE_INFINITY;

        
        if(!networkTerminatesAtMRCA())
        	return Double.NEGATIVE_INFINITY;
        
        return logHR;
    }

    
    public double RemoveAllEmptyNetworkSegments() {
    	double logHR = 0.0;
    	
        List<RecombinationNetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        List<RecombinationNetworkEdge> removableEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.isRecombination())
                .filter(e -> e.breakPoints.isEmpty())
                .filter(e -> e.parentNode.isCoalescence())
                .collect(Collectors.toList());
        
        int nrRemoved = 0;
        
        while (removableEdges.size()>0){
            nrRemoved++;

            logHR -= Math.log(1.0/(removableEdges.size()));
            int edgeInd = Randomizer.nextInt(removableEdges.size());     
            
        	logHR += removeEmptyRecombinationEdge(removableEdges.get(edgeInd));
        	
        	networkEdges = new ArrayList<>(network.getEdges());
        	
        	// In case and invalid network has been created
            if (logHR == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;

            
            removableEdges = networkEdges.stream()
                    .filter(e -> !e.isRootEdge())
                    .filter(e -> e.childNode.isRecombination())
                    .filter(e -> e.breakPoints.isEmpty())
                    .filter(e -> e.parentNode.isCoalescence())
                    .collect(Collectors.toList());            
        } 
        
        // probability of adding n empty edges in reverse move
        logHR += Math.log(Math.pow(lambda, nrRemoved)) -lambda -  Math.log(factorial(nrRemoved));
        
        if (!allEdgesAncestral()){
        	System.out.println(network);
        	throw new IllegalArgumentException("still has empty segments, should not happen ever!");        	
        }
        
        if(!networkTerminatesAtMRCA())
        	return Double.NEGATIVE_INFINITY;
        
        return logHR;
    }
    
    double removeEmptyRecombinationEdge(RecombinationNetworkEdge edgeToRemove) {
        double logHR = 0.0;
        

        double sourceTime = edgeToRemove.childNode.getHeight();
        RecombinationNetworkEdge sourceEdge = edgeToRemove.childNode.getChildEdges().get(0);
        RecombinationNetworkEdge destEdge = getSisterEdge(edgeToRemove);
        if (destEdge.childNode == edgeToRemove.childNode)
            destEdge = sourceEdge;
        double destTime = edgeToRemove.parentNode.getHeight();
        

        // Remove reassortment edge
        logHR += actuallyRemoveEmptyRecombinationEdge(edgeToRemove);

        if (logHR == Double.NEGATIVE_INFINITY)
            return Double.NEGATIVE_INFINITY;
        
        // HR contribution for reverse move
        Set<RecombinationNetworkEdge> finalNetworkEdges = network.getEdges();

        int nPossibleSourceEdges = (int)finalNetworkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .count();
        
        
        // change to -1 because after the move, 1 less can be a source edge
        logHR += Math.log(1.0/(double) nPossibleSourceEdges )
                + Math.log(1.0/sourceEdge.getLength());        

        List<RecombinationNetworkEdge> possibleDestEdges = finalNetworkEdges.stream()
                .collect(Collectors.toList());

        // works
        logHR += Math.log(1.0/finalNetworkEdges.size());
    

        double minDestTime = Math.max(destEdge.childNode.getHeight(), sourceTime);

        if (destEdge.isRootEdge()) {
        	// works
            logHR += -(1.0/emptyAlpha)*(destTime-minDestTime) + Math.log(1.0/emptyAlpha);
        } else {
            logHR += Math.log(1.0/(destEdge.parentNode.getHeight()-minDestTime));
            
        }

        return logHR;
    }


    double actuallyRemoveEmptyRecombinationEdge(RecombinationNetworkEdge edgeToRemove) {

        double logHR = 0.0;

        network.startEditing(this);

        RecombinationNetworkNode nodeToRemove = edgeToRemove.childNode;
        RecombinationNetworkEdge edgeToRemoveSpouse = getSpouseEdge(edgeToRemove);
        RecombinationNetworkNode edgeToRemoveSpouseParent = edgeToRemoveSpouse.parentNode;

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

        if (!networkTerminatesAtMRCA()){
            return Double.NEGATIVE_INFINITY;
        }


        return logHR;
    }
 


    protected boolean networkTerminatesAtMRCA() {
        List<RecombinationNetworkNode> sortedNodes = new ArrayList<>(network.getNodes());
        sortedNodes.sort(Comparator.comparingDouble(RecombinationNetworkNode::getHeight));
        List<RecombinationNetworkNode> sampleNodes = sortedNodes.stream().filter(RecombinationNetworkNode::isLeaf).collect(Collectors.toList());
        double maxSampleHeight = sampleNodes.get(sampleNodes.size()-1).getHeight();

        int lineages = 0;
        for (RecombinationNetworkNode node : sortedNodes) {
            switch(node.getChildEdges().size()) {
                case 2:
                    // Coalescence

                    lineages -= 1;
                    break;

                case 1:
                    // Recombination

                    if (lineages < 2 && node.getHeight() > maxSampleHeight)
                        return false;

                    lineages += 1;
                    break;

                case 0:
                    // Sample

                    lineages += 1;
                    break;
            }
        }

        return true;
    }
    
    /**
     * Check that each edge is ancestral to at least one segment.
     *
     * @return true if all edges are ancestral.
     */
    public boolean allEdgesAncestral() {
        Set<RecombinationNetworkNode> nodeList = networkInput.get().getNodes();
        for (RecombinationNetworkNode node : nodeList) {
            for (RecombinationNetworkEdge parentEdge : node.getParentEdges()) {
                if (parentEdge.breakPoints.isEmpty()) {
                    return false;
                }
            }
        }

        return true;
    }
    
//    /** 
//     * sanity checks whether the passingRange and loci on branches are compatible
//     * @return
//     */
//    public boolean allRangesMatch() {
//        Set<RecombinationNetworkNode> nodeList = networkInput.get().getNodes();
//        for (RecombinationNetworkNode node : nodeList) {
//            for (RecombinationNetworkEdge parentEdge : node.getParentEdges()) {
//                if (parentEdge.childNode.isRecombination()) {
//                    BreakPoints lc = parentEdge.breakPoints.copy();
//                    lc.andNot(parentEdge.passingRange);
//                    if (!lc.isEmpty()) {
//                    	System.out.println(parentEdge.breakPoints);
//                    	System.out.println(parentEdge.passingRange);
//                    	System.out.println(lc);
//                    	System.out.println(parentEdge.childNode.getHeight());
//                    	return false;
//                    }
//                }
//            }
//        }
//        return true;
//
//    }
    
    
    /** 
     * sanity checks that no breakpoints overlap
     * @return
     */
    public boolean allBreaksDiffer() {
        Set<RecombinationNetworkNode> nodeList = networkInput.get().getNodes();
        for (RecombinationNetworkNode node : nodeList) {
            for (RecombinationNetworkEdge parentEdge : node.getParentEdges()) {
            	for (int i = 0; i < parentEdge.breakPoints.breakPoints.size()-1;i++) {
            		if (parentEdge.breakPoints.breakPoints.get(i).getOverlap(parentEdge.breakPoints.breakPoints.get(i+1))!=null) {
//            			System.out.println(parentEdge.childNode.getHeight());
            			return false;
            		}
            	}

            }
        }
        return true;

    }
    
    /** 
     * sanity checks that no breakpoints overlap
     * @return
     */
    public boolean coalBPcheck() {
        List<RecombinationNetworkNode> nodeList = network.getNodes().stream()
                .filter(e -> e.isCoalescence())
                .collect(Collectors.toList());
        
        for (RecombinationNetworkNode node : nodeList) {
        	BreakPoints cp = node.getChildEdges().get(0).breakPoints.copy();
        	cp.or(node.getChildEdges().get(1).breakPoints);
        	
        	if (!cp.equals(node.getParentEdges().get(0).breakPoints)) {
        		System.out.println(node.getHeight());
        		return false;
        	}
        	
        }
        return true;

    }
    
    public boolean crecombBPcheck() {
        List<RecombinationNetworkNode> nodeList = network.getNodes().stream()
                .filter(e -> e.isRecombination())
                .collect(Collectors.toList());
        
        for (RecombinationNetworkNode node : nodeList) {
        	BreakPoints cp = node.getParentEdges().get(0).breakPoints.copy();
        	cp.or(node.getParentEdges().get(1).breakPoints);
        	
        	if (!cp.equals(node.getChildEdges().get(0).breakPoints)) {
        		System.out.println(node.getHeight());
        		return false;
        	}
        	
        }
        return true;

    }


    
//    public boolean splitsMakeSense() {
//        Set<RecombinationNetworkNode> nodeList = networkInput.get().getNodes();
//        for (RecombinationNetworkNode node : nodeList) {
//            for (RecombinationNetworkEdge parentEdge : node.getParentEdges()) {
//            	if (parentEdge.childNode.isRecombination()) {
//            		if (parentEdge.passingRange.getMin()!=0 && parentEdge.passingRange.getMax()!=network.totalLength-1) {
//            			return false;
//            		}
//            	}
//
//            }
//        }
//        return true;
//
//    }




    private int factorial(int k){
    	int f = 1;
    	for (int i = 2; i <= k; i++)
    		f*=k;
    	return f;
    }




    
}
