/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package recombination.annotator;

import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A rewrite of TreeAnnotator that outputs how often reassortment events happen on trunk branches vs. other branches 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class RecombinationAnnotator {
	/**
	 * performs all the removing things steps
	 * @param network
	 * @param segmentToRemove
	 */
	void pruneLociFromNetwork(RecombinationNetwork network, BreakPoints breakPointsToKeep){
    	
    	BreakPoints breakPointsToRemove = new BreakPoints();
    	if (!breakPointsToKeep.isEmpty()) {
    		breakPointsToRemove = new BreakPoints(network.totalLength); 
    		breakPointsToRemove.andNot(breakPointsToKeep);
    	}
    	
   		removeLoci(network, breakPointsToRemove);
    	// remove all empty edges in the segment
    	removeEmptyNetworkEdge(network);  
	}

	/**
	 * performs all the removing things steps
	 * @param network
	 * @param segmentToRemove
	 */
	void pruneNetwork(RecombinationNetwork network, BreakPoints breakPointsToKeep){
    	// remove all parts of the network that aren't informed by the genetic data
    	removeNonGeneticSegmentEdges(network);
    	
    	BreakPoints breakPointsToRemove = new BreakPoints();
    	if (!breakPointsToKeep.isEmpty()) {
    		breakPointsToRemove = new BreakPoints(network.totalLength); 
    		breakPointsToRemove.andNot(breakPointsToKeep);
    	}
    	
   		removeLoci(network, breakPointsToRemove);

    	// remove all loops
    	removeLoops(network);
    	// remove all empty edges in the segment
    	removeEmptyNetworkEdge(network);      
    	// remove all loops
    	removeLoops(network);

    	// remove all empty edges in the segment
    	removeEmptyNetworkEdge(network);  
	}
	
    /**
     * removes segments from network edges for which there is no
     * genetic information, i.e. segment that are above the segment tree root.
     * @param network
     */
    private void removeNonGeneticSegmentEdges(RecombinationNetwork network){
    	// remove segments from edges if they are "above" the segment tree root
    	removeLociFromEdge(network.getRootEdge(), new BreakPoints(network.totalLength));    	
    }
    
    /**
     * remove segment from edges and child edges until a coalescent event with both children
     * carrying the segment is reached
     * @param edge
     * @param segIdx
     */
    private void removeLociFromEdge(RecombinationNetworkEdge edge, BreakPoints bp){
    	if (bp.isEmpty())
    		return;
    	
    	// only use range represented on edge
    	BreakPoints new_bp = bp.copy();
    	new_bp.and(edge.breakPoints);   			

    	// remove the segment from the edge    	
    	edge.breakPoints.andNot(bp);
    	
    	// get all child segments
    	List<RecombinationNetworkEdge> childEdges = edge.childNode.getChildEdges();
    	if (childEdges.size()==1){
    		removeLociFromEdge(childEdges.get(0), new_bp);
    	}else if (childEdges.size()==2){
    		// check which loci ar on both children
        	BreakPoints overlap = childEdges.get(0).breakPoints.copy();
        	overlap.and(childEdges.get(1).breakPoints);
        	
        	new_bp.andNot(overlap);
        	
    		for (RecombinationNetworkEdge childEdge : childEdges)
    			removeLociFromEdge(childEdge, new_bp);
    		
    	}else{
    		throw new IllegalArgumentException("odd number of child edges");
    	}
    		
    }

    /**
     * removes all reticulation edges that start and end at the same place
     * @param network
     */
    private void removeLoops(RecombinationNetwork network){
    	List<RecombinationNetworkNode> reticulationNodes = network.getNodes().stream()
                .filter(e -> e.isRecombination())
                .filter(e -> e.getParentEdges().get(0).parentNode.equals(e.getParentEdges().get(1).parentNode))
                .filter(e -> !e.getParentEdges().get(1).breakPoints.isEmpty())
                .collect(Collectors.toList());
    	
    	// for each of these, check if the parents are the same node
    	while (!reticulationNodes.isEmpty()){
    		RecombinationNetworkNode node = reticulationNodes.get(0);
			// if this is the case, put all segment from 1 onto 0
			node.getParentEdges().get(0).breakPoints.or(node.getParentEdges().get(1).breakPoints);
			node.getParentEdges().get(1).breakPoints = new BreakPoints();
			
			removeEmptyNetworkEdge(network);
			reticulationNodes = network.getNodes().stream()
	                .filter(e -> e.isRecombination())
	                .filter(e -> e.getParentEdges().get(0).parentNode.ID==e.getParentEdges().get(1).parentNode.ID)
	                .filter(e -> !e.getParentEdges().get(1).breakPoints.isEmpty())
	                .collect(Collectors.toList());			
    	}
    }

    /**
     * removes segment with id segIdx from the network.
     * @param network
     * @param segIdx
     */
    private void removeLoci(RecombinationNetwork network, BreakPoints bp){
    	// get all networkNodes
    	Set<RecombinationNetworkEdge> networkEdges  = network.getEdges();
    	
    	// set carries segment nr segIdx to false for every node
    	for (RecombinationNetworkEdge edge : networkEdges){
    		edge.breakPoints.andNot(bp);
    	}
    }
    
    /**
     * removes all edges from the network that don't carry any segments
     * @param network
     */
    private void removeEmptyNetworkEdge(RecombinationNetwork network){
        List<RecombinationNetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        List<RecombinationNetworkEdge> removableEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.isRecombination())
                .filter(e -> e.breakPoints.isEmpty())
                .filter(e -> e.parentNode.isCoalescence())
                .collect(Collectors.toList());
        
        while (removableEdges.size()>0){
            int edgeInd = Randomizer.nextInt(removableEdges.size());     
            
        	removeEmptyReassortmentEdge(network, removableEdges.get(edgeInd));
        	
        	networkEdges = new ArrayList<>(network.getEdges());
            
            removableEdges = networkEdges.stream()
                    .filter(e -> !e.isRootEdge())
                    .filter(e -> e.childNode.isRecombination())
                    .filter(e -> e.breakPoints.isEmpty())
                    .filter(e -> e.parentNode.isCoalescence())
                    .collect(Collectors.toList());            
        } 

    }    

    private void removeEmptyReassortmentEdge(RecombinationNetwork network, RecombinationNetworkEdge edgeToRemove) {

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
    } 

    /**
     * Retrieve sister of given edge
     * @param childEdge child edge
     * @return sister of given child edge
     */
    private RecombinationNetworkEdge getSisterEdge(RecombinationNetworkEdge childEdge) {
        int idx = childEdge.parentNode.getChildEdges().indexOf(childEdge);
        int otherIdx = (idx + 1) % 2;

        return childEdge.parentNode.getChildEdges().get(otherIdx);
    }

    /**
     * Retrieve spouse of given edge
     * @param parentEdge parent edge
     * @return spouse of given parent edge
     */
    private RecombinationNetworkEdge getSpouseEdge(RecombinationNetworkEdge parentEdge) {
        int idx = parentEdge.childNode.getParentEdges().indexOf(parentEdge);
        int otherIdx = (idx + 1) % 2;

        return parentEdge.childNode.getParentEdges().get(otherIdx);
    }
    
}