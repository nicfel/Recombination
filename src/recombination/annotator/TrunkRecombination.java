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

import beast.core.util.Log;
import beast.evolution.tree.Tree;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * A rewrite of TreeAnnotator that outputs how often reassortment events happen on trunk branches vs. other branches 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class TrunkRecombination extends RecombinationAnnotator {

    private enum TrunkDefinition { MostRecentSample, TipDistance }
    
    List<RecombinationNetworkNode> allTrunkNodes;
    List<Double> leaveDistance;
    List<Boolean> isTrunkNode;

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile = new File("recombination_trunk.txt");
        double burninPercentage = 10.0;
        TrunkDefinition trunkDefinition = TrunkDefinition.MostRecentSample;
        double minTipDistance = 2.0;
        BreakPoints breakPoints = new BreakPoints();


        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
                    "Definition of the trunk: " + trunkDefinition + "\n" +
            		"minimal distance to a tip to be considered trunk\n" + 
                    "(ignored in MostRecentSample Case): " + minTipDistance;
        }
    }

    public TrunkRecombination(NetworkAnnotatorOptions options) throws IOException {

        // Display options:
        System.out.println(options + "\n");

        // Initialise reader
        RecombinationLogReader logReader = new RecombinationLogReader(options.inFile,
                options.burninPercentage);

        System.out.println(logReader.getNetworkCount() + " Networks in file.");

        System.out.println("The first " + logReader.getBurnin() +
                 " (" + options.burninPercentage + "%) ACGs will be discarded " +
                "to account for burnin.");

	      System.out.println("\nWriting output to " + options.outFile.getName()
	      + "...");
        
        // compute the pairwise reassortment distances 
        try (PrintStream ps = new PrintStream(options.outFile)) {
	        for (RecombinationNetwork network : logReader){	        	
//	        	pruneNetwork(network, options.breakPoints);
	        	if (options.trunkDefinition == TrunkDefinition.MostRecentSample)
	        		computeTrunkReassortment(network, ps);
	        	else
	        		computeTrunkReassortmentLeaveDist(network, ps, options.minTipDistance);
	        	
	        	ps.print("\n");
	        }
	        ps.close();
        }
        System.out.println("\nDone!");
    }
    
    /**
     * gets how many reticulation events happen on the trunk vs. not on the trunk
     * The trunk is defined as any edge of the network which is between a sample with height 0 and the root
     * @param network
     * @param ps
     */
    private void computeTrunkReassortment(RecombinationNetwork network, PrintStream ps){    	
//        throw new IllegalArgumentException("not implemented");
    	
    	
    	
    	
       
//        int leaf = 0;
//        int coal = 0;
//        
//        double leaf_length = 0;
//        double coal_length = 0;        
//        
//        for (RecombinationNetworkNode n :  network.getNodes().stream()
//                .filter(e -> e.isLeaf())
//                .collect(Collectors.toList())) {
//        	if (n.getParentEdges().get(0).getLength()<1) {
//	        	if (n.getParentEdges().get(0).parentNode.isRecombination()) {
//	        		leaf++;
//	        	}  
//	        	leaf_length += n.getParentEdges().get(0).getLength()* n.getParentEdges().get(0).breakPoints.getLength();
//        	}
//        }
//        
//        for (RecombinationNetworkNode n :  network.getNodes().stream()
//        		.filter(e -> !e.getParentEdges().get(0).isRootEdge())
//                .filter(e -> e.isCoalescence())
//                .collect(Collectors.toList())) {
//        	if (n.getParentEdges().get(0).parentNode.isRecombination()) {
//        		coal++;
//        	}  
//        	coal_length += n.getParentEdges().get(0).getLength()* n.getParentEdges().get(0).breakPoints.getLength();
//        }
//
//        
//        ps.print(coal + "\t" + (leaf) + "\t" + coal_length + "\t" + leaf_length);


    }
    
    private void getAllAncestralEdges(RecombinationNetworkNode node){
    	if (allTrunkNodes.indexOf(node)!=-1)
    		return;
    	
    	allTrunkNodes.add(node);
		
    	
    	for (RecombinationNetworkEdge parentEdge : node.getParentEdges()){
    		if (parentEdge.isRootEdge()){
    			return;
    		}else{
    			getAllAncestralEdges(parentEdge.parentNode);   			
    		}
			
		}
    }
    
   
    int onTrunkCount;
    double trunkLength;
    Map<Integer, BreakPoints> rootBreaks;
    Map<Integer, Double> restLength;

    /**
     * gets how many reticulation events happen on the trunk vs. not on the trunk
     * The trunk is define as any edge on the network that has descendents that are more than minTipDistance 
     * away from that node
     * @param network
     * @param ps
     * @param minTipDistance
     */
    private void computeTrunkReassortmentLeaveDist(RecombinationNetwork network, PrintStream ps, double minTipDistance){    	
        
    	rootBreaks = new HashMap<Integer, BreakPoints>();
    	traversalRoots(network.getRootEdge(), new BreakPoints(network.totalLength));
    	
    	BreakPoints bp = new BreakPoints();
    	for (Integer i : rootBreaks.keySet()) {
    		bp.or(rootBreaks.get(i));
    	}
    	
        // get the length of the network        
        List<RecombinationNetworkNode> allLeafs = network.getNodes().stream()
                .filter(e -> e.isLeaf())
                .collect(Collectors.toList());
        
          
        
        restLength = new HashMap<>();
		for (RecombinationNetworkNode l : allLeafs) {
			labelDistances(l.getParentEdges().get(0),
					0.0, minTipDistance, 
					l.getParentEdges().get(0).breakPoints.copy());
		}
		
		
		
        // get the length of the network        
        List<RecombinationNetworkEdge> allNetworkEdges = network.getEdges().stream()
        		.filter(e -> !e.isRootEdge())
                .collect(Collectors.toList());
		

        // get the length of the network        
        List<RecombinationNetworkEdge> allTrunkEdges = network.getEdges().stream()
        		.filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.visited)
                .collect(Collectors.toList());
		
		double trunkWeightedLength = 0;
        
		for (RecombinationNetworkEdge e : allTrunkEdges) {
			BreakPoints bpEdge = e.childNode.dirtyBreakPoints.andCopy(e.breakPoints);
			trunkWeightedLength += e.getLength() * bpEdge.getLength();
		}
		
		for (Integer e : restLength.keySet()) {
			for (RecombinationNetworkEdge edge : allNetworkEdges) {
				if (edge.ID==e && !edge.childNode.visited) {
					BreakPoints bpEdge = edge.childNode.dirtyBreakPoints.andCopy(edge.breakPoints);
//					BreakPoints bpEdge = edge.breakPoints.copy();
					trunkWeightedLength += restLength.get(edge.ID)*bpEdge.getLength();
//					trunkWeightedLength += restLength.get(edge.ID);
				}
			}
		}

				
//		System.out.println(network);
//		System.exit(0);
        
		List<RecombinationNetworkNode> trunkRecombinationNodes = network.getNodes().stream()
                .filter(e -> e.isRecombination())
                .filter(e -> e.visited)
                .collect(Collectors.toList());
		
		int onTrunk = 0;
        
		for (RecombinationNetworkNode n : trunkRecombinationNodes) {
			BreakPoints bp1 = n.dirtyBreakPoints.andCopy(n.getParentEdges().get(0).breakPoints);
			BreakPoints bp2 = n.dirtyBreakPoints.andCopy(n.getParentEdges().get(1).breakPoints);
//			BreakPoints bp1 = n.getParentEdges().get(0).breakPoints.copy();
//			BreakPoints bp2 = n.getParentEdges().get(1).breakPoints.copy();

			if (!bp1.isEmpty() && !bp2.isEmpty())
					onTrunk++;

		}
		
		
		
		restLength = new HashMap<>();
		for (RecombinationNetworkNode l : allLeafs) {
			l.visited=true;
			labelDistances(l.getParentEdges().get(0),
					0.0, 0.0, 
					l.getParentEdges().get(0).breakPoints.copy());
		}

		
        // get the length of the network        
        List<RecombinationNetworkEdge> allEdges = network.getEdges().stream()
        		.filter(e -> !e.isRootEdge())
                .filter(e -> e.childNode.visited)
                .collect(Collectors.toList());
		
		double totalWeightedLength = 0;
        
		for (RecombinationNetworkEdge e : allEdges) {
			BreakPoints bpEdge = e.childNode.dirtyBreakPoints.andCopy(e.breakPoints);
//			BreakPoints bpEdge = e.breakPoints.copy();
			totalWeightedLength += e.getLength()*bpEdge.getLength();
		}
        
		List<RecombinationNetworkNode> recombinationNodes = network.getNodes().stream()
                .filter(e -> e.isRecombination())
                .filter(e -> e.visited)
                .collect(Collectors.toList());
		
		int totalRecombination = 0;
        
		for (RecombinationNetworkNode n : recombinationNodes) {
			BreakPoints bp1 = n.dirtyBreakPoints.andCopy(n.getParentEdges().get(0).breakPoints);
			BreakPoints bp2 = n.dirtyBreakPoints.andCopy(n.getParentEdges().get(1).breakPoints);
			
//			BreakPoints bp1 = n.getParentEdges().get(0).breakPoints.copy();
//			BreakPoints bp2 = n.getParentEdges().get(1).breakPoints.copy();


			if (!bp1.isEmpty() && !bp2.isEmpty())
				totalRecombination++;

		}

        
        ps.print(onTrunk + "\t" + (totalRecombination-onTrunk) + "\t" + trunkWeightedLength + "\t" + (totalWeightedLength-trunkWeightedLength));

    }
    
    private void labelDistances(RecombinationNetworkEdge edge, double leafDist, double minLeaveDist, BreakPoints bp) {
    	    	
   	
    	edge.childNode.dirtyBreakPoints.or(bp);
    	
    	if (leafDist<=minLeaveDist && (leafDist+edge.getLength())>minLeaveDist) {
    		if (restLength.containsKey(edge.ID)) {
    			restLength.replace(edge.ID, Math.max(restLength.get(edge.ID), leafDist+edge.getLength()-minLeaveDist));
    		}else {
    			restLength.put(edge.ID, leafDist+edge.getLength()-minLeaveDist);
    		}
    	}
    	
    	leafDist += edge.getLength();
    	   		
    	
    	if (edge.parentNode.isRecombination()) {
    		
    		
    		if (leafDist > minLeaveDist)
    			edge.parentNode.visited=true;

    		BreakPoints bp1 = edge.parentNode.getParentEdges().get(0).breakPoints.andCopy(bp);
    		BreakPoints bp2 = edge.parentNode.getParentEdges().get(1).breakPoints.andCopy(bp);
    		
    		
    		if (!bp1.isEmpty())
    			labelDistances(edge.parentNode.getParentEdges().get(0), leafDist, minLeaveDist, bp1);
    		
    		if (!bp2.isEmpty()) 
    			labelDistances(edge.parentNode.getParentEdges().get(1), leafDist, minLeaveDist, bp2);
    		
    	}else {
    		// check overlap
    		BreakPoints overlap = bp.copy();
    		if (overlap.isEmpty())
    			return;

    		if (leafDist > minLeaveDist)
    			edge.parentNode.visited = true;
    		
    		//check if it is a local root
    		if (rootBreaks.containsKey(edge.parentNode.getParentEdges().get(0).ID)) {
    			overlap.andNot(rootBreaks.get(edge.parentNode.getParentEdges().get(0).ID));
    		}
    		if (overlap.isEmpty())
    			return;    		
    		
    		labelDistances(edge.parentNode.getParentEdges().get(0), leafDist, minLeaveDist, overlap);
    	}
    }    

    
    private void getTrunkRateUpwards(RecombinationNetworkEdge edge, double leafDist, double minLeaveDist, BreakPoints bp) {
    	// check if the parent Edge is trunk, if yes, go onwards
    	if (!edge.parentNode.dirtyBreakPoints.isEmpty()) {
        	if (leafDist < minLeaveDist) {
    	    	if ((leafDist+edge.getLength()) > minLeaveDist) {
    	    		trunkLength += (leafDist+edge.getLength() - minLeaveDist)*bp.getLength();
    	    	}    		 
        	}else {
        		trunkLength += edge.getLength()*bp.getLength();
        	}
    	}    	
    	
		// check overlap
		BreakPoints overlap = bp.andNotCopy(edge.parentNode.dirtyBreakPoints);
		if (overlap.isEmpty())
			return;

    	
    	
    	leafDist += edge.getLength();
    		
    	
    	if (edge.parentNode.isRecombination()) {
    		BreakPoints bp1 = edge.parentNode.getParentEdges().get(0).breakPoints.andCopy(bp);
    		BreakPoints bp2 = edge.parentNode.getParentEdges().get(1).breakPoints.andCopy(bp);
    		
    		if (!bp1.isEmpty() && !bp2.isEmpty())    		
	    		if (leafDist > minLeaveDist)
	    			onTrunkCount++;
    		
    		if (!bp1.isEmpty())
    			getTrunkRateUpwards(edge.parentNode.getParentEdges().get(0), 0.0, minLeaveDist, bp1);
    		
    		if (!bp2.isEmpty()) 
    			getTrunkRateUpwards(edge.parentNode.getParentEdges().get(1), 0.0, minLeaveDist, bp2);
    		
    	}else {
    		
    		edge.parentNode.dirtyBreakPoints.or(overlap);
    		
    		//check if it is a local root
    		if (rootBreaks.containsKey(edge.parentNode.getParentEdges().get(0).ID)) {
//    			System.out.println(rootBreaks.get(edge.parentNode.getParentEdges().get(0).ID));
    			overlap.andNot(rootBreaks.get(edge.parentNode.getParentEdges().get(0).ID));
    		}
    		if (overlap.isEmpty())
    			return;

    		
			getTrunkRateUpwards(edge.parentNode.getParentEdges().get(0), leafDist, minLeaveDist, overlap);
    	}
    }
    

    
    private void getTotalRateUpwards(RecombinationNetworkEdge edge, BreakPoints bp) {
		trunkLength += edge.getLength()*bp.getLength();
    		
    	
    	if (edge.parentNode.isRecombination()) {
    		BreakPoints bp1 = edge.parentNode.getParentEdges().get(0).breakPoints.andCopy(bp);
    		BreakPoints bp2 = edge.parentNode.getParentEdges().get(1).breakPoints.andCopy(bp);
    		
    		if (!bp1.isEmpty() && !bp2.isEmpty())    		
    			onTrunkCount++;
    		
    		if (!bp1.isEmpty())
    			getTotalRateUpwards(edge.parentNode.getParentEdges().get(0), bp1);
    		
    		if (!bp2.isEmpty()) 
    			getTotalRateUpwards(edge.parentNode.getParentEdges().get(1), bp2);
    		
    	}else {
    		// check overlap
    		BreakPoints overlap = bp.andNotCopy(edge.parentNode.dirtyBreakPoints);
    		if (overlap.isEmpty())
    			return;
    		
    		edge.parentNode.dirtyBreakPoints.or(overlap);
    		
    		//check if it is a local root
    		if (rootBreaks.containsKey(edge.parentNode.getParentEdges().get(0).ID)) {
//    			System.out.println(rootBreaks.get(edge.parentNode.getParentEdges().get(0).ID));
    			overlap.andNot(rootBreaks.get(edge.parentNode.getParentEdges().get(0).ID));
    		}
    		if (overlap.isEmpty())
    			return;

    		
    		getTotalRateUpwards(edge.parentNode.getParentEdges().get(0), overlap);
    	}
    }

    
    private void getAllAncestralEdgesLeaveDist(RecombinationNetworkNode node, double dist, double threshold){
    	int index = allTrunkNodes.indexOf(node);
    	if (index==-1){
    		allTrunkNodes.add(node);
    		leaveDistance.add(dist);
    	}else{
    		if (leaveDistance.get(index)>dist)
    			return;
    		else if (leaveDistance.get(index) > threshold)
    			return;
			else
    			leaveDistance.set(index, dist);
    	}
    	
    	for (RecombinationNetworkEdge parentEdge : node.getParentEdges()){
    		if (parentEdge.isRootEdge()){
    			return;
    		}else{
    			getAllAncestralEdgesLeaveDist(parentEdge.parentNode, dist+parentEdge.getLength(), threshold);   			
    		}			
		}
    }
        
    private void traversalRoots(RecombinationNetworkEdge edge, BreakPoints breakPoints) {
    	BreakPoints bp = breakPoints.copy();
    	bp.and(edge.breakPoints);
    	
    	if (bp.isEmpty())
    		return;

    	RecombinationNetworkNode node = edge.childNode;
    	if (node.isCoalescence()) {
    		//get which loci coalesced here
    		BreakPoints bp1 = node.getChildEdges().get(0).breakPoints.copy();
    		bp1.and(node.getChildEdges().get(1).breakPoints); 
    		bp1.and(bp);
    		

    		if (!bp1.isEmpty()) {
     			if (rootBreaks.containsKey(node.getParentEdges().get(0).ID)) {
    				rootBreaks.get(node.getParentEdges().get(0).ID).or(bp1);
    			}else {
    				rootBreaks.put(node.getParentEdges().get(0).ID, bp1);
    			}
    		}
    		
    		
    		// get which loci did not coalesce
    		BreakPoints bp2 = bp.copy();
    		bp2.andNot(bp1);
    		
    		traversalRoots(node.getChildEdges().get(0), bp2);
    		traversalRoots(node.getChildEdges().get(1), bp2);    				
    	}else if (node.isRecombination()) {
    		traversalRoots(node.getChildEdges().get(0), bp);
    	}else {
    		return;
    	}	

	}
    
	void getTrunkRate(RecombinationNetworkNode node, BreakPoints computeFor_BP, 
			int prev_edge_ID, BreakPoints prev_Pointer, double leafDist, double minLeafDist) {   
		

    	if (computeFor_BP.isEmpty())
    		return;   
    	      
        if (node.isRecombination()) {
        	for (RecombinationNetworkEdge edge : node.getParentEdges()) {
        		BreakPoints bp = computeFor_BP.copy();      		
        		bp.andPR(edge.passingRange); 
        		if (!bp.isEmpty()) {   			
        			
        			getTrunkRate(edge.parentNode, bp, prev_edge_ID, prev_Pointer, leafDist + edge.getLength(), minLeafDist);
        		}
        	}       	
        }else { 	
        	// make a copy of the BP's
        	BreakPoints computeFor = computeFor_BP.copy();
        	RecombinationNetworkEdge edge = node.getParentEdges().get(0);    
      	
        	// compute with breakpoints are "visibly" coalescing at this node
        	if (node.overlap==null) {
        		node.overlap = node.getChildEdges().get(0).breakPoints.andCopy(node.getChildEdges().get(1).breakPoints);
        	}
        	
        	// test if compute for is visibly coalescing here
    		BreakPoints cf_only = computeFor.andNotCopy(node.overlap);

    		if (!cf_only.isEmpty()) {    
    			getTrunkRate(edge.parentNode, cf_only, prev_edge_ID, prev_Pointer, leafDist, minLeafDist);
                // see "how" much is left of the compute for BP
    			computeFor.andNot(cf_only);
    		}
    		
    		if (computeFor.isEmpty())
    			return;
    		
    		
    		boolean exists = false;
    		
    		for (int i = 0; i < node.prevPointer.size(); i++) {
	    		if (node.prevPointer.get(i)==prev_edge_ID &&
	    				node.dummy2.get(i).equals(prev_Pointer)) {	
		    		node.dummy.get(i).or(computeFor);		    		
		    		exists = true;
	    		}
    		}
    		
    		if (!exists) {
	    		node.dummy.add(computeFor.copy());		    		
	    		node.prevPointer.add(prev_edge_ID);
	    		node.dummy2.add(prev_Pointer.copy());
    		}
    		   		

    		for (int i = 0; i < node.dummy.size();i++) {
    			if (node.prevPointer.get(i)!=prev_edge_ID) {    			
	        		// get the overlap
	        		if (node.dummy.get(i).overlapFast(computeFor)) {        			
		        		BreakPoints bp_here = node.dummy.get(i).andCopy(computeFor);
		        		computeFor.andNot(bp_here);
	                	// only pass on loci for which the root has not been reached yet.		                	
                		node.dummy3.or(bp_here);                
	        		} 
	        	}
    		}
    		
    		
    		
    		if (node.dummy3.equals(node.overlap)) {
        		for (int i = 0; i < node.dummy.size(); i++) {
        			for (int j = i + 1; j < node.dummy.size(); j++) {
	        			if (node.prevPointer.get(i)!=node.prevPointer.get(j)) {
	        				if (node.dummy.get(i).overlapFast(node.dummy.get(j))) {
		        				BreakPoints bp1 = node.dummy.get(i).andCopy(node.dummy.get(j));
	        		        	if (bp1.overlapFast(node.dirtyBreakPoints) || 
	        		        			node.getChildEdges().get(0).isDirty()==Tree.IS_FILTHY ||
	        		        			node.getChildEdges().get(1).isDirty()==Tree.IS_FILTHY ) {	  
	        		        		
	        	        		}else {    
	        	        		}	        		        		        		        	
	                        	if (!edge.isRootEdge()) {
	                        		BreakPoints rootBp = rootBreaks.get(edge.ID);
	                        		if (rootBp!=null) {
	                        			if (!bp1.overlapFast(rootBp)) {
	                        				getTrunkRate(edge.parentNode, bp1, edge.ID, bp1, leafDist, minLeafDist);	                        				
	                        			}
	                        		}else {
	                        			getTrunkRate(edge.parentNode, bp1, edge.ID, bp1, leafDist, minLeafDist);
	                        		}
	                        	}
	        				}
	        			}
        			}
        		}
    		}

    		
		}        
    }

    
    
    
   

    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(NetworkAnnotatorOptions options) {

        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("Reassortment Event Trunk Mapper");

        JLabel logFileLabel = new JLabel("Reassortment Network log file:");
        JLabel outFileLabel = new JLabel("Output file:");
        JLabel burninLabel = new JLabel("Burn-in percentage:");
        JLabel trunkDefinitionLabel = new JLabel("Trunk definition:");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");

        JTextField minTipDistance = new JTextField(20);
        minTipDistance.setEditable(true);
//        minTipDistance.setEnabled(false);        

        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.burninPercentage));
        burninSlider.setMajorTickSpacing(50);
        burninSlider.setMinorTickSpacing(10);
        burninSlider.setPaintTicks(true);
        burninSlider.setPaintLabels(true);
        burninSlider.setSnapToTicks(true);

        JComboBox<TrunkDefinition> heightMethodCombo = new JComboBox<>(TrunkDefinition.values());

//        JSlider thresholdSlider = new JSlider(JSlider.HORIZONTAL,
//                0, 100, (int)(options.convSupportThresh));
//        thresholdSlider.setMajorTickSpacing(50);
//        thresholdSlider.setMinorTickSpacing(10);
//        thresholdSlider.setPaintTicks(true);
//        thresholdSlider.setPaintLabels(true);
//        thresholdSlider.setSnapToTicks(true);

        Container cp = dialog.getContentPane();
        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
        cp.setLayout(boxLayout);

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);
        layout.setAutoCreateGaps(true);
        layout.setAutoCreateContainerGaps(true);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(outFileLabel)
                        .addComponent(burninLabel)
                        .addComponent(trunkDefinitionLabel))
//                        .addComponent(thresholdLabel)
//                        .addComponent(geneFlowCheckBox))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename)
                        .addComponent(burninSlider)
                        .addComponent(heightMethodCombo)
//                        .addComponent(thresholdSlider)
                        .addComponent(minTipDistance))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFileButton)
                        .addComponent(outFileButton))
                );

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(inFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(inFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(burninLabel)
                        .addComponent(burninSlider,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(trunkDefinitionLabel)
                        .addComponent(heightMethodCombo,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(minTipDistance))
                );

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
            options.trunkDefinition = (TrunkDefinition)heightMethodCombo.getSelectedItem();
            options.minTipDistance = Double.parseDouble(minTipDistance.getText());
            dialog.setVisible(false);
        });
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

        JFileChooser inFileChooser = new JFileChooser();
        inFileButton.addActionListener(e -> {
            inFileChooser.setDialogTitle("Select Reassortment Network log file to summarize");
            if (options.inFile == null)
                inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = inFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.inFile = inFileChooser.getSelectedFile();
                inFilename.setText(inFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
            outFileChooser.setDialogTitle("Select output file name.");
            if (options.inFile != null)
                outFileChooser.setCurrentDirectory(options.inFile);
            else
                outFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
            }
        });

        cp.add(buttonPanel);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        return !canceled[0];
    }

    /**
     * Prepare JFrame to which ACGAnnotator output streams will be
     * directed.
     */
    private static void setupGUIOutput() {

        JFrame frame = new JFrame();
        frame.setTitle("Reassortment Event Locator");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JTextArea textArea = new JTextArea(25, 80);
        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
        textArea.setEditable(false);
        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> System.exit(0));
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

        // Redirect streams to output window:
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                SwingUtilities.invokeLater(() -> {
                    if ((char)b == '\r') {
                        int from = textArea.getText().lastIndexOf("\n") + 1;
                        int to = textArea.getText().length();
                        textArea.replaceRange(null, from, to);
                    } else
                        textArea.append(String.valueOf((char) b));
                });
            }
        };

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));

        frame.pack();
        frame.setVisible(true);
    }

    public static String helpMessage =
            "TrunkReassortment - counts how many reassortment events happened on trunk and non-trunk nodes.\n"
                    + "\n"
                    + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-trunkDefinition {MostRecentSample, TipDistance} Choose trunk definition method.\n"
                    + "                         (default MostRecentSample)\n"
                    + "-burnin percentage       Choose _percentage_ of log to discard\n"
                    + "                         in order to remove burn-in period.\n"
                    + "                         (Default 10%)\n"
                    + "-minTipDistance     		minimum distance between internal network node\n"
                    + "                         and tip node such that the internal node is considered trunk.\n"
                    + "                         If not  specified, the trunk is any node between samples\n"
                    + "                         height=0 and the root.\n"
                    + "\n"
                    + "If no output file is specified, output is written to a file\n"
                    + "named 'reassortment_distances.txt'.";

    /**
     * Print usage info and exit.
     */
    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    /**
     * Display error, print usage and exit with error.
     */
    public static void printUsageAndError(String errMsg) {
        System.err.println(errMsg);
        System.err.println(helpMessage);
        System.exit(1);
    }

    /**
     * Retrieve TrunkReassortment options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, NetworkAnnotatorOptions options) {
        int i=0;
        while (args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

                case "-burnin":
                    if (args.length<=i+1)
                        printUsageAndError("-burnin must be followed by a number (percent)");

                    try {
                        options.burninPercentage = Double.parseDouble(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
                    }

                    if (options.burninPercentage<0 || options.burninPercentage>100) {
                        printUsageAndError("Burnin percentage must be >= 0 and < 100.");
                    }

                    i += 1;
                    break;
                case "-trunkDefinition":
                    if (args.length<=i+1) {
                        printUsageAndError("-trunkDefinition must be either mostRecentSample or minTipDistance.");
                    }

                    try {
                    	if (args[i + 1].equals("mostRecentSample"))
                    		options.trunkDefinition = TrunkDefinition.MostRecentSample;
                    	else if (args[i + 1].equals("minTipDistance"))
                    		options.trunkDefinition = TrunkDefinition.TipDistance;
                    	else
                    		throw new NumberFormatException();

                    } catch (NumberFormatException e) {
                        printUsageAndError("trunkDefinition must be either mostRecentSample or minTipDistance.");
                    }

                    i += 1;
                    break;


                case "-minTipDistance":
                    if (args.length<=i+1) {
                        printUsageAndError("-minTipDistance must be followed by a number.");
                    }

                    try {
                        options.minTipDistance =
                                Double.parseDouble(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("minTipDistance must be a positive number. ");
                     }

                    i += 1;
                    break;
                    
                case "-subsetRange":
                    if (args.length<=i+1) {
                        printUsageAndError("-subsetRange must be a range in the format of 0-100.");
                    }

                    try {
                    	String[] argarray = args[i + 1].split(",");
                    	List<Integer> bp_list = new ArrayList<>();
                    	for (int j = 0; j < argarray.length; j++) {
                    		String[] tmp = argarray[j].split("-");
                    		bp_list.add(Integer.parseInt(tmp[0]));
                    		bp_list.add(Integer.parseInt(tmp[1]));
                    	}
                		options.breakPoints.init(bp_list);
                    } catch (NumberFormatException e) {
                        printUsageAndError("removeSegments must be an array of integers separated by commas if more than one");
                     }

                    i += 1;
                    break;
                    
                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }

        if (i >= args.length)
            printUsageAndError("No input file specified.");
        else
            options.inFile = new File(args[i]);

        if (i+1<args.length)
            options.outFile = new File(args[i+1]);
    }

    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
    	NetworkAnnotatorOptions options = new NetworkAnnotatorOptions();

        if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
                Log.warning.println("Error setting cross-platform look and feel.");
            }

            try {
                SwingUtilities.invokeAndWait(() -> {
                    if (!getOptionsGUI(options))
                        System.exit(0);

                    setupGUIOutput();
                });
            } catch (InterruptedException | InvocationTargetException e) {
                e.printStackTrace();
            }


        } else {
            getCLIOptions(args, options);
        }

        // Run ACGAnnotator
        try {
            new TrunkRecombination(options);

        } catch (Exception e) {
            if (args.length == 0) {
                JOptionPane.showMessageDialog(null, e.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                System.err.println("Error: " + e.getMessage());
                e.printStackTrace();
                System.err.println();
                System.err.println(helpMessage);
            }

            System.exit(1);
        }
    }
}