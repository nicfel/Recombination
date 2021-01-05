package recombination.annotator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import javax.swing.plaf.synth.SynthSeparatorUI;

import beast.evolution.tree.Tree;
import beast.math.statistic.DiscreteStatistics;
import coalre.networkannotator.NetworkCladeSystem.BitSetArray;
import coalre.networkannotator.NetworkCladeSystem.DummyClade;
import coalre.networkannotator.NetworkCladeSystem.ReassortmentClade;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

/**
 * extracted from TreeAnnotator
 */
public class RecombinationNetworkCladeSystem {

//    protected Map<Double, DummyClade> newCoalescentCladeMap;
//    protected List<Map<Double, DummyClade>> coalescentCladeMap;
//    protected Map<BitSetArray, ReassortmentClade> cladeMap = new HashMap<>();
    
    protected List<DummyClade> cladeMap = new ArrayList<>();
   
    public List<String> leafNodeMap;
    
    protected Map<BitSet, DummyClade> leafCladeMap = new HashMap<>();
    
    
    private Map<Integer, List<BreakPoints>> visitedBP;
    private Map<Integer, List<BitSet>> visitedBits;
    private Map<Integer, Double> nodeProbs;
    
    private Map<Integer, List<BreakPoints>> nodeBP;
    private Map<Integer, List<BitSet>> nodeBits;
    private Map<Integer, List<Double>> nodeHeights;

    
    protected boolean started;
    
    public RecombinationNetworkCladeSystem() { 
    }
    
    /**
     * adds all leaf labels to a list
     * @param leafNodes
     * @param nrSegments
     */
    public void setLeafLabels(List<RecombinationNetworkNode> leafNodes, int totalLength){
    	leafNodeMap = new ArrayList<String>();
    	for (RecombinationNetworkNode leaf : leafNodes)
    		leafNodeMap.add(leaf.getTaxonLabel());    	
    }


    /**
     * adds all the clades in the tree
     */
    public void add(RecombinationNetwork network, boolean includeTips) { 	
    	
        // Recurse over the tree and add all the clades (or increment their
        // frequency if already present). The root clade is added too (for
        // annotation purposes).   
		visitedBP = new HashMap<>(); 
		visitedBits = new HashMap<>(); 
		
		for (RecombinationNetworkNode l : network.getLeafNodes())
			addCladesUpwards(l, l.getParentEdges().get(0).breakPoints, new BitSet());
    } 
    
    private void addCladesUpwards(RecombinationNetworkNode node, BreakPoints computeFor_BP, BitSet bits) {
    	if (computeFor_BP.isEmpty())
    		return;    	
    	
    	if (node.isLeaf()) {
            int index = getTaxonIndex(node);
   			bits.set(index);
   			BreakPoints computeFor = computeFor_BP.copy();
   			
			addCoalescentClade(bits, node, computeFor_BP);	       			

   			addCladesUpwards(node.getParentEdges().get(0).parentNode, computeFor, bits);
        }else if (node.isRecombination()) {
        	for (RecombinationNetworkEdge edge : node.getParentEdges()) {
       			BreakPoints computeFor = computeFor_BP.copy();
       			computeFor.and(edge.breakPoints.copy());
       			BitSet biton = (BitSet) bits.clone();
       			addCladesUpwards(edge.parentNode, computeFor, biton);
        	}        	
        }else {
   			BreakPoints computeFor = computeFor_BP.copy();
   			RecombinationNetworkEdge edge = node.getParentEdges().get(0);
   			
        	// compute with breakpoints are "visibly" coalescing at this node
        	BreakPoints overlap = node.getChildEdges().get(0).breakPoints.copy();        	
        	overlap.and(node.getChildEdges().get(1).breakPoints);	
        	
        	// test if compute for is visibly coalescing here
    		BreakPoints cf_only = computeFor.copy();
    		cf_only.andNot(overlap);
    		
    		if (!cf_only.isEmpty()) {    			
                if (!edge.isRootEdge()) {        
           			BitSet biton = (BitSet) bits.clone();
           			addCladesUpwards(edge.parentNode, cf_only, biton);
            	}
                // see "how" much is left of the compute for BP
    			computeFor.andNot(cf_only);
    		}

    		BreakPoints bp_in = computeFor.copy();
    		if (visitedBP.get(node.ID)!=null) {
	    		for (int i = 0; i < visitedBP.get(node.ID).size();i++) {
	        		BreakPoints bp_here = visitedBP.get(node.ID).get(i).copy();
	        		// get the overlap
	        		bp_here.and(bp_in);       		
	        		if (!bp_here.isEmpty()) {
	           			BitSet biton = (BitSet) bits.clone();
	           			biton.or(visitedBits.get(node.ID).get(i));
	        			addCoalescentClade(biton, node, bp_here);	        			
		        		computeFor.andNot(bp_here);
		                if (!edge.isRootEdge()) {
		                	bp_here.and(edge.breakPoints);

		           			addCladesUpwards(edge.parentNode, bp_here, biton);
		                }
	        		}  
	        	}
    		}else{
	    		visitedBP.put(node.ID, new ArrayList<>());
	    		visitedBits.put(node.ID, new ArrayList<>());
    		}
    		visitedBP.get(node.ID).add(computeFor.copy());
    		visitedBits.get(node.ID).add((BitSet) bits.clone());
        }
    	
    	
    }
    
    private void cladeProbsUpwards(RecombinationNetworkNode node, BreakPoints computeFor_BP, BitSet bits, int totalLength) {
    	if (computeFor_BP.isEmpty())
    		return;
    	
    	if (node.isLeaf()) {
            int index = getTaxonIndex(node);
   			bits.set(index);
   			BreakPoints computeFor = computeFor_BP.copy();
   			cladeProbsUpwards(node.getParentEdges().get(0).parentNode, computeFor, bits, totalLength);
        }else if (node.isRecombination()) {
        	for (RecombinationNetworkEdge edge : node.getParentEdges()) {
       			BreakPoints computeFor = computeFor_BP.copy();
       			computeFor.and(edge.breakPoints.copy());
       			BitSet biton = (BitSet) bits.clone();

       			cladeProbsUpwards(edge.parentNode, computeFor, biton, totalLength);
        	}        	
        }else {
   			BreakPoints computeFor = computeFor_BP.copy();
   			RecombinationNetworkEdge edge = node.getParentEdges().get(0);
   			
        	// compute with breakpoints are "visibly" coalescing at this node
        	BreakPoints overlap = node.getChildEdges().get(0).breakPoints.copy();        	
        	overlap.and(node.getChildEdges().get(1).breakPoints);	
        	
        	// test if compute for is visibly coalescing here
    		BreakPoints cf_only = computeFor.copy();
    		cf_only.andNot(overlap);
    		
    		if (!cf_only.isEmpty()) {    			
                if (!edge.isRootEdge()) {  
           			BitSet biton = (BitSet) bits.clone();
                	cladeProbsUpwards(edge.parentNode, cf_only, biton, totalLength);
            	}
                // see "how" much is left of the compute for BP
    			computeFor.andNot(cf_only);
    		}

    		BreakPoints bp_in = computeFor.copy();
    		if (visitedBP.get(node.ID)!=null) {
	    		for (int i = 0; i < visitedBP.get(node.ID).size();i++) {
	        		BreakPoints bp_here = visitedBP.get(node.ID).get(i).copy();
	        		// get the overlap
	        		bp_here.and(bp_in);       		
	        		if (!bp_here.isEmpty()) {
	           			BitSet biton = (BitSet) bits.clone();
	           			biton.or(visitedBits.get(node.ID).get(i));
	        			getCoalescentCladeProb(biton, node, bp_here, totalLength);	        			
		        		computeFor.andNot(bp_here);
		                if (!edge.isRootEdge()) {
		                	bp_here.and(edge.breakPoints);
		                	cladeProbsUpwards(edge.parentNode, bp_here, biton, totalLength);
		                }
	        		}  
	        	}
    		}else{
	    		visitedBP.put(node.ID, new ArrayList<>());
	    		visitedBits.put(node.ID, new ArrayList<>());
    		}
    		visitedBP.get(node.ID).add(computeFor.copy());
    		visitedBits.get(node.ID).add((BitSet) bits.clone());
		}
    }
    
    /**
     * adds a new dummy clade to a coalescent node in the network
     * @param bits
     * @param segment
     **/    	
    private void addCoalescentClade(BitSet bits, RecombinationNetworkNode node, BreakPoints bp) {
    	DummyClade clade = new DummyClade(bits, bp);    	    	
    	int i = cladeMap.indexOf(clade);

        if (i != -1) {
        	cladeMap.get(i).bp.add(bp.copy());
        	cladeMap.get(i).attributeValues.add(null);
        }else{
        	clade.attributeValues.add(null);
        	cladeMap.add(clade);
        }
        
    }    
    
    /**
     * adds a new dummy clade to a coalescent node in the network
     * @param bits
     * @param totalNumber 
     * @param segment
     **/    	
    private void getCoalescentCladeProb(BitSet bits, RecombinationNetworkNode node, BreakPoints bp, int totalLength) {   	
    	double weighted_support = 0;
    	boolean found = false;
    	for (DummyClade clade : cladeMap) {
    		if (clade.bits.equals(bits)) {
    			for (BreakPoints bp_clade : clade.bp) {
        			BreakPoints bp_in = bp.copy();
        			bp_in.and(bp_clade);
        			weighted_support += bp_in.getGeneticLength()/totalLength;
        			found = true;

    			}
    			
    	    	if (weighted_support==0) {
    	    		System.out.println(clade.bp);
    	    		System.out.println(clade);
    	    		System.exit(0);
    	    	}

    		}
    	}
    	
    	if (!found) {
    		System.out.println(node.getHeight());
    		System.out.println(bp);
    		System.out.println(bits);
    		throw new IllegalArgumentException("clade not found");
    	}
    	
    	if (nodeProbs.get(node.ID)==null)
    		nodeProbs.put(node.ID, weighted_support);
		else
    		nodeProbs.replace(node.ID, nodeProbs.get(node.ID)+ weighted_support);
    }    

    
    /**
     * adds new dummy clades to the new clade map
     * @param bits
     * @param segment
     * @param height
     */
    protected void addReassortmentClade(BitSet bits, RecombinationNetworkNode node) {
//        DummyClade clade = newReassortmentCladeMap.get(segment).get(height);
//        if (clade == null) {
//            clade = new DummyClade(bits, height, isLeft);
//            newReassortmentCladeMap.get(segment).put(height, clade);
//        }else{
//        	throw new IllegalArgumentException("reassortment clade should never be found");
//        }
    }   
       
    public void collectAttributes(RecombinationNetwork network, Set<String> attributeNames, boolean includeTips) {
    	
		visitedBP = new HashMap<>(); 
		visitedBits = new HashMap<>(); 
		
		for (RecombinationNetworkNode l : network.getLeafNodes())
        	collectAttributes(l, attributeNames, l.getParentEdges().get(0).breakPoints, includeTips, new BitSet());
    }
    
    private void collectAttributes(RecombinationNetworkNode node, Set<String> attributeNames, BreakPoints computeFor_BP, boolean includeTips, BitSet bits) {
    	if (computeFor_BP.isEmpty())
    		return;    	
    	
    	if (node.isLeaf()) {
            int index = getTaxonIndex(node);
   			bits.set(index);
   			BreakPoints computeFor = computeFor_BP.copy();
   			
   			if (includeTips)
    			collectCoalescentClade(bits, node, computeFor_BP, attributeNames);	        			

   			
   			collectAttributes(node.getParentEdges().get(0).parentNode, attributeNames, computeFor, includeTips, bits);
        }else if (node.isRecombination()) {
        	for (RecombinationNetworkEdge edge : node.getParentEdges()) {
       			BreakPoints computeFor = computeFor_BP.copy();
       			computeFor.and(edge.breakPoints.copy());
       			BitSet biton = (BitSet) bits.clone();
       			collectAttributes(edge.parentNode, attributeNames, computeFor, includeTips, biton);
        	}        	
        }else {
   			BreakPoints computeFor = computeFor_BP.copy();
   			RecombinationNetworkEdge edge = node.getParentEdges().get(0);   			
   			
        	// compute with breakpoints are "visibly" coalescing at this node
        	BreakPoints overlap = node.getChildEdges().get(0).breakPoints.copy();        	
        	overlap.and(node.getChildEdges().get(1).breakPoints);	
        	
        	// test if compute for is visibly coalescing here
    		BreakPoints cf_only = computeFor.copy();
    		cf_only.andNot(overlap);
    		
    		if (!cf_only.isEmpty()) {    			
                if (!edge.isRootEdge()) {        
           			BitSet biton = (BitSet) bits.clone();
           			collectAttributes(edge.parentNode, attributeNames, cf_only, includeTips, biton);
            	}
                // see "how" much is left of the compute for BP
    			computeFor.andNot(cf_only);
    		}

    		BreakPoints bp_in = computeFor.copy();
    		if (visitedBP.get(node.ID)!=null) {
	    		for (int i = 0; i < visitedBP.get(node.ID).size();i++) {
	        		BreakPoints bp_here = visitedBP.get(node.ID).get(i).copy();
	        		// get the overlap
	        		bp_here.and(bp_in);       		
	        		if (!bp_here.isEmpty()) {
	           			BitSet biton = (BitSet) bits.clone();
	           			biton.or(visitedBits.get(node.ID).get(i));
	        			collectCoalescentClade(biton, node, bp_here, attributeNames);	        			
		        		computeFor.andNot(bp_here);
		                if (!edge.isRootEdge()) {
		                	bp_here.and(edge.breakPoints);
		                	collectAttributes(edge.parentNode, attributeNames, bp_here, includeTips, biton);
		                }
	        		}  
	        	}
    		}else{
	    		visitedBP.put(node.ID, new ArrayList<>());
	    		visitedBits.put(node.ID, new ArrayList<>());
    		}
    		visitedBP.get(node.ID).add(computeFor.copy());
    		visitedBits.get(node.ID).add((BitSet) bits.clone());
        }
    }
    
    private void collectCoalescentClade(BitSet bits, RecombinationNetworkNode node, BreakPoints bp, Set<String> attributeNames) {   	
    	DummyClade clade = new DummyClade(bits, bp);
    	    	
    	int index = cladeMap.indexOf(clade);  	
    	
        if (index != -1) {
            if (clade.attributeValues == null) {
                clade.attributeValues = new ArrayList<>();
            }
            
            int i = 0;
            Object[] values = new Object[attributeNames.size()];
            for (String attributeName : attributeNames) {

                Object value;
                switch (attributeName) {
                    case "height":
                        value = node.getHeight();
                        break;
//                    case "length":
//                        value = getBranchLength(node);
//                        break;
                    default:
                    	throw new IllegalArgumentException("Summary not implemented for values other than height and posterior");
//                        value = node.getMetaData(attributeName);
//                        if (value instanceof String && ((String) value).startsWith("\"")) {
//                            value = ((String) value).replaceAll("\"", "");
//                        }
//                        break;
                }
                values[i] = value;

                i++;
            }
            cladeMap.get(index).attributeValues.add(values);
            cladeMap.get(index).bp.add(bp.copy());
        }       
    }    
         
    public void summarizeAttributes(RecombinationNetwork network, Set<String> attributeNames, int heightSummary, int nrNetworks, boolean onTarget) {
    	
		visitedBP = new HashMap<>(); 
		visitedBits = new HashMap<>(); 
		
		nodeBP = new HashMap<>(); 
		nodeBits = new HashMap<>(); 
		nodeHeights = new HashMap<>();
		
		for (RecombinationNetworkNode l : network.getLeafNodes())
			summarizeAttributes(l, attributeNames, l.getParentEdges().get(0).breakPoints, new BitSet(), heightSummary, nrNetworks, onTarget);
		
		for (RecombinationNetworkEdge e : network.getEdges())
			e.visited=false;

		summarizeCoalescentNodes(network.getRootEdge(), attributeNames, heightSummary, nrNetworks, onTarget);
    }
        
    private void summarizeAttributes(RecombinationNetworkNode node, Set<String> attributeNames, 
    		BreakPoints computeFor_BP, BitSet bits, int heightSummary, int nrNetworks, boolean onTarget) {
    	
    	if (computeFor_BP.isEmpty())
    		return;    	
    	
    	if (node.isLeaf()) {
            int index = getTaxonIndex(node);
   			bits.set(index);
   			BreakPoints computeFor = computeFor_BP.copy();
   			
   			if (nodeBP.get(node.ID)==null) {
   				nodeBP.put(node.ID, new ArrayList<>());
   				nodeBits.put(node.ID, new ArrayList<>());
   			}
			nodeBP.get(node.ID).add(computeFor_BP.copy());
			nodeBits.get(node.ID).add((BitSet) bits.clone());   	
   			
   			summarizeAttributes(node.getParentEdges().get(0).parentNode, attributeNames, computeFor, bits, heightSummary, nrNetworks, onTarget);
        }else if (node.isRecombination()) {
        	for (RecombinationNetworkEdge edge : node.getParentEdges()) {
       			BreakPoints computeFor = computeFor_BP.copy();
       			computeFor.and(edge.breakPoints.copy());
       			BitSet biton = (BitSet) bits.clone();
       			summarizeAttributes(edge.parentNode, attributeNames, computeFor, biton, heightSummary, nrNetworks, onTarget);
        	}        	
        }else {
   			BreakPoints computeFor = computeFor_BP.copy();
   			RecombinationNetworkEdge edge = node.getParentEdges().get(0);
   			
   			
        	// compute with breakpoints are "visibly" coalescing at this node
        	BreakPoints overlap = node.getChildEdges().get(0).breakPoints.copy();        	
        	overlap.and(node.getChildEdges().get(1).breakPoints);	
        	
        	// test if compute for is visibly coalescing here
    		BreakPoints cf_only = computeFor.copy();
    		cf_only.andNot(overlap);
    		
    		if (!cf_only.isEmpty()) {    			
                if (!edge.isRootEdge()) {        
           			BitSet biton = (BitSet) bits.clone();
           			summarizeAttributes(edge.parentNode, attributeNames, cf_only, biton, heightSummary, nrNetworks, onTarget);
            	}
                // see "how" much is left of the compute for BP
    			computeFor.andNot(cf_only);
    		}

    		BreakPoints bp_in = computeFor.copy();
    		if (visitedBP.get(node.ID)!=null) {
	    		for (int i = 0; i < visitedBP.get(node.ID).size();i++) {
	        		BreakPoints bp_here = visitedBP.get(node.ID).get(i).copy();
	        		// get the overlap
	        		bp_here.and(bp_in);       		
	        		if (!bp_here.isEmpty()) {
	           			BitSet biton = (BitSet) bits.clone();
	           			biton.or(visitedBits.get(node.ID).get(i));
	           			
	           			if (nodeBP.get(node.ID)==null) {
	           				nodeBP.put(node.ID, new ArrayList<>());
	           				nodeBits.put(node.ID, new ArrayList<>());
	           				nodeHeights.put(node.ID, new ArrayList<>());
	           			}
	        			nodeBP.get(node.ID).add(computeFor_BP.copy());
	        			nodeBits.get(node.ID).add((BitSet) biton.clone());   			

	        			computeFor.andNot(bp_here);
		                if (!edge.isRootEdge()) {
		                	bp_here.and(edge.breakPoints);
		                	summarizeAttributes(edge.parentNode, attributeNames, bp_here, biton, heightSummary, nrNetworks, onTarget);
		                }
	        		}  
	        	}
    		}else{
	    		visitedBP.put(node.ID, new ArrayList<>());
	    		visitedBits.put(node.ID, new ArrayList<>());
    		}
    		visitedBP.get(node.ID).add(computeFor.copy());
    		visitedBits.get(node.ID).add((BitSet) bits.clone());
        }
	}
    
    
    private void summarizeCoalescentNodes(RecombinationNetworkEdge edge, Set<String> attributeNames, 
    		int heightSummary, int nrNetworks, boolean onTarget) {

    	if (edge.visited)
    		return;
    	
    	edge.visited=true;
    	
    	RecombinationNetworkNode node = edge.childNode;
    	
    	if (node.isLeaf()) {  			  			
        	summarizeClade(node, attributeNames, heightSummary, nrNetworks, onTarget);
        }else if (node.isRecombination()) {
        	for (RecombinationNetworkEdge e : node.getChildEdges()) {
        		summarizeCoalescentNodes(e, attributeNames, heightSummary, nrNetworks, onTarget);
        	}        	
        }else {
        	summarizeClade(node, attributeNames, heightSummary, nrNetworks, onTarget);
        	for (RecombinationNetworkEdge e : node.getChildEdges()) {
        		summarizeCoalescentNodes(e, attributeNames, heightSummary, nrNetworks, onTarget);
        	}        	
        }
	}
    

    
    private void summarizeClade(RecombinationNetworkNode node, Set<String> attributeNames, int heightSummary, int nrNetworks, boolean onTarget) {   	
    	List<BreakPoints>  bp_list = nodeBP.get(node.ID);
    	List<BitSet>  bits_list = nodeBits.get(node.ID);
    	
    	if (bp_list==null) {
    		return;
    	}
    	
		// keeps track of the height of the target network,
		Double targetHeight = null;
		
		int k = 0;
		
        for (String attributeName : attributeNames) {
            switch (attributeName) {
            	case "height":
            		node.setMetaData("");
            		
            		double totalLength = 0;
            		double consideredLength = 0;
            		
            		List<HeightWeights> heights = new ArrayList<>();
            		
            		for (int i = 0; i < bp_list.size(); i++) {
            			consideredLength += bp_list.get(i).getGeneticLength();
            					
            			DummyClade indexClade = new DummyClade(bits_list.get(i), bp_list.get(i));
            			int index = cladeMap.indexOf(indexClade);            			
            			
            			if (index==-1)
            				throw new IllegalArgumentException("didn't find clade");
            			
            			int j = 0;
            			          			
                		for (BreakPoints bp_here : cladeMap.get(index).bp) {
                			BreakPoints bp_tmp = bp_here.copy();
                			bp_tmp.and(bp_list.get(i));
                			if (!bp_tmp.isEmpty()) {
                				if (cladeMap.get(index).attributeValues.get(j)!=null) {
                					heights.add(new HeightWeights((Double) cladeMap.get(index).attributeValues.get(j)[k], bp_tmp.getGeneticLength()));
	                    			totalLength += bp_tmp.getGeneticLength();
                				}
                			}
                			j++;
                		}  
            		}
            		
            		
            		if(onTarget)
            			targetHeight = node.getHeight();     
            		
            		
		            
            		double posterior = (double) totalLength / (double) (nrNetworks*consideredLength);
            		
        		            		
            		double[] quantiles = getQuantiles(heights, totalLength);
            		
            		
            		if(!onTarget){
	            		if (heights.size()>0){
							if (heightSummary==1){
								node.setHeight(quantiles[3]);
							}else if (heightSummary==2){
							    node.setHeight(quantiles[1]);
							}
	            		}
            		}

            		double minHPD,maxHPD;
            		if (heights.size()>0){
	                    minHPD = quantiles[0];
	                    maxHPD = quantiles[2];
            		}else{
            			minHPD = node.getHeight();
            			maxHPD = node.getHeight();
            		}	            		
            		
            		if (targetHeight!=null){
            			node.setMetaData(",posterior=" + posterior +
            					",targetHeight=" + targetHeight +
            					",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
		            		node.getMetaData() + "");

            		}else{
            			node.setMetaData(",posterior=" + posterior +
            					",height_95%_HPD={" + minHPD + "," + maxHPD + "}" + 
		            		node.getMetaData() + "");
            		}

                   
                case "length":
                    break;
                default:
                	throw new IllegalArgumentException("");
            }
            k++;
        }
    	
        
    }   
    
    private double[] getQuantiles(List<HeightWeights> heights, double totalLength) {
    	double[] quantiles = new double[4];
    	// sort data
    	HeightsComparator comparator = new HeightsComparator();
        // sort the array    	
    	Collections.sort(heights, comparator);
    	double currweights = 0.0;
    	double mean = 0.0;
    	boolean lower=false,median=false,upper=false;
    	for (HeightWeights height : heights) {
    		mean += height.height*height.weight;
    		currweights += height.weight;
    		if (!lower && currweights>0.025*totalLength) {
    			quantiles[0] = height.height;
    			lower=true;
    		}
    		if (!median && currweights>0.5*totalLength) {
    			quantiles[1] = height.height;
    			median=true;
    		}
    		if (!upper && currweights>0.975*totalLength) {
    			quantiles[2] = height.height;
    			upper=true;
    		}
    	}    	
    	
    	quantiles[3] = mean/currweights;
    	return quantiles;
    }

    public double getLogCladeCredibility(RecombinationNetwork network, int totalNumber){
		visitedBP = new HashMap<>(); 
		visitedBits = new HashMap<>();
		nodeProbs = new HashMap<>();
    	
		for (RecombinationNetworkNode l : network.getLeafNodes())
			cladeProbsUpwards(l, l.getParentEdges().get(0).breakPoints, new BitSet(), network.totalLength);
    			
    	double logP = 0.0;
    	for (Integer id : nodeProbs.keySet()) {
    		logP += Math.log(nodeProbs.get(id)/totalNumber);
    	}    		
		return logP;
    }
    
    /**
     * get the index of a leaf node
     */ 
    protected int getTaxonIndex(RecombinationNetworkNode leaf){
    	return leafNodeMap.indexOf(leaf.getTaxonLabel());
    }
        
    /**
     * clade needed to temporarily store reassortment clades before putting them all together
     *
     */
    public class DummyClade {
    	public DummyClade(BitSet bits, BreakPoints bp) {
	        this.bits = (BitSet) bits.clone();
            count = 1;
            this.bp = new ArrayList<>();
            attributeValues = new ArrayList<>();
            this.bp.add(bp.copy());
    	}
        
        public List<Object[]> getAttributeValues() {
            return attributeValues;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final DummyClade clade = (DummyClade) o;
            
            return !(bits != null ? !bits.equals(clade.bits) : clade.bits != null);

        }

        @Override
        public int hashCode() {
            return (bits != null ? bits.hashCode() : 0);
        }

        @Override
        public String toString() {
            return  count + " " + bp;
       }
        
        BitSet bits;
        List<BreakPoints> bp;
        int count;
        List<Object[]> attributeValues = null;
    }
   
    public class HeightWeights {
    	public HeightWeights(double height, double weight) {
    		this.height = height;
    		this.weight = weight;
    	}
    	
    	public String toString() {
    		return "" + height;
    	}

    	double height;
    	double weight;
    } 
    
    
    /**
     * SiteComparator is used for ordering the sites,
     * which makes it easy to identify patterns.
     */
    class HeightsComparator implements Comparator<HeightWeights> {
        @Override
		public int compare(HeightWeights o1, HeightWeights o2) {
            if (o1.height > o2.height) {
                return 1;
            }
            if (o1.height < o2.height) {
                return -1;
            }
            
            return 0;
        }
    } // class SiteComparator


    
    
}