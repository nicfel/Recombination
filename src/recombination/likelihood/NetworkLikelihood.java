/*
* File TreeLikelihood.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/


package recombination.likelihood;

import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;
import recombination.util.Partials;

import java.util.*;
import java.util.stream.Collectors;

@Description("Calculates the probability of sequence data on a beast.tree given a site and substitution model using " +
        "a variant of the 'peeling algorithm'. For details, see" +
        "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class NetworkLikelihood extends GenericTreeLikelihood {

    final public Input<RecombinationNetwork> networkInput = new Input<>("recombinationNetwork", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    final public Input<Boolean> m_useAmbiguities = new Input<>("useAmbiguities", "flag to indicate that sites containing ambiguous states should be handled instead of ignored (the default)", false);
    final public Input<Boolean> m_useTipLikelihoods = new Input<>("useTipLikelihoods", "flag to indicate that partial likelihoods are provided at the tips", false);
    final public Input<String> implementationInput = new Input<>("implementation", "name of class that implements this treelikelihood potentially more efficiently. "
    		+ "This class will be tried first, with the TreeLikelihood as fallback implementation. "
    		+ "When multi-threading, multiple objects can be created.", "beast.evolution.likelihood.BeagleTreeLikelihood");
    
    public static enum Scaling {none, always, _default};
    final public Input<Scaling> scaling = new Input<>("scaling", "type of scaling to use, one of " + Arrays.toString(Scaling.values()) + ". If not specified, the -beagle_scaling flag is used.", Scaling._default, Scaling.values());
    

    /**
     * calculation engine *
     */
    protected NetworkLikelihoodCore likelihoodCore;
//    protected BeagleTreeLikelihood beagle;

    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected SubstitutionModel substitutionModel;
    protected SiteModel.Base m_siteModel;
    protected BranchRateModel.Base branchRateModel;
    boolean[] computeForPatterns;

    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;


//    protected double[] m_fRootPartials;

    protected double[] probabilities;

    protected int matrixSizeExtended;
    protected int matrixSize;
    
    private double[] mat1;
    private double[] mat2;

    /**
     * flag to indicate ascertainment correction should be applied *
     */
    protected boolean useAscertainedSitePatterns = false;

    /**
     * dealing with proportion of site being invariant *
     */
    double proportionInvariant = 0;
    List<Integer> constantPattern = null;
    
    HashMap<Integer, BreakPoints> rootBreaks;
    
    Partials rootPartials;
    
    /**
     * Dummy node to deal with subs models requiring nodes
     */
    Node dummyNode;
    
    Map<Integer, List<BreakPoints>> passOnPointer;
    Map<Integer, List<BreakPoints>> storedPassOnPointer;
    
    Map<Integer, List<BreakPoints>> passOnRange;
    Map<Integer, List<BreakPoints>> storedPassOnRange;

    
    Map<Integer, List<Integer>> passOnEdge;
    Map<Integer, List<Integer>> storedPassOnEdge;

    Map<Integer, Double> nodeHeight;
    Map<Integer, Double> storedNodeHeight;
    
    Map<Integer, Double> nodeLogP;
    Map<Integer, Double> storedNodeLogP;
    
    List<String> taxaLabels;
    
    public NetworkLikelihood() {
    	treeInput.setRule(Validate.OPTIONAL);
    }



    @Override
    public void initAndValidate() {
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        
        int stateCount = dataInput.get().getMaxStateCount();
        
        if (stateCount == 4) {
            likelihoodCore = new BeerNetworkLikelihood4();
        } else {
        	throw new IllegalArgumentException("likelihood calculation currently only implemented for 4 states");
//            likelihoodCore = new BeerNetworkLikelihoodCore(stateCount);
        }

        String className = getClass().getSimpleName();

        Alignment alignment = dataInput.get();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        initCore();
        
        passOnPointer = new HashMap<>();
        passOnRange = new HashMap<>();
        passOnEdge = new HashMap<>();
        nodeLogP = new HashMap<>();
        
        rootPartials = new Partials(2, 1);
        
        matrixSizeExtended = (stateCount + 1) * (stateCount + 1);
        matrixSize = stateCount * stateCount;
        probabilities = new double[(stateCount + 1) * (stateCount + 1)];
        Arrays.fill(probabilities, 1.0);
        
        mat1 = new double[m_siteModel.getCategoryCount()*matrixSize];
        mat2 = new double[m_siteModel.getCategoryCount()*matrixSize];


        if (dataInput.get().isAscertained) {
            useAscertainedSitePatterns = true;
        }
    }


    protected void initCore() {
        likelihoodCore.initialize(
                dataInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities.get(),
                networkInput.get().getLeafNodes().size()
        );
        
        initPartials(networkInput.get().getEdges());
    	
        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(networkInput.get(), dataInput.get().getPatternCount());
        } else {
            setStates(networkInput.get(), dataInput.get().getPatternCount());
        }
        hasDirt = Tree.IS_FILTHY;
    }
    
    private void initPartials(Set<RecombinationNetworkEdge> edges) {
        
        // init partials
    	List<RecombinationNetworkEdge> edgesList = edges.stream()
                .collect(Collectors.toList());
    	
    	int partialLength = dataInput.get().getPatternCount() * dataInput.get().getDataType().getStateCount();
    	for (RecombinationNetworkEdge e : edgesList) {
    		likelihoodCore.initPartials(e.ID, partialLength);
    		if (hasDirt==Tree.IS_FILTHY)
    			likelihoodCore.cleanPartialsNode(e.ID);
    	}

    }

    /**
     * This method samples the sequences based on the tree and site model.
     */
    @Override
	public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
    }

    /**
     * set leaf states in likelihood core *
     */
    protected void setStates(RecombinationNetwork network, int patternCount) {
    	
        List<RecombinationNetworkNode> networkNodes = new ArrayList<>(network.getNodes());    	
    	List<RecombinationNetworkNode> leafs = networkNodes.stream()
                .filter(e -> e.isLeaf())
                .collect(Collectors.toList());
    	
    	taxaLabels = new ArrayList<>();
    	int j=0;
    	
    	Alignment data = dataInput.get();
    	for (RecombinationNetworkNode l : leafs) {            
    		taxaLabels.add(l.getTaxonLabel());
            int i;
            int[] states = new int[patternCount];
            int taxonIndex = getTaxonIndex(l.getTaxonLabel(), data);
            for (i = 0; i < patternCount; i++) {
                int code = data.getPattern(taxonIndex, i);
                int[] statesForCode = data.getDataType().getStatesForCode(code);
                if (statesForCode.length==1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.
            }
            likelihoodCore.setStates(j, states);
            j++;
    	}
    }
    
    

    /**
     *
     * @param taxon the taxon name as a string
     * @param data the alignment
     * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
     *         or -1 if the taxon is not in the alignment.
     */
    private int getTaxonIndex(String taxon, Alignment data) {
        int taxonIndex = data.getTaxonIndex(taxon);
        if (taxonIndex == -1) {
        	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
            }
            if (taxonIndex == -1) {
            	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
            }
        }
        return taxonIndex;
	}

	/**
     * set leaf partials in likelihood core *
     */
    protected void setPartials(RecombinationNetwork network, int patternCount) {
        List<RecombinationNetworkNode> networkNodes = new ArrayList<>(network.getNodes());    	
    	List<RecombinationNetworkNode> leafs = networkNodes.stream()
                .filter(e -> e.isLeaf())
                .collect(Collectors.toList());
    	
        Alignment data = dataInput.get();
    	for (RecombinationNetworkNode l : leafs) {
            int states = data.getDataType().getStateCount();
            double[] partials = new double[patternCount * states];
            int k = 0;
            int taxonIndex = getTaxonIndex(l.getTaxonLabel(), data);
            for (int patternIndex_ = 0; patternIndex_ < patternCount; patternIndex_++) {                
                double[] tipLikelihoods = data.getTipLikelihoods(taxonIndex,patternIndex_);
                if (tipLikelihoods != null) {
                	for (int state = 0; state < states; state++) {
                		partials[k++] = tipLikelihoods[state];
                	}
                }
                else {
                	int stateCount = data.getPattern(taxonIndex, patternIndex_);
	                boolean[] stateSet = data.getStateSet(stateCount);
	                for (int state = 0; state < states; state++) {
	                	 partials[k++] = (stateSet[state] ? 1.0 : 0.0);                
	                }
                }
            }    
//            System.arraycopy(partials, 0, l.partials, 0, partials.length);
        }
    }

//    // for testing
//    public double[] getRootPartials() {
//        return m_fRootPartials.clone();
//    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;
        
    
    @Override
    public double calculateLogP() {
    	   	
        final RecombinationNetwork network = networkInput.get();
        if (network.resume) {
        	initCore();
        	network.resume=false;
        }
        

        List<RecombinationNetworkNode> nodes = network.getNodes().stream().filter(e -> !e.isLeaf()).collect(Collectors.toList());
        Set<RecombinationNetworkEdge> edges = network.getEdges();
        		
        List<Integer> edgeIDs = new ArrayList<>();
        
    	nodeHeight = new HashMap<>();
        
        for (RecombinationNetworkEdge edge : edges.stream().collect(Collectors.toList())) {
        	edge.visited=false;
        	edgeIDs.add(edge.ID);
        	if (edge.childNode.isLeaf())
        		nodeHeight.put(taxaLabels.indexOf(edge.childNode.getTaxonLabel()), edge.childNode.getHeight());
        	else
        		nodeHeight.put(edge.ID, edge.childNode.getHeight());
        }
        
        likelihoodCore.cleanPartials(edgeIDs);
        likelihoodCore.cleanMatrix(edgeIDs);
    	
      	// init partials that have not yet been initialized
    	initPartials(edges);
    	    	
    	// set dummy nodes
    	for (RecombinationNetworkNode n : nodes) {
    		if (!n.isLeaf()) {
	    		n.dummy = new ArrayList<>();
	    		n.dummy2 = new ArrayList<>();
	    		n.dummy3 = new BreakPoints();
	    		n.dummy4 = new ArrayList<>();
	    		n.edgePointer = new ArrayList<>();
	    		n.prevPointer = new ArrayList<>();
	    		n.overlap = null;
	    		n.visited = false;
    		}
    	}

    	setDirty(network, edges);
    	// check where the roots of local trees are
    	rootBreaks = new HashMap<>();    
    	
    	List<Integer> dirtyRoot = new ArrayList<>();    	
    	traversalRoots(network.getRootEdge(), new BreakPoints(network.totalLength), dirtyRoot);

    	if (hasDirt==Tree.IS_FILTHY)
    		likelihoodCore.debug=true;
    	else
    		likelihoodCore.debug=false;
    	
    	try {
//    		if (hasDirt == Tree.IS_FILTHY) {
		    	for (RecombinationNetworkNode n : network.getNodes().stream().filter(e -> e.isLeaf()).collect(Collectors.toList())) {
		    		getCoalChildren(n.getParentEdges().get(0).parentNode, n.getParentEdges().get(0).breakPoints, taxaLabels.indexOf(n.getTaxonLabel()), n.getParentEdges().get(0).breakPoints);
		    	}
//    		}else {
//    			List<RecombinationNetworkEdge> es = edges.stream()
//    					.filter(e -> !e.isRootEdge())
//    					.filter(e -> !e.parentNode.dirtyBreakPoints.isEmpty())
//						.filter(e -> e.childNode.dirtyBreakPoints.isEmpty())
//    					.collect(Collectors.toList());
//    			    			
//    	        for (RecombinationNetworkEdge e : es) {
//    	        	e.visited=true;
//    	        }
//    			
//    	        for (RecombinationNetworkEdge e : es) {
//    	        	if (e.childNode.isLeaf()) {
//    		    		getCoalChildren(e.parentNode,e.breakPoints, e.ID, e.breakPoints);
//    	        	}else {
//        	        	if (passOnPointer.get(e.ID)!=null) { //start of calculations is "above" local root and doesn't affect likelihood
//			        		for (int i = 0; i < passOnPointer.get(e.ID).size(); i++) {
//			        			if (nodeHeight.get(passOnEdge.get(e.ID).get(i))==null) {
//			        				System.out.println(e.ID);
//			        				System.out.println(passOnEdge.get(e.ID));
//			        				System.exit(0);
//			        			}
//		    					getCoalChildren(e.parentNode, passOnRange.get(e.ID).get(i), passOnEdge.get(e.ID).get(i), passOnPointer.get(e.ID).get(i));   
//			        		}
//        	        	}
//    	        	}
//    	        }
//    		}

    		calcLogP(dirtyRoot);
        }catch (ArithmeticException e) {
        	return Double.NEGATIVE_INFINITY;
        }
    	
    	// TODO investigate why it was NaN
    	if (logP==Double.NaN) {
    		logP =  Double.NEGATIVE_INFINITY;
    		return logP;
    	}

        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
        	System.err.println("scaling not implementated and logP is negative Inf");
            return logP;
        }
        return logP;
    }
    
	int iiiiiii=0;
	int jjjjjjj=0;


    private void traversalRoots(RecombinationNetworkEdge edge, BreakPoints breakPoints, List<Integer> dirtyRoot) {
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
    			if (!node.dirtyBreakPoints.isEmpty() || 
    					hasDirt==Tree.IS_FILTHY ||
    					node.getChildEdges().get(0).isDirty()==Tree.IS_FILTHY ||
    					node.getChildEdges().get(1).isDirty()==Tree.IS_FILTHY) {
         			if (!dirtyRoot.contains(node.getParentEdges().get(0).ID)) {
         				dirtyRoot.add(node.getParentEdges().get(0).ID);
         			}

    			}
     			if (rootBreaks.containsKey(node.getParentEdges().get(0).ID)) {
    				rootBreaks.get(node.getParentEdges().get(0).ID).or(bp1);
    			}else {
    				rootBreaks.put(node.getParentEdges().get(0).ID, bp1);
    			}
    		}
    		
    		
    		// get which loci did not coalesce
    		BreakPoints bp2 = bp.copy();
    		bp2.andNot(bp1);
    		
    		traversalRoots(node.getChildEdges().get(0), bp2, dirtyRoot);
    		traversalRoots(node.getChildEdges().get(1), bp2, dirtyRoot);    				
    	}else if (node.isRecombination()) {
    		traversalRoots(node.getChildEdges().get(0), bp, dirtyRoot);
    	}else {
    		return;
    	}	

	}
   
    void calcLogP(List<Integer> dirtyRoot) {    	
    	logP=0;
        // No parent this is the root of the beast.tree -
        // calculate the pattern likelihoods
        final double[] frequencies = //m_pFreqs.get().
                substitutionModel.getFrequencies();

        final double[] proportions = m_siteModel.getCategoryProportions(dummyNode);

        for (Integer key : rootBreaks.keySet()) {
        	if (dirtyRoot.contains(key)) {
        		double nlp = likelihoodCore.integratePartials(proportions, frequencies, dataInput.get(), key, rootBreaks.get(key));
        		if (nodeLogP.containsKey(key)) {
        			nodeLogP.replace(key, nlp);
        		}else {
        			nodeLogP.put(key, nlp);
        		}
        		logP += nlp;
        	}else {
        		logP += nodeLogP.get(key);
        	}
        }
    }    
    
    private void computeForPatterns(BreakPoints computeFor) {
    	computeForPatterns = new boolean[dataInput.get().getPatternCount()];
    	for (int i=0; i < computeFor.size();i++) {
    		for (int j=computeFor.getRange(i).from; j<=computeFor.getRange(i).to; j++)
    			computeForPatterns[dataInput.get().getPatternIndex(j)] = true;
    	}
	}

	void getCoalChildren(RecombinationNetworkNode node, BreakPoints computeFor_BP, 
			int prev_edge_ID, BreakPoints prev_Pointer) {   
		

    	if (computeFor_BP.isEmpty())
    		return;   
      
        if (node.isRecombination()) {
        	for (RecombinationNetworkEdge edge : node.getParentEdges()) {
        		BreakPoints bp = computeFor_BP.copy();      		
        		bp.andPR(edge.passingRange); 
        		if (!bp.isEmpty()) {        		
	        		getCoalChildren(edge.parentNode, bp, prev_edge_ID, prev_Pointer);
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
                getCoalChildren(edge.parentNode, cf_only, prev_edge_ID, prev_Pointer);
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
	        		        		
	        		                if (m_siteModel.integrateAcrossCategories()) {
	        		                	computeForPatterns(bp1);		        		                	
	        		                	mat1 = getLengthMatrix(node.getHeight() - nodeHeight.get(node.prevPointer.get(i)));
	        		                	mat2 = getLengthMatrix(node.getHeight() - nodeHeight.get(node.prevPointer.get(j)));    
	        		                	
	        		        			try{
	        		        				likelihoodCore.calculatePartials(node.prevPointer.get(i), node.prevPointer.get(j),	        		        			
	        		        					edge.ID, bp1, node.dummy2.get(i), node.dummy2.get(j), 
	        		        					computeForPatterns, mat1, mat2);
	        		        			}catch (Exception e){
	        		        				System.out.println(e);
	        		        				System.out.println(networkInput.get());
	        		        				System.out.println(node.getHeight());
	        		        				System.out.println(bp1);
	        		        				System.out.println(node.dummy2);
	        		        				System.out.println(node.dummy);
	        		        				System.out.println(node.prevPointer);
	        		        				System.out.println(i + " " + j);
	        		        				System.exit(0);
	        		        			}
	        		                } else {
	        		                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
	        		                }	        		                
	        	        		}else {    
	        	    				// check label overlap
	        	                    likelihoodCore.reassignLabels(edge.ID, bp1, node.dirtyBreakPoints);
	        	        		}	        		        		        		        	
	        		        	
//                            	updateEdgeInfo(edge, bp1, edge.ID, bp1);
	                        	if (!edge.isRootEdge()) {
	                        		BreakPoints rootBp = rootBreaks.get(edge.ID);
	                        		if (rootBp!=null) {
	                        			if (!bp1.overlapFast(rootBp)) {
		                        			getCoalChildren(edge.parentNode, bp1, edge.ID, bp1);	                        				
	                        			}
	                        		}else {
	                        			getCoalChildren(edge.parentNode, bp1, edge.ID, bp1);
	                        		}
	                        	}
	        				}
	        			}
        			}
        		}
    		}

    		
		}        
    }
	
	private void updateEdgeInfo(RecombinationNetworkEdge edge, BreakPoints computeFor, int ID, BreakPoints prevPointer) {
		if (!passOnRange.containsKey(edge.ID)) {
    		passOnRange.put(edge.ID, new ArrayList<>());
    		passOnPointer.put(edge.ID, new ArrayList<>());
    		passOnEdge.put(edge.ID, new ArrayList<>());
    		
	    	passOnRange.get(edge.ID).add(computeFor.copy());
	    	passOnEdge.get(edge.ID).add(ID);
	    	passOnPointer.get(edge.ID).add(prevPointer.copy());
	    	
	    	edge.visited=true;	    
	    	
		}else {
	    	if (!edge.visited) {
	    		passOnRange.replace(edge.ID, new ArrayList<>());
	    		passOnPointer.replace(edge.ID, new ArrayList<>());
	    		passOnEdge.replace(edge.ID, new ArrayList<>());
	    	}	    	
    	
	    	passOnRange.get(edge.ID).add(computeFor.copy());
	    	passOnEdge.get(edge.ID).add(ID);
	    	passOnPointer.get(edge.ID).add(prevPointer.copy());
	    	
	    	edge.visited=true; 	
		}
		

	}
  	
	void setDirty(RecombinationNetwork network, Set<RecombinationNetworkEdge> edges) {
		if (hasDirt==Tree.IS_FILTHY) {
	    	// check which edges and break points need recomputation
	    	for (RecombinationNetworkEdge e : edges.stream().collect(Collectors.toList())) {
	    		e.makeDirty(Tree.IS_FILTHY);
	    	}
	    	return;
		}			
		
    	// check which edges and break points need recomputation
    	for (RecombinationNetworkEdge e : edges.stream().filter(e -> e.isDirty()==Tree.IS_FILTHY).collect(Collectors.toList())) {
    		upwardsTraversalDirtyEdges(e);
    	}
	}	
   
	void upwardsTraversalDirtyEdges(RecombinationNetworkEdge edge) {
    	if (edge.isRootEdge())        	
        	return;   
    	
   		edge.parentNode.dirtyBreakPoints = new BreakPoints(networkInput.get().totalLength);
    	edge.makeDirty(Tree.IS_FILTHY);

//    	
    	for (RecombinationNetworkEdge e : edge.parentNode.getParentEdges()) {  
    		if (e.isDirty()!=Tree.IS_FILTHY)
    			upwardsTraversalDirtyEdges(e);        	
    	}       	

    }
	
	private double[] getLengthMatrix(double length) {
		double[] mat = new double[mat1.length];
	    final double branchRate = branchRateModel.getRateForBranch(dummyNode);
	    for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
	    	
	        final double jointBranchRate = m_siteModel.getRateForCategory(i, dummyNode) * branchRate;
	        substitutionModel.getTransitionProbabilities(dummyNode, length, 0, jointBranchRate, probabilities);
	        System.arraycopy(probabilities, 0, mat, i * matrixSize, matrixSize);
	    }	   
	    return mat;
	}
   
//    /* return copy of pattern log likelihoods for each of the patterns in the alignment */
//	public double [] getPatternLogLikelihoods() {
////		if (beagle != null) {
////			return beagle.getPatternLogLikelihoods();
////		}
//		return patternLogLikelihoods.clone();
//	} // getPatternLogLikelihoods

    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
//        if (beagle != null) {
//            return beagle.requiresRecalculation();
//        }
        hasDirt = Tree.IS_CLEAN;

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        
        return networkInput.get().somethingIsDirty();
    }

    @Override
    public void store() {
//		if (iiiiii>300000) {
//	        System.out.println(passOnPointer.containsKey(1040018857));
//	        System.out.println("store");
//		}
//        if (beagle != null) {
//            beagle.store();
//            super.store();
//            return;
//        }
        if (likelihoodCore != null) {
            likelihoodCore.store();
        }
        
        
        storedNodeHeight = new HashMap<>();
        storedNodeLogP = new HashMap<>();
        for (Integer key : nodeHeight.keySet()) {
        	storedNodeHeight.put(key, nodeHeight.get(key));
        	if (nodeLogP.containsKey(key)) {
        		storedNodeLogP.put(key, nodeLogP.get(key));
        	}
        }

        
//        storedPassOnPointer = new HashMap<>();
//        storedPassOnRange = new HashMap<>();
//        storedPassOnEdge = new HashMap<>();
//        
//        for (Integer key : passOnPointer.keySet()) {
//        	if (storedNodeHeight.containsKey(key)) {
//	        	storedPassOnPointer.put(key, new ArrayList<>());
//	        	storedPassOnRange.put(key, new ArrayList<>());
//	        	storedPassOnEdge.put(key, new ArrayList<>(passOnEdge.get(key)));
//	        	
//	        	for (BreakPoints bp : passOnPointer.get(key))
//	        		storedPassOnPointer.get(key).add(bp.copy());
//	        	
//	        	for (BreakPoints bp : passOnRange.get(key))
//	        		storedPassOnRange.get(key).add(bp.copy());
//        	}
//       }
        
//        System.out.println(passOnEdge.get(-487675077));
//        System.out.println(nodeHeight.get(1593237398));
//        System.out.println(storedPassOnEdge.get(1591136059));
        

        super.store();
        
//        System.out.println("store");
//        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
    }

    @Override
    public void restore() {
//		if (iiiiii>300000) 
//			System.out.println("restore");
//        if (beagle != null) {
//            beagle.restore();
//            super.restore();
//            return;
//        }
        if (likelihoodCore != null) {
            likelihoodCore.restore();
//            networkInput.get().restore();
        }
        
//        Map<Integer, List<BreakPoints>> tmp = passOnPointer;
//        passOnPointer = storedPassOnPointer;
//        storedPassOnPointer = tmp;
//
//        Map<Integer, List<BreakPoints>> tmp2 = passOnRange;
//        passOnRange = storedPassOnRange;
//        storedPassOnRange = tmp2;
//
//        Map<Integer, List<Integer>> tmp3 = passOnEdge;
//        passOnEdge = storedPassOnEdge;
//        storedPassOnEdge = tmp3;
//        
        Map<Integer, Double> tmp4 = nodeHeight;
        nodeHeight = storedNodeHeight;
        storedNodeHeight = tmp4;
//        
        Map<Integer, Double> tmp5 = nodeLogP;
        nodeLogP = storedNodeLogP;
        storedNodeLogP = tmp5;

        
        super.restore();
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    @Override
	public List<String> getArguments() {
        return Collections.singletonList(dataInput.get().getID());
    }

    /**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
    @Override
	public List<String> getConditions() {
        return m_siteModel.getConditions();
    }

} // class TreeLikelihood
