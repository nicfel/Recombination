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
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.util.*;
import java.util.stream.Collectors;

@Description("Calculates the probability of sequence data on a beast.tree given a site and substitution model using " +
        "a variant of the 'peeling algorithm'. For details, see" +
        "Felsenstein, Joseph (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. J Mol Evol 17 (6): 368-376.")
public class NetworkLikelihood extends GenericNetworkLikelihood {

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

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    /**
     * memory allocation for likelihoods for each of the patterns *
     */
    protected double[] patternLogLikelihoods;
    /**
     * memory allocation for the root partials *
     */
    protected double[] m_fRootPartials;
    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    protected double[] probabilities;

    protected int matrixSizeExtended;
    protected int matrixSize;

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
    
    /**
     * Dummy node to deal with subs models requiring nodes
     */
    Node dummyNode;

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
            likelihoodCore = new BeerNetworkLikelihoodCore(stateCount);
        }

        String className = getClass().getSimpleName();

        Alignment alignment = dataInput.get();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
//        if (proportionInvariant > 0) {
//            calcConstantPatternIndices(patterns, stateCount);
//        }

        initCore();
        
        patternLogLikelihoods = new double[dataInput.get().getSiteCount()];
        m_fRootPartials = new double[dataInput.get().getSiteCount() * stateCount];
        matrixSizeExtended = (stateCount + 1) * (stateCount + 1);
        matrixSize = stateCount * stateCount;
        probabilities = new double[(stateCount + 1) * (stateCount + 1)];
        Arrays.fill(probabilities, 1.0);

        if (dataInput.get().isAscertained) {
            useAscertainedSitePatterns = true;
        }
    }


//    /**
//     * Determine indices of m_fRootProbabilities that need to be updates
//     * // due to sites being invariant. If none of the sites are invariant,
//     * // the 'site invariant' category does not contribute anything to the
//     * // root probability. If the site IS invariant for a certain character,
//     * // taking ambiguities in account, there is a contribution of 1 from
//     * // the 'site invariant' category.
//     */
//    void calcConstantPatternIndices(final int patterns, final int stateCount) {
//        constantPattern = new ArrayList<>();
//        for (int i = 0; i < patterns; i++) {
//            final int[] pattern = dataInput.get().getPattern(i);
//            final boolean[] isInvariant = new boolean[stateCount];
//            Arrays.fill(isInvariant, true);
//            for (final int state : pattern) {
//                final boolean[] isStateSet = dataInput.get().getStateSet(state);
//                if (m_useAmbiguities.get() || !dataInput.get().getDataType().isAmbiguousCode(state)) {
//                    for (int k = 0; k < stateCount; k++) {
//                        isInvariant[k] &= isStateSet[k];
//                    }
//                }
//            }
//            for (int k = 0; k < stateCount; k++) {
//                if (isInvariant[k]) {
//                    constantPattern.add(i * stateCount + k);
//                }
//            }
//        }
//    }
//
    protected void initCore() {
        likelihoodCore.initialize(
                dataInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true, m_useAmbiguities.get(),
                networkInput.get().getLeafNodes().size()
        );
        
        initPartials();
    	


        if (m_useAmbiguities.get() || m_useTipLikelihoods.get()) {
            setPartials(networkInput.get(), dataInput.get().getPatternCount());
        } else {
            setStates(networkInput.get(), dataInput.get().getPatternCount());
        }
        hasDirt = Tree.IS_FILTHY;
    }
    
    private void initPartials() {
        
        // init partials
        List<RecombinationNetworkNode> networkNodes = new ArrayList<>(networkInput.get().getNodes());    	
    	List<RecombinationNetworkNode> nodes = networkNodes.stream()
    			.filter(e -> !e.isLeaf())
                .collect(Collectors.toList());
    	
    	int partialLength = dataInput.get().getPatternCount() * dataInput.get().getDataType().getStateCount();
    	for (RecombinationNetworkNode n : nodes) {
    		likelihoodCore.initPartials(n, partialLength);
    		if (hasDirt==Tree.IS_FILTHY)
    			likelihoodCore.cleanPartialsNode(n);
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
    	
    	Alignment data = dataInput.get();
    	for (RecombinationNetworkNode l : leafs) {            
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
            likelihoodCore.setStates(l, states);
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

    // for testing
    public double[] getRootPartials() {
        return m_fRootPartials.clone();
    }

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
        
        List<RecombinationNetworkNode> nodes = network.getNodes().stream().filter(e -> !e.isLeaf()).collect(Collectors.toList());
        List<RecombinationNetworkEdge> edges = network.getEdges().stream().collect(Collectors.toList());
        		
        List<Integer> nodeIDs = new ArrayList<>();
        List<Integer> edgeIDs = new ArrayList<>();
        
        
        for (RecombinationNetworkNode node : nodes) 
        	nodeIDs.add(node.ID);
        for (RecombinationNetworkEdge edge : edges) 
        	edgeIDs.add(edge.ID);
        
        likelihoodCore.cleanPartials(nodeIDs);
        likelihoodCore.cleanMatrix(edgeIDs);

        for (RecombinationNetworkEdge e : edges) 
    		e.visited = false;        
    	
      	// init partials that have not yet been initialized
    	initPartials();
    	
    	// set dummy nodes
    	for (RecombinationNetworkNode n : network.getNodes().stream().filter(e -> !e.isLeaf()).collect(Collectors.toList())) {
    		n.dummy = new ArrayList<>();
    		n.dummy2 = new ArrayList<>();
    		n.prevLength = new ArrayList<>();
    		n.prevPointer = new ArrayList<>();
    	}
    	
    	// check where the roots of local trees are
    	rootBreaks = new HashMap<>();    	
    	traversalRoots(network.getRootEdge(), new BreakPoints(network.totalLength));
    	
    	
    	setDirty(network);
    	
//    	System.out.println(network);
//    	System.exit(0);

    	if (hasDirt==Tree.IS_FILTHY)
    		likelihoodCore.debug=true;
    	else
    		likelihoodCore.debug=false;

    	try {
	    	for (RecombinationNetworkNode n : network.getNodes().stream().filter(e -> e.isLeaf()).collect(Collectors.toList())) {
	    		upwardsTraversalBP(n, n.getParentEdges().get(0).breakPoints, null, n.getParentEdges().get(0).breakPoints, 0.0);
	    	}

    		calcLogP(network.getRootEdge());

        }catch (ArithmeticException e) {
        	return Double.NEGATIVE_INFINITY;
        }

        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
        	System.err.println("scaling not implementated and logP is negative Inf");
            return logP;
        }
        return logP;
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
    			if (rootBreaks.containsKey(node.ID)) {
    				rootBreaks.get(node.ID).or(bp1);
    			}else {
    				rootBreaks.put(node.ID, bp1);
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


	void calcLogP() {
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
            }
        }
    }
    
    void calcLogP(RecombinationNetworkEdge edge) {
    	
        if (edge.isRootEdge()) {
            // No parent this is the root of the beast.tree -
            // calculate the pattern likelihoods
            final double[] frequencies = //m_pFreqs.get().
                    substitutionModel.getFrequencies();

            final double[] proportions = m_siteModel.getCategoryProportions(dummyNode);
            likelihoodCore.integratePartials(proportions, m_fRootPartials, dataInput.get(), rootBreaks);
            

            if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                proportionInvariant = m_siteModel.getProportionInvariant();
                // some portion of sites is invariant, so adjust root partials for this
                for (final int i : constantPattern) {
                    m_fRootPartials[i] += proportionInvariant;
                }
            }
            likelihoodCore.calculateLogLikelihoods(m_fRootPartials, frequencies, patternLogLikelihoods);
        }    	
    	
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < patternLogLikelihoods.length; i++) {
                logP += patternLogLikelihoods[i];
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

	void upwardsTraversalBP(RecombinationNetworkNode node, BreakPoints computeFor_BP, 
			RecombinationNetworkEdge prev_edge, BreakPoints prev_Pointer, double old_time) {    
		
    	if (computeFor_BP.isEmpty())
    		return;    	
      
        if (node.isLeaf()) {        	
        	upwardsTraversalBP(node.getParentEdges().get(0).parentNode, computeFor_BP, node.getParentEdges().get(0), computeFor_BP, old_time + node.getParentEdges().get(0).getLength());
        }else if (node.isRecombination()) {
        	BreakPoints computeFor = computeFor_BP.copy();
        	computeFor.and(node.dirtyBreakPoints); 
       	
        	for (RecombinationNetworkEdge edge : node.getParentEdges()) {
        		BreakPoints bp = computeFor_BP.copy();
        		bp.and(edge.breakPoints);                
            	upwardsTraversalBP(edge.parentNode, bp, prev_edge, prev_Pointer, old_time + edge.getLength());
        	}
        	
        }else {
        	// make a copy of the BP's
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
                	upwardsTraversalBP(edge.parentNode, cf_only, prev_edge, prev_Pointer, old_time + edge.getLength());
            	}
                // see "how" much is left of the compute for BP
    			computeFor.andNot(cf_only);
    		}

    		BreakPoints bp_in = computeFor.copy();
    		for (int i = 0; i < node.dummy.size();i++) {
        		BreakPoints bp_here = node.dummy.get(i).copy();
        		// get the overlap
        		bp_here.and(bp_in);       		
        		if (!bp_here.isEmpty()) {     
//                	BreakPoints cp = bp_here.copy();
//                	cp.and(node.dirtyBreakPoints);        
                	
                	
//	        		if (cp.equals(bp_here)) {
                	if (node.getChildEdges().get(0).isDirty()==Tree.IS_FILTHY ||
                			node.getChildEdges().get(1).isDirty()==Tree.IS_FILTHY) {
	        			
		                if (m_siteModel.integrateAcrossCategories()) {
//		                	System.out.println(node.getHeight());
		                	computeForPatterns(bp_here);		  
		        			double[] mat1 = getLengthMatrix(old_time);
		        			double[] mat2 = getLengthMatrix(node.prevLength.get(i));     	
	        			
		        			likelihoodCore.calculatePartials(prev_edge, node.prevPointer.get(i), 
		        					node, bp_here, prev_Pointer, node.dummy2.get(i), 
		        					computeForPatterns, mat1, mat2);
		                } else {
		                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
		                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
		                }
	        		}else {
	    				// check label overlap
	                    likelihoodCore.reassignLabels(node, bp_here, node.dirtyBreakPoints);
	        		}
	        		
	        		computeFor.andNot(bp_here);
	                if (!edge.isRootEdge()) {
	                	BreakPoints passOn = bp_here.copy();	                	
	                	// only pass on loci for which the root has not been reached yet.
	                	if (rootBreaks.containsKey(node.ID)) {
	                		passOn.andNot(rootBreaks.get(node.ID));
	                	}
	                	upwardsTraversalBP(edge.parentNode, passOn, edge, passOn, edge.getLength());
	                }

        		}  
        	}
    		// save the info for when the node is passed again.
    		node.dummy2.add(prev_Pointer.copy());
    		node.dummy.add(computeFor);
    		node.prevLength.add(old_time);
    		node.prevPointer.add(prev_edge);
		}        
    }
    
	
	void setDirty(RecombinationNetwork network) {
		if (hasDirt==Tree.IS_FILTHY) {
	    	// check which edges and break points need recomputation
	    	for (RecombinationNetworkEdge e : network.getEdges().stream().collect(Collectors.toList())) {
	    		e.makeDirty(Tree.IS_FILTHY);
	    	}    	

	    	for (RecombinationNetworkNode n : network.getNodes().stream().filter(e -> e.isCoalescence()).collect(Collectors.toList())) {
       			n.dirtyBreakPoints = new BreakPoints(networkInput.get().totalLength);
	    	}
	    	return;
		}			
		
    	// check which edges and break points need recomputation
    	for (RecombinationNetworkEdge e : network.getEdges().stream().filter(e -> e.isDirty()==Tree.IS_FILTHY).collect(Collectors.toList())) {
    		upwardsTraversalDirtyEdges(e);
    	}    	

	}	
   
	void upwardsTraversalDirtyEdges(RecombinationNetworkEdge edge) {
    	if (edge.isRootEdge())        	
        	return;    

    	edge.parentNode.dirtyBreakPoints = new BreakPoints(networkInput.get().totalLength);
    	edge.makeDirty(Tree.IS_FILTHY);

    	for (RecombinationNetworkEdge e : edge.parentNode.getParentEdges()) {  
    		if (e.isDirty()!=Tree.IS_FILTHY)
    			upwardsTraversalDirtyEdges(e);        	
    	}       	

    }
	
	private double[] getLengthMatrix(double length) {
        double[] newmat = new double[m_siteModel.getCategoryCount()*matrixSize];

	    final double branchRate = branchRateModel.getRateForBranch(dummyNode);
	    for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
	    	
	        final double jointBranchRate = m_siteModel.getRateForCategory(i, dummyNode) * branchRate;
	        substitutionModel.getTransitionProbabilities(dummyNode, length, 0, jointBranchRate, probabilities);
	        System.arraycopy(probabilities, 0, newmat, i * matrixSize, matrixSize);
	    }	    
	    return newmat;
	}


    
	private int updateEdgeMatrix(RecombinationNetworkEdge edge) {
		int updated = 0;
		if (edge.isDirty()==Tree.IS_FILTHY || hasDirt!=Tree.IS_CLEAN) {
	    	if (!edge.visited) {
	            final double branchRate = branchRateModel.getRateForBranch(dummyNode);
	            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
	            	
	                final double jointBranchRate = m_siteModel.getRateForCategory(i, dummyNode) * branchRate;
	                substitutionModel.getTransitionProbabilities(dummyNode, edge.parentNode.getHeight(), edge.childNode.getHeight(), jointBranchRate, probabilities);
	                likelihoodCore.setEdgeMatrix(edge, i, probabilities);
	            }
	            edge.visited = true;
	    	}	    	
	    	updated = 2;
		}else if (edge.isDirty()==Tree.IS_DIRTY) {
			updated = 1;
		}
	    return updated;
	}
	
	private int updateCompute(int newCompute, int compute) {
		return newCompute > compute ? newCompute : compute;
	}


    /* return copy of pattern log likelihoods for each of the patterns in the alignment */
	public double [] getPatternLogLikelihoods() {
//		if (beagle != null) {
//			return beagle.getPatternLogLikelihoods();
//		}
		return patternLogLikelihoods.clone();
	} // getPatternLogLikelihoods

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
//    	System.out.println("store");
//        if (beagle != null) {
//            beagle.store();
//            super.store();
//            return;
//        }
        if (likelihoodCore != null) {
            likelihoodCore.store();
//            networkInput.get().storeLikelihoods();
        }
        super.store();
//        System.out.println("store");
//        System.arraycopy(m_branchLengths, 0, storedBranchLengths, 0, m_branchLengths.length);
    }

    @Override
    public void restore() {
//    	System.out.println("restore");
//        if (beagle != null) {
//            beagle.restore();
//            super.restore();
//            return;
//        }
        if (likelihoodCore != null) {
            likelihoodCore.restore();
//            networkInput.get().restore();
        }
        super.restore();
//        double[] tmp = m_branchLengths;
//        m_branchLengths = storedBranchLengths;
//        storedBranchLengths = tmp;
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
