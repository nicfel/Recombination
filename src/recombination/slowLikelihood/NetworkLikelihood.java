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


package recombination.slowLikelihood;

import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import recombination.alignment.RecombinationAlignment;
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
    
    /**
     * Dummy node to deal with subs models requiring nodes
     */
    Node dummyNode;

    @Override
    public void initAndValidate() {
//        // sanity check: alignment should have same #taxa as tree
//        if (dataInput.get().getTaxonCount() != networkInput.get().getLeafNodeCount()) {
//            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
//        }
////        beagle = null;
////        beagle = new BeagleTreeLikelihood();
////        try {
////	        beagle.initByName(
////                    "data", dataInput.get(), "tree", treeInput.get(), "siteModel", siteModelInput.get(),
////                    "branchRateModel", branchRateModelInput.get(), "useAmbiguities", m_useAmbiguities.get(), 
////                    "useTipLikelihoods", m_useTipLikelihoods.get(),"scaling", scaling.get().toString());
////	        if (beagle.beagle != null) {
////	            //a Beagle instance was found, so we use it
////	            return;
////	        }
////        } catch (Exception e) {
////			// ignore
////		}
////        // No Beagle instance was found, so we use the good old java likelihood core
//        beagle = null;

//        int nodeCount = networkInput.get().getNodeCount();
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
        
//        m_branchLengths = new double[nodeCount];
//        storedBranchLengths = new double[nodeCount];

        int stateCount = dataInput.get().getMaxStateCount();
        int patterns = dataInput.get().getPatternCount();
        
//        if (stateCount == 4) {
//            likelihoodCore = new BeerNetworkLikelihoodCore4();
//        } else {
            likelihoodCore = new BeerNetworkLikelihoodCore(stateCount);
//        }

        String className = getClass().getSimpleName();

        RecombinationAlignment alignment = dataInput.get();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
//        if (proportionInvariant > 0) {
//            calcConstantPatternIndices(patterns, stateCount);
//        }

        initCore();
        
        patternLogLikelihoods = new double[patterns];
        m_fRootPartials = new double[patterns * stateCount];
        matrixSize = (stateCount + 1) * (stateCount + 1);
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
                true, m_useAmbiguities.get()
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
                .collect(Collectors.toList());
    	
    	int partialLength = dataInput.get().getPatternCount() * dataInput.get().getDataType().getStateCount();
    	for (RecombinationNetworkNode n : nodes) 
    		likelihoodCore.initPartials(n, partialLength);

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
    	
    	for (RecombinationNetworkNode l : leafs) {
            RecombinationAlignment data = dataInput.get();
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
    private int getTaxonIndex(String taxon, RecombinationAlignment data) {
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
    	
    	for (RecombinationNetworkNode l : leafs) {
            RecombinationAlignment data = dataInput.get();
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
    	
//    	System.out.println("clears");
    	

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
    	
    	    	
//    	System.out.println(network);
    	// init partials that have not yet been initialized
    	initPartials();
    	
    	// set dummy nodes
    	for (RecombinationNetworkNode n : network.getNodes().stream().filter(e -> !e.isLeaf()).collect(Collectors.toList())) {
    		n.dummy = new ArrayList<>();
    		n.computeOnwards = new ArrayList<>();
    	}


    	try {
        	for (RecombinationNetworkNode n : network.getNodes().stream().filter(e -> e.isLeaf()).collect(Collectors.toList())) {
        		upwardsTraversal(n, n.getParentEdges().get(0).breakPoints, false, null);
        	}
    		calcLogP(network.getRootEdge());

        }catch (ArithmeticException e) {
        	return Double.NEGATIVE_INFINITY;
        }
//    	System.out.println("   ");
//    	System.exit(0);

        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
//        	System.exit(0);
//            m_nScale = 0;
//            m_fScale *= 1.01;
//            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
//            likelihoodCore.setUseScaling(m_fScale);
//            likelihoodCore.unstore();
//            hasDirt = Tree.IS_FILTHY;
//            traverse(network.getRootEdge(), network.getRootEdge().breakPoints);
//            calcLogP();
            return logP;
        }
        return logP;
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
            likelihoodCore.integratePartials(edge, proportions, m_fRootPartials);

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
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
            }
        }

    }

    
    void upwardsTraversal(RecombinationNetworkNode node, BreakPoints computeFor_BP, boolean compute, RecombinationNetworkEdge prev_edge) {
    	BreakPoints computeFor = computeFor_BP.copy();
    	
    	if (computeFor.isEmpty())
    		return;
    	
      
        if (node.isLeaf()) {        	
        	RecombinationNetworkEdge edge = node.getParentEdges().get(0);
        	compute = updateEdgeMatrix(edge) ? true : compute;
        	upwardsTraversal(node.getParentEdges().get(0).parentNode, computeFor, compute, edge);
        }else if (node.isRecombination()) {
        	
        	for (RecombinationNetworkEdge edge : node.getParentEdges()) {
        		BreakPoints bp = computeFor.copy();
        		bp.and(edge.breakPoints);    
        		// update partials for this node and for the breakpoints
        		if (compute) {
	        		if (m_siteModel.integrateAcrossCategories()) {
	                		likelihoodCore.calculatePartialsRecombination(prev_edge, node, bp);
	                } else {
	                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
	                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
	                }
        		}
                
            	compute = updateEdgeMatrix(edge) ? true : compute;
            	upwardsTraversal(edge.parentNode, bp, compute, edge);
        	}
        	
        }else {
//        	if (node.getHeight()==2.169092342538857)
//        		System.out.println(node.dummy);
        	for (RecombinationNetworkEdge edge : node.getChildEdges())
        		if (edge.isDirty()!=Tree.IS_CLEAN)
        			compute=true;
        	
        	boolean edgeUpdated = false;
        	RecombinationNetworkEdge edge = node.getParentEdges().get(0);
            if (!edge.isRootEdge()) {
            	edgeUpdated = updateEdgeMatrix(edge);
            }
        	
        	// traverse the computefor part that is only on one edge onwards
        	BreakPoints overlap = node.getChildEdges().get(0).breakPoints.copy();
        	overlap.and(node.getChildEdges().get(1).breakPoints);
        	
        	
    		BreakPoints cf_only = computeFor.copy();
    		cf_only.andNot(overlap);
    		if (!cf_only.isEmpty()) {
    			if (compute) {
//    				likelihoodCore.debug = true;
//    			}
	                if (m_siteModel.integrateAcrossCategories()) {
	                    likelihoodCore.calculatePartialsRecombination(prev_edge, node, cf_only);
	                } else {
	                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
	                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
	                }       
    			}
                likelihoodCore.debug = false;
                if (!edge.isRootEdge()) {
                	compute = edgeUpdated ? true : compute;
                	upwardsTraversal(edge.parentNode, cf_only, compute, edge);
            	}

    			computeFor.andNot(cf_only);
    		}
    		
//        	compute = true;

    		BreakPoints bp_in = computeFor.copy();
    		for (int i = 0; i < node.dummy.size();i++) {
        		BreakPoints bp_here = node.dummy.get(i).copy();
        		// get the overlap
        		bp_here.and(bp_in);
        		if (!bp_here.isEmpty()) {	  

        			BreakPoints compute1 = bp_here.copy();
        			BreakPoints compute2 = bp_here.copy();
	        		
	        		compute1.and(node.getChildEdges().get(0).breakPoints);
	        		compute2.and(node.getChildEdges().get(1).breakPoints);
	        		
	        		
                	compute = node.computeOnwards.get(i) ? true : compute;

	        		
	        		if (compute) {
		                if (m_siteModel.integrateAcrossCategories()) {
		                    likelihoodCore.calculatePartials(node.getChildEdges().get(0), node.getChildEdges().get(1), node, bp_here, compute1, compute2);
		                } else {
		                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
		                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
		                }
	        		}
//	        		}else {	        			
	            			
//	        	    	if (hasDirt==Tree.IS_CLEAN && !compute && !node.isLeaf())
//	        	    		System.out.println(node.getHeight() + " " + bp_here);

//	        			System.out.println(compute1 + " " + compute2 + " " + node.computeOnwards.get(i) + " " + compute);
//	        		}  
	        		
	        		computeFor.andNot(bp_here);

	        		

	                if (!edge.isRootEdge()) {
	                	compute = edgeUpdated ? true : compute;
	                	upwardsTraversal(edge.parentNode, bp_here, compute, edge);
	                }

        		}     		
        	}
    		node.dummy.add(computeFor);
    		node.computeOnwards.add(compute);
        }
        
    }
    
	private boolean updateEdgeMatrix(RecombinationNetworkEdge edge) {
		boolean updated = false;
		if (edge.isDirty()!=Tree.IS_CLEAN || hasDirt!=Tree.IS_CLEAN) {
	    	if (!edge.visited) {
	            final double branchRate = branchRateModel.getRateForBranch(dummyNode);
	            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
	            	
	                final double jointBranchRate = m_siteModel.getRateForCategory(i, dummyNode) * branchRate;
	                substitutionModel.getTransitionProbabilities(dummyNode, edge.parentNode.getHeight(), edge.childNode.getHeight(), jointBranchRate, probabilities);    	
	                likelihoodCore.setEdgeMatrix(edge, 0, probabilities);
	            }
	            edge.visited = true;
	    	}
	    	updated = true;
		}
	    return updated;
	}


//    /* Assumes there IS a branch rate model as opposed to traverse() */
//    int traverse(final RecombinationNetworkEdge edge, BreakPoints computeFor) {
//    	
//        int update = (edge.isDirty() | (hasDirt));
//        
//        
//        if (computeFor.isEmpty())
//        	return update;
//        
//        if (!edge.isRootEdge() && (update > -1)) {
//        	if (!edge.visited) {
//	            final double branchRate = branchRateModel.getRateForBranch(dummyNode);
//	//            final double branchTime = edge.getLength() * branchRate;
//	            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
//	            	
//	                final double jointBranchRate = m_siteModel.getRateForCategory(i, dummyNode) * branchRate;
//	                substitutionModel.getTransitionProbabilities(dummyNode, edge.parentNode.getHeight(), edge.childNode.getHeight(), jointBranchRate, probabilities);
//	//            	if (edge.childNode.getHeight()==0 )
//	//            		System.out.println("update " + branchRate + " " + jointBranchRate);
//	
//	//                if (edge.isDirty()==Tree.IS_CLEAN)
//	//                	if((edge.matrixList[0] -probabilities[0])!=0.0)
//	//                		System.exit(0);
//	
//	                edge.matrixList = new double[probabilities.length];
//	        		System.arraycopy(probabilities, 0, edge.matrixList,0, probabilities.length);
//	            }
//	            edge.visited = true;
//        	}
//
//            update |= Tree.IS_DIRTY;
//        }   
////        if (!edge.isRootEdge())
////	        if (edge.childNode.getHeight()<32 && edge.childNode.getHeight()>31)
////	        	System.out.println(edge.childNode.getHeight() + " " + edge.childNode.partials[0]);
//
//        // If the node is internal, update the partial likelihoods.
//        if (edge.childNode.isCoalescence()) {
//        	BreakPoints compute1 = computeFor.copy();
//        	BreakPoints compute2 = computeFor.copy();
//        	
//        	
//        	
//            // Traverse down the two child nodes
//            final RecombinationNetworkEdge child1 = edge.childNode.getChildEdges().get(0); //Two children
//        	compute1.and(child1.breakPoints);
//        	final int update1 = traverse(child1, compute1);
//
//            final RecombinationNetworkEdge child2 = edge.childNode.getChildEdges().get(1); //Two children
//            compute2.and(child2.breakPoints);
//        	final int update2 = traverse(child2, compute2);
//        	        	
//            // If either child node was updated then update this node too
//            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {
//            	
//                update |= (update1 | update2);
//                
//                if (m_siteModel.integrateAcrossCategories()) {
//                    likelihoodCore.calculatePartials(child1, child2, edge.childNode, computeFor, compute1, compute2);
////                    if (!edge.isRootEdge())
////            	        if (edge.childNode.getHeight()<32 && edge.childNode.getHeight()>31)
////            	        	System.out.println(edge.childNode.getHeight() + " " + edge.childNode.partials[0]);
//
//                } else {
//                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
//                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
//                }
//
//                if (edge.isRootEdge()) {
//                    // No parent this is the root of the beast.tree -
//                    // calculate the pattern likelihoods
//                    final double[] frequencies = //m_pFreqs.get().
//                            substitutionModel.getFrequencies();
//
//                    final double[] proportions = m_siteModel.getCategoryProportions(dummyNode);
//                    likelihoodCore.integratePartials(edge, proportions, m_fRootPartials);
//
//                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
//                        proportionInvariant = m_siteModel.getProportionInvariant();
//                        // some portion of sites is invariant, so adjust root partials for this
//                        for (final int i : constantPattern) {
//                            m_fRootPartials[i] += proportionInvariant;
//                        }
//                    }
//                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, frequencies, patternLogLikelihoods);
//                }
//            }
//        }else if (edge.childNode.isRecombination()) {
//        	BreakPoints compute = computeFor.copy();
//
//            final RecombinationNetworkEdge child1 = edge.childNode.getChildEdges().get(0); //Two children
//            compute.and(child1.breakPoints);
//            final int update1 = traverse(child1, compute);
//            update |= update1;
//            if (update1 != Tree.IS_CLEAN) {
//                if (m_siteModel.integrateAcrossCategories()) {
//                    likelihoodCore.calculatePartialsRecombination(child1, edge.childNode, compute);
//                } else {
//                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
//                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
//                }
//            }
//        }
//        return update;
//    }

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
