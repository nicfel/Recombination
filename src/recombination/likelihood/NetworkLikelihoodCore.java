/*
* File LikelihoodCore.java
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

import java.util.HashMap;
import java.util.List;

import beast.base.evolution.alignment.Alignment;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

/**
 * The likelihood core is the class that performs the calculations
 * in the peeling algorithm (see Felsenstein, Joseph (1981).
 * Evolutionary trees from DNA sequences: a maximum likelihood approach.
 * J Mol Evol 17 (6): 368-376.). It does this by calculating the partial
 * results for all sites, node by node. The order in which the nodes
 * are visited is controlled by the TreeLikelihood. T
 * <p/>
 * In order to reuse computations of previous likelihood calculations,
 * a current state, and a stored state are maintained. Again, the TreeLikelihood
 * controls when to update from current to stored and vice versa. So, in
 * LikelihoodCore implementations, duplicates need to be kept for all partials.
 * Also, a set of indices to indicate which of the data is stored state and which
 * is current state is commonly the most efficient way to sort out which is which.
 */

abstract public class NetworkLikelihoodCore {

    /**
     * reserve memory for partials, indices and other
     * data structures required by the core *
     */
    abstract public void initialize(int patternCount, int matrixCount, boolean integrateCategories, boolean useAmbiguities, int nrOfNodes);

    /**
     * clean up after last likelihood calculation, if at all required *
     */
    @Override
	abstract public void finalize() throws java.lang.Throwable;


    /**
     * flag to indicate whether scaling should be used in the
     * likelihood calculation. Scaling can help in dealing with
     * numeric issues (underflow).
     */
    boolean m_bUseScaling = false;
    
    boolean debug = false;

    abstract public void setUseScaling(double scale);

    public boolean getUseScaling() {
        return m_bUseScaling;
    }

    /**
     * return the cumulative scaling effect. Should be zero if no scaling is used *
     */
    abstract public double getLogScalingFactor(int patternIndex_);

    /**
     * Calculate partials for node node3, with children node1 and node2Index.
     * NB Depending on whether the child nodes contain states or partials, the
     * calculation differs-*
     */
    abstract public void calculatePartials(Integer ID1, Integer ID2, Integer ID3, 
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
    		boolean[] computeForPatterns, double[] matrices1, double[] matrices2);
    
    //abstract public void calculatePartials(int node1, int node2Index, int node3, int[] matrixMap);

    /**
     * integrate partials over categories (if any). *
     */
    abstract public double integratePartials(double[] proportions, double[] frequencies, Alignment data, Integer i, BreakPoints rootBreaks);

    /**
     * calculate log likelihoods at the root of the tree,
     * using frequencies as root node distribution.
     * outLogLikelihoods contains the resulting probabilities for each of
     * the patterns *
     */
    public void processStack() {
    }

    abstract public void setEdgeMatrix(RecombinationNetworkEdge edge, int i, double[] matrix);
    
    abstract public void setStates(int j, int[] states);

    abstract public void initPartials(Integer iD, int length);

    abstract public void cleanMatrix(List<Integer> edges);

    abstract public void cleanPartials(List<Integer> nodes);

    /**
     * store current state *
     */
    abstract public void store();

    abstract public void rapidStore();
    
    /**
     * reset current state to stored state, only used when switching from non-scaled to scaled or vice versa *
     */
    abstract public void unstore();

    /**
     * restore state *
     */
    abstract public void restore();
//    /** do internal diagnosics, and suggest an alternative core if appropriate **/ 
//    abstract LikelihoodCore feelsGood();

	protected abstract void cleanPartialsNode(Integer ID);

	protected abstract void checkLabels(Integer ID, BreakPoints computeFor);

	protected abstract int reassignLabels(Integer ID, BreakPoints computeFor, BreakPoints dirtyEdges);
	
	protected abstract void dummy(Integer ID, BreakPoints bp);

	
}
