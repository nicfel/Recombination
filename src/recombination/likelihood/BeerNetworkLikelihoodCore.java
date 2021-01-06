package recombination.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import beast.core.Logger;
import beast.evolution.alignment.Alignment;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;
import recombination.util.Partials;

/**
 * standard likelihood core, uses no caching *
 */
public class BeerNetworkLikelihoodCore extends NetworkLikelihoodCore {
    protected int nrOfStates;
    protected int nrOfNodes;
    protected int nrOfPatterns;
    protected int partialsSize;
    protected int matrixSize;
    protected int nrOfMatrices;

    protected boolean integrateCategories;


    protected int[] currentMatrixIndex;
    protected int[] storedMatrixIndex;
    protected int[] currentPartialsIndex;
    protected int[] storedPartialsIndex;
    
    protected Partials partials;

    protected HashMap<Integer, double[]> matrix;
    protected HashMap<Integer, double[]> storedMatrix;
    
    
    public HashMap<Integer, int[]> states;

    protected boolean useScaling = false;

    protected double[][][] scalingFactors;

    private double scalingThreshold = 1.0E-100;
    double SCALE = 2;

    public BeerNetworkLikelihoodCore(int nrOfStates) {
        this.nrOfStates = nrOfStates;
        matrix = new HashMap<>();
        states = new HashMap<>();        
    } // c'tor


    /**
     * Calculates partial likelihoods at a node when both children have exactly known states (e.g. for leaves).
     */
    protected void calculateStatesStatesPruning(Integer ID1, Integer ID2, Integer ID3,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
    		boolean[] computeForPatterns, double[] matrices1, double[] matrices2) {
            }
    
    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */
    protected void calculateStatesPartialsPruning(Integer ID1, Integer ID2, Integer ID3,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
    		boolean[] computeForPatterns, double[] matrices1, double[] matrices2) { 
    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected void calculatePartialsPartialsPruning(Integer ID1, Integer ID2, Integer ID3,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
    		boolean[] computeForPatterns, double[] matrices1, double[] matrices2) {
    }

    
    /**
     * Integrates partials across categories.
     *
     * @param inPartials  the array of partials to be integrated
     * @param proportions the proportions of sites in each category
     * @param outPartials an array into which the partials will go
     */
    @Override
	protected double calculateIntegratePartials(double[] proportions, double[] frequencies, Alignment data, HashMap<Integer, BreakPoints> rootBreaks) {
    }

    
   
    /**
     * Calculates pattern log likelihoods at a node.
     *
     * @param partials          the partials used to calculate the likelihoods
     * @param frequencies       an array of state frequencies
     * @param outLogLikelihoods an array into which the likelihoods will go
     */
    @Override
	public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
        int v = 0;
        for (int k = 0; k < outLogLikelihoods.length; k++) {
            double sum = 0.0;
            for (int i = 0; i < nrOfStates; i++) {

                sum += frequencies[i] * partials[v];
                v++;
            }
            outLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
        }
    }


    /**
     * initializes partial likelihood arrays.
     *
     * @param nodeCount           the number of nodes in the tree
     * @param patternCount        the number of patterns
     * @param matrixCount         the number of matrices (i.e., number of categories)
     * @param integrateCategories whether sites are being integrated over all matrices
     */
    @Override
	public void initialize(int patternCount, int matrixCount, boolean integrateCategories, boolean useAmbiguities, int nrOfNodes) {
        this.nrOfPatterns = patternCount;
        this.nrOfMatrices = matrixCount;
        this.integrateCategories = integrateCategories;
        matrixSize = nrOfStates * nrOfStates;
        partials = new Partials(nrOfNodes+50,nrOfMatrices*nrOfPatterns*nrOfStates);
    }

    /**
     * cleans up and deallocates arrays.
     */
    @Override
	public void finalize() throws java.lang.Throwable {
        nrOfPatterns = 0;
        scalingFactors = null;
    }

    @Override
    public void setUseScaling(double scale) {
        useScaling = (scale != 1.0);

        if (useScaling) {
            scalingFactors = new double[2][nrOfNodes][nrOfPatterns];
        }
    }


    /**
     * Calculates partial likelihoods at a node.
     *
     * @param nodeIndex1 the 'child 1' node
     * @param nodeIndex2 the 'child 2' node
     * @param nodeIndex3 the 'parent' node
     */
    @Override
	public void calculatePartials(Integer ID1, Integer ID2, Integer ID3, 
			BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
			boolean[] computeForPatterns, double[] matrices1, double[] matrices2) {

        if (states.containsKey(ID1)) {
            if (states.containsKey(ID2)) {
                calculateStatesStatesPruning(ID1,ID2,ID3,computeFor,compute1,compute2,computeForPatterns, matrices1, matrices2);
            } else {
                calculateStatesPartialsPruning(ID1,ID2,ID3,computeFor,compute1,compute2,computeForPatterns, matrices1, matrices2);
            }
        } else {
            if (states.containsKey(ID2)) {
                calculateStatesPartialsPruning(ID2,ID1,ID3,computeFor,compute2,compute1,computeForPatterns, matrices2, matrices1);
            } else {
                calculatePartialsPartialsPruning(ID1,ID2,ID3,computeFor,compute1,compute2,computeForPatterns, matrices1, matrices2);
            }
        }
        
        ensureLables(ID3, computeFor);

    }
    
	@Override
	public double integratePartials(double[] proportions, double[] frequencies, Alignment data, HashMap<Integer, BreakPoints> rootBreaks) {
        double logP = 0;
        int[] nrInPattern;
        double[] outPartials = new double[nrOfStates*nrOfPatterns];
        
        for (Integer i : rootBreaks.keySet()) {
        	for (BreakPoints bp : partials.getBreaks(i)) {
        		if (!bp.isEmpty()) {
	        		BreakPoints bp1 = bp.copy();
	        		bp1.and(rootBreaks.get(i));
	        		
	        		if (!bp1.isEmpty()) {
	        			
	        	        nrInPattern = new int[nrOfPatterns];      
	        			computeInPatterns(nrInPattern, data, bp1);
	        			
	                    double[] inPartials = partials.getPartials(i, bp);                    
	                    
	                    int u = 0;
	                    int v = 0;
	                    for (int k = 0; k < nrOfPatterns; k++) {
	                    	if (nrInPattern[k]!=0) {
		                        for (int s = 0; s < nrOfStates; s++) {
	
		                            outPartials[u] = inPartials[v] * proportions[0];
		                            u++;
		                            v++;
		                        }
	                    	}else {
	                    		u+=nrOfStates;
	                    		v+=nrOfStates;
	                    	}
	                    }


	                    for (int l = 1; l < nrOfMatrices; l++) {
	                        u = 0;
	                        for (int k = 0; k < nrOfPatterns; k++) {
	                        	if (nrInPattern[k]!=0) {
		                            for (int s = 0; s < nrOfStates; s++) {
	
		                                outPartials[u] += inPartials[v] * proportions[l];
		                                u++;
		                                v++;
		                            }
	                        	}else {
		                    		u+=nrOfStates;
		                    		v+=nrOfStates;
	                        	}
	                        }
	                    }
	                    u = 0;
                        for (int k = 0; k < nrOfPatterns; k++) {
                        	if (nrInPattern[k]!=0) {
                        		double sum = 0;
                                for (int s = 0; s < nrOfStates; s++) {
                                    sum += frequencies[s] * outPartials[u];
                                    u++;
                                }
                                logP += Math.log(sum)*nrInPattern[k];
                        	}else {
                        		u+=nrOfStates;
                        	}
                        }	                    
	        		}   	        		
        		}         		
        	}        	
        }  
        return logP;
    }
	
	
    private void computeInPatterns(int[] nrInPattern, Alignment data, BreakPoints computeFor) {
        for (int j = 0; j < computeFor.size();j++) {
        	for (int m = computeFor.getRange(j).from; m <= computeFor.getRange(j).to; m++) {
        		nrInPattern[data.getPatternIndex(m)]++;
        	}			
		}
	}



    /**
     * Scale the partials at a given node. This uses a scaling suggested by Ziheng Yang in
     * Yang (2000) J. Mol. Evol. 51: 423-432
     * <p/>
     * This function looks over the partial likelihoods for each state at each pattern
     * and finds the largest. If this is less than the scalingThreshold (currently set
     * to 1E-40) then it rescales the partials for that pattern by dividing by this number
     * (i.e., normalizing to between 0, 1). It then stores the log of this scaling.
     * This is called for every internal node after the partials are calculated so provides
     * most of the performance hit. Ziheng suggests only doing this on a proportion of nodes
     * but this sounded like a headache to organize (and he doesn't use the threshold idea
     * which improves the performance quite a bit).
     *
     * @param nodeIndex
     */
    protected void scalePartials(int nodeIndex) {
        int u = 0;

        for (int i = 0; i < nrOfPatterns; i++) {

            double scaleFactor = 0.0;
            int v = u;
            for (int k = 0; k < nrOfMatrices; k++) {
                for (int j = 0; j < nrOfStates; j++) {
                    if (partials[currentPartialsIndex[nodeIndex]][nodeIndex][v] > scaleFactor) {
                        scaleFactor = partials[currentPartialsIndex[nodeIndex]][nodeIndex][v];
                    }
                    v++;
                }
                v += (nrOfPatterns - 1) * nrOfStates;
            }

            if (scaleFactor < scalingThreshold) {

                v = u;
                for (int k = 0; k < nrOfMatrices; k++) {
                    for (int j = 0; j < nrOfStates; j++) {
                        partials[currentPartialsIndex[nodeIndex]][nodeIndex][v] /= scaleFactor;
                        v++;
                    }
                    v += (nrOfPatterns - 1) * nrOfStates;
                }
                scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][i] = Math.log(scaleFactor);

            } else {
                scalingFactors[currentPartialsIndex[nodeIndex]][nodeIndex][i] = 0.0;
            }
            u += nrOfStates;


        }
    }

    
    
    /**
     * This function returns the scaling factor for that pattern by summing over
     * the log scalings used at each node. If scaling is off then this just returns
     * a 0.
     *
     * @return the log scaling factor
     */
    @Override
	public double getLogScalingFactor(int patternIndex_) {
//    	if (m_bUseScaling) {
//    		return -(m_nNodeCount/2) * Math.log(SCALE);
//    	} else {
//    		return 0;
//    	}        
        double logScalingFactor = 0.0;
        if (useScaling) {
            for (int i = 0; i < nrOfNodes; i++) {
                logScalingFactor += scalingFactors[currentPartialsIndex[i]][i][patternIndex_];
            }
        }
        return logScalingFactor;
    }


    /**
    /**
     * Store current state
     */
    @Override
    public void restore() {
        // Rather than copying the stored stuff back, just swap the pointers...
        HashMap<Integer, double[]> tmp1 = matrix;
        matrix = storedMatrix;
        storedMatrix = tmp1;
        
        partials.restore();

    }

    @Override
	public void unstore() {
//        System.arraycopy(storedMatrixIndex, 0, currentMatrixIndex, 0, nrOfNodes);
//        System.arraycopy(storedPartialsIndex, 0, currentPartialsIndex, 0, nrOfNodes);
    }

    /**
     * Restore the stored state
     */
    @Override
    public void store() {
    	storeMatrix();
    	partials.store();
    }
    
    private void storeMatrix() {
    	storedMatrix = new HashMap<>();
    	for (Integer e : matrix.keySet()) {
    		double[] oldmat = matrix.get(e);
            double[] newmat = new double[oldmat.length];
    		System.arraycopy(oldmat, 0, newmat, 0, oldmat.length);
    		storedMatrix.put(e, newmat);
    	}
    }
    
    
   	@Override
	public void setEdgeMatrix(RecombinationNetworkEdge edge, int matrixIndex, double[] matrix) {
		if (this.matrix.containsKey(edge.ID)) {
	        double[] newmat = this.matrix.get(edge.ID);
	        System.arraycopy(matrix, 0, newmat,
	        		matrixIndex * matrixSize, matrixSize);
		}else {
	        double[] newmat = new double[nrOfMatrices*matrixSize];
	        System.arraycopy(matrix, 0, newmat,
	        		matrixIndex * matrixSize, matrixSize);
			this.matrix.put(edge.ID, newmat);
		}		
	}
	
	@Override
    public void initPartials(Integer ID, int length) {
		partials.addNewNode(ID);
	}
	
	@Override
	public void setStates(RecombinationNetworkEdge edge, int[] states) {
        int[] newstates = new int[states.length];
		System.arraycopy(states, 0, newstates, 0, states.length);
		this.states.put(edge.ID, newstates);
	}


	@Override
	public void cleanMatrix(List<Integer> edgesIDs) {
		List<Integer> remove = new ArrayList<>();
		for (Integer e : matrix.keySet()) {
			if (!edgesIDs.contains(e))
				remove.add(e);
		}
		for (Integer e : remove)
			matrix.remove(e);

	}


	@Override
	public void cleanPartials(List<Integer> nodeIDs) {
				
		List<Integer> remove = new ArrayList<>();
		for (Integer n : partials.keySet()) {
			if (!nodeIDs.contains(n))
				remove.add(n);
		}
		for (Integer n : remove)
			partials.remove(n);
	}


	@Override
	public void rapidStore() {
		// TODO Auto-generated method stub
		
	}


	@Override
	protected void cleanPartialsNode(Integer ID) {
		partials.removeBreaks(ID);
	}

	protected void ensureLables(Integer ID, BreakPoints computeFor) {

		List<BreakPoints> breaks = partials.getBreaks(ID);
		
		for (BreakPoints bp : breaks) {
			if (!bp.isEmpty()) {
				if (!computeFor.equals(bp)) {
					BreakPoints cp = computeFor.copy();
					cp.and(bp);
					if (!cp.isEmpty()) {
						BreakPoints cp2 = bp.copy();
						cp2.andNot(computeFor);
						partials.replaceBreaks(ID, bp, cp2);
					}				
				}
			}
		}	
		
	}
	
	@Override
	protected void checkLabels(Integer ID, BreakPoints computeFor) {

		List<BreakPoints> breaks = partials.getBreaks(ID);
		
        if (breaks.contains(computeFor))
			return;

		for (BreakPoints bp : breaks) {
			if (!bp.isEmpty()) {
				if (!computeFor.equals(bp)) {
					BreakPoints cp = computeFor.copy();
					cp.andNot(bp);
					if (cp.isEmpty()) {
						partials.replaceBreaks(ID, bp, computeFor.copy());
	
					}			
				}
			}
		}
		ensureLables(ID, computeFor);

		return;

	}
	
	@Override
	protected int reassignLabels(Integer ID, BreakPoints computeFor, BreakPoints dirtyEdges) {

		List<BreakPoints> breaks = partials.getBreaks(ID);
		
		if (breaks.contains(computeFor))
			return 0;

		BreakPoints cp = computeFor.copy();
		cp.andNot(dirtyEdges);
		
		if (cp.isEmpty()) {
			throw new IllegalArgumentException("caching issue with dirty edges");
		}

		
		for (BreakPoints bp : breaks) {
			if (!bp.isEmpty()) {
				if (!computeFor.equals(bp)) {
					BreakPoints cp_tmp = cp.copy();
					cp_tmp.andNot(bp);
					if (cp_tmp.isEmpty()) {
				        partials.replaceBreaks(ID, bp, computeFor.copy());
					}			
				}
			}
		}
		ensureLables(ID, computeFor);

		return 0;
	}


	@Override
	protected void dummy(Integer ID, BreakPoints bp) {
		List<BreakPoints> breaks = partials.getBreaks(ID);
		System.out.println(breaks);		
	}

	

} // class BeerLikelihoodCore
