package recombination.likelihood;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import beast.base.inference.Logger;
import beast.base.evolution.alignment.Alignment;
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

        if (ID1>=0) {
            if (ID2>=0) {
                calculateStatesStatesPruning(ID1,ID2,ID3,computeFor,compute1,compute2,computeForPatterns, matrices1, matrices2);
            } else {
                calculateStatesPartialsPruning(ID1,ID2,ID3,computeFor,compute1,compute2,computeForPatterns, matrices1, matrices2);
            }
        } else {
            if (ID2>=0) {
                calculateStatesPartialsPruning(ID2,ID1,ID3,computeFor,compute2,compute1,computeForPatterns, matrices2, matrices1);
            } else {
                calculatePartialsPartialsPruning(ID1,ID2,ID3,computeFor,compute1,compute2,computeForPatterns, matrices1, matrices2);
            }
        }        
        ensureLables(ID3, computeFor);

    }
    
	@Override
	public double integratePartials(double[] proportions, double[] frequencies, Alignment data, Integer i, BreakPoints rootBreaks) {
        double logP = 0;
        int[] nrInPattern;
        double[] outPartials = new double[nrOfStates*nrOfPatterns];
        
    	for (BreakPoints bp : partials.getBreaks(i)) {
    		if (!bp.isEmpty()) {
        		BreakPoints bp1 = bp.copy();
        		bp1.and(rootBreaks);
        		
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
//        // Rather than copying the stored stuff back, just swap the pointers...
//        HashMap<Integer, double[]> tmp1 = matrix;
//        matrix = storedMatrix;
//        storedMatrix = tmp1;
        
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
//    	storeMatrix();
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
	public void setStates(int j, int[] states) {
        int[] newstates = new int[states.length];
		System.arraycopy(states, 0, newstates, 0, states.length);
		this.states.put(j, newstates);
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
					if (computeFor.overlapFast(bp)) {
						BreakPoints cp2 = bp.andNotCopy(computeFor);
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
		
		for (BreakPoints bp : breaks) {
			if (!bp.isEmpty()) {
				if (!computeFor.equals(bp)) {
					BreakPoints cp_tmp = computeFor.copy();
					cp_tmp.andNot(bp);
					if (cp_tmp.isEmpty()) {
				        partials.replaceBreaks(ID, bp, computeFor.copy());
					}			
				}
			}
		}

		return 0;
	}


	@Override
	protected void dummy(Integer ID, BreakPoints bp) {
		List<BreakPoints> breaks = partials.getBreaks(ID);
		System.out.println(breaks);		
	}

	

} // class BeerLikelihoodCore
