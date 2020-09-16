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
    
    protected Partials partialsNew;

    protected HashMap<Integer, double[]> matrix;
    protected HashMap<Integer, double[]> storedMatrix;
    
    protected List<Integer> touched;

    
    public HashMap<String, int[]> states;

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
    protected void calculateStatesStatesPruning(RecombinationNetworkEdge edge1, RecombinationNetworkEdge edge2, RecombinationNetworkNode node,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
    		boolean[] computeForPatterns, double[] matrices1, double[] matrices2) {
        
    	
        // compute the breakpoints that are on both edges or only on either edge
        BreakPoints joint = compute1.copy();
        BreakPoints e1 = compute1.copy();
        BreakPoints e2 = compute2.copy();
        
        joint.and(compute2);
        e1.andNot(compute2);
        e2.andNot(compute1);      
        
        double[] mat1 = matrix.get(edge1.ID);
        double[] mat2 = matrix.get(edge2.ID);
        
        double[] partials_parent = partials.get(node.ID).get(computeFor);
        int[] states_child1 = states.get(edge1.childNode.getTaxonLabel());
        int[] states_child2 = states.get(edge2.childNode.getTaxonLabel());
        
        if (partials_parent==null) {
        	 partials.get(node.ID).put(computeFor, new double[nrOfMatrices*nrOfPatterns*nrOfStates]);
        	 partials_parent = partials.get(node.ID).get(computeFor);
        }
        	
        
        
       int v = 0;
       for (int l = 0; l < nrOfMatrices; l++) {
            for (int k = 0; k < nrOfPatterns; k++) {
                int state1 = states_child1[k];
                int state2 = states_child2[k];

                int w = l * matrixSize;

                if (state1 < nrOfStates && state2 < nrOfStates) {

                    for (int i = 0; i < nrOfStates; i++) {

                    	partials_parent[v] = mat1[w + state1] * mat2[w + state2];

                        v++;
                        w += nrOfStates;
                    }

                } else if (state1 < nrOfStates) {
                    // child 2 has a gap or unknown state so treat it as unknown

                    for (int i = 0; i < nrOfStates; i++) {

                    	partials_parent[v] = mat1[w + state1];

                        v++;
                        w += nrOfStates;
                    }
                } else if (state2 < nrOfStates) {
                    // child 2 has a gap or unknown state so treat it as unknown

                    for (int i = 0; i < nrOfStates; i++) {

                    	partials_parent[v] = mat2[w + state2];

                        v++;
                        w += nrOfStates;
                    }
                } else {
                    // both children have a gap or unknown state so set partials to 1

                    for (int j = 0; j < nrOfStates; j++) {
                    	partials_parent[v] = 1.0;
                        v++;
                    }
                }
            }
        }
    }
    
    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */
    protected void calculateStatesPartialsPruning(RecombinationNetworkEdge edge1, RecombinationNetworkEdge edge2, RecombinationNetworkNode node,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
    		boolean[] computeForPatterns, double[] matrices1, double[] matrices2) { 
        
        double[] mat1 = this.matrix.get(edge1.ID);
        double[] mat2 = this.matrix.get(edge2.ID);
        
        double[] partials_parent = this.partials.get(node.ID).get(computeFor);
        int[] states_child1 = this.states.get(edge1.childNode.getTaxonLabel());
        double[] partials_child2 = this.partials.get(edge2.childNode.ID).get(compute2);
        
        
        if (partials_parent==null) {
       	 	partials.get(node.ID).put(computeFor, new double[nrOfMatrices*nrOfPatterns*nrOfStates]);
       	 	partials_parent = partials.get(node.ID).get(computeFor);
        }
        
        
        double sum, tmp;

        int u = 0;
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {
            for (int k = 0; k < nrOfPatterns; k++) {

                int state1 = states_child1[k];

                int w = l * matrixSize;

                if (state1 < nrOfStates) {


                    for (int i = 0; i < nrOfStates; i++) {

                        tmp = mat1[w + state1];

                        sum = 0.0;
                        for (int j = 0; j < nrOfStates; j++) {
                            sum += mat2[w] * partials_child2[v + j];
                            w++;
                        }

                        partials_parent[u] = tmp * sum;
                        u++;
                    }

                    v += nrOfStates;
                } else {
                    // Child 1 has a gap or unknown state so don't use it

                    for (int i = 0; i < nrOfStates; i++) {

                        sum = 0.0;
                        for (int j = 0; j < nrOfStates; j++) {
                            sum += mat2[w] * partials_child2[v + j];
                            w++;
                        }

                        partials_parent[u] = sum;
                        u++;
                    }

                    v += nrOfStates;
                }
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected void calculatePartialsPartialsPruning(RecombinationNetworkEdge edge1, RecombinationNetworkEdge edge2, RecombinationNetworkNode node,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
    		boolean[] computeForPatterns, double[] matrices1, double[] matrices2) {
    	
    	
        
        double[] mat1 = this.matrix.get(edge1.ID);
        double[] mat2 = this.matrix.get(edge2.ID);
        
        double[] partials_parent = this.partials.get(node.ID).get(computeFor);
        double[] partials_child1 = this.partials.get(edge1.childNode.ID).get(compute1);
        double[] partials_child2 = this.partials.get(edge2.childNode.ID).get(compute2);
               
        if (partials_parent==null) {
          	 partials.get(node.ID).put(computeFor, new double[nrOfMatrices*nrOfPatterns*nrOfStates]);
          	 partials_parent = partials.get(node.ID).get(computeFor);
        }
        
        double sum1, sum2;
        

        int u = 0;
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {

                int w = l * matrixSize;

                for (int i = 0; i < nrOfStates; i++) {

                    sum1 = sum2 = 0.0;

                    for (int j = 0; j < nrOfStates; j++) {
                        sum1 += mat1[w] 
                        		* partials_child1[v + j];
                        sum2 += mat2[w] 
                        		* partials_child2[v + j];

                        w++;
                    }
                    partials_parent[u] = sum1 * sum2;
                    u++;
                }
                v += nrOfStates;
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node when both children have exactly known states (e.g. for leaves).
     */
    protected void calculateStatesPruning(BreakPoints carries, RecombinationNetworkEdge edge, RecombinationNetworkNode node, boolean[] computeForPatterns) {
    	
        double[] mat = this.matrix.get(edge.ID);
        
        double[] partials_parent = this.partials.get(node.ID).get(carries);
        int[] states_child = states.get(edge.childNode.getTaxonLabel());
        
        if (partials_parent==null) {
         	 partials.get(node.ID).put(carries, new double[nrOfMatrices*nrOfPatterns*nrOfStates]);
         	 partials_parent = partials.get(node.ID).get(carries);
        }
        
        
        int v = 0;
        for (int l = 0; l < nrOfMatrices; l++) {
             for (int k = 0; k < nrOfPatterns; k++) {
                 int state1 = states_child[k];

                 int w = l * matrixSize;

                 if (state1 < nrOfStates) {
                     // child 2 has a gap or unknown state so treat it as unknown

                     for (int i = 0; i < nrOfStates; i++) {

                     	partials_parent[v] = mat[w + state1];

                         v++;
                         w += nrOfStates;
                     }
                 } else {
                     // both children have a gap or unknown state so set partials to 1

                     for (int j = 0; j < nrOfStates; j++) {
                     	partials_parent[v] = 1.0;
                         v++;
                     }
                 }
             }
         }
        


    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected void calculatePartialsPruning(BreakPoints carries, RecombinationNetworkEdge edge, RecombinationNetworkNode node, BreakPoints oldEdgePointer, boolean[] computeForPatterns) {

    	
        double[] mat = matrix.get(edge.ID);
        double[] partials_parent = partials.get(node.ID).get(carries);
        double[] partials_child = partials.get(edge.childNode.ID).get(oldEdgePointer);
        
        
        if (partials_parent==null) {
         	 partials.get(node.ID).put(carries, new double[nrOfMatrices*nrOfPatterns*nrOfStates]);
         	 partials_parent = partials.get(node.ID).get(carries);
        }
        
        double sum1;

        int u = 0;
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {

                int w = l * matrixSize;

                for (int i = 0; i < nrOfStates; i++) {

                    sum1 = 0.0;

                    for (int j = 0; j < nrOfStates; j++) {
                        sum1 += mat[w] * partials_child[v + j];
                        w++;
                    }

                    partials_parent[u] = sum1;
                    u++;
                }
                v += nrOfStates;
            }
        }
    }
    
    /**
     * Integrates partials across categories.
     *
     * @param inPartials  the array of partials to be integrated
     * @param proportions the proportions of sites in each category
     * @param outPartials an array into which the partials will go
     */
    @Override
	protected void calculateIntegratePartials(double[] proportions, double[] outPartials, Alignment data, HashMap<Integer, BreakPoints> rootBreaks) {
		int[] offset = new int[nrOfMatrices];
        for (int l = 1; l < nrOfMatrices; l++) {
        	offset[l] = l*nrOfStates*nrOfPatterns;
        }       
        
        
        for (Integer i : rootBreaks.keySet()) {
        	for (BreakPoints bp : partialsNew.getBreaks(i)) {
        		if (!bp.isEmpty()) {
	        		BreakPoints bp1 = bp.copy();
	        		bp1.and(rootBreaks.get(i));
	        		if (!bp1.isEmpty()) {
	                    double[] inPartials = partialsNew.getPartials(i, bp);
	                    for (int j = 0; j < bp.size();j++) {
		            		for (int m = bp1.getRange(j).from; m <= bp1.getRange(j).to; m++) {
		            			int stateStart = m*nrOfStates;
		            			int patternStart = data.getPatternIndex(m)*nrOfStates;
		    		            int u = stateStart;
		    		            int v = patternStart;
		    		            
		    	                for (int k = 0; k < nrOfStates; k++) {
		    	                	outPartials[u] = inPartials[v] * proportions[0];
		    	                    u++;
		    	                    v++;
		    	                }		
		    		
		    		            for (int l = 1; l < nrOfMatrices; l++) {
		    		                u = stateStart;	
		    		                v = patternStart+offset[l];
		    	                    for (int k = 0; k < nrOfStates; k++) {
		
		    	                    	outPartials[u] += inPartials[v] * proportions[l];
		    	                        u++;
		    	                        v++;
		    	                    }
		    		                
		    		            }
		    	        	}
	                    }
	        		}     
        		}        		
        	}        	
        }

//        for (BreakPoints bp : rootBreaks){
//            double[] inPartials = partialsNew.getPartials(node.ID, bp);
//            for (int j = 0; j< bp.size();j++) {
//
//        		for (int m = bp.getRange(j).from; m <= bp.getRange(j).to; m++) {
//        			int stateStart = m*nrOfStates;
//        			int patternStart = data.getPatternIndex(m)*nrOfStates;
//		            int u = stateStart;
//		            int v = patternStart;
//		            
//	                for (int i = 0; i < nrOfStates; i++) {
//	                	outPartials[u] = inPartials[v] * proportions[0];
//	                    u++;
//	                    v++;
//	                }		
//		
//		            for (int l = 1; l < nrOfMatrices; l++) {
//		                u = stateStart;	
//		                v = patternStart+offset[l];
//	                    for (int i = 0; i < nrOfStates; i++) {
//
//	                    	outPartials[u] += inPartials[v] * proportions[l];
//	                        u++;
//	                        v++;
//	                    }
//		                
//		            }
//	        	}
//        	}
//        }
        
        
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
        partialsNew = new Partials(nrOfNodes+50,nrOfMatrices*nrOfPatterns*nrOfStates);
        
        
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
	public void calculatePartials(RecombinationNetworkEdge edge1, RecombinationNetworkEdge edge2, RecombinationNetworkNode node, 
			BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
			boolean[] computeForPatterns, double[] matrices1, double[] matrices2) {
    	

        if (states.containsKey(edge1.childNode.getTaxonLabel())) {
            if (states.containsKey(edge2.childNode.getTaxonLabel())) {
                calculateStatesStatesPruning(edge1,edge2,node,computeFor,compute1,compute2,computeForPatterns, matrices1, matrices2);
            } else {
                calculateStatesPartialsPruning(edge1,edge2,node,computeFor,compute1,compute2,computeForPatterns, matrices1, matrices2);
            }
        } else {
            if (states.containsKey(edge2.childNode.getTaxonLabel())) {
                calculateStatesPartialsPruning(edge2,edge1,node,computeFor,compute2,compute1,computeForPatterns, matrices2, matrices1);
            } else {
                calculatePartialsPartialsPruning(edge1,edge2,node,computeFor,compute1,compute2,computeForPatterns, matrices1, matrices2);
            }
        }
        ensureLables(node, computeFor);
        
        if (!touched.contains(node.ID)) {
        	touched.add(node.ID);
        }
    }
    
    @Override
	public void calculatePartialsRecombination(RecombinationNetworkEdge edge, RecombinationNetworkNode node, BreakPoints compute1, BreakPoints computeFor, boolean[] computeForPatterns) {
    	// check if this computation may have already been solved
        if (states.containsKey(edge.childNode.getTaxonLabel())) {
        	calculateStatesPruning(compute1, edge, node,computeForPatterns);
        } else {
        	calculatePartialsPruning(compute1, edge, node, computeFor,computeForPatterns);
        }
        
        ensureLables(node, compute1);
        
        
        if (!touched.contains(node.ID)) {
        	touched.add(node.ID);
        }
        
    }


	@Override
	public void integratePartials(double[] proportions, double[] outPartials, Alignment data, HashMap<Integer, BreakPoints> rootBreaks) {
        calculateIntegratePartials(proportions, outPartials, data, rootBreaks);
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
        
        partialsNew.restore();

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
    	partialsNew.store();
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
    public void initPartials(RecombinationNetworkNode node, int length) {
		partialsNew.addNewNode(node.ID);
	}
	
	@Override
	public void setStates(RecombinationNetworkNode node, int[] states) {
        int[] newstates = new int[states.length];
		System.arraycopy(states, 0, newstates, 0, states.length);
		this.states.put(node.getTaxonLabel(), newstates);
	}


	@Override
	public void cleanMatrix(List<Integer> edgesIDs) {
		touched = new ArrayList<>();
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
		for (Integer n : partialsNew.keySet()) {
			if (!nodeIDs.contains(n))
				remove.add(n);
		}
		for (Integer n : remove)
			partialsNew.remove(n);
	}


	@Override
	public void rapidStore() {
		// TODO Auto-generated method stub
		
	}


	@Override
	protected void cleanPartialsNode(RecombinationNetworkNode node) {
		partialsNew.removeBreaks(node.ID);
	}

	protected void ensureLables(RecombinationNetworkNode node, BreakPoints computeFor) {
		List<BreakPoints> breaks = partialsNew.getBreaks(node.ID);
		


		for (BreakPoints bp : breaks) {
			if (!bp.isEmpty()) {
				if (!computeFor.equals(bp)) {
					BreakPoints cp = computeFor.copy();
					cp.and(bp);
					if (!cp.isEmpty()) {
				        if (!touched.contains(node.ID)) {
				        	touched.add(node.ID);
				        }
						BreakPoints cp2 = bp.copy();
						cp2.andNot(computeFor);
						partialsNew.replaceBreaks(node.ID, bp, cp2);
					}				
				}
			}
		}	
		
	}
	
	@Override
	protected void checkLabels(RecombinationNetworkNode node, BreakPoints computeFor) {

		List<BreakPoints> breaks = partialsNew.getBreaks(node.ID);
		
        if (breaks.contains(computeFor))
			return;

		for (BreakPoints bp : breaks) {
			if (!bp.isEmpty()) {
				if (!computeFor.equals(bp)) {
					BreakPoints cp = computeFor.copy();
					cp.andNot(bp);
					if (cp.isEmpty()) {
				        if (!touched.contains(node.ID)) {
				        	touched.add(node.ID);
				        }
						partialsNew.replaceBreaks(node.ID, bp, computeFor.copy());
	
					}			
				}
			}
		}
		ensureLables(node, computeFor);

		return;

	}
	
	@Override
	protected int reassignLabels(RecombinationNetworkNode node, BreakPoints computeFor, BreakPoints dirtyEdges) {
		List<BreakPoints> breaks = partialsNew.getBreaks(node.ID);
		
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
				        if (!touched.contains(node.ID)) {
				        	touched.add(node.ID);
				        }
				        partialsNew.replaceBreaks(node.ID, bp, computeFor.copy());
					}			
				}
			}
		}
		ensureLables(node, computeFor);

		return 0;
	}

	

} // class BeerLikelihoodCore
