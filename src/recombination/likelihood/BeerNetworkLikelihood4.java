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

/**
 * standard likelihood core, uses no caching *
 */
public class BeerNetworkLikelihood4 extends BeerNetworkLikelihoodCore {
    public BeerNetworkLikelihood4() {
    	super(4);
    }


    /**
     * Calculates partial likelihoods at a node when both children have exactly known states (e.g. for leaves).
     */
    protected void calculateStatesStatesPruning(Integer ID1, Integer ID2, Integer ID3,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
    		boolean[] computeForPatterns, double[] matrices1, double[] matrices2) {
        
        double[] partials3 = partials.getPartialsOperation(ID3, computeFor);
        int[] stateIndex1 = states.get(ID1);
        int[] stateIndex2 = states.get(ID2);

        
        int v=0;
        	
        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {
            	if (computeForPatterns[k]) {
	
	                int state1 = stateIndex1[k];
	                int state2 = stateIndex2[k];
	
	                int w = l * matrixSize;
	
	                if (state1 < 4 && state2 < 4) {
                    
	                    partials3[v] = matrices1[w + state1] * matrices2[w + state2];
	                    v++;
	                    w += 4;
	                    partials3[v] = matrices1[w + state1] * matrices2[w + state2];
	                    v++;
	                    w += 4;
	                    partials3[v] = matrices1[w + state1] * matrices2[w + state2];
	                    v++;
	                    w += 4;
	                    partials3[v] = matrices1[w + state1] * matrices2[w + state2];
	                    v++;
	                    w += 4;
	
	                } else if (state1 < 4) {
	                    // child 2 has a gap or unknown state so don't use it
	
	                    partials3[v] = matrices1[w + state1];
	                    v++;
	                    w += 4;
	                    partials3[v] = matrices1[w + state1];
	                    v++;
	                    w += 4;
	                    partials3[v] = matrices1[w + state1];
	                    v++;
	                    w += 4;
	                    partials3[v] = matrices1[w + state1];
	                    v++;
	                    w += 4;
	
	                } else if (state2 < 4) {
	                    // child 2 has a gap or unknown state so don't use it
	                    partials3[v] = matrices2[w + state2];
	                    v++;
	                    w += 4;
	                    partials3[v] = matrices2[w + state2];
	                    v++;
	                    w += 4;
	                    partials3[v] = matrices2[w + state2];
	                    v++;
	                    w += 4;
	                    partials3[v] = matrices2[w + state2];
	                    v++;
	                    w += 4;
	
	                } else {
	                    // both children have a gap or unknown state so set partials to 1
	                    partials3[v] = 1.0;
	                    v++;
	                    partials3[v] = 1.0;
	                    v++;
	                    partials3[v] = 1.0;
	                    v++;
	                    partials3[v] = 1.0;
	                    v++;
	                }
	            }else {
	            	 v+=4;
	            }
            }
        }      
        
//        System.out.println(node.getHeight() + "  fdslkhkjdlfskljfds");
    }
    
    /**
     * Calculates partial likelihoods at a node when one child has states and one has partials.
     */
    protected void calculateStatesPartialsPruning(Integer ID1, Integer ID2, Integer ID3,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
    		boolean[] computeForPatterns, double[] matrices1, double[] matrices2) { 
        
    	double sum;

        int u = 0;
        int v = 0;      
        
        double[] partials3 = partials.getPartialsOperation(ID3, computeFor);
        int[] states_child1 = this.states.get(ID1);
        double[] partials2 = partials.getPartials(ID2, compute2);
        
       
        for (int l = 0; l < nrOfMatrices; l++) {
            for (int k = 0; k < nrOfPatterns; k++) {
            	if (computeForPatterns[k]) {
	
	                int state1 = states_child1[k];
	
	                int w = l * matrixSize;
	
	                if (state1 < 4) {	
	                    sum = matrices2[w] * partials2[v];
	                    sum += matrices2[w + 1] * partials2[v + 1];
	                    sum += matrices2[w + 2] * partials2[v + 2];
	                    sum += matrices2[w + 3] * partials2[v + 3];
	                    partials3[u] = matrices1[w + state1] * sum;
	                    u++;
	
	                    sum = matrices2[w + 4] * partials2[v];
	                    sum += matrices2[w + 5] * partials2[v + 1];
	                    sum += matrices2[w + 6] * partials2[v + 2];
	                    sum += matrices2[w + 7] * partials2[v + 3];
	                    partials3[u] = matrices1[w + 4 + state1] * sum;
	                    u++;
	
	                    sum = matrices2[w + 8] * partials2[v];
	                    sum += matrices2[w + 9] * partials2[v + 1];
	                    sum += matrices2[w + 10] * partials2[v + 2];
	                    sum += matrices2[w + 11] * partials2[v + 3];
	                    partials3[u] = matrices1[w + 8 + state1] * sum;
	                    u++;
	
	                    sum = matrices2[w + 12] * partials2[v];
	                    sum += matrices2[w + 13] * partials2[v + 1];
	                    sum += matrices2[w + 14] * partials2[v + 2];
	                    sum += matrices2[w + 15] * partials2[v + 3];
	                    partials3[u] = matrices1[w + 12 + state1] * sum;
	                    u++;
	
	                    v += 4;
	
	                } else {
	                    // Child 1 has a gap or unknown state so don't use it
	
	
	                    sum = matrices2[w] * partials2[v];
	                    sum += matrices2[w + 1] * partials2[v + 1];
	                    sum += matrices2[w + 2] * partials2[v + 2];
	                    sum += matrices2[w + 3] * partials2[v + 3];
	                    partials3[u] = sum;
	                    u++;
	
	                    sum = matrices2[w + 4] * partials2[v];
	                    sum += matrices2[w + 5] * partials2[v + 1];
	                    sum += matrices2[w + 6] * partials2[v + 2];
	                    sum += matrices2[w + 7] * partials2[v + 3];
	                    partials3[u] = sum;
	                    u++;
	
	                    sum = matrices2[w + 8] * partials2[v];
	                    sum += matrices2[w + 9] * partials2[v + 1];
	                    sum += matrices2[w + 10] * partials2[v + 2];
	                    sum += matrices2[w + 11] * partials2[v + 3];
	                    partials3[u] = sum;
	                    u++;
	
	                    sum = matrices2[w + 12] * partials2[v];
	                    sum += matrices2[w + 13] * partials2[v + 1];
	                    sum += matrices2[w + 14] * partials2[v + 2];
	                    sum += matrices2[w + 15] * partials2[v + 3];
	                    partials3[u] = sum;
	                    u++;
	
	                    v += 4;
	                }
            	}else {
            		partials3[v] = Double.NaN;
            		v+=4;
            		u+=4;
                }
            }
        }
        

    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected void calculatePartialsPartialsPruning(Integer ID1, Integer ID2, Integer ID3,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, 
    		boolean[] computeForPatterns, double[] matrices1, double[] matrices2) {
    	
        double[] partials3 = partials.getPartialsOperation(ID3, computeFor);
        
        
        double[] partials1 = partials.getPartials(ID1, compute1);
        double[] partials2 = partials.getPartials(ID2, compute2);
               

        int u = 0;
        int v = 0;
        double sum1,sum2;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {
            	if (computeForPatterns[k]) {
	
	                int w = l * matrixSize;
	
	                sum1 = matrices1[w] * partials1[v];
	                sum2 = matrices2[w] * partials2[v];
	                sum1 += matrices1[w + 1] * partials1[v + 1];
	                sum2 += matrices2[w + 1] * partials2[v + 1];
	                sum1 += matrices1[w + 2] * partials1[v + 2];
	                sum2 += matrices2[w + 2] * partials2[v + 2];
	                sum1 += matrices1[w + 3] * partials1[v + 3];
	                sum2 += matrices2[w + 3] * partials2[v + 3];
	                partials3[u] = sum1 * sum2;
	                u++;
	
	                sum1 = matrices1[w + 4] * partials1[v];
	                sum2 = matrices2[w + 4] * partials2[v];
	                sum1 += matrices1[w + 5] * partials1[v + 1];
	                sum2 += matrices2[w + 5] * partials2[v + 1];
	                sum1 += matrices1[w + 6] * partials1[v + 2];
	                sum2 += matrices2[w + 6] * partials2[v + 2];
	                sum1 += matrices1[w + 7] * partials1[v + 3];
	                sum2 += matrices2[w + 7] * partials2[v + 3];
	                partials3[u] = sum1 * sum2;
	                u++;
	
	                sum1 = matrices1[w + 8] * partials1[v];
	                sum2 = matrices2[w + 8] * partials2[v];
	                sum1 += matrices1[w + 9] * partials1[v + 1];
	                sum2 += matrices2[w + 9] * partials2[v + 1];
	                sum1 += matrices1[w + 10] * partials1[v + 2];
	                sum2 += matrices2[w + 10] * partials2[v + 2];
	                sum1 += matrices1[w + 11] * partials1[v + 3];
	                sum2 += matrices2[w + 11] * partials2[v + 3];
	                partials3[u] = sum1 * sum2;
	                u++;
	
	                sum1 = matrices1[w + 12] * partials1[v];
	                sum2 = matrices2[w + 12] * partials2[v];
	                sum1 += matrices1[w + 13] * partials1[v + 1];
	                sum2 += matrices2[w + 13] * partials2[v + 1];
	                sum1 += matrices1[w + 14] * partials1[v + 2];
	                sum2 += matrices2[w + 14] * partials2[v + 2];
	                sum1 += matrices1[w + 15] * partials1[v + 3];
	                sum2 += matrices2[w + 15] * partials2[v + 3];
	                partials3[u] = sum1 * sum2;
	                u++;
	
	                v += 4;
	            }else {
	            	v+=4;
               	 	u+=4;
               }
            }
        }


    }
    
	@Override
	public double integratePartials(double[] proportions, double[] frequencies, Alignment data, HashMap<Integer, BreakPoints> rootBreaks) {
		
        double logP = 0;
        int[] nrInPattern;
        double[] outPartials = new double[4*nrOfPatterns];
        
        
        
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
		                            outPartials[u+1] = inPartials[v+1] * proportions[0];
		                            outPartials[u+2] = inPartials[v+2] * proportions[0];
		                            outPartials[u+3] = inPartials[v+3] * proportions[0];
		                        }
	                    	}
                    		u+=4;
                    		v+=4;
	                    }


	                    for (int l = 1; l < nrOfMatrices; l++) {
	                        u = 0;
	                        for (int k = 0; k < nrOfPatterns; k++) {
	                        	if (nrInPattern[k]!=0) {
	                                outPartials[u] += inPartials[v] * proportions[l];
	                                outPartials[u+1] += inPartials[v+1] * proportions[l];
	                                outPartials[u+2] += inPartials[v+2] * proportions[l];
	                                outPartials[u+3] += inPartials[v+3] * proportions[l];
	                        	}
	                    		u+=4;
	                    		v+=4;

	                        }
	                    }
	                    u = 0;
                        for (int k = 0; k < nrOfPatterns; k++) {
                        	if (nrInPattern[k]!=0) {
                        		double sum = frequencies[0] * outPartials[u];
                        		sum+= frequencies[1] * outPartials[u+1];
                        		sum+= frequencies[2] * outPartials[u+2];
                        		sum+= frequencies[3] * outPartials[u+3];
                                logP += Math.log(sum)*nrInPattern[k];
                        	
                        	}
                        	u+=4;
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


    
} // class BeerLikelihoodCore
