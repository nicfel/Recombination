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
    protected void calculateStatesStatesPruning(RecombinationNetworkEdge edge1, RecombinationNetworkEdge edge2, RecombinationNetworkNode node,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, boolean[] computeForPatterns) {
        
    	        
        double[] matrices1 = matrix.get(edge1.ID);
        double[] matrices2 = matrix.get(edge2.ID);
        
        double[] partials3 = partialsNew.getPartialsOperation(node.ID, computeFor);
        int[] stateIndex1 = states.get(edge1.childNode.getTaxonLabel());
        int[] stateIndex2 = states.get(edge2.childNode.getTaxonLabel());
        

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
    protected void calculateStatesPartialsPruning(RecombinationNetworkEdge edge1, RecombinationNetworkEdge edge2, RecombinationNetworkNode node,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, boolean[] computeForPatterns) { 
        
    	double sum;

        int u = 0;
        int v = 0;
        
        double[] matrices1 = this.matrix.get(edge1.ID);
        double[] matrices2 = this.matrix.get(edge2.ID);
        
        double[] partials3 = partialsNew.getPartialsOperation(node.ID, computeFor);
        int[] states_child1 = this.states.get(edge1.childNode.getTaxonLabel());
        double[] partials2 = partialsNew.getPartials(edge2.childNode.ID, compute2);
       
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
    protected void calculatePartialsPartialsPruning(RecombinationNetworkEdge edge1, RecombinationNetworkEdge edge2, RecombinationNetworkNode node,
    		BreakPoints computeFor, BreakPoints compute1, BreakPoints compute2, boolean[] computeForPatterns) {
    	
        double[] matrices1 = this.matrix.get(edge1.ID);
        double[] matrices2 = this.matrix.get(edge2.ID);
                
        double[] partials3 = partialsNew.getPartialsOperation(node.ID, computeFor);
        double[] partials1 = partialsNew.getPartials(edge1.childNode.ID, compute1);
        double[] partials2 = partialsNew.getPartials(edge2.childNode.ID, compute2);
               

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
            		partials3[v] = Double.NaN;
	              	 v+=4;
	               	 u+=4;
               }
            }
        }


    }

    /**
     * Calculates partial likelihoods at a node when both children have exactly known states (e.g. for leaves).
     */
    protected void calculateStatesPruning(BreakPoints carries, RecombinationNetworkEdge edge, RecombinationNetworkNode node, boolean[] computeForPatterns) {
    	
        double[] matrices1 = this.matrix.get(edge.ID);
        
        double[] partials3 = partialsNew.getPartialsOperation(node.ID, carries);
        int[] states_child = states.get(edge.childNode.getTaxonLabel());
        
        
        int v = 0;
        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {
            	if (computeForPatterns[k]){
	
	                int state1 = states_child[k];
	
	                int w = l * matrixSize;
	
	                if (state1 < 4) {
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
            		partials3[v] = Double.NaN;
            		v+=4;
                  	
            	}
            }
        }
    }

    /**
     * Calculates partial likelihoods at a node when both children have partials.
     */
    protected void calculatePartialsPruning(BreakPoints carries, RecombinationNetworkEdge edge, RecombinationNetworkNode node, BreakPoints oldEdgePointer, boolean[] computeForPatterns) {
    	
        double[] partials3 = partialsNew.getPartialsOperation(node.ID, carries);
        double[] matrices1 = matrix.get(edge.ID);
        double[] partials1 = partialsNew.getPartials(edge.childNode.ID, oldEdgePointer);

        double sum1;

        int u = 0;
        int v = 0;

        for (int l = 0; l < nrOfMatrices; l++) {

            for (int k = 0; k < nrOfPatterns; k++) {
	            	if (computeForPatterns[k]) {
	
	                int w = l * matrixSize;
	
	                sum1 = matrices1[w] * partials1[v];
	                sum1 += matrices1[w + 1] * partials1[v + 1];
	                sum1 += matrices1[w + 2] * partials1[v + 2];
	                sum1 += matrices1[w + 3] * partials1[v + 3];
	                partials3[u] = sum1;
	                u++;
	
	                sum1 = matrices1[w + 4] * partials1[v];
	                sum1 += matrices1[w + 5] * partials1[v + 1];
	                sum1 += matrices1[w + 6] * partials1[v + 2];
	                sum1 += matrices1[w + 7] * partials1[v + 3];
	                partials3[u] = sum1;
	                u++;
	
	                sum1 = matrices1[w + 8] * partials1[v];
	                sum1 += matrices1[w + 9] * partials1[v + 1];
	                sum1 += matrices1[w + 10] * partials1[v + 2];
	                sum1 += matrices1[w + 11] * partials1[v + 3];
	                partials3[u] = sum1;
	                u++;
	
	                sum1 = matrices1[w + 12] * partials1[v];
	                sum1 += matrices1[w + 13] * partials1[v + 1];
	                sum1 += matrices1[w + 14] * partials1[v + 2];
	                sum1 += matrices1[w + 15] * partials1[v + 3];
	                partials3[u] = sum1;
	                u++;
	
	                v += 4;
	            }else {
            		partials3[v] = Double.NaN;
	              	 v+=4;
	               	 u+=4;
               }
            }
        }
//        if (error)
//        	System.out.println(node.getHeight() + " " + carries);


    }

} // class BeerLikelihoodCore
