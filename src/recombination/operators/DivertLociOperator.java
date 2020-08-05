package recombination.operators;

import beast.core.Input;
import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class DivertLociOperator extends EmptyEdgesRecombinationNetworkOperator {

    public Input<Double> scaleFactorInput = new Input<>(
            "scaleFactor",
            "Scale factor tuning parameter.",
            1000.0);

    double lambdaDiversion;
    public int totalLength;
    
    public void initAndValidate() {
    	lambdaDiversion = scaleFactorInput.get();    	
    	super.initAndValidate();
    	totalLength = network.totalLength;
    }
	
    @Override
    public double networkProposal() {
    	    	
        double logHR = 0.0;
        
        List<RecombinationNetworkEdge> sourceEdges = network.getEdges().stream()
                .filter(e -> e.childNode.isRecombination())
                .filter(e -> !e.breakPoints.isEmpty())
                .collect(Collectors.toList());

        if (sourceEdges.isEmpty())
            return Double.NEGATIVE_INFINITY;

        logHR -= Math.log(1.0/sourceEdges.size());

        RecombinationNetworkEdge sourceEdge = sourceEdges.get(Randomizer.nextInt(sourceEdges.size()));
        RecombinationNetworkEdge destEdge = getSpouseEdge(sourceEdge);

        network.startEditing(this);
        
        System.out.println(network);
        System.out.println(sourceEdge.childNode.getHeight());
        BreakPoints newRange = getUniformRange(sourceEdge);
        
        
        if (newRange.getLength()>sourceEdge.passingRange.getLength()) {
            BreakPoints lociToDivert = destEdge.breakPoints.copy();
            BreakPoints rangeToDivert = destEdge.passingRange.copy();
            
	        if (lociToDivert.isEmpty()) {
		        destEdge.passingRange.andNot(rangeToDivert);
		    	sourceEdge.passingRange.or(rangeToDivert);
		    	return logHR;
	        }
            
        	rangeToDivert.and(newRange);
	    	lociToDivert.and(rangeToDivert);

	        if (lociToDivert.isEmpty()) {
		        destEdge.passingRange.andNot(rangeToDivert);
		    	sourceEdge.passingRange.or(rangeToDivert);
		    	return logHR;
	        }else {			
		        logHR -= addLociToAncestors(sourceEdge, lociToDivert);
		        logHR += removeLociFromAncestors(destEdge, lociToDivert);
		        
		        destEdge.passingRange.andNot(rangeToDivert);
		    	sourceEdge.passingRange.or(rangeToDivert);
	        }

        }else {
            BreakPoints lociToDivert = sourceEdge.breakPoints.copy();
        	BreakPoints rangeToDivert = sourceEdge.passingRange.copy();
        	rangeToDivert.andNot(newRange);
	    	lociToDivert.and(rangeToDivert);
	    	
	    	
	        if (lociToDivert.isEmpty()) {
		    	sourceEdge.passingRange.andNot(rangeToDivert);
		    	if (destEdge.passingRange==null)
		    		destEdge.passingRange = rangeToDivert;
		    	else
		    		destEdge.passingRange.or(rangeToDivert);
	        }else {		
		        logHR -= addLociToAncestors(destEdge, lociToDivert);
		        logHR += removeLociFromAncestors(sourceEdge, lociToDivert);
		        
		    	sourceEdge.passingRange.andNot(rangeToDivert);
		    	if (destEdge.passingRange==null)
		    		destEdge.passingRange = rangeToDivert;
		    	else
		    		destEdge.passingRange.or(rangeToDivert);
	        }

        }
        


        int reverseSourceEdgeCount = (int)(network.getEdges().stream()
                .filter(e -> e.childNode.isRecombination())
                .filter(e -> !e.breakPoints.isEmpty())
                .count());

        logHR += Math.log(1.0/reverseSourceEdgeCount);
        
        return logHR;
    }


    /**
     * Remove segments from this edge and ancestors.
     *
     * @param edge edge at which to start removal
     * @param segsToRemove segments to remove from edge and ancestors
     * @return log probability of reverse operation
     */
    public double removeLociFromAncestors(RecombinationNetworkEdge edge, BreakPoints rangeToRemove) {
        double logP = 0.0;

        rangeToRemove = rangeToRemove.copy();
        
        rangeToRemove.and(edge.breakPoints);
        
        if (rangeToRemove.isEmpty())
            return logP;
        
        edge.breakPoints.andNot(rangeToRemove.copy());       
                       
        if (edge.isRootEdge())
            return logP;

        if (edge.parentNode.isRecombination()) {
        	
        	// get the breakpoints before the operation
        	int old_pr1 = edge.parentNode.getParentEdges().get(0).passingRange.getMax();
        	int old_pr2 = edge.parentNode.getParentEdges().get(1).passingRange.getMax();
        	
            logP += removeLociFromAncestors(edge.parentNode.getParentEdges().get(0), rangeToRemove);
            logP += removeLociFromAncestors(edge.parentNode.getParentEdges().get(1), rangeToRemove);
                      
    		int min1 = edge.parentNode.getParentEdges().get(0).breakPoints.getMin();
    		int max1 = edge.parentNode.getParentEdges().get(0).breakPoints.getMax();
    		
    		int min2 = edge.parentNode.getParentEdges().get(1).breakPoints.getMin();
    		int max2 = edge.parentNode.getParentEdges().get(1).breakPoints.getMax();
    		
            if (min1==-1 && min2==-1) {
    			// both passing ranges are null
            	logP += Math.log(0.5) + Math.log(1.0/(totalLength-1));
            }else if (min1==-1){
            	if (max2==totalLength-1) {
                	logP += Math.log(1.0/(min2));
            	}else if (min2==0) {
            		logP += Math.log(1.0/(totalLength-max2-1));
            	}else {
                	logP += Math.log(0.5);
                	if (old_pr2>=max2) {
                		logP += Math.log(1.0/(totalLength-max2-1));
                	}else if (old_pr2<min2) {
                		logP += Math.log(1.0/(min2));
                	}else {
                		throw new IllegalArgumentException("should not happen, or at least is not accounted for");
                	}
            	}
            }else if(min2==-1) {
            	if (max1==totalLength-1) {
                	logP += Math.log(1.0/(min1+1));
            	}else if (min1==0) {
            		logP += Math.log(1.0/(totalLength-max1-1));
            	}else {
                	logP += Math.log(0.5);
                	if (old_pr1>=max1) {
                		logP += Math.log(1.0/(totalLength-max1-1));
                	}else if (old_pr1<min1) {
                		logP += Math.log(1.0/(min1));
                	}else {
                		throw new IllegalArgumentException("should not happen, or at least is not accounted for");
                	}
            	}
            }else {
            	double diff = Math.max(min2-max1, min1-max2); 
        		logP += Math.log(1.0/(diff));        		
            }
            
            
        } else {
        	rangeToRemove.andNot(getSisterEdge(edge).breakPoints);
            logP += removeLociFromAncestors(edge.parentNode.getParentEdges().get(0), rangeToRemove);
        }
        return logP;
    }

    /**
     * Add segments to this edge and ancestors.
     *
     * @param edge edge at which to start addition
     * @param segsToAdd segments to add to the edge and ancestors
     * @return log probability of operation
     */
    public double addLociToAncestors(RecombinationNetworkEdge edge, BreakPoints rangeToAdd) {
        double logP = 0.0;

        rangeToAdd = rangeToAdd.copy();
               
        if (rangeToAdd.isEmpty())
            return logP;

        rangeToAdd.andNot(edge.breakPoints);

        if (rangeToAdd.isEmpty())
            return logP;        

        edge.breakPoints.or(rangeToAdd);

        if (edge.isRootEdge())
            return logP;        
        
        if (edge.parentNode.isRecombination()) {        	
        	// resample the passing Range between the boundries given by the left and right breakpoints
        	logP += resamplePassingRange(edge.parentNode.getParentEdges().get(0),
        			edge.parentNode.getParentEdges().get(1));  
        	
            BreakPoints rangeToAddLeft = rangeToAdd.copy();
            BreakPoints rangeToAddRight = rangeToAdd.copy();
                        
            rangeToAddLeft.and(edge.parentNode.getParentEdges().get(0).passingRange);
            rangeToAddRight.and(edge.parentNode.getParentEdges().get(1).passingRange);
                        
            logP += addLociToAncestors(edge.parentNode.getParentEdges().get(0), rangeToAddLeft);
            logP += addLociToAncestors(edge.parentNode.getParentEdges().get(1), rangeToAddRight);
        } else {
            logP += addLociToAncestors(edge.parentNode.getParentEdges().get(0), rangeToAdd);
        }

        return logP;
    }
        

    
    public BreakPoints getNewRangeToDivert(RecombinationNetworkEdge sourceEdge) {
    	
    	BreakPoints rangeToDivert = new BreakPoints();
    	
    	
    	int newBreakPoint = Randomizer.nextInt((int) sourceEdge.breakPoints.getLength()) + sourceEdge.breakPoints.getMin();

    	
    	if (Randomizer.nextBoolean()) {
    		rangeToDivert = new BreakPoints(0,newBreakPoint);
    	}else {
    		rangeToDivert = new BreakPoints(newBreakPoint, totalLength-1);
    	}    	    	    	    	    	
    	return rangeToDivert;
    }
    
    
    public BreakPoints getUniformRange(RecombinationNetworkEdge sourceEdge) {
    	
    	BreakPoints newRange = new BreakPoints();
    	
 
    	int newBreakPoint = Randomizer.nextInt(totalLength);

    	
    	if (sourceEdge.passingRange.getMin()==0) {
    		newRange = new BreakPoints(0,newBreakPoint);
    	}else {
    		newRange = new BreakPoints(newBreakPoint, totalLength-1);
    	}    	    	    	    	    	
    	return newRange;
    }
    
   
    private double resamplePassingRange(RecombinationNetworkEdge edge1,
			RecombinationNetworkEdge edge2) {
    	
    	double logHR = 0.0;

		int min1 = edge1.breakPoints.getMin();
		int max1 = edge1.breakPoints.getMax();
		
		int min2 = edge2.breakPoints.getMin();
		int max2 = edge2.breakPoints.getMax();
    	
		if (min1==-1 && min2==-1) {
			// both passing ranges are null
			logHR += Math.log(0.5) + Math.log(1.0/(totalLength-1));
			if (Randomizer.nextBoolean()) {
				int start = Randomizer.nextInt(totalLength-1);
				edge1.setPassingRange(start+1, totalLength-1);
				edge2.setPassingRange(0,start);
			}else {
				int start = Randomizer.nextInt(totalLength-1);
				edge1.setPassingRange(0, start);
				edge2.setPassingRange(start+1, totalLength-1);
			}	
		}else if(min1==-1) {
			// 1 is null
			logHR += sampleNewPassingRange(edge2,edge1,min2,max2);
		}else if(min2==-1) {
			// 2 is null
			logHR += sampleNewPassingRange(edge1,edge2,min1,max1);
		}else {
			// resample between the two edges
			if (max1>max2) {
				int newBreakPoint = Randomizer.nextInt(min1-max2)+max2;
				edge1.setPassingRange(newBreakPoint+1, totalLength-1);
				edge2.setPassingRange(0, newBreakPoint);
				logHR += Math.log(1.0/(min1-max2));
			}else {
				int newBreakPoint = Randomizer.nextInt(min2-max1)+max1;
				edge2.setPassingRange(newBreakPoint+1, totalLength-1);
				edge1.setPassingRange(0, newBreakPoint);
				logHR += Math.log(1.0/(min2-max1));
			}
		}   	
		return logHR;
	}
    
    
    private double sampleNewPassingRange(RecombinationNetworkEdge edge1,
			RecombinationNetworkEdge edge2, int min, int max) {
    	double logHR = 0.0;
    	if (min==0 && max==totalLength-1) {
			// passing range 2 is null
			edge1.setPassingRange(0, totalLength-1);
			edge2.passingRange=null;
    	}else if (min>0 && max<totalLength-1) {
    		logHR += Math.log(0.5);
			if (Randomizer.nextBoolean()) {		
				int diff = totalLength-max-1;
				logHR += Math.log(1.0/diff);				
				int start = Randomizer.nextInt(diff)+max+1;
				edge2.setPassingRange(start+1, totalLength-1);
				edge1.setPassingRange(0,start);
			}else {
				logHR += Math.log(1.0/min);
				int end = Randomizer.nextInt(min);
				edge2.setPassingRange(0, end);
				edge1.setPassingRange(end+1,totalLength-1);
			}			
		}else if (min>0){
			logHR += Math.log(1.0/min);

			int end = Randomizer.nextInt(min);
			edge2.setPassingRange(0, end);
			edge1.setPassingRange(end+1,totalLength-1);

		}else {
			int diff = totalLength-max-1;
			logHR += Math.log(1.0/diff);

			int start = Randomizer.nextInt(diff)+max+1;
			edge2.setPassingRange(start+1, totalLength-1);
			edge1.setPassingRange(0,start);
		}   

    	return logHR;
    }
}
