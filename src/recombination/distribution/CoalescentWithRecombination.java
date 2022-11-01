package recombination.distribution;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import recombination.statistics.RecombinationNetworkStatsLogger;

import java.util.List;


/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a reassortment network using under" +
        " the framework of Mueller (2018).")
public class CoalescentWithRecombination extends RecombinationNetworkDistribution {
	
	public Input<Function> recombinationRateInput = new Input<>(
	        "recombinationRate",
            "recombination rate (per lineage per unit time)",
            Input.Validate.REQUIRED);
	
	public Input<RealParameter> relativeRecombinationRateInput = new Input<>(
	        "relativeRecombinationRate",
            "relative recombination rate (per lineage per unit time) for a specific part of the genome",
            Input.Validate.OPTIONAL);

	public Input<PopulationFunction> populationFunctionInput = new Input<>(
	        "populationModel",
            "Population model.",
            Input.Validate.REQUIRED);
	
	public Input<Boolean> conditionOnCoalescentEventsInput = new Input<>(
	        "conditionOnCoalescentEvents",
            "if true, only coalescent events are allowed after the .",
            true);	
	
	public Input<Double> maxHeightRatioInput = new Input<>(
	        "maxHeightRatio",
            "if specified, above the ratio, only coalescent events are allowed.", Double.POSITIVE_INFINITY);	

	public Input<Double> redFactorInput = new Input<>(
	        "redFactor",
            "by how much the recombination rate should be reduced after reaching the maxHeightRatio.", 0.2);	




    public PopulationFunction populationFunction;
    public Function recombinationRate;
    public RealParameter relativeRecombinationRate;
    public RecombinationNetworkIntervals intervals;
    double[] recRates;
    public double redFactor;
        
    @Override
    public void initAndValidate(){
        populationFunction = populationFunctionInput.get();
        recombinationRate = recombinationRateInput.get();
        intervals = networkIntervalsInput.get();
        redFactor = redFactorInput.get();
                
        
        if (intervals.hasBreakPoints) {
        	if (relativeRecombinationRateInput.get()==null)
        		throw new IllegalArgumentException("relative recombination rates have to be defined");        	
       	
        	relativeRecombinationRate = relativeRecombinationRateInput.get();
        	relativeRecombinationRate.setDimension(intervals.recBP.length);
        	recRates = new double[intervals.recBP.length];
        }else {
        	recRates = new double[1];
        }
    }
    
    public double calculateLogP() {
    	logP = 0;  	    	
    	
    	// Calculate tree intervals
    	List<RecombinationNetworkEvent> networkEventList = intervals.getRecombinationNetworkEventList();
    	
    	if (recRates.length>1)
    		for (int i = 0; i < recRates.length; i++)
    			recRates[i] = Math.exp(relativeRecombinationRate.getArrayValue(i)) * recombinationRate.getArrayValue();
		else
			recRates[0] = recombinationRate.getArrayValue();
    	

    	// get the mrca of all loci trees
    	double lociMRCA = conditionOnCoalescentEventsInput.get() ? RecombinationNetworkStatsLogger.getMaxLociMRCA(intervals.recombinationNetworkInput.get()) : Double.POSITIVE_INFINITY;

    	RecombinationNetworkEvent prevEvent = null;
    	
    	for (RecombinationNetworkEvent event : networkEventList) {
        	if (prevEvent != null)
        		logP += intervalContribution(prevEvent, event, lociMRCA);

        	switch (event.type) {
				case COALESCENCE:
					logP += coalesce(event);
					break;

				case SAMPLE:
					break;

				case RECOMBINATION:
					logP += recombination(event, lociMRCA);
					break;
			}

       		if (logP==Double.NEGATIVE_INFINITY)
       			break;
       		

        	prevEvent = event;
        } 

    	return logP;
    }
    
	private double recombination(RecombinationNetworkEvent event, double lociMRCA) {
		if (conditionOnCoalescentEventsInput.get()) {
	        if (event.time<=(lociMRCA*maxHeightRatioInput.get())) {
	    		int i = 0; 
	    		while (!intervals.recBP[i].contains(event.breakPoint)) {
	    			i++;
	    		}        	
	    		return Math.log(recRates[i]);
	        }else {
	    		int i = 0; 
	    		while (!intervals.recBP[i].contains(event.breakPoint)) {
	    			i++;
	    		}        	
	    		return Math.log(redFactor*recRates[i]);
	        }        		
		}else{
    		int i = 0; 
    		while (!intervals.recBP[i].contains(event.breakPoint)) {
    			i++;
    		}        	
    		return Math.log(recRates[i]);
		}
	}

	private double coalesce(RecombinationNetworkEvent event) {
		return Math.log(1.0/populationFunction.getPopSize(event.time));
	}

	private double intervalContribution(RecombinationNetworkEvent prevEvent, RecombinationNetworkEvent event, double lociMRCA) {
        double result = 0.0;
		if (conditionOnCoalescentEventsInput.get()) {
	        if (event.time<(lociMRCA*maxHeightRatioInput.get())) {   
	    		for (int i = 0; i < intervals.recBP.length; i++) {
	                result += -recRates[i] * prevEvent.totalRecombinationObsProb[i]
	                        * (event.time - prevEvent.time);
	            }  
	        }else if (prevEvent.time<(lociMRCA*maxHeightRatioInput.get())) {
	        	
	    		for (int i = 0; i < intervals.recBP.length; i++) {
	                result += -recRates[i] * prevEvent.totalRecombinationObsProb[i]
	                        * (lociMRCA*maxHeightRatioInput.get() - prevEvent.time);
	            }  
        		for (int i = 0; i < intervals.recBP.length; i++) {        			
                    result += -redFactor*recRates[i] * prevEvent.totalRecombinationObsProb[i]
                            * (event.time - lociMRCA*maxHeightRatioInput.get());

                }	        
	        }else {
        		for (int i = 0; i < intervals.recBP.length; i++) {        			
                    result += -redFactor*recRates[i] * prevEvent.totalRecombinationObsProb[i]
                            * (event.time - prevEvent.time);

                }        		
	        }
	   }else {
    		for (int i =0; i < intervals.recBP.length; i++) {
                result += -recRates[i] * prevEvent.totalRecombinationObsProb[i]
                        * (event.time - prevEvent.time);
            }       	
		}
		
        result += -0.5*prevEvent.lineages*(prevEvent.lineages-1)
                * populationFunction.getIntegral(prevEvent.time, event.time);
		
		return result;
	}	
	
//    @Override
//    protected boolean requiresRecalculation() {
//    	if (networkIntervalsInput.get().isDirtyCalculation())
//    		return true;
//    	
//        if (((CalculationNode) populationFunctionInput.get()).isDirtyCalculation())
//			return true;
//                
//    	if (relativeRecombinationRateInput.get()==null)
//            if (relativeRecombinationRateInput.get().isDirtyCalculation())
//            	return true;        
//
//    	if (((RealParameter) recombinationRateInput.get()).isDirtyCalculation())
//			return true;
//
//    	return super.requiresRecalculation();       
//    }
//	
}
