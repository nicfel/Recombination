package recombination.distribution;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.evolution.tree.coalescent.PopulationFunction;
import coalre.distribution.NetworkDistribution;
import coalre.distribution.NetworkEvent;
import coalre.distribution.NetworkIntervals;

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

	public Input<PopulationFunction> populationFunctionInput = new Input<>(
	        "populationModel",
            "Population model.",
            Input.Validate.REQUIRED);

    private PopulationFunction populationFunction;
    private Function recombinationRate;
    private RecombinationNetworkIntervals intervals;
    

    @Override
    public void initAndValidate(){
        populationFunction = populationFunctionInput.get();
        recombinationRate = recombinationRateInput.get();
        intervals = networkIntervalsInput.get();
    }

    public double calculateLogP() {
    	logP = 0;

    	// Calculate tree intervals
    	List<RecombinationNetworkEvent> networkEventList = intervals.getRecombinationNetworkEventList();

    	RecombinationNetworkEvent prevEvent = null;

    	for (RecombinationNetworkEvent event : networkEventList) {
        	if (prevEvent != null)
        		logP += intervalContribution(prevEvent, event);

        	switch (event.type) {
				case COALESCENCE:
					logP += coalesce(event);
					break;

				case SAMPLE:
					break;

				case RECOMBINATION:
					logP += recombination(event);
					break;
			}

       		if (logP==Double.NEGATIVE_INFINITY)
       			break;

        	prevEvent = event;
        }
    	
		return logP;
    }
    
	private double recombination(RecombinationNetworkEvent event) {
        return Math.log(recombinationRate.getArrayValue() * event.lociToSort);
	}

	private double coalesce(RecombinationNetworkEvent event) {
		return Math.log(1.0/populationFunction.getPopSize(event.time));
	}

	private double intervalContribution(RecombinationNetworkEvent prevEvent, RecombinationNetworkEvent event) {
        double result = 0.0;
        
        result += -recombinationRate.getArrayValue() * prevEvent.totalRecombinationObsProb
                * (event.time - prevEvent.time);          
        
		result += -0.5*prevEvent.lineages*(prevEvent.lineages-1)
                * populationFunction.getIntegral(prevEvent.time, event.time);
		
		return result;
	}
}
