package recombination.operators;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Package;
import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class GibbsOperatorAboveLociRoots extends RecombinationNetworkOperator {

    public Input<RealParameter> recombinationRateInput = new Input<>("recombinationRate",
            "Rate of recombination (per lineage per unit time)", Validate.REQUIRED);

    public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
            "Population model to use.", Validate.REQUIRED);

    
    private PopulationFunction populationFunction;
    private RealParameter recombinationRate;

    @Override
    public void initAndValidate() {    	
        populationFunction = populationFunctionInput.get();
        recombinationRate = recombinationRateInput.get();
    	
        super.initAndValidate();
    }

    @Override
    public double networkProposal() {
    	return resimulate();
    	
    }

    double resimulate() {
    	network.startEditing(this);
    	// get the place where to cut
    	double maxHeight = getMaxLociMRCA(network.getRootEdge());

    	// get all network edges 
        List<RecombinationNetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

        // keep only those that coexist at the time of maxHeight
        List<RecombinationNetworkEdge> startingEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.parentNode.getHeight()>maxHeight)
                .filter(e -> e.childNode.getHeight()<=maxHeight)
               .collect(Collectors.toList());
        
        
        if (startingEdges.size()==0)
        	return Double.NEGATIVE_INFINITY;        
                
       // simulate the rest of the network starting from mxHeight
        double currentTime = maxHeight;
        double timeUntilNextSample = Double.POSITIVE_INFINITY;
        do {

            // get the current propensities
            int k = startingEdges.size();

            double currentTransformedTime = populationFunction.getIntensity(currentTime);
            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
            double timeToNextCoal = populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;

            double timeToNextReass = k>=1 ? Randomizer.nextExponential(k*recombinationRate.getValue()*(networkInput.get().totalLength-1)) : Double.POSITIVE_INFINITY;

            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextReass);
            if (timeUntilNextEvent < timeUntilNextSample) {
                currentTime += timeUntilNextEvent;
                if (timeUntilNextEvent == timeToNextCoal)
                    coalesce(currentTime, startingEdges);
                else
                    recombine(currentTime, startingEdges);
            }

        }
        while (startingEdges.size() > 1);
        
        network.setRootEdge(startingEdges.get(0));
                
        return Double.POSITIVE_INFINITY;
    }

    double getMaxLociMRCA(RecombinationNetworkEdge edge){
    	RecombinationNetworkNode node = edge.childNode;
    	if (node.isCoalescence()) {
    		BreakPoints bp1 = node.getChildEdges().get(0).breakPoints.copy();
    		bp1.and(node.getChildEdges().get(1).breakPoints);
    		if (!bp1.isEmpty()) {
    			return node.getHeight();
    		}
    		return Math.max(getMaxLociMRCA(node.getChildEdges().get(0)), getMaxLociMRCA(node.getChildEdges().get(1)));
    	}else if (node.isRecombination()) {
    		return getMaxLociMRCA(node.getChildEdges().get(0));
    	}else {
    		return -1.0;
    	}    	
    }
    
    private void coalesce(double coalescentTime, List<RecombinationNetworkEdge> extantLineages) {
        // Sample the pair of lineages that are coalescing:
    	RecombinationNetworkEdge lineage1 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
    	RecombinationNetworkEdge lineage2;
        do {
            lineage2 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
        } while (lineage1 == lineage2);
        
        // Create coalescent node
        RecombinationNetworkNode coalescentNode = new RecombinationNetworkNode();
        coalescentNode.setHeight(coalescentTime)
                .addChildEdge(lineage1)
                .addChildEdge(lineage2);
        lineage1.parentNode = coalescentNode;
        lineage2.parentNode = coalescentNode;       

        
        // Merge segment flags:
        BreakPoints breakPoints = lineage1.breakPoints.copy();
        breakPoints.or(lineage2.breakPoints);
        


        // Create new lineage
        RecombinationNetworkEdge lineage = new RecombinationNetworkEdge(null, coalescentNode, breakPoints);
        coalescentNode.addParentEdge(lineage);

        extantLineages.remove(lineage1);
        extantLineages.remove(lineage2);
        extantLineages.add(lineage);
        
        

    }

    private void recombine(double reassortmentTime, List<RecombinationNetworkEdge> extantLineages) {
    	RecombinationNetworkEdge lineage = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
    	    	
    	int breakpoint = Randomizer.nextInt(network.totalLength-1);
    	
   	
    	// check if this breakpoint on this lineage would lead to a recombination event that can be observed
    	if (!lineage.breakPoints.withinLimits(breakpoint)) {
    		return;
    	}    	
    	
    	lineage.breakPoints.computeLeftAndRight(breakpoint);
    	
        // Create reassortment node
        RecombinationNetworkNode node = new RecombinationNetworkNode();
        node.setHeight(reassortmentTime).addChildEdge(lineage);

        // Create reassortment lineages
        RecombinationNetworkEdge leftLineage = new RecombinationNetworkEdge(null, node, lineage.breakPoints.getLeft());
        RecombinationNetworkEdge rightLineage = new RecombinationNetworkEdge(null, node, lineage.breakPoints.getRight());
            
        // add the breakPoints to the edges
        leftLineage.setPassingRange(0, breakpoint);
        rightLineage.setPassingRange(breakpoint+1, network.totalLength-1);
        
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

        extantLineages.remove(lineage);
        extantLineages.add(leftLineage);
        extantLineages.add(rightLineage);        
    }
    
}
