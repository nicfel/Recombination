package recombination.operators;

import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Package;
import beast.util.Randomizer;
import recombination.distribution.CoalescentWithRecombination;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;
import recombination.statistics.RecombinationNetworkStatsLogger;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class GibbsOperatorAboveLociRoots extends RecombinationNetworkOperator {

    public Input<CoalescentWithRecombination> coalescentDistrInput = new Input<>("coalescentWithRecombination",
            "Mean of exponential used for choosing root attachment times.",
            Input.Validate.REQUIRED);
    
    
    private CoalescentWithRecombination coalDistr;

    @Override
    public void initAndValidate() {    	
    	coalDistr = coalescentDistrInput.get();
        super.initAndValidate();
    }

    @Override
    public double networkProposal() {
    	return resimulate();
    	
    }

    double resimulate() {
    	network.startEditing(this);
    	// get the place where to cut
    	double maxHeight = RecombinationNetworkStatsLogger.getMaxLociMRCA(network);
    	
    	if (maxHeight == network.getRootEdge().childNode.getHeight())
    		return Double.NEGATIVE_INFINITY;

    	// get all network edges 
        List<RecombinationNetworkEdge> networkEdges = new ArrayList<>(network.getEdges());
        
        for (RecombinationNetworkEdge e : networkEdges) 
    		e.visited = false;        


        // keep only those that coexist at the time of maxHeight
        List<RecombinationNetworkEdge> startingEdges = networkEdges.stream()
                .filter(e -> !e.isRootEdge())
                .filter(e -> e.parentNode.getHeight()>maxHeight)
                .filter(e -> e.childNode.getHeight()<=maxHeight)
               .collect(Collectors.toList());
        

        if (startingEdges.size()==0) {
        	System.out.println(maxHeight);
        	System.out.println(network);
        	throw new IllegalArgumentException("should not arrive here");
        }
        
        
        
                
        // simulate the rest of the network starting from mxHeight
        double currentTime = maxHeight;        
        
        double[] recRates = new double[coalDistr.intervals.recBP.length];
        
    	if (recRates.length>1)
    		for (int i = 0; i < recRates.length; i++)
    			recRates[i] = Math.exp(coalDistr.relativeRecombinationRate.getArrayValue(i)) * coalDistr.recombinationRate.getArrayValue()*(coalDistr.intervals.recBP[i].getLength()-1);
		else
			recRates[0] = coalDistr.recombinationRate.getArrayValue()*(networkInput.get().totalLength-1);

        double sumRates = 0;
        
        for (int i = 0; i < recRates.length;i++)
        	sumRates+=recRates[i];
        
        double recChangeTime = maxHeight*coalDistr.maxHeightRatioInput.get(); 
        

        
        do {

            // get the current propensities
            int k = startingEdges.size();
            
            double currentTransformedTime = coalDistr.populationFunction.getIntensity(currentTime);
            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
            double timeToNextCoal = coalDistr.populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;
            
            double timeToNextReass = k>=1 ? Randomizer.nextExponential(k * sumRates) : Double.POSITIVE_INFINITY;

            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextReass);
            if ((timeUntilNextEvent+currentTime)>recChangeTime) {
            	            	
	            currentTime = recChangeTime;
	            sumRates *= coalDistr.redFactor;   

	            for (int i = 0; i < recRates.length;i++)
	            	recRates[i] *= coalDistr.redFactor;
	            
	            recChangeTime = Double.POSITIVE_INFINITY;

            }else {
	            currentTime += timeUntilNextEvent;
	            if (timeUntilNextEvent == timeToNextCoal) {
	                coalesce(currentTime, startingEdges);
	            }else {
            		recombine(currentTime, startingEdges, recRates, sumRates);
	            }
            }

        }
        while (startingEdges.size() > 1);
        
        network.setRootEdge(startingEdges.get(0));
        
        return Double.POSITIVE_INFINITY;
    }

    
    private void coalesce(double coalescentTime, List<RecombinationNetworkEdge> extantLineages) {
        // Sample the pair of lineages that are coalescing:
    	RecombinationNetworkEdge lineage1 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
    	RecombinationNetworkEdge lineage2;
        do {
            lineage2 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
        } while (lineage1 == lineage2);
        
        // Create coalescent node
        RecombinationNetworkNode coalescentNode = new RecombinationNetworkNode(network.nodeEdgeIDs);
        coalescentNode.setHeight(coalescentTime)
                .addChildEdge(lineage1)
                .addChildEdge(lineage2);
        lineage1.parentNode = coalescentNode;
        lineage2.parentNode = coalescentNode;       

        
        // Merge segment flags:
        BreakPoints breakPoints = lineage1.breakPoints.copy();
        breakPoints.or(lineage2.breakPoints);

        // Create new lineage
        RecombinationNetworkEdge lineage = new RecombinationNetworkEdge(null, coalescentNode, breakPoints, null, network.nodeEdgeIDs);
        coalescentNode.addParentEdge(lineage);

        extantLineages.remove(lineage1);
        extantLineages.remove(lineage2);
        extantLineages.add(lineage);
        
        

    }

    private void recombine(double reassortmentTime, List<RecombinationNetworkEdge> extantLineages, double[] recRates, double sumRates) {
    	RecombinationNetworkEdge lineage = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
    	
    	// sample on which section the breakpoint will have occured.
    	double cumsum = recRates[0]/sumRates; 
		double rand = Randomizer.nextDouble();
		int section = 0;
		while (section < recRates.length-1) {
			if (rand<=cumsum)
				break;
			section++;
			cumsum += recRates[section]/sumRates;
		}    	
    	int breakpoint = Randomizer.nextInt(coalDistr.intervals.recBP[section].getLengthInt()-1) + coalDistr.intervals.recBP[section].getMin();
    	
    	// check if this breakpoint on this lineage would lead to a recombination event that can be observed
    	if (!lineage.breakPoints.withinLimits(breakpoint)) {
    		return;
    	}    	
    	
    	lineage.breakPoints.computeLeftAndRight(breakpoint);
    	
        // Create reassortment node
        RecombinationNetworkNode node = new RecombinationNetworkNode(network.nodeEdgeIDs);
        node.setHeight(reassortmentTime).addChildEdge(lineage);

        // Create reassortment lineages
        RecombinationNetworkEdge leftLineage = new RecombinationNetworkEdge(null, node, lineage.breakPoints.getLeft(), new BreakPoints(0,breakpoint), network.nodeEdgeIDs);
        RecombinationNetworkEdge rightLineage = new RecombinationNetworkEdge(null, node, lineage.breakPoints.getRight(), new BreakPoints(breakpoint+1, network.totalLength-1), network.nodeEdgeIDs);
            
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
