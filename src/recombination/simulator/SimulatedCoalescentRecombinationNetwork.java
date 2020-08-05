package recombination.simulator;

import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;

public class SimulatedCoalescentRecombinationNetwork extends RecombinationNetwork {

    public Input<RealParameter> recombinationRateInput = new Input<>("recombinationRate",
            "Rate of recombination (per lineage per unit time)", Validate.REQUIRED);

    public Input<Function> binomialProbInput = new Input<>("binomialProb",
            "Probability parameter in binomial reassortment distribution.");

    public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
            "Population model to use.", Validate.REQUIRED);

    public Input<TraitSet> traitSetInput = new Input<>("traitSet",
            "Trait set used to assign leaf ages.");

    public Input<TaxonSet> taxonSetInput = new Input<>("taxonSet",
            "Taxon set used to define leaves");

    public Input<String> fileNameInput = new Input<>("fileName",
            "Name of file to write simulated network to.");
    
    public Input<Integer> totalLengthInput = new Input<>("totalLength",
            "total length of the alignment", Validate.REQUIRED);


    private PopulationFunction populationFunction;
    private RealParameter recombinationRate;
    private Function binomialProb;

    public void initAndValidate() {

        populationFunction = populationFunctionInput.get();
        binomialProb = binomialProbInput.get();
        totalLength = totalLengthInput.get();
        
        recombinationRate = recombinationRateInput.get();

        if (totalLength < 2) {
            throw new IllegalArgumentException("the length of the alignment has to be at least 2");
        }

        // Set up sample nodes:

        List<RecombinationNetworkNode> sampleNodes = new ArrayList<>();
    

        TaxonSet taxonSet = null;
        if (traitSetInput.get() != null)
            taxonSet = traitSetInput.get().taxaInput.get();
        else if (taxonSetInput.get() != null)
            taxonSet = taxonSetInput.get();
        else
            throw new IllegalArgumentException("Taxon set must be specified " +
                    "using either taxonSet, traitSet or provided by a segmentTree input.");

        TraitSet traitSet = null;
        if (traitSetInput.get() != null)
            traitSet = traitSetInput.get();

        for (int taxonIndex=0; taxonIndex<taxonSet.getTaxonCount(); taxonIndex++) {
            String taxonName = taxonSet.getTaxonId(taxonIndex);

            RecombinationNetworkNode sampleNode = new RecombinationNetworkNode();
            sampleNode.setTaxonLabel(taxonName);
            sampleNode.setTaxonIndex(taxonIndex);

            if (traitSet != null)
                sampleNode.setHeight(traitSet.getValue(taxonName));
            else
                sampleNode.setHeight(0.0);

            sampleNodes.add(sampleNode);
        }

        // Perform network simulation:
        simulateRecombinationNetwork(sampleNodes);

        // Write simulated network to file if requested
        if (fileNameInput.get() != null) {
            try (PrintStream ps = new PrintStream(fileNameInput.get())) {

                ps.println(toString());

            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error writing to output file '"
                        + fileNameInput.get() + "'.");
            }
        }

        super.initAndValidate();
    }

    private double getBinomialProb() {
        return binomialProb != null
                ? binomialProb.getArrayValue()
                : 0.5;
    }

    /**
     * Simulate network under coalescent with reassortment model.
     * @param sampleNodes network nodes corresponding to samples.
     */
    public void simulateRecombinationNetwork(List<RecombinationNetworkNode> sampleNodes) {

        List<RecombinationNetworkNode> remainingSampleNodes = new ArrayList<>(sampleNodes);
        List<RecombinationNetworkEdge> extantLineages = new ArrayList<>();

        remainingSampleNodes.sort(Comparator.comparingDouble(RecombinationNetworkNode::getHeight));

        double currentTime = 0;
        double timeUntilNextSample;
        do {
            // get the timing of the next sampling event
            if (!remainingSampleNodes.isEmpty()) {
                timeUntilNextSample = remainingSampleNodes.get(0).getHeight() - currentTime;
            } else {
                timeUntilNextSample = Double.POSITIVE_INFINITY;
            }

            // get the current propensities
            int k = extantLineages.size();

            double currentTransformedTime = populationFunction.getIntensity(currentTime);
            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
            double timeToNextCoal = populationFunction.getInverseIntensity(
                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;

            double timeToNextReass = k>=1 ? Randomizer.nextExponential(k*recombinationRate.getValue()*(totalLength-1)) : Double.POSITIVE_INFINITY;
            
            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextReass);
            if (timeUntilNextEvent < timeUntilNextSample) {
                currentTime += timeUntilNextEvent;
                if (timeUntilNextEvent == timeToNextCoal) {
                    coalesce(currentTime, extantLineages);
                }else {
                    recombine(currentTime, extantLineages);
                }
            } else {
                currentTime += timeUntilNextSample;
                sample(remainingSampleNodes, extantLineages);
            }

        }
        while (extantLineages.size() > 1 || !remainingSampleNodes.isEmpty());

        setRootEdge(extantLineages.get(0));
        this.totalLength = totalLengthInput.get();
    }

    private void sample(List<RecombinationNetworkNode> remainingSampleNodes, List<RecombinationNetworkEdge> extantLineages) {
        // sample the network node
    	RecombinationNetworkNode n = remainingSampleNodes.get(0);

        // Create corresponding lineage
        BreakPoints breakPoints = new BreakPoints(totalLength);        
        
        RecombinationNetworkEdge lineage = new RecombinationNetworkEdge(null, n, breakPoints);
        
        extantLineages.add(lineage);
        n.addParentEdge(lineage);

        remainingSampleNodes.remove(0);
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
    	    	
    	int breakpoint = Randomizer.nextInt(totalLength-1);
    	
    	// check if this breakpoint on this lineage would lead to a recombination event that can be observed
    	if (!lineage.breakPoints.withinLimits(breakpoint))
    		return;
    	
    	
    	lineage.breakPoints.computeLeftAndRight(breakpoint);
    	
        // Create reassortment node
        RecombinationNetworkNode node = new RecombinationNetworkNode();
        node.setHeight(reassortmentTime).addChildEdge(lineage);

        // Create reassortment lineages
        RecombinationNetworkEdge leftLineage = new RecombinationNetworkEdge(null, node, lineage.breakPoints.getLeft());
        RecombinationNetworkEdge rightLineage = new RecombinationNetworkEdge(null, node, lineage.breakPoints.getRight());
            
        // add the breakPoints to the edges
        leftLineage.setPassingRange(0, breakpoint);
        rightLineage.setPassingRange(breakpoint+1, totalLength-1);
        
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

        extantLineages.remove(lineage);
        extantLineages.add(leftLineage);
        extantLineages.add(rightLineage);        
    }

}
