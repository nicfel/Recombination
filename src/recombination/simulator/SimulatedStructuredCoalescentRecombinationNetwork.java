package recombination.simulator;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;
import recombination.util.NodeEdgeID;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class SimulatedStructuredCoalescentRecombinationNetwork extends RecombinationNetwork {

    public Input<RealParameter> recombinationRateInput = new Input<>("recombinationRate",
            "Rate of recombination (per lineage per unit time)", Validate.REQUIRED);

    public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
            "Population model to use.", Validate.REQUIRED);

    public Input<TraitSet> traitSetInput = new Input<>("traitSet",
            "Trait set used to assign leaf ages.");

    public Input<TaxonSet> taxonSetInput = new Input<>("taxonSet",
            "Taxon set used to define leaves");
    
    public Input<Tree> treeInput = new Input<>("tree",
            "tree from where to get taxa etc.");


    public Input<String> fileNameInput = new Input<>("fileName",
            "Name of file to write simulated network to.");
    
    public Input<Integer> totalLengthInput = new Input<>("totalLength",
            "total length of the alignment");
    
    public Input<Alignment> dataInput = new Input<>("data",
            "total length of the alignment");
    
    public Input<Boolean> IgnoreUnknownInput = new Input<>("IgnoreUnknown",
            "if true, N and gap positions are completely ignored", false);
    
    public Input<Boolean> conditionCoalescenceInput = new Input<>("conditionCoalescence",
            "if true, only coalescent events can happen after all loci have reached their mrca", false);

	public Input<String> recombinationRatesChangePointsInput = new Input<>(
	        "recombinationRatesChangePoints",
            "if true, only coalescent events are allowed after the .",
            Input.Validate.OPTIONAL);
	
	public Input<RealParameter> relativeRecombinationRateInput = new Input<>(
	        "relativeRecombinationRate",
            "relative recombination rate (per lineage per unit time) for a specific part of the genome",
            Input.Validate.OPTIONAL);

    private PopulationFunction populationFunction;
    private RealParameter recombinationRate;
    public BreakPoints[] recBP;


    public void initAndValidate() {
    	nodeEdgeIDs = new NodeEdgeID();

        populationFunction = populationFunctionInput.get();
        if (dataInput.get()!=null)
        	totalLength = dataInput.get().getSiteCount();
        else
        	totalLength = totalLengthInput.get();
        
        recombinationRate = recombinationRateInput.get();

        if (totalLength < 2) {
            throw new IllegalArgumentException("the length of the alignment has to be at least 2");
        }

        // Set up sample nodes:

        List<RecombinationNetworkNode> sampleNodes = new ArrayList<>();
    

        TaxonSet taxonSet = null;
        
        if (traitSetInput.get() != null) {
            taxonSet = traitSetInput.get().taxaInput.get();
        }else if (taxonSetInput.get() != null) {
            taxonSet = taxonSetInput.get();
    	}else if (treeInput.get()!=null) {
            taxonSet = treeInput.get().m_taxonset.get();
        }
        else {
            throw new IllegalArgumentException("Taxon set must be specified " +
                    "using either taxonSet, traitSet or provided by a tree input.");
        }

        TraitSet traitSet = null;
        if (traitSetInput.get() != null) {
            traitSet = traitSetInput.get();
        }else if (treeInput.get()!=null && treeInput.get().hasDateTrait()) {	
        	traitSet = treeInput.get().getDateTrait();
        }

        for (int taxonIndex=0; taxonIndex<taxonSet.getTaxonCount(); taxonIndex++) {
            String taxonName = taxonSet.getTaxonId(taxonIndex);

            RecombinationNetworkNode sampleNode = new RecombinationNetworkNode(nodeEdgeIDs);
            sampleNode.setTaxonLabel(taxonName);
            sampleNode.setTaxonIndex(taxonIndex);

            if (traitSet != null)
                sampleNode.setHeight(traitSet.getValue(taxonName));
            else
                sampleNode.setHeight(0.0);

            sampleNodes.add(sampleNode);
        }
        
        if (recombinationRatesChangePointsInput.get()!=null) {
        	
        	String[] tmp = recombinationRatesChangePointsInput.get().trim().split("\\s+");
        	recBP = new BreakPoints[tmp.length+1];
        	List<Integer> bpList = new ArrayList<>();
        	bpList.add(0);
        	for (int i = 0; i < tmp.length; i++) {            	
            	bpList.add(Integer.parseInt(tmp[i]));
            	bpList.add(Integer.parseInt(tmp[i]));
        		if (Integer.parseInt(tmp[i])>totalLength) {
            		throw new IllegalArgumentException("recombination rate break point value is larger than the network size");
        		}
        	}   
        	bpList.add(totalLength-1);
        	
        	for (int i = 0; i < recBP.length; i++) 
        		recBP[i] = new BreakPoints(bpList.get(2*i), bpList.get(2*i+1));
        }else {
        	recBP = new BreakPoints[1];
    		recBP[0] = new BreakPoints(0, totalLength-1);
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

    /**
     * Simulate network under coalescent with reassortment model.
     * @param sampleNodes network nodes corresponding to samples.
     */
    public void simulateRecombinationNetwork(List<RecombinationNetworkNode> sampleNodes) {

        List<RecombinationNetworkNode> remainingSampleNodes = new ArrayList<>(sampleNodes);
        List<RecombinationNetworkEdge> extantLineages = new ArrayList<>();

        remainingSampleNodes.sort(Comparator.comparingDouble(RecombinationNetworkNode::getHeight));
        
        
        double[] recRates = new double[recBP.length];
        
    	if (recRates.length>1)
    		for (int i = 0; i < recRates.length; i++)
    			recRates[i] = Math.exp(relativeRecombinationRateInput.get().getArrayValue(i)) * recombinationRate.getArrayValue()*(recBP[i].getLength()-1);
		else
			recRates[0] = recombinationRate.getArrayValue()*(recBP[0].getLength()-1);

        double sumRates = 0;
        
        for (int i = 0; i < recRates.length;i++)
        	sumRates+=recRates[i];
        
//        double recChangeTime = maxHeight*coalDistr.maxHeightRatioInput.get(); 


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

//            double totObsProb = 0;
//    		for (int i = 0; i < extantLineages.size(); i++)
//    			totObsProb += extantLineages.get(i).breakPoints.getLength()-1;
    			
            double timeToNextReass = k>=1 ? Randomizer.nextExponential(k*sumRates) : Double.POSITIVE_INFINITY;
            
            boolean allowRecomb = true;
            if (timeUntilNextSample == Double.POSITIVE_INFINITY && conditionCoalescenceInput.get()) {
            	allowRecomb = getOverlap(new BreakPoints(), extantLineages);
            }
            
            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextReass);
            if (timeUntilNextEvent < timeUntilNextSample) {
                currentTime += timeUntilNextEvent;
                if (timeUntilNextEvent == timeToNextCoal) {
                    coalesce(currentTime, extantLineages);
                }else if (allowRecomb){
                    recombine(currentTime, extantLineages, recRates, sumRates);
                }
            } else {
                currentTime += timeUntilNextSample;
                sample(remainingSampleNodes, extantLineages);
            }

        }
        while (extantLineages.size() > 1 || !remainingSampleNodes.isEmpty());

        setRootEdge(extantLineages.get(0));
    }

    private boolean getOverlap(BreakPoints breakPoints, List<RecombinationNetworkEdge> extantLineages) {
		for (RecombinationNetworkEdge edge : extantLineages) {
			BreakPoints cp = breakPoints.copy();
			cp.and(edge.breakPoints);
			if (!cp.isEmpty()) {
				return true;				
			}
			breakPoints.or(edge.breakPoints);
		}
		return false;
	}

	private void sample(List<RecombinationNetworkNode> remainingSampleNodes, List<RecombinationNetworkEdge> extantLineages) {
        // sample the network node
    	RecombinationNetworkNode n = remainingSampleNodes.get(0);

        // Create corresponding lineage
        BreakPoints breakPoints;
        if (IgnoreUnknownInput.get() && dataInput.get()!=null)
        	breakPoints= getGappedBreakPoints(dataInput.get(), n);
        else
        	breakPoints= new BreakPoints(totalLength);

        
        RecombinationNetworkEdge lineage = new RecombinationNetworkEdge(null, n, breakPoints, null, nodeEdgeIDs);
        
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
        RecombinationNetworkNode coalescentNode = new RecombinationNetworkNode(nodeEdgeIDs);
        coalescentNode.setHeight(coalescentTime)
                .addChildEdge(lineage1)
                .addChildEdge(lineage2);
        lineage1.parentNode = coalescentNode;
        lineage2.parentNode = coalescentNode;       

        
        // Merge segment flags:
        BreakPoints breakPoints = lineage1.breakPoints.copy();
        breakPoints.or(lineage2.breakPoints);


        // Create new lineage
        RecombinationNetworkEdge lineage = new RecombinationNetworkEdge(null, coalescentNode, breakPoints, null, nodeEdgeIDs);
        coalescentNode.addParentEdge(lineage);

        extantLineages.remove(lineage1);
        extantLineages.remove(lineage2);
        extantLineages.add(lineage);       
    }

    private void recombine(double reassortmentTime, List<RecombinationNetworkEdge> extantLineages) {
    	RecombinationNetworkEdge lineage = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
    	    	
    	int breakpoint = Randomizer.nextInt(totalLength-1);
    	
    	// check if this breakpoint on this lineage would lead to a recombination event that can be observed
    	if (!lineage.breakPoints.withinLimits(breakpoint)) {
    		return;
    	}    	
    	
    	lineage.breakPoints.computeLeftAndRight(breakpoint);
    	
        // Create reassortment node
        RecombinationNetworkNode node = new RecombinationNetworkNode(nodeEdgeIDs);
        node.setHeight(reassortmentTime).addChildEdge(lineage);

        // Create reassortment lineages
        RecombinationNetworkEdge leftLineage = new RecombinationNetworkEdge(null, node, lineage.breakPoints.getLeft(), new BreakPoints(0,breakpoint), nodeEdgeIDs);
        RecombinationNetworkEdge rightLineage = new RecombinationNetworkEdge(null, node, lineage.breakPoints.getRight(), new BreakPoints(breakpoint+1, totalLength-1), nodeEdgeIDs);
            
        // add the breakPoints to the edges
        leftLineage.setPassingRange(0, breakpoint);
        rightLineage.setPassingRange(breakpoint+1, totalLength-1);
        
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

        extantLineages.remove(lineage);
        extantLineages.add(leftLineage);
        extantLineages.add(rightLineage);        
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
    	
    	int breakpoint = Randomizer.nextInt(recBP[section].getLengthInt()-1) + recBP[section].getMin();
    	
    	// check if this breakpoint on this lineage would lead to a recombination event that can be observed
    	if (!lineage.breakPoints.withinLimits(breakpoint)) {
    		return;
    	}    	
    	
    	lineage.breakPoints.computeLeftAndRight(breakpoint);
    	
        // Create reassortment node
        RecombinationNetworkNode node = new RecombinationNetworkNode(nodeEdgeIDs);
        node.setHeight(reassortmentTime).addChildEdge(lineage);

        // Create reassortment lineages
        RecombinationNetworkEdge leftLineage = new RecombinationNetworkEdge(null, node, lineage.breakPoints.getLeft(), new BreakPoints(0,breakpoint), nodeEdgeIDs);
        RecombinationNetworkEdge rightLineage = new RecombinationNetworkEdge(null, node, lineage.breakPoints.getRight(), new BreakPoints(breakpoint+1, totalLength-1), nodeEdgeIDs);
            
        // add the breakPoints to the edges
        leftLineage.setPassingRange(0, breakpoint);
        rightLineage.setPassingRange(breakpoint+1, totalLength-1);
        
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

        extantLineages.remove(lineage);
        extantLineages.add(leftLineage);
        extantLineages.add(rightLineage);        
    }
    


    private BreakPoints getGappedBreakPoints(Alignment data, RecombinationNetworkNode n) {
        int taxonIndex = data.getTaxonIndex(n.getTaxonLabel());
        if (taxonIndex == -1) {
        	if (n.getTaxonLabel().startsWith("'") || n.getTaxonLabel().startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(n.getTaxonLabel().substring(1, n.getTaxonLabel().length() - 1));
            }
            if (taxonIndex == -1) {
            	throw new RuntimeException("Could not find sequence " + n.getTaxonLabel() + " in the alignment");
            }
        }
        
        int code, states;
        int[] statesForCode;
        boolean on = false;
        List<Integer> bp_list = new ArrayList<>();
    	for (int i = 0; i < data.getSiteCount() ;i++) {
            code = data.getPattern(taxonIndex, dataInput.get().getPatternIndex(i));
            
            statesForCode = data.getDataType().getStatesForCode(code);
            if (statesForCode.length==1)
                states = statesForCode[0];
            else
                states = code; // Causes ambiguous states to be ignored.

            if (i==0 && states<data.getMaxStateCount()) {
            	bp_list.add(i);
            	on = true;
            }
            
            if (on && states>=data.getMaxStateCount()) {
            	bp_list.add(i-1);
            	on = false;
            }      
            
            if (!on && states<data.getMaxStateCount()) {
            	bp_list.add(i);
            	on = true;
            }           
            	
    	}
    	
    	if (on) {
        	bp_list.add(data.getSiteCount()-1);
    	}    	
    	
    	
    	BreakPoints bp = new BreakPoints();
    	bp.init(bp_list);
		return bp;
	}

    
}
