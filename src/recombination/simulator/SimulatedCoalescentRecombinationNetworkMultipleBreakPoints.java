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
import recombination.network.BreakPoints.Range;
import recombination.util.NodeEdgeID;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class SimulatedCoalescentRecombinationNetworkMultipleBreakPoints extends RecombinationNetwork {

    public Input<RealParameter> recombinationRateInput = new Input<>("recombinationRate",
            "Rate of recombination (per lineage per unit time)", Validate.REQUIRED);

    public Input<RealParameter> NeInput = new Input<>("Ne",
            "Ne for the different locations.", Validate.REQUIRED);
    
    public Input<RealParameter> migrationRatesInput = new Input<>("migrationRates",
            "migration rates input for the different locations", Validate.REQUIRED);
    
    public Input<TraitSet> typeTraitInput = new Input<>("typeTrait", "Type trait set. ", Validate.REQUIRED);

    
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
	
    public Input<Double> lambdaInput = new Input<>("lambda",
            "lambda for Poisson distribution of number of breakpoints", 2.0);


    private RealParameter Ne;
    private RealParameter migrationRates;
    private RealParameter recombinationRate;
    public BreakPoints[] recBP;
    private int states;


    public void initAndValidate() {
    	nodeEdgeIDs = new NodeEdgeID();

    	Ne = NeInput.get();
    	migrationRates = migrationRatesInput.get();
    	
    	states = Ne.getDimension();
    	
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
        List<RecombinationNetworkEdge>[] extantLineages = new List[states];
        for (int i = 0; i < states; i++)
        	extantLineages[i] = new ArrayList<>();

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
        int sumlins = 0;

        do {
            // get the timing of the next sampling event
            if (!remainingSampleNodes.isEmpty()) {
                timeUntilNextSample = remainingSampleNodes.get(0).getHeight() - currentTime;
            } else {
                timeUntilNextSample = Double.POSITIVE_INFINITY;
            }

            // get the current propensities
            int[] k = new int[states];

            int locNextCoal = -1;
            double timeToNextCoal = Double.POSITIVE_INFINITY;
            
            for (int i = 0; i < states; i++) {
            	k[i] = extantLineages[i].size();
            	double coaltime = Randomizer.nextExponential(k[i] * (k[i] - 1) / Ne.getArrayValue(i));
            	if (coaltime < timeToNextCoal) {
            		locNextCoal = i;
            		timeToNextCoal = coaltime;
            	}
            }         
            
            int locNextMig = -1;
            int fromMig=-1, toMig = -1;
            double timeToNextMig = Double.POSITIVE_INFINITY;
            int c = 0;
            for (int a = 0; a < states; a++) {
            	for (int b=0; b < states; b++) {
            		if (a!=b) {
		            	double migtime = Randomizer.nextExponential(k[a] * migrationRates.getArrayValue(c));
		            	if (migtime < timeToNextMig) {
		            		locNextMig = c;
		            		timeToNextMig = migtime;
		            		fromMig=a;
		            		toMig=b;		            				
		            	}
		            	c++;
            		}
            	}
            }        
            
            int locNextRecomb = -1;
            double timeToNextRecomb = Double.POSITIVE_INFINITY;
            for (int i = 0; i < states; i++) {
            	double rectime = Randomizer.nextExponential(k[i] * sumRates);
            	if (rectime < locNextRecomb) {
            		locNextRecomb = i;
            		timeToNextRecomb = rectime;
            	}
            }         

            
            boolean allowRecomb = true;
            if (timeUntilNextSample == Double.POSITIVE_INFINITY && conditionCoalescenceInput.get()) {
            	allowRecomb = getOverlap(new BreakPoints(), extantLineages);
            }
            
            // next event time
            double timeUntilNextEvent = Math.min(timeToNextCoal, Math.min(timeToNextRecomb, timeToNextMig));
            
            
            if (timeUntilNextEvent < timeUntilNextSample ) {
                currentTime += timeUntilNextEvent;
                if (timeUntilNextEvent == timeToNextCoal) {
                    coalesce(currentTime, extantLineages[locNextCoal]);
                }else if (timeUntilNextEvent == timeToNextMig) {
                	migrate(currentTime, extantLineages, fromMig, toMig);
                }else if (allowRecomb){
                    recombine(currentTime, extantLineages[locNextRecomb], recRates, sumRates);
                }
            } else {
                currentTime += timeUntilNextSample;
                sample(remainingSampleNodes, extantLineages);
            }
            
            sumlins = 0;
            for (int i = 0; i < states; i++)
            	sumlins+=extantLineages[i].size();
            		

        }
        while (sumlins > 1 || !remainingSampleNodes.isEmpty());

        for (int i = 0; i < states; i++) {
        	if(extantLineages[i].size()>0) {
                setRootEdge(extantLineages[i].get(0));               

        	}
        }

    }

    private boolean getOverlap(BreakPoints breakPoints, List<RecombinationNetworkEdge>[] extantLineages) {
    	for (int i = 0; i < states; i++) {
			for (RecombinationNetworkEdge edge : extantLineages[i]) {
				BreakPoints cp = breakPoints.copy();
				cp.and(edge.breakPoints);
				if (!cp.isEmpty()) {
					return true;				
				}
				breakPoints.or(edge.breakPoints);
			}
    	}
		return false;
	}

	private void sample(List<RecombinationNetworkNode> remainingSampleNodes, List<RecombinationNetworkEdge>[] extantLineages) {
        // sample the network node
    	RecombinationNetworkNode n = remainingSampleNodes.get(0);

        // Create corresponding lineage
        BreakPoints breakPoints;
        if (IgnoreUnknownInput.get() && dataInput.get()!=null)
        	breakPoints= getGappedBreakPoints(dataInput.get(), n);
        else
        	breakPoints= new BreakPoints(totalLength);

        
        RecombinationNetworkEdge lineage = new RecombinationNetworkEdge(null, n, breakPoints, null, nodeEdgeIDs);
        
        extantLineages[n.getTypeIndex()].add(lineage);
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
        coalescentNode.setTypeIndex(lineage1.childNode.getTypeIndex());
        
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
	
	private void migrate(double coalescentTime, List<RecombinationNetworkEdge>[] extantLineages, int fromMig, int toMig) {
        // Sample the pair of lineages that are coalescing:
    	RecombinationNetworkEdge lineage1 = extantLineages[fromMig].get(Randomizer.nextInt(extantLineages[fromMig].size()));

        // Create coalescent node
        RecombinationNetworkNode migrationNode = new RecombinationNetworkNode(nodeEdgeIDs);
        migrationNode.setHeight(coalescentTime)
                .addChildEdge(lineage1);
        
        migrationNode.setTypeIndex(toMig);
        
        lineage1.parentNode = migrationNode;

        
        // Merge segment flags:
        BreakPoints breakPoints = lineage1.breakPoints.copy();

        // Create new lineage
        RecombinationNetworkEdge lineage = new RecombinationNetworkEdge(null, migrationNode, breakPoints, null, nodeEdgeIDs);
        migrationNode.addParentEdge(lineage);

        extantLineages[fromMig].remove(lineage1);
        extantLineages[toMig].add(lineage);       
    }

    private void recombine(double reassortmentTime, List<RecombinationNetworkEdge> extantLineages, double[] recRates, double sumRates) {
    	RecombinationNetworkEdge lineage = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
    	
    	long nrBreakpoints = Randomizer.nextPoisson(lambdaInput.get());
    	
    	if (nrBreakpoints==0)
    		return;
    	
    	List<Integer> breaks = new ArrayList<>();
    	
    	for (int i =0; i < nrBreakpoints; i++) {
    	
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
	    	
			breaks.add(Randomizer.nextInt(recBP[section].getLengthInt()-1) + recBP[section].getMin());

    	}
    	
    	Collections.sort(breaks);
    	
    	List<Integer> bp = new ArrayList<>();
    	
    	bp.add(0);
    	int c=0;
    	for (Integer breakpoint : breaks) {
    		if ( c % 2==0)
    			bp.add(breakpoint);
			else
    			bp.add(breakpoint+1);
    		c++;
    	}
    	    	    	
    	if (bp.size()%2!=0)
    		bp.add(totalLength-1);
    	
    	// compute the passing ranges
    	BreakPoints pr1 = new BreakPoints();
    	pr1.init(bp);
    	
    	BreakPoints pr2 = new BreakPoints(0,totalLength-1);
    	pr2.andNot(pr1);
    	
//    	System.out.println();
//    	System.out.println(lineage.breakPoints);
//    	System.out.println(pr1);
//    	System.out.println(pr2);
//    	System.exit(0);
    	
    	BreakPoints left = lineage.breakPoints.andCopy(pr1);
		BreakPoints right = lineage.breakPoints.andCopy(pr2);
		
//		System.out.println(left);
//		System.out.println(right);
		
		if (left.isEmpty() || right.isEmpty())
			return;

    	
        // Create reassortment node
        RecombinationNetworkNode node = new RecombinationNetworkNode(nodeEdgeIDs);
        node.setHeight(reassortmentTime).addChildEdge(lineage);

        // Create reassortment lineages
        RecombinationNetworkEdge leftLineage = new RecombinationNetworkEdge(null, node, left, pr1, nodeEdgeIDs);
        RecombinationNetworkEdge rightLineage = new RecombinationNetworkEdge(null, node, right, pr2, nodeEdgeIDs);
            
        // add the breakPoints to the edges
//        leftLineage.setPassingRange(0, pr1);
//        rightLineage.setPassingRange(breakpoint+1, totalLength-1);
        
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
