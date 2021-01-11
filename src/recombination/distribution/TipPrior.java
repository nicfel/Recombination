package recombination.distribution;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.math.distributions.ParametricDistribution;
import coalre.network.Network;
import coalre.network.NetworkNode;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkNode;


@Description("Prior over set of taxa, useful for defining monophyletic constraints and "
        + "distributions over MRCA times or (sets of) tips of trees")
public class TipPrior extends Distribution {
    public final Input<RecombinationNetwork> networkInput = new Input<>("network", "the network containing the taxon set", Validate.REQUIRED);
    public final Input<List<TaxonSet>> taxonsetInput = new Input<>("taxonset",
            "set of taxa for which prior information is available", new ArrayList<>());
    public final Input<Boolean> isMonophyleticInput = new Input<>("monophyletic",
            "whether the taxon set is monophyletic (forms a clade without other taxa) or nor. Default is false.", false);
    public final Input<List<ParametricDistribution>> distInput = new Input<>("distr",
            "distribution used to calculate prior over MRCA time, "
                    + "e.g. normal, beta, gamma. If not specified, monophyletic must be true", new ArrayList<>());
    public Input<Double> dateOffsetInput = new Input<>("dateOffset", 
    		"keeps track of how much the dates have change", Validate.REQUIRED);


    /**
     * shadow members *
     */
    ParametricDistribution dist;
    RecombinationNetwork network;
    
    // number of taxa in taxon set
    int nrOfTaxa = -1;
    // array of flags to indicate which taxa are in the set
    Set<String> isInTaxaSet = new LinkedHashSet<>();

    // array of indices of taxa
    int[] taxonIndex;
    // stores time to be calculated
    double MRCATime = -1;
    double storedMRCATime = -1;
    
    boolean initialised = false;
    
    RecombinationNetworkNode operatingNode;
    Double dateOffset;
    
    List<String> tipNames;

    @Override
    public void initAndValidate() {
//        dist = distInput.get();
        network = networkInput.get();
        
        getTipNames();
    }



    // A lightweight version for finding the most recent common ancestor of a group of taxa.
    // return the node-ref of the MRCA.

    @Override
    public double calculateLogP() {
    	if (tipNames==null || tipNames.size()==0) {
    		getTipNames();
    	}
        logP = 0;
        int i = 0;
    	for (final RecombinationNetworkNode taxon : network.getLeafNodes()) {
    		int index = tipNames.indexOf(taxon.getTaxonLabel());
    		if (index!=-1) {
                logP += distInput.get().get(index).logDensity(dateOffsetInput.get() - taxon.getHeight());
                i++;
    		}
    	}
    	
    	if (i!=tipNames.size()){
    		getTipNames();
    		logP = 0;
        	for (final RecombinationNetworkNode taxon : network.getLeafNodes()) {
        		int index = tipNames.indexOf(taxon.getTaxonLabel());
        		if (index!=-1) {
                    logP += distInput.get().get(index).logDensity(dateOffsetInput.get() - taxon.getHeight());
        		}
        	}
    	}
    	    	
        return logP;
    }
    
    private void getTipNames() {
    	tipNames = new ArrayList<>();
        for (TaxonSet taxon : taxonsetInput.get()) {
        	tipNames.add(taxon.getTaxonId(0));
        }
    }
    

    @Override
    public void store() {
        super.store();
    }

    @Override
    public void restore() {
        super.restore();
    }

    @Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation();
    }


    /**
     * Loggable interface implementation follows *
     */
    @Override
    public void init(final PrintStream out) {
    	if (tipNames==null || tipNames.size()==0) {
    		getTipNames();
    	}

    	for (String s : tipNames)
    		out.print(s + ".height" + "\t");
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	if (tipNames==null || tipNames.size()==0) {
    		getTipNames();
    	}

    	double[] heights = new double[tipNames.size()];
    	for (final RecombinationNetworkNode taxon : network.getLeafNodes()) {
    		int index = tipNames.indexOf(taxon.getTaxonLabel());
    		if (index!=-1) {
    			heights[index] = dateOffsetInput.get() - taxon.getHeight();
    		}
    	}
    	for (int i = 0; i < heights.length; i++)
    		out.print(heights[i] + "\t");
    }

    @Override
    public void close(final PrintStream out) {
        // nothing to do
    }

    /**
     * Valuable interface implementation follows, first dimension is log likelihood, second the time *
     */
    @Override
    public int getDimension() {
        return 2;
    }

    @Override
    public double getArrayValue() {
    	if (Double.isNaN(logP)) {
    		try {
    			calculateLogP();
    		}catch (Exception e) {
    			logP  = Double.NaN;
    		}
    	}
        return logP;
    }

    @Override
    public void sample(final State state, final Random random) {
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }
}