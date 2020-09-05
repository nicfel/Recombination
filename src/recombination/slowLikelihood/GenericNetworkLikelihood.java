package recombination.slowLikelihood;


import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.tree.TreeInterface;
import recombination.alignment.RecombinationAlignment;
import recombination.network.RecombinationNetwork;





@Description("Generic tree likelihood for an alignment given a generic SiteModel, " +
		"a beast tree and a branch rate model")
// Use this as base class to define any non-standard TreeLikelihood.
// Override Distribution.calculatLogP() to make this class functional.
//
// TODO: This could contain a generic traverse() method that takes dirty trees in account.
//
public class GenericNetworkLikelihood extends Distribution {
    
    final public Input<RecombinationAlignment> dataInput = new Input<>("recombinationData", "sequence data for the beast.tree", Validate.REQUIRED);

    final public Input<RecombinationNetwork> networkInput = new Input<>("recombinationNetwork", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
    
    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");

    
    
	@Override
	public List<String> getArguments() {return null;}

	@Override
	public List<String> getConditions() {return null;}

	@Override
	public void sample(State state, Random random) {}

}
