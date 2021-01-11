package recombination.util;

import beast.evolution.likelihood.ThreadedTreeLikelihood;

public class DummyTreeLikelihood extends ThreadedTreeLikelihood {
	
    @Override
    public void initAndValidate() {
    }
    
    @Override
    public double calculateLogP() {
		return 0.0;
    }
	
    @Override
    protected boolean requiresRecalculation() {
		return false;
    }
    
    @Override
    public void store() {
    }

    @Override
    public void restore() {
    }



}
