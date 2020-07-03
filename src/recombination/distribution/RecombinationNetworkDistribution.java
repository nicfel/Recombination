package recombination.distribution;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import coalre.distribution.NetworkIntervals;
import beast.core.State;

import java.util.List;
import java.util.Random;

public class RecombinationNetworkDistribution extends Distribution {
	
    public Input<RecombinationNetworkIntervals> networkIntervalsInput = new Input<>("networkIntervals",
            "Structured Intervals for a phylogenetic beast tree", Validate.REQUIRED);

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }

}
