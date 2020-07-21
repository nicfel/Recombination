package recombination.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.math.Binomial;
import beast.util.Randomizer;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

public abstract class RecombinationNetworkOperator extends Operator {

    public Input<RecombinationNetwork> networkInput = new Input<>("network",
            "Network on which to operate",
            Input.Validate.REQUIRED);

    public RecombinationNetwork network;
    List<Tree> segmentTrees;

    @Override
    public void initAndValidate() {
        network = networkInput.get();
    }

    public double proposal() {
        double logHR = networkProposal();

        return logHR;
    }

    /**
     * Propose a new network state.  The network state will be
     * used to update the segment trees once the network proposal
     * is complete.
     *
     * @return log of HR for proposal.
     */
    protected abstract double networkProposal();

    /**
     * Retrieve sister of given edge
     * @param childEdge child edge
     * @return sister of given child edge
     */
    protected RecombinationNetworkEdge getSisterEdge(RecombinationNetworkEdge childEdge) {
        int idx = childEdge.parentNode.getChildEdges().indexOf(childEdge);
        int otherIdx = (idx + 1) % 2;

        return childEdge.parentNode.getChildEdges().get(otherIdx);
    }

    /**
     * Retrieve spouse of given edge
     * @param parentEdge parent edge
     * @return spouse of given parent edge
     */
    protected RecombinationNetworkEdge getSpouseEdge(RecombinationNetworkEdge parentEdge) {
        int idx = parentEdge.childNode.getParentEdges().indexOf(parentEdge);
        int otherIdx = (idx + 1) % 2;

        return parentEdge.childNode.getParentEdges().get(otherIdx);
    }
    
    public void replace(final RecombinationNetworkNode node, final RecombinationNetworkEdge child, final RecombinationNetworkEdge replacement) {
    	node.removeChildEdge(child);
    	node.addChildEdge(replacement);
    }

    /**
     * Will be used in the next version of BEAST to prevent date trait cloning
     * from breaking the BEAuti model.
     */
    public boolean notCloneable() {
        return true;
    }
}
