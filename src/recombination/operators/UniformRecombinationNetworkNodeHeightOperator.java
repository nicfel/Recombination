package recombination.operators;

import beast.base.util.Randomizer;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.util.List;
import java.util.stream.Collectors;

public class UniformRecombinationNetworkNodeHeightOperator extends RecombinationNetworkOperator {

    @Override
    public double networkProposal() {
        List<RecombinationNetworkNode> networkNodes = network.getNodes().stream()
                .filter(n -> !n.isLeaf())
                .filter(n -> !n.getParentEdges().get(0).isRootEdge())
                .collect(Collectors.toList());

        if (networkNodes.isEmpty())
            return Double.NEGATIVE_INFINITY;

        RecombinationNetworkNode node = networkNodes.get(Randomizer.nextInt(networkNodes.size()));

        double maxHeight = Double.POSITIVE_INFINITY;
        for (RecombinationNetworkEdge parentEdge : node.getParentEdges())
            if (parentEdge.parentNode.getHeight() < maxHeight)
                maxHeight = parentEdge.parentNode.getHeight();

        double minHeight = Double.NEGATIVE_INFINITY;
        for (RecombinationNetworkEdge childEdge : node.getChildEdges())
            if (childEdge.childNode.getHeight() > minHeight)
                minHeight = childEdge.childNode.getHeight();

        network.startEditing(this);
        node.setHeightFilty(minHeight + Randomizer.nextDouble()*(maxHeight-minHeight));

        return 0.0;
    }

}
