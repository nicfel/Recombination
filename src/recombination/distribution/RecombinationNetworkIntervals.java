package recombination.distribution;



import beast.core.CalculationNode;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import recombination.network.RecombinationNetwork;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Nicola Felix Mueller
 */
public class RecombinationNetworkIntervals extends CalculationNode {
    public Input<RecombinationNetwork> recombinationNetworkInput = new Input<>("recombinationNetwork",
            "recombinationNetwork for which to calculate the intervals", Validate.REQUIRED);

    public Input<RealParameter> binomialProbInput = new Input<>("binomialProb",
            "Probability of a given segment choosing a particular parent.");

    private RecombinationNetwork recombinationNetwork;

    private List<RecombinationNetworkEvent> recombinationNetworkEventList, storedRecombinationNetworkEventList;

    public boolean eventListDirty = true;

    @Override
    public void initAndValidate() {
        recombinationNetwork = recombinationNetworkInput.get();

        storedRecombinationNetworkEventList = new ArrayList<>();
    }

    List<RecombinationNetworkEvent> getRecombinationNetworkEventList() {
        update();

        return recombinationNetworkEventList;
    }

    public double getBinomialProb() {
        return binomialProbInput.get() != null
                ? binomialProbInput.get().getArrayValue()
                : 0.5;
    }

    void update() {
        if (!eventListDirty)
            return;

        recombinationNetworkEventList = recombinationNetwork.getNodes().stream().map(n -> {
            RecombinationNetworkEvent event = new RecombinationNetworkEvent();
            event.time = n.getHeight();
            event.node = n;
            switch(n.getChildCount()) {
                case 0:
                    event.type = RecombinationNetworkEvent.RecombinationNetworkEventType.SAMPLE;
                    break;

                case 1:
                    event.type = RecombinationNetworkEvent.RecombinationNetworkEventType.RECOMBINATION;
                    break;

                case 2:
                    event.type = RecombinationNetworkEvent.RecombinationNetworkEventType.COALESCENCE;
                    break;

                default:
                    throw new RuntimeException("RecombinationNetwork node has illegal number of children.");
            }
            return event;
        }).sorted(Comparator.comparingDouble(e -> e.time)).collect(Collectors.toList());

        int lineages = 0;
        double totalReassortmentObsProb = 0;

        for (RecombinationNetworkEvent event : recombinationNetworkEventList) {
            switch(event.type) {
                case SAMPLE:
                    lineages += 1;
                    totalReassortmentObsProb += event.node.getParentEdges().get(0).getRecombinationObsProb(getBinomialProb());
                    break;

                case RECOMBINATION:
                    lineages += 1;
                    totalReassortmentObsProb -= event.node.getChildEdges().get(0).getRecombinationObsProb(getBinomialProb());
                    totalReassortmentObsProb += event.node.getParentEdges().get(0).getRecombinationObsProb(getBinomialProb());
                    totalReassortmentObsProb += event.node.getParentEdges().get(1).getRecombinationObsProb(getBinomialProb());

//                    event.segsToSort = event.node.getChildEdges().get(0).hasSegments.cardinality();
//                    event.segsSortedLeft = event.node.getParentEdges().get(0).hasSegments.cardinality();
                    break;

                case COALESCENCE:
                    lineages -= 1;
                    totalReassortmentObsProb -= event.node.getChildEdges().get(0).getRecombinationObsProb(getBinomialProb());
                    totalReassortmentObsProb -= event.node.getChildEdges().get(1).getRecombinationObsProb(getBinomialProb());
                    totalReassortmentObsProb += event.node.getParentEdges().get(0).getRecombinationObsProb(getBinomialProb());
                    break;
            }

            event.lineages = lineages;
            event.totalReassortmentObsProb = totalReassortmentObsProb;
        }

        eventListDirty = false;
    }

    @Override
    protected boolean requiresRecalculation() {
        eventListDirty = true;

        return true;
    }

    @Override
    protected void restore() {
        List<RecombinationNetworkEvent> tmp = recombinationNetworkEventList;
        recombinationNetworkEventList = storedRecombinationNetworkEventList;
        storedRecombinationNetworkEventList = tmp;

        super.restore();
    }

    @Override
    protected void store() {
        storedRecombinationNetworkEventList.clear();
        storedRecombinationNetworkEventList.addAll(recombinationNetworkEventList);

        super.store();
    }
}