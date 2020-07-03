package recombination.distribution;

import recombination.network.RecombinationNetworkNode;

import java.util.BitSet;

public class RecombinationNetworkEvent {
    public enum RecombinationNetworkEventType { SAMPLE, COALESCENCE, RECOMBINATION }

    public RecombinationNetworkEventType type;
    public double time;

    /**
     * Number of segments on a reassorting lineage.
     */
    int segsToSort;

    /**
     * Number of segments sent to the first parent.
     */
    int segsSortedLeft;

    public int lineages;
    public double totalReassortmentObsProb;

    /**
     * Only used when setting up event list.
     * May not point to a compatible node at other times.
     */
    public RecombinationNetworkNode node;
}
