package recombination.annotator;

import recombination.network.RecombinationNetwork;

/**
 * @author Tim Vaughan <tgvaughan@gmail.com>
 */
public interface NetworkLogReader extends Iterable<RecombinationNetwork> {

    int getNetworkCount();
    int getCorrectedNetworkCount();
}
