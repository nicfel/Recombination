package recombination.statistics;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Tree;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

public class RecombinationNetworkStatsLogger extends BEASTObject implements Loggable {


    public Input<RecombinationNetwork> networkInput = new Input<>("recombinationNetwork",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);
    
    public Input<Boolean> logObservableInput = new Input<>("logObservable",
            "if true, it only logs events that happened more recent than the loci mrca.",
            false);



    RecombinationNetwork network;
    boolean logObservable = false;

    public RecombinationNetworkStatsLogger() { }

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        logObservable = logObservableInput.get(); 
    }

    @Override
    public void init(PrintStream out) {

        String prefix = network.getID() == null ? "networkStat." : network.getID() + ".";

        if (logObservable)
	        out.print(prefix + "obsHeight\t" +
	                prefix + "obsTotalLength\t" +
	                prefix + "obsReassortmentNodeCount\t");
        else
        	out.print(prefix + "height\t" +
	                prefix + "totalLength\t" +
	                prefix + "recombinationNodeCount\t");

    }

    @Override
    public void log(long sample, PrintStream out) {
        if (logObservable){
    		double rootHeight = getMaxLociMRCA(network);
			
        	out.print(rootHeight + "\t" +
	                getTotalEdgeLength(network, rootHeight) + "\t" +
	                getRecombinationCount(network, rootHeight) + "\t");

        	
        }else{   
        	out.print(getTotalHeight(network) + "\t" +
	                getTotalEdgeLength(network) + "\t" +
	                getRecombinationCount(network) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {

    }

    public static int getRecombinationCount(RecombinationNetwork network) {
        return (int)network.getNodes().stream().filter(RecombinationNetworkNode::isRecombination).count();
    }

    public static double getTotalEdgeLength(RecombinationNetwork network) {
        return network.getEdges().stream().filter(e -> !e.isRootEdge()).
                map(RecombinationNetworkEdge::getLength).reduce((l1, l2) -> l1+l2).get();
    }

    public static double getTotalHeight(RecombinationNetwork network) {
        return network.getRootEdge().childNode.getHeight();
    }
    
    public static double getMaxLociMRCA(RecombinationNetwork network){
        for (RecombinationNetworkEdge e : network.getEdges()) 
    		e.visited = false;     
        
        return getMaxLociMRCATraverse(network.getRootEdge());
    }

    
    public static double getMaxLociMRCATraverse(RecombinationNetworkEdge edge){
    	if (edge.visited)
    		return -1;
    	
    	edge.visited=true;
    	
    	RecombinationNetworkNode node = edge.childNode;
    	if (node.isCoalescence()) {
    		BreakPoints bp1 = node.getChildEdges().get(0).breakPoints.copy();
    		bp1.and(node.getChildEdges().get(1).breakPoints);
    		if (!bp1.isEmpty()) {
    			return node.getHeight();
    		}
    		return Math.max(getMaxLociMRCATraverse(node.getChildEdges().get(0)), getMaxLociMRCATraverse(node.getChildEdges().get(1)));
    	}else if (node.isRecombination()) {
    		return getMaxLociMRCATraverse(node.getChildEdges().get(0));
    	}else {
    		return -1.0;
    	}    	
    }
    
    public static int getRecombinationCount(RecombinationNetwork network, double rootHeight) {
        return (int)network.getNodes().stream()
        		.filter(e -> e.isRecombination())
				.filter(e -> e.getHeight()<rootHeight)
        		.count();
    }

    public static double getTotalEdgeLength(RecombinationNetwork network, double rootHeight) {
        return network.getEdges().stream()
        		.filter(e -> !e.isRootEdge())
        		.filter(e -> e.parentNode.getHeight()<rootHeight).
                map(RecombinationNetworkEdge::getLength).reduce((l1, l2) -> l1+l2).get();
    }
}
