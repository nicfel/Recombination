package recombination.statistics;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Tree;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

public class LocusStatsLogger extends BEASTObject implements Loggable {


    public Input<RecombinationNetwork> networkInput = new Input<>("recombinationNetwork",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);
    
    public Input<Integer> locusInput = new Input<>("locus",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);



    RecombinationNetwork network;
    boolean logObservable = false;
    int locus;

    public LocusStatsLogger() { }

    @Override
    public void initAndValidate() {
        network = networkInput.get();
        locus = locusInput.get();
    }

    @Override
    public void init(PrintStream out) {

        String prefix = network.getID() == null ? "" : network.getID() + ".";
        
//    	out.print(prefix + "locus." + locus + ".height\t" +
//                prefix + "locus." + locus + ".length\t");
    	out.print(           prefix + "locus." + locus + ".length\t");


    }

    @Override
    public void log(long sample, PrintStream out) {
        	
//    	out.print(getLocusHeight(network, locus) + "\t" +
//    			getLocusLength(network, locus) + "\t");
    	out.print(getLocusHeight(network, locus) + "\t");

        
    }

    @Override
    public void close(PrintStream out) {

    }

    public static int getRecombinationCount(RecombinationNetwork network) {
        return (int)network.getNodes().stream().filter(RecombinationNetworkNode::isRecombination).count();
    }

    public static double getLocusLength(RecombinationNetwork network, int locus) {
        return getLocusLength(network.getRootEdge(), locus, false);
    }

    public static double getLocusHeight(RecombinationNetwork network, int locus) {
        return getLocusRoot(network.getRootEdge(), locus);
    }
    
    private static double getLocusRoot(RecombinationNetworkEdge edge, int locus) {
    	if (edge.childNode.isCoalescence()) {
    		if (edge.childNode.getChildEdges().get(0).breakPoints.contains(locus) && 
    				edge.childNode.getChildEdges().get(1).breakPoints.contains(locus))
    			return edge.childNode.getHeight(); 
    	}
    	
    	for (RecombinationNetworkEdge e : edge.childNode.getChildEdges()) {
    		if (e.breakPoints.contains(locus))
    			return getLocusRoot(e, locus);
    	}
    	
    	return -1;
    }
    
    private static double getLocusLength(RecombinationNetworkEdge edge, int locus, boolean started) {
    	double length = 0.0;
    	
    	if (started)
    		length = edge.getLength();
    	
    	if (edge.childNode.isCoalescence()) {
    		if (edge.childNode.getChildEdges().get(0).breakPoints.contains(locus) && 
    				edge.childNode.getChildEdges().get(1).breakPoints.contains(locus))
    			started = true;    		 
    	}
    	
    	for (RecombinationNetworkEdge e : edge.childNode.getChildEdges()) {
    		if (e.breakPoints.contains(locus))
    			length += getLocusLength(e, locus, started);
    	}
    	
    	return length;
    }
    

    
    
    
}
