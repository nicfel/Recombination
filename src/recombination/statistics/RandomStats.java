package recombination.statistics;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
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

public class RandomStats extends BEASTObject implements Loggable {


    public Input<RecombinationNetwork> networkInput = new Input<>("recombinationNetwork",
            "Network for which to log statistics.",
            Input.Validate.REQUIRED);


    RecombinationNetwork network;
    boolean logObservable = false;

    public RandomStats() { }

    @Override
    public void initAndValidate() {
        network = networkInput.get();
    }

    @Override
    public void init(PrintStream out) {

        String prefix = network.getID() == null ? "networkStat." : network.getID() + ".";

        if (logObservable) {
        	throw new IllegalArgumentException("conditioning on only observed nodes not implemented yet");
        }else {
        	out.print(prefix + "mean\t" +
	                prefix + "std\t");
        	out.print(prefix + "0\t" +
	                prefix + "1\t");
        }

    }
    

    @Override
    public void log(long sample, PrintStream out) {
        int[] lala = new int[2];
        if (logObservable){
        	throw new IllegalArgumentException("conditioning on only observed nodes not implemented yet");
        	
        }else{   
            List<RecombinationNetworkNode> possibleSourceEdges = network.getNodes().stream()
            .filter(e -> e.isRecombination())
            .collect(Collectors.toList());
    
		    double mean = 0;
		    double var = 0;
		    List<Double> values = new ArrayList<>();
		    for (RecombinationNetworkNode n : possibleSourceEdges) {
	    		int max = Math.min(n.getParentEdges().get(0).breakPoints.getMax(), n.getParentEdges().get(1).breakPoints.getMax());
				int min = Math.max(n.getParentEdges().get(0).breakPoints.getMin(), n.getParentEdges().get(1).breakPoints.getMin());
				values.add((double)Randomizer.nextInt(min-max)+max+1);
				mean += values.get(values.size()-1);
		    }
		    mean /= possibleSourceEdges.size();
		    
		    for (Double i : values) {
		    	var += Math.pow(i -mean, 2);
		    }
		    
		    var /= possibleSourceEdges.size();


        	out.print(mean + "\t" +
        			var + "\t");
		    for (RecombinationNetworkNode n : possibleSourceEdges) {
		    	if (n.getChildEdges().get(0).breakPoints.contains(0))
		    		lala[0]++;
		    	if (n.getChildEdges().get(0).breakPoints.contains(1))
		    		lala[1]++;

		    }
        	out.print(lala[0] + "\t" +
        			lala[1] + "\t");


        	
        }
    }

    @Override
    public void close(PrintStream out) {

    }
}
