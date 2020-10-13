package recombination.statistics;

import java.util.ArrayList;
import java.util.List;

import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;

public class DotConverter {
	
	public static List<String> getDotFormat(RecombinationNetwork network) {
		List<String> dotString = new ArrayList<>();
		List<Integer> nodeNames = new ArrayList<>();
		
		dotString.add("digraph network {\n");
		for (RecombinationNetworkEdge edge : network.getEdges()) {
			if (!edge.isRootEdge()) {			
			
				int indexP = nodeNames.indexOf(edge.parentNode.ID);			
				if (indexP==-1) {
					nodeNames.add(edge.parentNode.ID);
					indexP = nodeNames.size()-1;
				}
				
				if (!edge.childNode.isLeaf()) {
					int indexC = nodeNames.indexOf(edge.childNode.ID);
					
					
					if (indexC==-1) {
						nodeNames.add(edge.childNode.ID);
						indexC = nodeNames.size()-1;
					}	
					
					dotString.add("\t" + indexP + " -> " + indexC + "\n");
				}else {
					dotString.add("\t" + indexP + " -> " + edge.childNode.getTaxonLabel() + "\n");				
				}
			}
			
		}
		dotString.add("}\n");

		
		return dotString;
	}

}
