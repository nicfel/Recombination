package recombination.util;

import java.util.ArrayList;
import java.util.List;

import beast.util.Randomizer;

public class NodeEdgeID {
	private List<Integer> edgeids = new ArrayList<>();
	private List<Integer> nodeids = new ArrayList<>();

	private List<Integer> storedEdgeids = new ArrayList<>();
	private List<Integer> storedNodeids = new ArrayList<>();
		
	public int getNewNodeID(){
		int newID = Randomizer.nextInt();
		while (nodeids.contains(newID))
			newID = Randomizer.nextInt();
		
		nodeids.add(newID);
		return newID;
	}
	
	public int getNewEdgeID(){
		int newID = Randomizer.nextInt();
		while (edgeids.contains(newID))
			newID = Randomizer.nextInt();
		
		edgeids.add(newID);
		return newID;
	}
	
	public void store() {
		storedEdgeids = new ArrayList<>(edgeids);
		storedNodeids = new ArrayList<>(nodeids);		
	}
	
	public void restore() {
		List<Integer> tmp1 = storedEdgeids;
		storedEdgeids = edgeids;
		edgeids = tmp1;
				
				
		List<Integer> tmp2 = storedNodeids;
		storedNodeids = nodeids;
		nodeids = tmp2;
	}
	
	public void purgeNodeIDs(List<Integer> id){
		for (int i = nodeids.size()-1; i!=0; i--) {
			if (!id.contains(nodeids.get(i)))
				nodeids.remove(i);
		}
	}
	
	public void purgeEdgeIDs(List<Integer> id){
		for (int i = edgeids.size()-1; i!=0; i--) {
			if (!id.contains(edgeids.get(i)))
				edgeids.remove(i);
		}
	}


	

}
