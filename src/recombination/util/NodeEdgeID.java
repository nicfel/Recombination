package recombination.util;

import java.util.ArrayList;
import java.util.List;

import beast.util.Randomizer;

public class NodeEdgeID {
	private static List<Integer> edgeids = new ArrayList<>();
	private static List<Integer> nodeids = new ArrayList<>();

	private static List<Integer> storedEdgeids = new ArrayList<>();
	private static List<Integer> storedNodeids = new ArrayList<>();
		
	public static int getNewNodeID(){
		int newID = Randomizer.nextInt();
		while (nodeids.contains(newID))
			newID = Randomizer.nextInt();
		
		nodeids.add(newID);
		return newID;
	}
	
	public static int getNewEdgeID(){
		int newID = Randomizer.nextInt();
		while (edgeids.contains(newID))
			newID = Randomizer.nextInt();
		
		edgeids.add(newID);
		return newID;
	}

	
//	public static void removeNodeID(int id){
//		nodeids.remove(nodeids.indexOf(id));
//	}
//	
//	public static void removeEdgeID(int id){
//		edgeids.remove(edgeids.indexOf(id));
//	}
	
	public static void store() {
		storedEdgeids = new ArrayList<>(edgeids);
		storedNodeids = new ArrayList<>(nodeids);		
	}
	
	public static void restore() {
		List<Integer> tmp1 = storedEdgeids;
		storedEdgeids = edgeids;
		edgeids = tmp1;
				
				
		List<Integer> tmp2 = storedNodeids;
		storedNodeids = nodeids;
		nodeids = tmp2;
	}
	
	public static void purgeNodeIDs(List<Integer> id){
		for (int i = nodeids.size()-1; i!=0; i--) {
			if (!id.contains(nodeids.get(i)))
				nodeids.remove(i);
		}
	}
	
	public static void purgeEdgeIDs(List<Integer> id){
		for (int i = edgeids.size()-1; i!=0; i--) {
			if (!id.contains(edgeids.get(i)))
				edgeids.remove(i);
		}
	}


	

}
