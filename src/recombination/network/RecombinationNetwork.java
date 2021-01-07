package recombination.network;

import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.network.parser.NetworkBaseVisitor;
import coalre.network.parser.NetworkLexer;
import coalre.network.parser.NetworkParser;
import recombination.alignment.RecombinationAlignment;
import recombination.util.NodeEdgeID;

import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;
import org.fest.swing.util.Range;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class RecombinationNetwork extends StateNode {

    protected RecombinationNetworkEdge rootEdge;

    protected RecombinationNetworkEdge storedRootEdge;
        
    public Integer totalLength = null;
    
    public boolean resume = false;
    
    public NodeEdgeID nodeEdgeIDs;
    
    // keeps track of whether the last operator was the gibbs operator, meaning the likelihood does not requrie updating;
    public boolean wasGibbs;
    
    
    public RecombinationNetwork() {
    }

    public RecombinationNetwork(RecombinationNetworkEdge rootEdge) {
        this.rootEdge = rootEdge;
    }

    public RecombinationNetwork(String newickString) {
        fromExtendedNewick(newickString);
    }

    public RecombinationNetwork(String newickString, TaxonSet taxonSet) {
        fromExtendedNewick(newickString);

        for (RecombinationNetworkNode leafNode : getLeafNodes())
            leafNode.setTaxonIndex(taxonSet.getTaxonIndex(leafNode.getTaxonLabel()));
    }
    

    @Override
    public void initAndValidate() { }

    /**
     * @return the root edge of the network
     */
	public RecombinationNetworkEdge getRootEdge() {
        return rootEdge;
    }

    /**
     * Set the root edge of the network.
     *
     * @param rootEdge the new root edge.
     */
    public void setRootEdge(RecombinationNetworkEdge rootEdge) {
        this.rootEdge = rootEdge;
    }

    /**
     * @return set of node objects comprising network
     */
    public Set<RecombinationNetworkNode> getNodes() {
        Set<RecombinationNetworkNode> recombinationNetworkNodeSet = new HashSet<>();

        getNodesRecurse(rootEdge, recombinationNetworkNodeSet);

        return recombinationNetworkNodeSet;
    }

    private void getNodesRecurse(RecombinationNetworkEdge lineage, Set<RecombinationNetworkNode> recombinationNetworkNodeSet) {

        if (recombinationNetworkNodeSet.contains(lineage.childNode))
            return;

        recombinationNetworkNodeSet.add(lineage.childNode);

        for (RecombinationNetworkEdge childLineage : lineage.childNode.getChildEdges())
            getNodesRecurse(childLineage, recombinationNetworkNodeSet);
    }

    /**
     * @return set of leaf nodes in recombinationNetwork
     */
    public Set<RecombinationNetworkNode> getLeafNodes() {
        return getNodes().stream().filter(RecombinationNetworkNode::isLeaf).collect(Collectors.toSet());
    }

    /**
     * @return set of internal nodes in recombinationNetwork
     */
    public Set<RecombinationNetworkNode> getInternalNodes() {
        return getNodes().stream().filter(n -> !n.isLeaf()).collect(Collectors.toSet());
    }

    /**
     * @return set of edge objects comprising recombinationNetwork
     */
    public Set<RecombinationNetworkEdge> getEdges() {
        Set<RecombinationNetworkEdge> recombinationNetworkEdgeSet = new HashSet<>();

        getEdgesRecurse(rootEdge, recombinationNetworkEdgeSet);

        return recombinationNetworkEdgeSet;
    }


    private void getEdgesRecurse(RecombinationNetworkEdge edge, Set<RecombinationNetworkEdge> recombinationNetworkEdgeSet) {

        if (recombinationNetworkEdgeSet.contains(edge))
            return;

        recombinationNetworkEdgeSet.add(edge);
        for (RecombinationNetworkEdge childEdge : edge.childNode.getChildEdges())
            getEdgesRecurse(childEdge, recombinationNetworkEdgeSet);
    }

    /**
     * @return Extended Newick representation of recombinationNetwork
     */
    public String getExtendedNewick() {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), false, new BreakPoints(totalLength)) + ";";
    }
    
    public String getExtendedNewick(BreakPoints followBreaks) {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), false, followBreaks) + ";";
    }


    /**
     * @return Extended Newick representation of recombinationNetwork, with
     *         segment presence annotation.
     */
    public String getExtendedNewickVerbose() {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), true, new BreakPoints(totalLength)) + ";";
    }
    
    public String getExtendedNewickVerbose(BreakPoints followBreaks) {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), true, followBreaks) + ";";
    }


    private String getExtendedNewick(RecombinationNetworkEdge currentEdge, List<RecombinationNetworkNode> seenReassortmentNodes, 
    		List<Boolean> isTraverseEdge, boolean verbose, BreakPoints followBreaks) {
        StringBuilder result = new StringBuilder();

        boolean traverse = true;
        int hybridID = -1;
        boolean printMetaData = true;
        if (currentEdge.childNode.isRecombination()) {
        	
            hybridID = seenReassortmentNodes.indexOf(currentEdge.childNode);
            if (hybridID<0) {
	        	List<RecombinationNetworkEdge> parentEdges = currentEdge.childNode.getParentEdges();
	        	
	        	// get the other edge
	        	RecombinationNetworkEdge otherEdge;
	        	if (parentEdges.get(0).equals(currentEdge))
	        		otherEdge = parentEdges.get(1);
	        	else
	        		otherEdge = parentEdges.get(0);
	        	
	        	BreakPoints bp1 = currentEdge.breakPoints.copy();
	        	BreakPoints bp2 = otherEdge.breakPoints.copy();
	        	bp1.and(followBreaks);
	        	bp2.and(followBreaks);
	        	
	        	// check which edge is the main edge
	        	if (bp1.getGeneticLength() < bp2.getGeneticLength()){
	                traverse = false;
	                seenReassortmentNodes.add(currentEdge.childNode);
	                isTraverseEdge.add(true);
	                hybridID = seenReassortmentNodes.size()-1;
	        	}else{
	                seenReassortmentNodes.add(otherEdge.childNode);
	                isTraverseEdge.add(false);
	                hybridID = seenReassortmentNodes.size()-1;	        		
	        	}
            }else{
            	traverse = isTraverseEdge.get(hybridID);
            }        
        }

        if (traverse && !currentEdge.childNode.isLeaf()) {
            result.append("(");

            boolean isFirst = true;
            for (RecombinationNetworkEdge childEdge : currentEdge.childNode.getChildEdges()) {
                if (isFirst)
                    isFirst = false;
                else
                    result.append(",");

                result.append(getExtendedNewick(childEdge, seenReassortmentNodes, isTraverseEdge, verbose, followBreaks));
            }

            result.append(")");
        }

        if (currentEdge.childNode.getTaxonLabel() != null)
            result.append(currentEdge.childNode.getTaxonLabel());

        if (hybridID>=0) {
            result.append("#H").append(hybridID);
        }

        result.append("[&");
        
        result.append("loci={").append(currentEdge.breakPoints);
        result.append("}");
        if (currentEdge.passingRange!=null && !currentEdge.passingRange.isEmpty()) {
	        result.append(",pr={").append(currentEdge.passingRange);
	        result.append("}");
        }
//        result.append(",id=\"").append(currentEdge.ID).append("\"");


        
//        for (int i =0; i < totalLength; i++)
//            result.append(",").append(i).append("=").append(currentEdge.breakPoints.contains(i));

        	
//        result.append(",dirty=").append(currentEdge.isDirty);
//        result.append(",dirtyEdges={").append(currentEdge.childNode.dirtyBreakPoints);
//        result.append("}");
//        if (currentEdge.childNode.dummy4!=null) {
//	        if (currentEdge.childNode.dummy4.size()>0) {
//		        result.append(",dummy=\"");
//		        
//		        for (int i = 0; i < currentEdge.childNode.dummy4.size()-1; i++)
//		        	result.append(currentEdge.childNode.dummy4.get(i)).append(":");
//		        
//	        	result.append(currentEdge.childNode.dummy4.get(currentEdge.childNode.dummy4.size()-1));
//		        
//		        result.append("\"");
//	        }
//        }


        

        result.append(",length=").append(currentEdge.breakPoints.getGeneticLengthInt());
        
        if (currentEdge.childNode.getTypeLabel() != null) 
        		result.append(",state=").append(currentEdge.childNode.getTypeLabel());
        
        if (currentEdge.childNode.metaDataString != null) 
			result.append(currentEdge.childNode.getMetaData());

        result.append("]");

        if (currentEdge.parentNode != null)
            result.append(":").append(currentEdge.parentNode.getHeight() - currentEdge.childNode.getHeight());
        else
            result.append(":0.0");

        return result.toString();
    }

    public void fromExtendedNewick(String newickStr) {
    	nodeEdgeIDs = new NodeEdgeID();

        CharStream inputStream = CharStreams.fromString(newickStr);
        NetworkLexer lexer = new NetworkLexer(inputStream);
        CommonTokenStream tokenStream = new CommonTokenStream(lexer);
        NetworkParser parser = new NetworkParser(tokenStream);
        ParseTree tree = parser.network();

        RecombinationNetworkBuilderVisitor builder = new RecombinationNetworkBuilderVisitor();
        rootEdge = builder.visit(tree);

        List<RecombinationNetworkNode> nodes = new ArrayList<>(getNodes());
        totalLength = -1;
        for (RecombinationNetworkNode n : nodes) {
        	if (n.isLeaf()) {
        		totalLength = Math.max(totalLength, n.getParentEdges().get(0).breakPoints.getMax()+1);
        	}
        } 
        for (RecombinationNetworkNode n : nodes) {
        	if (n.isRecombination()) {
        		// set passing range
        		int min1 = n.getParentEdges().get(0).breakPoints.getMin();
        		int max1 = n.getParentEdges().get(0).breakPoints.getMax();
        		int min2 = n.getParentEdges().get(1).breakPoints.getMin();
        		int max2 = n.getParentEdges().get(1).breakPoints.getMax();
        		
        		if (min1>max2) {
        			int newBp = Randomizer.nextInt(min1-max2)+max2;
        			n.getParentEdges().get(0).setPassingRange(newBp+1, totalLength-1);
        			n.getParentEdges().get(1).setPassingRange(0, newBp);
        		}else if (min2>max1){
        			int newBp = Randomizer.nextInt(min2-max1)+max1;
        			n.getParentEdges().get(1).setPassingRange(newBp+1, totalLength-1);
        			n.getParentEdges().get(0).setPassingRange(0, newBp);
        		}else {
        			int newBp = Randomizer.nextInt(totalLength);
        			n.getParentEdges().get(1).setPassingRange(newBp+1, totalLength-1);
        			n.getParentEdges().get(0).setPassingRange(0, newBp);
        		}
        	}
        }  
    }

    /** StateNode implementation: **/

    @Override
	public String toString() {
        return getExtendedNewick();
    }

    @Override
    public void setEverythingDirty(boolean isDirty) {
    	for (RecombinationNetworkEdge e : getEdges().stream().collect(Collectors.toList())) {
            if (!isDirty) {
            	e.isDirty = Tree.IS_CLEAN;
            }else {
            	e.isDirty = Tree.IS_FILTHY;
            }
    	}
    	if (!isDirty) {
        	for (RecombinationNetworkNode n : getNodes().stream().collect(Collectors.toList())) {
            	n.dirtyBreakPoints = new BreakPoints();
        	}
    	}
        setSomethingIsDirty(isDirty);
    }

    @Override
    public StateNode copy() {
        return new RecombinationNetwork(rootEdge.getCopy());
    }

    @Override
    public void assignTo(StateNode other) {
        RecombinationNetwork otherRecombinationNetwork = (RecombinationNetwork) other;

        otherRecombinationNetwork.rootEdge = rootEdge;
        otherRecombinationNetwork.storedRootEdge = null;
        otherRecombinationNetwork.totalLength = totalLength;
        nodeEdgeIDs = new NodeEdgeID();
    }

    @Override
    public void assignFrom(StateNode other) {
        assignFromFragile(other);
        setID(other.getID());
    }

    @Override
    public void assignFromFragile(StateNode other) {
        RecombinationNetwork otherRecombinationNetwork = (RecombinationNetwork) other;

        // Save taxon indices.
        Map<String, Integer> taxonToIndexMap = new HashMap<>();
        getLeafNodes().forEach(n -> taxonToIndexMap.put(n.getTaxonLabel(), n.getTaxonIndex()));

        rootEdge = otherRecombinationNetwork.rootEdge;
        storedRootEdge = null;
        totalLength = otherRecombinationNetwork.totalLength;

        // Restore taxon indices
        getLeafNodes().forEach(n -> n.setTaxonIndex(taxonToIndexMap.get(n.getTaxonLabel())));
        resume = true;
    }

    @Override
    public void fromXML(org.w3c.dom.Node node) {
        fromExtendedNewick(node.getTextContent().replaceAll("&amp;", "&"));
    }

    @Override
    public int scale(double scale) {
        return 0;
    }

    @Override
    protected void store() {
//    	System.out.println("store");
        storedRootEdge = rootEdge.getCopy();
        
        List<Integer> nodeIDs = new ArrayList<>();
        List<Integer> edgeIDs = new ArrayList<>();
        for (RecombinationNetworkNode node : getNodes().stream().filter(e -> !e.isLeaf()).collect(Collectors.toList())) 
        	nodeIDs.add(node.ID);
        for (RecombinationNetworkEdge edge : getEdges()) 
        	edgeIDs.add(edge.ID);
        
        nodeEdgeIDs.purgeNodeIDs(nodeIDs);
        nodeEdgeIDs.purgeEdgeIDs(edgeIDs);
        nodeEdgeIDs.store();
    }
    
   
    @Override
    public void restore() {
//    	System.out.println("restore");

        RecombinationNetworkEdge tmp = storedRootEdge;
        storedRootEdge = rootEdge;
        rootEdge = tmp;
        hasStartedEditing = false;
        
    	for (RecombinationNetworkEdge e : getEdges().stream().collect(Collectors.toList()))
           	e.isDirty = Tree.IS_CLEAN;
    	
    	// restore id mapping
    	nodeEdgeIDs.restore();
    }

    @Override
    public int getDimension() {
        return 0;
    }

    @Override
    public double getArrayValue() {
        return 0;
    }

    @Override
    public double getArrayValue(int dim) {
        return 0;
    }

    /** Loggable implementation: **/

    @Override
    public void init(PrintStream out) {
        out.println("#nexus");
        out.println("begin trees;");
    }

    @Override
    public void close(PrintStream out) {
        out.println("end trees;");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.println("tree STATE_" + sample + " = " + getExtendedNewick());
    }


    /**
     * Visitor class used to build recombinationNetwork from parser-generated AST.
     */
    class RecombinationNetworkBuilderVisitor extends NetworkBaseVisitor<RecombinationNetworkEdge> {

        Map<Integer, RecombinationNetworkNode> seenHybrids;
        Map<RecombinationNetworkEdge, Double> edgeLengths;

        private void convertEdgeLengthsToNodeHeights() {

        }

        double getMaxRootToLeafTime(RecombinationNetworkNode node, Set<RecombinationNetworkNode> seenNodes) {

            if (seenNodes.contains(node))
                return 0.0;

            seenNodes.add(node);

            double maxTime = 0.0;
            for (RecombinationNetworkEdge childEdge : node.getChildEdges()) {
                RecombinationNetworkNode childNode = childEdge.childNode;
                childNode.setHeight(node.getHeight()-edgeLengths.get(childEdge));

                double thisTime = edgeLengths.get(childEdge) +
                        getMaxRootToLeafTime(childNode, seenNodes);
                if (thisTime > maxTime)
                    maxTime = thisTime;
            }

            return maxTime;
        }

        void shiftNodeHeights(double maxRTLT, RecombinationNetworkNode node, Set<RecombinationNetworkNode> seenNodes) {
            if (seenNodes.contains(node))
                return;

            seenNodes.add(node);

            node.setHeight(node.getHeight() + maxRTLT);

            for (RecombinationNetworkEdge childEdge : node.getChildEdges())
                shiftNodeHeights(maxRTLT, childEdge.childNode, seenNodes);
        }

        @Override
        public RecombinationNetworkEdge visitNetwork(NetworkParser.NetworkContext ctx) {
            seenHybrids = new HashMap<>();
            edgeLengths = new HashMap<>();

            RecombinationNetworkEdge rootEdge = visit(ctx.node());

            Set<RecombinationNetworkNode> seenNodes = new HashSet<>();
            RecombinationNetworkNode rootNode = rootEdge.childNode;
            rootNode.setHeight(0.0);
            double maxRTLT = getMaxRootToLeafTime(rootNode, seenNodes);

            seenNodes.clear();
            shiftNodeHeights(maxRTLT, rootEdge.childNode, seenNodes);

            return rootEdge;
        }

        private String removeQuotes(String str) {

            String[] quoteChars = {"\"", "'"};

            for (String quoteChar : quoteChars) {
                if (str.startsWith(quoteChar) && str.endsWith(quoteChar) && str.length() >= 2)
                    str = str.substring(1, str.length() - 1);
            }

            return str;
        }

        @Override
        public RecombinationNetworkEdge visitNode(NetworkParser.NodeContext ctx) {

            visit(ctx.post());

            RecombinationNetworkNode node;

            if (ctx.post().hybrid() != null) {
                int hybridID = Integer.valueOf(ctx.post().hybrid().id.getText());

                if (seenHybrids.containsKey(hybridID)) {
                    node = seenHybrids.get(hybridID);
                } else {
                    node = new RecombinationNetworkNode(nodeEdgeIDs);
                    seenHybrids.put(hybridID, node);
                }
            } else {
                node = new RecombinationNetworkNode(nodeEdgeIDs);
            }

            if (ctx.post().label() != null)
                node.setTaxonLabel(removeQuotes(ctx.post().label().getText()));

            for (NetworkParser.NodeContext childNodeCtx : ctx.node()) {
                RecombinationNetworkEdge childEdge = visit(childNodeCtx);
                childEdge.parentNode = node;
                node.addChildEdge(childEdge);
            }

            boolean lociProcessed = false;
           
            List<Integer> breakPoints = new ArrayList<Integer>();
            List<Integer> splitPoints = new ArrayList<Integer>();
            if (ctx.post().meta() != null
                    && ctx.post().meta().attrib() != null) {

                for (NetworkParser.AttribContext attribCtx : ctx.post().meta().attrib()) {
                    if (!removeQuotes(attribCtx.attribKey.getText()).equals("loci"))
                        continue;

                    if (attribCtx.attribValue().vector() == null)
                        continue;

                    for (NetworkParser.AttribValueContext attribValueCtx : attribCtx.attribValue().vector().attribValue()) {
                    	if (attribValueCtx.getText().contains("null")){
                    	}else {
	                    	String[] attr = attribValueCtx.getText().split("-");
	                    	breakPoints.add(Integer.valueOf(attr[0]));
	                    	breakPoints.add(Integer.valueOf(attr[1]));
                    	}
                    }

                    lociProcessed = true;
                    break;
                }

            }

            if (!lociProcessed) {
                throw new RuntimeException("Loci attribute missing/malformed " +
                        "for edge in input recombinationNetwork string.");
            }
            
            BreakPoints bp = new BreakPoints();
            bp.init(breakPoints);
            RecombinationNetworkEdge edge = new RecombinationNetworkEdge(null, node, bp, null, nodeEdgeIDs);
            node.addParentEdge(edge);

            if (ctx.post().length == null) {
                throw new RuntimeException("Edge missing length in input " +
                        "recombinationNetwork string.");
            }

            edgeLengths.put(edge, Double.valueOf(ctx.post().length.getText()));

            return edge;
        }
    }

    /**
     * Will be used in the next version of BEAST to prevent date trait cloning
     * from breaking the BEAuti model.
     */
    public boolean notCloneable() {
        return true;
    }
    
    /**
     * Update segment tree below networkNode to match tree implied by network.
     *
     * @param networkNode node below which segment tree is updated
     * @param segmentIdx  index of segment tree
     * @param cladeNodes  map containing links from existing tree clades to tree nodes
     * @param nodeBin     list containing node objects available for new tree clades
     * @return clade corresponding to networkNode
     */
    public Node getLocusChildren(RecombinationNetworkNode networkNode, Integer locus) {

        List<Node> children = new ArrayList<>();
        
        for (RecombinationNetworkEdge childEdge : networkNode.getChildEdges()) {
            if (childEdge.breakPoints.contains(locus)) {
            	children.add(getLocusChildren(childEdge.childNode, locus));
            }
        }
        if (children.size()==2) {
        	Node treeNode = new Node();
        	treeNode.setHeight(networkNode.getHeight());
        	treeNode.addChild(children.get(0));
        	treeNode.addChild(children.get(1));     
        	return treeNode;
        }else if (children.size()==1) {
        	return children.get(0);
        }else {
        	Node treeNode = new Node();
        	treeNode.setHeight(networkNode.getHeight());
        	treeNode.setID(networkNode.getTaxonLabel());
        	return treeNode;
        }       
    }    

    @Override
    public boolean somethingIsDirty() {    	
        return this.hasStartedEditing;
    }

    
}