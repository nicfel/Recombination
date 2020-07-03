package recombination.network;

import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import coalre.network.parser.NetworkBaseVisitor;
import coalre.network.parser.NetworkLexer;
import coalre.network.parser.NetworkParser;
import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;

import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

public class RecombinationNetwork extends StateNode {

    protected RecombinationNetworkEdge rootEdge;

    protected RecombinationNetworkEdge storedRootEdge;

    protected Integer segmentCount = null;
    
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
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), false, Integer.MAX_VALUE) + ";";
    }
    
    public String getExtendedNewick(int followSegment) {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), false, followSegment) + ";";
    }


    /**
     * @return Extended Newick representation of recombinationNetwork, with
     *         segment presence annotation.
     */
    public String getExtendedNewickVerbose() {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), true, Integer.MAX_VALUE) + ";";
    }
    
    public String getExtendedNewickVerbose(int followSegment) {
        return getExtendedNewick(rootEdge, new ArrayList<>(), new ArrayList<>(), true, followSegment) + ";";
    }


    private String getExtendedNewick(RecombinationNetworkEdge currentEdge, List<RecombinationNetworkNode> seenReassortmentNodes, 
    		List<Boolean> isTraverseEdge, boolean verbose, int followSegment) {
        StringBuilder result = new StringBuilder();

        boolean traverse = true;
        int hybridID = -1;
        boolean printMetaData = true;
        if (currentEdge.childNode.isReassortment()) {
        	
            hybridID = seenReassortmentNodes.indexOf(currentEdge.childNode);
            if (hybridID<0) {
	        	List<RecombinationNetworkEdge> parentEdges = currentEdge.childNode.getParentEdges();
	        	
	        	// get the other edge
	        	RecombinationNetworkEdge otherEdge;
	        	if (parentEdges.get(0).equals(currentEdge))
	        		otherEdge = parentEdges.get(1);
	        	else
	        		otherEdge = parentEdges.get(0);
	        	
	        	// check which edge is the main edge
	        	if (otherEdge.hasSegments.get(followSegment)){
	                traverse = false;
	                seenReassortmentNodes.add(currentEdge.childNode);
	                isTraverseEdge.add(true);
	                hybridID = seenReassortmentNodes.size()-1;
	        	}else if(currentEdge.hasSegments.get(followSegment)){
	                seenReassortmentNodes.add(otherEdge.childNode);
	                isTraverseEdge.add(false);
	                hybridID = seenReassortmentNodes.size()-1;	        		
	        	}else if (currentEdge.hasSegments.cardinality() < otherEdge.hasSegments.cardinality()){
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

                result.append(getExtendedNewick(childEdge, seenReassortmentNodes, isTraverseEdge, verbose, followSegment));
            }

            result.append(")");
        }

        if (currentEdge.childNode.getTaxonLabel() != null)
            result.append(currentEdge.childNode.getTaxonLabel());

        if (hybridID>=0) {
            result.append("#H").append(hybridID);
        }

        result.append("[&");
        result.append("segments=").append(currentEdge.hasSegments);
        if (verbose) {
            for (int segIdx=0; segIdx<getSegmentCount(); segIdx++) {
                result.append(",seg").append(segIdx).append("=")
                        .append(currentEdge.hasSegments.get(segIdx));
            }
        }
        result.append(",segsCarried=").append(currentEdge.hasSegments.cardinality());
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

        CharStream inputStream = CharStreams.fromString(newickStr);
        NetworkLexer lexer = new NetworkLexer(inputStream);
        CommonTokenStream tokenStream = new CommonTokenStream(lexer);
        NetworkParser parser = new NetworkParser(tokenStream);
        ParseTree tree = parser.network();

        RecombinationNetworkBuilderVisitor builder = new RecombinationNetworkBuilderVisitor();
        rootEdge = builder.visit(tree);

        List<RecombinationNetworkNode> leafNodes = new ArrayList<>(getLeafNodes());
    }

    /** StateNode implementation: **/

    @Override
	public String toString() {
        return getExtendedNewick();
    }

    @Override
    public void setEverythingDirty(boolean isDirty) {
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
        otherRecombinationNetwork.segmentCount = null;
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
        segmentCount = null;

        // Restore taxon indices
        getLeafNodes().forEach(n -> n.setTaxonIndex(taxonToIndexMap.get(n.getTaxonLabel())));
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
        storedRootEdge = rootEdge.getCopy();
    }

    @Override
    public void restore() {
        RecombinationNetworkEdge tmp = storedRootEdge;
        storedRootEdge = rootEdge;
        rootEdge = tmp;
        hasStartedEditing = false;
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


//    /**
//     * Visitor class used to build recombinationNetwork from parser-generated AST.
//     */
//    class RecombinationNetworkBuilderVisitor extends RecombinationNetworkBaseVisitor<RecombinationNetworkEdge> {
//
//        Map<Integer, RecombinationNetworkNode> seenHybrids;
//        Map<RecombinationNetworkEdge, Double> edgeLengths;
//
//        private void convertEdgeLengthsToNodeHeights() {
//
//        }
//
//        double getMaxRootToLeafTime(RecombinationNetworkNode node, Set<RecombinationNetworkNode> seenNodes) {
//
//            if (seenNodes.contains(node))
//                return 0.0;
//
//            seenNodes.add(node);
//
//            double maxTime = 0.0;
//            for (RecombinationNetworkEdge childEdge : node.getChildEdges()) {
//                RecombinationNetworkNode childNode = childEdge.childNode;
//                childNode.setHeight(node.getHeight()-edgeLengths.get(childEdge));
//
//                double thisTime = edgeLengths.get(childEdge) +
//                        getMaxRootToLeafTime(childNode, seenNodes);
//                if (thisTime > maxTime)
//                    maxTime = thisTime;
//            }
//
//            return maxTime;
//        }
//
//        void shiftNodeHeights(double maxRTLT, RecombinationNetworkNode node, Set<RecombinationNetworkNode> seenNodes) {
//            if (seenNodes.contains(node))
//                return;
//
//            seenNodes.add(node);
//
//            node.setHeight(node.getHeight() + maxRTLT);
//
//            for (RecombinationNetworkEdge childEdge : node.getChildEdges())
//                shiftNodeHeights(maxRTLT, childEdge.childNode, seenNodes);
//        }
//
//        @Override
//        public RecombinationNetworkEdge visitRecombinationNetwork(RecombinationNetworkParser.RecombinationNetworkContext ctx) {
//            seenHybrids = new HashMap<>();
//            edgeLengths = new HashMap<>();
//
//            RecombinationNetworkEdge rootEdge = visit(ctx.node());
//
//            Set<RecombinationNetworkNode> seenNodes = new HashSet<>();
//            RecombinationNetworkNode rootNode = rootEdge.childNode;
//            rootNode.setHeight(0.0);
//            double maxRTLT = getMaxRootToLeafTime(rootNode, seenNodes);
//
//            seenNodes.clear();
//            shiftNodeHeights(maxRTLT, rootEdge.childNode, seenNodes);
//
//            return rootEdge;
//        }
//
//        private String removeQuotes(String str) {
//
//            String[] quoteChars = {"\"", "'"};
//
//            for (String quoteChar : quoteChars) {
//                if (str.startsWith(quoteChar) && str.endsWith(quoteChar) && str.length() >= 2)
//                    str = str.substring(1, str.length() - 1);
//            }
//
//            return str;
//        }
//
//        @Override
//        public RecombinationNetworkEdge visitNode(RecombinationNetworkParser.NodeContext ctx) {
//
//            visit(ctx.post());
//
//            RecombinationNetworkNode node;
//
//            if (ctx.post().hybrid() != null) {
//                int hybridID = Integer.valueOf(ctx.post().hybrid().id.getText());
//
//                if (seenHybrids.containsKey(hybridID)) {
//                    node = seenHybrids.get(hybridID);
//                } else {
//                    node = new RecombinationNetworkNode();
//                    seenHybrids.put(hybridID, node);
//                }
//            } else {
//                node = new RecombinationNetworkNode();
//            }
//
//            if (ctx.post().label() != null)
//                node.setTaxonLabel(removeQuotes(ctx.post().label().getText()));
//
//            for (RecombinationNetworkParser.NodeContext childNodeCtx : ctx.node()) {
//                RecombinationNetworkEdge childEdge = visit(childNodeCtx);
//                childEdge.parentNode = node;
//                node.addChildEdge(childEdge);
//            }
//
//            boolean segmentsProcessed = false;
//            BitSet hasSegments = new BitSet();
//            if (ctx.post().meta() != null
//                    && ctx.post().meta().attrib() != null) {
//
//                for (RecombinationNetworkParser.AttribContext attribCtx : ctx.post().meta().attrib()) {
//                    if (!removeQuotes(attribCtx.attribKey.getText()).equals("segments"))
//                        continue;
//
//                    if (attribCtx.attribValue().vector() == null)
//                        continue;
//
//                    for (RecombinationNetworkParser.AttribValueContext attribValueCtx : attribCtx.attribValue().vector().attribValue())
//                        hasSegments.set(Integer.valueOf(attribValueCtx.getText()));
//
//                    segmentsProcessed = true;
//                    break;
//                }
//
//            }
//
//            if (!segmentsProcessed) {
//                throw new RuntimeException("Segment attribute missing/malformed " +
//                        "for edge in input recombinationNetwork string.");
//            }
//
//            RecombinationNetworkEdge edge = new RecombinationNetworkEdge(null, node, hasSegments);
//            node.addParentEdge(edge);
//
//            if (ctx.post().length == null) {
//                throw new RuntimeException("Edge missing length in input " +
//                        "recombinationNetwork string.");
//            }
//
//            edgeLengths.put(edge, Double.valueOf(ctx.post().length.getText()));
//
//            return edge;
//        }
//    }
//
    /**
     * Will be used in the next version of BEAST to prevent date trait cloning
     * from breaking the BEAuti model.
     */
    public boolean notCloneable() {
        return true;
    }
}