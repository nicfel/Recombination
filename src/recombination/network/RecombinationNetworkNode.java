package recombination.network;


import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import recombination.network.BreakPoints.Range;
import recombination.util.NodeEdgeID;

public class RecombinationNetworkNode {

    /**
     * Taxon corresponding to this node (usually only used for leaves)
     */
    protected String taxonLabel;

    /**
     * Index of taxon corresponding to this node (usually only used for leaves)
     */
    protected int taxonIndex = -1;
    
    /**
     * Label of type corresponding to this node (used only for structured models)
     */
    protected String typeLabel;
    
    /**
     * Index of type corresponding to this node (used only for structured models)
     */
    protected int typeIndex;

    /**
     * height of this node.
     */
    protected double height = Double.MAX_VALUE;
        
    /**
     * Arbitrarily labeled metadata on this node. Needed for network summary only
     */
    protected String metaDataString;
    
    /**
     * Dummy breakpoints used for the the likelihood calculations
     */
    public List<BreakPoints> dummy;
    public List<BreakPoints> dummy2;
    public BreakPoints dummy3;
    public List<BreakPoints> dummy4;
    public List<int[]> edgePointer;
    public List<Double> prevLength;
    public List<Integer> prevPointer;
    public BreakPoints dirtyBreakPoints;
    public BreakPoints overlap;
    public boolean visited;
    
    

    
    
    public Integer ID;
    


    List<RecombinationNetworkEdge> children = new ArrayList<>();
    List<RecombinationNetworkEdge> parents = new ArrayList<>();

    public RecombinationNetworkNode(NodeEdgeID nodeEdgeIDs) {
    	ID = nodeEdgeIDs.getNewNodeID();
    	dirtyBreakPoints = new BreakPoints();
    }
    
    public RecombinationNetworkNode(int id) {   	
    	ID = id;
    	dirtyBreakPoints = new BreakPoints();
    	
    }


    public double getHeight() {
        return height;
    }

    public RecombinationNetworkNode setHeight(final double height) {
        this.height = height;        
        for (RecombinationNetworkEdge e : getChildEdges())
        	e.makeDirty(Tree.IS_DIRTY);
        return this;
    }   
    
    public RecombinationNetworkNode setHeightFilty(final double height) {
        this.height = height;        
        for (RecombinationNetworkEdge e : getChildEdges())
        	e.makeDirty(Tree.IS_FILTHY);
        return this;
    }   
    
    public RecombinationNetworkNode setFilty() {
        for (RecombinationNetworkEdge e : getChildEdges())
        	e.makeDirty(Tree.IS_FILTHY);
        return this;
    }
    

    


    /**
     * @return parent node, or null if this is root *
     */
    public List<RecombinationNetworkEdge> getParentEdges() {
        return parents;
    }

    public RecombinationNetworkNode addParentEdge(RecombinationNetworkEdge newParentEdge) {
        parents.add(newParentEdge);
        newParentEdge.childNode = this;
        newParentEdge.makeDirty(Tree.IS_DIRTY);
//        makeAllDirty();
        return this;
    }
    
    public RecombinationNetworkNode removeParentEdge(RecombinationNetworkEdge parentEdge) {
        parents.remove(parentEdge);
        parentEdge.childNode = null;
        parentEdge.makeDirty(Tree.IS_DIRTY);

//        makeAllDirty();
        return this;
    }
    
    public List<RecombinationNetworkEdge> getChildEdges() {
        return children;
    }

    public RecombinationNetworkNode addChildEdge(RecombinationNetworkEdge newChildEdge) {
        children.add(newChildEdge);
        newChildEdge.parentNode = this;
        
        newChildEdge.makeDirty(Tree.IS_DIRTY);
        
//        makeAllDirty();
        return this;
    }
    

    public RecombinationNetworkNode removeChildEdge(RecombinationNetworkEdge childEdge) {
        children.remove(childEdge);
        childEdge.parentNode = null;
        childEdge.makeDirty(Tree.IS_DIRTY);
        return this;
    }
    
    public RecombinationNetworkEdge otherParentEdge(RecombinationNetworkEdge edge) {
    	return parents.get(0).ID==edge.ID ? parents.get(1):parents.get(0);
    }
    
    public RecombinationNetworkEdge otherChildEdge(RecombinationNetworkEdge edge) {
    	return children.get(0).ID==edge.ID ? children.get(1):children.get(0);
    }

      

    /**
     * @return true iff current node is a leaf node.
     */
    public boolean isLeaf() {
        return children.size() == 0;
    }

    /**
     * @return true iff current node is a reassortment node.
     */
    public boolean isRecombination() {
        return parents.size()>1;
    }

    /**
     * @return true iff current node is a coalescence node.
     */
    public boolean isCoalescence() {
        return children.size() > 1;
    }

    public int getChildCount() {
        return children.size();
    }

    public int getParentCount() {
        return parents.size();
    }

    /**
     * @return the label of the taxon (if any) corresponding to this node
     */
    public String getTaxonLabel() {
        return taxonLabel;
    }

    /**
     * Sets the label of the taxon corresponding to this node.
     *
     * @param taxonLabel the new taxon label
     */
    public void setTaxonLabel(String taxonLabel) {
        this.taxonLabel = taxonLabel;
    }


    /**
     * @return the index of the taxon (if any) corresponding to this node
     */
    public int getTaxonIndex() {
        return taxonIndex;
    }

    /**
     * Sets the index of the taxon corresponding to this node.
     *
     * @param taxonIndex new taxon index
     */
    public void setTaxonIndex(int taxonIndex) {
        this.taxonIndex = taxonIndex;
    }
    
    /**
     * @return the label of the type (if any) corresponding to this node
     */
    public String getTypeLabel() {
        return typeLabel;
    }
    
    /**
     * Sets the type corresponding to this node.
     *
     * @param stateLabel new state label
     */
    public void setTypeLabel(String typeLabel) {
        this.typeLabel = typeLabel;
    }
    
    /**
     * @return the label of the state (if any) corresponding to this node
     */
    public int getTypeIndex() {
        return typeIndex;
    }
    
    /**
     * Sets the state corresponding to this node.
     *
     * @param stateLabel new state label
     */
    public void setTypeIndex(int typeIndex) {
        this.typeIndex = typeIndex;
    }  
    

    /**
     * Set a meta data value
     * @param pattern
     * @param value
     */
    public void setMetaData(String metaDataString) {
        this.metaDataString = metaDataString;
    }
    
    public String getMetaData() {
        return metaDataString;
    }   

}
