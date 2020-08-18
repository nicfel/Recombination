package recombination.operators;

import beast.core.Description;
import beast.core.Input;
import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;


/**
 * Implements the subnet slide move. General workflow:
 * 1. Choose an edge to move and a child it will carry
 * 2. Make a copy of subnet starting with this child edge
 * 3. Attach a new coppy to the new parent with randomly drawn height
 * 4. Rearrange segments
 * 5. Delete subnet starting at the child in the old position
 */
@Description("Moves the height of an internal node along the branch. " +
        "If it moves up, it can exceed the root and become a new root. " +
        "If it moves down, it may need to make a choice which branch to " +
        "slide down into.")
public class SubRecombinationNetworkSlide extends DivertLociOperator {

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta", true);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by network-height/log2(n-taxa).", -1.0);
    // shadows size
    double size;
    private double limit;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

        size = sizeInput.get();
        limit = limitInput.get();
	}

	@Override
	public double networkProposal() {

		double logHR = 0.0;
		network.startEditing(this);
		// 1. Choose a random edge, avoiding root
		List<RecombinationNetworkEdge> networkEdges = new ArrayList<>(network.getEdges());
		
		final List<RecombinationNetworkEdge> possibleEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> !e.parentNode.isRecombination())
				.collect(Collectors.toList());
		
		final int possibleNodes = possibleEdges.size();
        if (possibleNodes < 1) {
            return Double.NEGATIVE_INFINITY;
        }
        logHR -= Math.log(1.0/possibleNodes);        

		final RecombinationNetworkEdge iEdge = possibleEdges.get(Randomizer.nextInt(possibleNodes));
		final RecombinationNetworkNode i = iEdge.childNode;
		final RecombinationNetworkNode ip = iEdge.parentNode;
		final RecombinationNetworkEdge ipEdge = ip.getParentEdges().get(0); //ip not a reassortment, has only one parent
		final RecombinationNetworkEdge jEdge = getSisterEdge(iEdge);
		final RecombinationNetworkNode j = jEdge.childNode;
		final RecombinationNetworkNode pip = ipEdge.parentNode; 
		RecombinationNetworkEdge newChildEdge = null;
		
		
        // 2. choose a delta to move
        final double delta = getDelta();
        final double oldHeight = ip.getHeight();
        final double newHeight = oldHeight + delta;
        

        
        // 3. if the move is up
        if (delta > 0) {

        	// 3.1 if the topology will change
            if (pip != null &&  pip.getHeight() < newHeight) {

            	RecombinationNetworkNode newParent = pip;
                newChildEdge = ipEdge;
                while (newParent.getHeight() < newHeight) {
                	
                    if (newParent.isRecombination()) {
                    	if (Randomizer.nextBoolean()) {
                    		newChildEdge = newParent.getParentEdges().get(0);
                    		newParent = newChildEdge.parentNode;
                    	}
                    	else {
                    		newChildEdge = newParent.getParentEdges().get(1);
                    		newParent = newChildEdge.parentNode;
                    	}
                    	logHR -= Math.log(0.5);
                    } else {
                        // if not reassortment, has only one parent
                        newChildEdge = newParent.getParentEdges().get(0);
                        newParent = newChildEdge.parentNode;
                    }

                    if (newParent == null) break;
                }
                // the moved node 'p' would become a child of 'newParent'

                // 3.1.1 if creating a new root
                if (newChildEdge.isRootEdge()) {

                	RecombinationNetworkNode destNode = new RecombinationNetworkNode();
                	destNode.setHeight(newHeight);
                	
                	RecombinationNetworkEdge newIpEdge = new RecombinationNetworkEdge();
                	newIpEdge.breakPoints = newChildEdge.breakPoints.copy();
                	destNode.addParentEdge(newIpEdge);
                	destNode.addChildEdge(newChildEdge);
                	
                	RecombinationNetworkEdge iEdgeCopy = new RecombinationNetworkEdge();
                	iEdgeCopy.breakPoints = iEdge.breakPoints.copy();
                	i.addParentEdge(iEdgeCopy);
                	destNode.addChildEdge(iEdgeCopy);
                	
                	
                	BreakPoints segsToAdd = iEdgeCopy.breakPoints.copy();
                	segsToAdd.andNot(getSisterEdge(iEdgeCopy).breakPoints);
                	logHR -= addLociToAncestors(newIpEdge, segsToAdd);
                	
                	BreakPoints segsToRemove = iEdge.breakPoints.copy();
                	segsToRemove.andNot(getSisterEdge(iEdge).breakPoints);
                	logHR += removeLociFromAncestors(ipEdge, segsToRemove);
                	
                	
                	pip.removeChildEdge(ipEdge);
                	ip.removeChildEdge(jEdge);
                	pip.addChildEdge(jEdge);
                	
                	
                	i.removeParentEdge(iEdge);
                	ip.removeChildEdge(iEdge);
                	ip.removeParentEdge(ipEdge);
                	
                	
                	newIpEdge.parentNode = null;
                	network.setRootEdge(newIpEdge);
                	

                }
                // 3.1.2 no new root
                else {
                	RecombinationNetworkNode destNode = new RecombinationNetworkNode();
                	destNode.setHeight(newHeight);
                	
                	RecombinationNetworkEdge newIpEdge = new RecombinationNetworkEdge();
                	newIpEdge.breakPoints = newChildEdge.breakPoints.copy();
                	newParent.removeChildEdge(newChildEdge);
                	newParent.addChildEdge(newIpEdge);
                	destNode.addParentEdge(newIpEdge);
                	destNode.addChildEdge(newChildEdge);
                	
                	RecombinationNetworkEdge iEdgeCopy = new RecombinationNetworkEdge();
                	iEdgeCopy.breakPoints = iEdge.breakPoints.copy();
                	i.addParentEdge(iEdgeCopy);
                	destNode.addChildEdge(iEdgeCopy);
                	
          	
                	BreakPoints segsToAdd = iEdgeCopy.breakPoints.copy();
                	segsToAdd.andNot(getSisterEdge(iEdgeCopy).breakPoints);
                	logHR -= addLociToAncestors(newIpEdge, segsToAdd);

                	
                	BreakPoints segsToRemove = iEdge.breakPoints.copy();
                	segsToRemove.andNot(getSisterEdge(iEdge).breakPoints);
                	logHR += removeLociFromAncestors(ipEdge, segsToRemove);
                	
                	pip.removeChildEdge(ipEdge);
                	ip.removeChildEdge(jEdge);
                	pip.addChildEdge(jEdge);
                	
                	
                	i.removeParentEdge(iEdge);
                	ip.removeChildEdge(iEdge);
                	ip.removeParentEdge(ipEdge);
                }
                
                // 3.1.3 count the hypothetical sources of this destination.
                final int possibleSources = intersectingEdges(newChildEdge, oldHeight, null);

                logHR += Math.log(1.0/possibleSources);

            } else {
                ip.setHeight(newHeight);
            }
        }
        
        // 4 if we are sliding the subnet down.
        else {
        	if (i.getHeight() > newHeight) {
        		return Double.NEGATIVE_INFINITY;
        	}
  	
        	
            // 4.1 will the move change the topology
            if (j.getHeight() > newHeight) {	

            	final List<RecombinationNetworkEdge> possibleChildEdges = new ArrayList<>();
            	final int possibleDestinations = intersectingEdges(jEdge, newHeight, possibleChildEdges);
            	
                // if no valid destinations then return a failure
                if (possibleChildEdges.size() == 0) {
                    return Double.NEGATIVE_INFINITY;
                }
                
                logHR -= Math.log(1.0/possibleDestinations); 
                
                // pick a random parent/child destination edge uniformly from options
                final int childEdgeIndex = Randomizer.nextInt(possibleChildEdges.size());
                newChildEdge = possibleChildEdges.get(childEdgeIndex);
                RecombinationNetworkNode newParent = newChildEdge.parentNode;
            	
                
                // 4.1.1 if p was root
                if (ipEdge.isRootEdge()) {
                	RecombinationNetworkNode destNode = new RecombinationNetworkNode();
                	destNode.setHeight(newHeight);
                	
                	RecombinationNetworkEdge newIpEdge = new RecombinationNetworkEdge();
                	newIpEdge.breakPoints = newChildEdge.breakPoints.copy();
                	newParent.removeChildEdge(newChildEdge);
                	newParent.addChildEdge(newIpEdge);
                	destNode.addParentEdge(newIpEdge);
                	destNode.addChildEdge(newChildEdge);
                	
                	
                	RecombinationNetworkEdge iEdgeCopy = new RecombinationNetworkEdge();
                	iEdgeCopy.breakPoints = iEdge.breakPoints.copy();
                	i.addParentEdge(iEdgeCopy);
                	destNode.addChildEdge(iEdgeCopy);
                	
                	BreakPoints segsToAdd = iEdgeCopy.breakPoints.copy();
                	segsToAdd.andNot(getSisterEdge(iEdgeCopy).breakPoints);
                	logHR -= addLociToAncestors(newIpEdge, segsToAdd);
                	
                	BreakPoints segsToRemove = iEdge.breakPoints.copy();
                	segsToRemove.andNot(getSisterEdge(iEdge).breakPoints);
                	logHR += removeLociFromAncestors(ipEdge, segsToRemove);
                	
                	ip.removeChildEdge(jEdge);
                	
                	i.removeParentEdge(iEdge);
                	ip.removeChildEdge(iEdge);
                	ip.removeParentEdge(ipEdge);
                	
                	jEdge.parentNode = null;
                	network.setRootEdge(jEdge);

                } else {
                	RecombinationNetworkNode destNode = new RecombinationNetworkNode();
                	destNode.setHeight(newHeight);
                	
                	RecombinationNetworkEdge newIpEdge = new RecombinationNetworkEdge();
                	newIpEdge.breakPoints = newChildEdge.breakPoints.copy();
                	newParent.removeChildEdge(newChildEdge);
                	newParent.addChildEdge(newIpEdge);
                	destNode.addParentEdge(newIpEdge);
                	destNode.addChildEdge(newChildEdge);

                	RecombinationNetworkEdge iEdgeCopy = new RecombinationNetworkEdge();
                	iEdgeCopy.breakPoints = iEdge.breakPoints.copy();
                	i.addParentEdge(iEdgeCopy);
                	destNode.addChildEdge(iEdgeCopy);
                	
                	BreakPoints segsToAdd = iEdgeCopy.breakPoints.copy();
                	segsToAdd.andNot(getSisterEdge(iEdgeCopy).breakPoints);
                	logHR -= addLociToAncestors(newIpEdge, segsToAdd);
                	
                	BreakPoints segsToRemove = iEdge.breakPoints.copy();
                	segsToRemove.andNot(getSisterEdge(iEdge).breakPoints);
                	logHR += removeLociFromAncestors(ipEdge, segsToRemove);
                	
                	pip.removeChildEdge(ipEdge);
                	ip.removeChildEdge(jEdge);
                	pip.addChildEdge(jEdge);
                	
                	
                	i.removeParentEdge(iEdge);
                	ip.removeChildEdge(iEdge);

                }
                
                ip.setHeight(newHeight);
            	
            } else {
            	
                ip.setHeight(newHeight);
            }
        }
        
        if(!networkTerminatesAtMRCA())
        	return Double.NEGATIVE_INFINITY;
        

        
		List<RecombinationNetworkEdge> networkEdgesAfter = new ArrayList<>(network.getEdges());
		
		final List<RecombinationNetworkEdge> possibleEdgesAfter = networkEdgesAfter.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> !e.parentNode.isRecombination())
				.collect(Collectors.toList());
		
		final int possibleNodesAfter = possibleEdgesAfter.size();
        if (possibleNodesAfter < 1) {
            return Double.NEGATIVE_INFINITY;
        }
        
        logHR += Math.log(1.0/possibleNodesAfter);  

		return logHR;
	}
	
	 
    public double checkAncestralSegments(RecombinationNetworkEdge edge) {
    	double logHR = 0.0;
    	
    	RecombinationNetworkNode parentNode = edge.parentNode;
    	if (parentNode == null) {
    		return logHR;
    	}
    	
    	if (parentNode.getParentCount() > 1) {
    		RecombinationNetworkEdge parentEdge = parentNode.getParentEdges().get(0);
    		RecombinationNetworkEdge spouseEdge = parentNode.getParentEdges().get(1);
    		
    		BreakPoints segsToAdd = edge.breakPoints.copy();
    		segsToAdd.andNot(parentEdge.breakPoints);
    		segsToAdd.andNot(spouseEdge.breakPoints);
    		
    		if (!segsToAdd.isEmpty()) {
    			if (Randomizer.nextBoolean()) {
    				logHR -= addLociToAncestors(parentEdge, segsToAdd);
    			} else {
    				logHR -= addLociToAncestors(spouseEdge, segsToAdd);
    			}
    		}
    		
    		BreakPoints segsToRemoveParent = parentEdge.breakPoints.copy();
    		segsToRemoveParent.andNot(edge.breakPoints);
    		if (!segsToRemoveParent.isEmpty()) {
    			logHR += removeLociFromAncestors(parentEdge, segsToRemoveParent);
    		}
    		
    		BreakPoints segsToRemoveSpouse= spouseEdge.breakPoints.copy();
    		segsToRemoveSpouse.andNot(edge.breakPoints);
    		if (!segsToRemoveSpouse.isEmpty()) {
    			logHR += removeLociFromAncestors(spouseEdge, segsToRemoveSpouse);
    		}
	
    	} else {
    		RecombinationNetworkEdge parentEdge = parentNode.getParentEdges().get(0);
    		BreakPoints segsToAdd = edge.breakPoints.copy();
    		segsToAdd.andNot(parentEdge.breakPoints);
    		
    		if (!segsToAdd.isEmpty()) {
    			logHR -= addLociToAncestors(parentEdge, segsToAdd);
    		}
    		
    		BreakPoints segsToRemoveParent = parentEdge.breakPoints.copy();
    		segsToRemoveParent.andNot(edge.breakPoints);
    		if (!segsToRemoveParent.isEmpty()) {
    			logHR += removeLociFromAncestors(parentEdge, segsToRemoveParent);
    		}
    	}
    	return logHR;
    }
	
	



    private double getDelta() {
        if (!gaussianInput.get()) {
            return (Randomizer.nextDouble() * size) - (size / 2.0);
        } else {
            return Randomizer.nextGaussian() * size;
        }
    }
    
    private int intersectingEdges(RecombinationNetworkEdge edge, double height, List<RecombinationNetworkEdge> directChildEdges) {
        final RecombinationNetworkNode node = edge.childNode;
    	final RecombinationNetworkNode parent = edge.parentNode;
    	
        if (parent.getHeight() < height) return 0;


        if (node.getHeight() < height) {
            if (directChildEdges != null) directChildEdges.add(edge);
            return 1;
        }

        if (node.getChildCount() > 1 && node.getChildEdges().get(0).childNode == node.getChildEdges().get(1).childNode) {
        	int count = 0;
        	if (node.getChildEdges().get(0).childNode.getHeight() < height) {
        		if (directChildEdges != null) {
        			directChildEdges.add(node.getChildEdges().get(0));
        			directChildEdges.add(node.getChildEdges().get(1));
        		}
        		return 2;
        	}

        	count += intersectingEdges(node.getChildEdges().get(0).childNode.getChildEdges().get(0), height, directChildEdges);        	
        	return count;
        }
        if (node.isLeaf()) {
            // TODO: verify that this makes sense
            return 0;
        } else {
            int count = intersectingEdges(node.getChildEdges().get(0), height, directChildEdges);
            if (node.getChildCount() > 1) count += 
            		intersectingEdges(node.getChildEdges().get(1), height, directChildEdges);          
            return count;
        }
    }  
    
    /**
     * Simple (but probably too expensive) check for a kind of invalid network
     * which can result from an edge deletion operation: one in which the
     * network posesses nontrivial structure above the MRCA. (I.e. the MRCA
     * is not the root.)
     *
     * @return true if the network terminates at the true MRCA. (I.e. is valid.)
     */
    protected boolean networkTerminatesAtMRCA() {
        List<RecombinationNetworkNode> sortedNodes = new ArrayList<>(network.getNodes());
        sortedNodes.sort(Comparator.comparingDouble(RecombinationNetworkNode::getHeight));
        List<RecombinationNetworkNode> sampleNodes = sortedNodes.stream().filter(RecombinationNetworkNode::isLeaf).collect(Collectors.toList());
        double maxSampleHeight = sampleNodes.get(sampleNodes.size()-1).getHeight();

        int lineages = 0;
        for (RecombinationNetworkNode node : sortedNodes) {
            switch(node.getChildEdges().size()) {
                case 2:
                    // Coalescence

                    lineages -= 1;
                    break;

                case 1:
                    // Reassortment

                    if (lineages < 2 && node.getHeight() > maxSampleHeight)
                        return false;

                    lineages += 1;
                    break;

                case 0:
                    // Sample

                    lineages += 1;
                    break;
            }
        }

        return true;
    }
    
    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(size);
            final double f = Math.exp(delta);
            if( limit > 0 ) {
                final RecombinationNetwork network = networkInput.get();
                final double h = network.getRootEdge().childNode.getHeight();
                final double k = Math.log(network.getLeafNodes().size()) / Math.log(2);
                final double lim = (h / k) * limit;
                if( f <= lim ) {
                    size = f;
                }
            } else {
               size = f;
            }
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return size;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        size = value;
    }
    
    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        final double newDelta = size * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try decreasing size to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try increasing size to about " + formatter.format(newDelta);
        } else return "";
    }
}
