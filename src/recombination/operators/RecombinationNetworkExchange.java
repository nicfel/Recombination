package recombination.operators;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;


import beast.core.Input;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.operators.DivertSegmentOperator;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

public class RecombinationNetworkExchange extends DivertLociOperator {
	final public Input<Boolean> isNarrowInput = new Input<>("isNarrow",
			"if true (default) a narrow exchange is performed, "
					+ "otherwise a wide exchange", true);


	private boolean isNarrow;

	@Override
	public void initAndValidate() {
		super.initAndValidate();

		isNarrow = isNarrowInput.get();
	}

	@Override
	public double networkProposal() {

		double logHR;
		network.startEditing(this);

		if (isNarrow) {
			logHR = narrow(network);
		} else {
			logHR = wide(network);
		}
		return logHR;
	}
	
	private boolean hasCoalescenceKid(final RecombinationNetworkNode n) {
		return (n.getChildEdges().get(0).childNode.isCoalescence() ||
				n.getChildEdges().get(1).childNode.isCoalescence());
	}	


	/**
	 * Perform equivalent of narrow tree exchange on a network.
	 *
	 * @param	network
	 * @return	log of Hastings Ratio, or Double.NEGATIVE_INFINITY
	 * 			if proposal should not be accepted
	 */
	public double narrow(final RecombinationNetwork network) {
		
		double logHR = 0.0;

		List<RecombinationNetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

		final List<RecombinationNetworkEdge> possibleGrandParentEdges = networkEdges.stream()
				.filter(e -> e.childNode.isCoalescence())
				.filter(e -> hasCoalescenceKid(e.childNode))
				.collect(Collectors.toList());
		
		
		final int possibleGrandParents = possibleGrandParentEdges.size();
        if (possibleGrandParents < 1) {
            return Double.NEGATIVE_INFINITY;
        }
        logHR -= Math.log(1.0/possibleGrandParents);        

		final RecombinationNetworkEdge grandParentEdge = possibleGrandParentEdges.
				get(Randomizer.nextInt(possibleGrandParents));
		final RecombinationNetworkNode grandParent = grandParentEdge.childNode;

		final List<RecombinationNetworkEdge> possibleParentEdges = grandParent.getChildEdges();
		RecombinationNetworkEdge parentEdge = possibleParentEdges.get(0);
		RecombinationNetworkEdge auntEdge = possibleParentEdges.get(1);
		
		RecombinationNetworkNode parent = parentEdge.childNode;
		RecombinationNetworkNode aunt = auntEdge.childNode;
		
		
		if (parent.getHeight() < aunt.getHeight()) {
			auntEdge = possibleParentEdges.get(0);
			parentEdge = possibleParentEdges.get(1);

			parent = parentEdge.childNode;
			aunt = auntEdge.childNode;
		}

		if( !parent.isCoalescence()) {
			return Double.NEGATIVE_INFINITY;
		}

		
		final int childId = Randomizer.nextInt(parent.getChildEdges().size());
		final RecombinationNetworkEdge childEdge = parent.getChildEdges().get(childId);

		
		logHR += exchangeEdges(childEdge, auntEdge, parent, grandParent);
		
		networkEdges = new ArrayList<>(network.getEdges());
		
		final List<RecombinationNetworkEdge> possibleGrandParentEdgesAfter = networkEdges.stream()
				.filter(e -> e.childNode.isCoalescence())
				.filter(e -> hasCoalescenceKid(e.childNode))
				.collect(Collectors.toList());
		
		final int possibleGrandParentsAfter = possibleGrandParentEdgesAfter.size();

		logHR += Math.log(1.0/possibleGrandParentsAfter);
		
		
		return logHR;
	}

	/**
	 * Perform equivalent of wide tree exchange on a network.
	 *
	 * @param	network
	 * @return	log of Hastings Ratio, or Double.NEGATIVE_INFINITY
	 * 			if proposal should not be accepted
	 */
	public double wide(final RecombinationNetwork network) {
		double logHR = 0.0;

		List<RecombinationNetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

		final List<RecombinationNetworkEdge> possibleEdges = networkEdges.stream()
				.filter(e -> !e.isRootEdge())
				.filter(e -> !e.parentNode.isRecombination())
				.collect(Collectors.toList());
		
		final int nPossibleEdges = possibleEdges.size();
		logHR -= Math.log(1.0/(double)nPossibleEdges);

		final RecombinationNetworkEdge iEdge = possibleEdges.
				get(Randomizer.nextInt(possibleEdges.size()));
		final RecombinationNetworkNode i = iEdge.childNode;

		RecombinationNetworkEdge jEdge = iEdge;

		while(jEdge == iEdge) {
			jEdge = possibleEdges.
					get(Randomizer.nextInt(possibleEdges.size()));
		}
		final RecombinationNetworkNode j = jEdge.childNode;

		final RecombinationNetworkNode p = iEdge.parentNode;
		final RecombinationNetworkNode jP = jEdge.parentNode;


		if ((p != jP) && (i !=jP) && (j != p)
				&& (j.getHeight() < p.getHeight())
				&& (i.getHeight() < jP.getHeight())) {


			logHR += exchangeEdges(iEdge, jEdge, p, jP);
			
			
			networkEdges = new ArrayList<>(network.getEdges());
			
			final List<RecombinationNetworkEdge> possibleEdgesAfter = networkEdges.stream()
					.filter(e -> !e.isRootEdge())
					.filter(e -> !e.parentNode.isRecombination())
					.collect(Collectors.toList());
			
			final int nPossibleEdgesAfter = possibleEdgesAfter.size();
			logHR += Math.log(1.0/(double)nPossibleEdgesAfter);
	        
			return logHR;
		}
		else {
			return Double.NEGATIVE_INFINITY;
		}
	}


	/* exchange sub-nets whose root are i and j */
	protected double exchangeEdges(RecombinationNetworkEdge iEdge, RecombinationNetworkEdge jEdge,
			RecombinationNetworkNode p, RecombinationNetworkNode jP) {
		double logHR = 0.0;
		
		final RecombinationNetworkEdge pEdge = p.getParentEdges().get(0);
		final RecombinationNetworkEdge jPEdge = jP.getParentEdges().get(0);
		
		// After the exchange we want to add segments to the new ancestors
		// and remove from the old. Have to be careful not to remove segments
		// of siblings.
		final BreakPoints iSegs = iEdge.breakPoints;
		final BreakPoints jSegs = jEdge.breakPoints;
		
		final BreakPoints iSegsToRemove = iSegs.copy();
		iSegsToRemove.andNot(getSisterEdge(iEdge).breakPoints);
		iSegsToRemove.andNot(jSegs);

		final BreakPoints jSegsToRemove = jSegs.copy();
		jSegsToRemove.andNot(getSisterEdge(jEdge).breakPoints);
		jSegsToRemove.andNot(iSegs);
		
		final BreakPoints iSegsToAdd = iSegs.copy();
		iSegsToAdd.andNot(jSegs);
		
		final BreakPoints jSegsToAdd = jSegs.copy();
		jSegsToAdd.andNot(iSegs);
		
		p.removeChildEdge(iEdge);
		jP.removeChildEdge(jEdge);
		p.addChildEdge(jEdge);
		jP.addChildEdge(iEdge);
		
		logHR += removeLociFromAncestors(jPEdge, jSegsToRemove);
		logHR += removeLociFromAncestors(pEdge, iSegsToRemove);
		
		logHR -= addLociToAncestors(jPEdge, iSegsToAdd);
		logHR -= addLociToAncestors(pEdge, jSegsToAdd);
		
		return logHR;
	}

}
