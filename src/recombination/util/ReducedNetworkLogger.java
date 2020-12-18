/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com> adapted by Nicola Mueller
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package recombination.util;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.datatype.DataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.util.BEASTClassLoader;
import beast.util.PackageManager;
import beast.util.Randomizer;
import feast.nexus.CharactersBlock;
import feast.nexus.NexusBuilder;
import feast.nexus.TaxaBlock;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Tim Vaughan, adapted by Nicola Mueller for recombination networks
 */
@Description("A more flexible alignment simulator. Doesn't require " +
        "pre-specification of number of taxa.")
public class ReducedNetworkLogger extends BEASTObject implements Loggable {
    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

    final public Input<RecombinationNetwork> networkInput = new Input<>("recombinationNetwork", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

    
	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}


    
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
    	RecombinationNetwork network = getPrunedNetwork(networkInput.get());
        out.println("tree STATE_" + sample + " = " + network.getExtendedNewick());
    }

	public RecombinationNetwork getPrunedNetwork(RecombinationNetwork recombinationNetwork) {
		RecombinationNetwork network = new RecombinationNetwork(recombinationNetwork.getRootEdge().getCopy());
		network.totalLength = recombinationNetwork.totalLength;
		for (RecombinationNetworkNode leaf : network.getLeafNodes()) {
			BreakPoints pruneBP = new BreakPoints(network.totalLength);
			pruneBP.andNot(getGappedBreakPoints(dataInput.get(), leaf));
			prune(leaf.getParentEdges().get(0), pruneBP);
		}
		return network;
	}
	
	public RecombinationNetwork getPrunedNetwork(RecombinationNetwork recombinationNetwork, Alignment data) {
		RecombinationNetwork network = new RecombinationNetwork(recombinationNetwork.getRootEdge().getCopy());
		network.totalLength = recombinationNetwork.totalLength;
		for (RecombinationNetworkNode leaf : network.getLeafNodes()) {
			BreakPoints pruneBP = new BreakPoints(network.totalLength);
			pruneBP.andNot(getGappedBreakPoints(data, leaf));
			prune(leaf.getParentEdges().get(0), pruneBP);
		}
		return network;
	}

	
    private void prune(RecombinationNetworkEdge edge, BreakPoints pruneBP) {
		if (pruneBP.isEmpty())
			return;
		
		edge.breakPoints.andNot(pruneBP);
		
		if (edge.isRootEdge())
			return;
		
		RecombinationNetworkNode n = edge.parentNode;
		
		if (n.isCoalescence()) {
			// check if loop
			if (n.getChildEdges().get(0).childNode.ID==n.getChildEdges().get(1).childNode.ID) {
				prune(n.getParentEdges().get(0), pruneBP);
			}else {				
				RecombinationNetworkEdge otherEdge = edge.childNode.ID == n.getChildEdges().get(0).childNode.ID ? n.getChildEdges().get(1) : n.getChildEdges().get(0);
				BreakPoints bp = pruneBP.copy();
				bp.andNot(otherEdge.breakPoints);
				
				prune(n.getParentEdges().get(0), bp);
			}			
		}else {
			for (RecombinationNetworkEdge e : n.getParentEdges()) {
				BreakPoints bp = pruneBP.copy();
				bp.and(e.passingRange);
				prune(e, bp);				
			}
		}		
		
	}



	private BreakPoints getGappedBreakPoints(Alignment data, RecombinationNetworkNode n) {
        int taxonIndex = data.getTaxonIndex(n.getTaxonLabel());
        if (taxonIndex == -1) {
        	if (n.getTaxonLabel().startsWith("'") || n.getTaxonLabel().startsWith("\"")) {
                taxonIndex = data.getTaxonIndex(n.getTaxonLabel().substring(1, n.getTaxonLabel().length() - 1));
            }
            if (taxonIndex == -1) {
            	throw new RuntimeException("Could not find sequence " + n.getTaxonLabel() + " in the alignment");
            }
        }
        
        int code, states;
        int[] statesForCode;
        boolean on = false;
        List<Integer> bp_list = new ArrayList<>();
        
        
    	for (int i = 0; i < data.getSiteCount() ;i++) {
            code = data.getPattern(taxonIndex, data.getPatternIndex(i));
            
            statesForCode = data.getDataType().getStatesForCode(code);
            if (statesForCode.length==1)
                states = statesForCode[0];
            else
                states = code; // Causes ambiguous states to be ignored.

            if (i==0 && states<data.getMaxStateCount()) {
            	bp_list.add(i);
            	on = true;
            }
            
            if (on && states>=data.getMaxStateCount()) {
            	bp_list.add(i-1);
            	on = false;
            }      
            
            if (!on && states<data.getMaxStateCount()) {
            	bp_list.add(i);
            	on = true;
            }           
            	
    	}
    	
    	if (on) {
        	bp_list.add(data.getSiteCount()-1);
    	}    	
    	
    	
    	BreakPoints bp = new BreakPoints();
    	bp.init(bp_list);
		return bp;
	}



}