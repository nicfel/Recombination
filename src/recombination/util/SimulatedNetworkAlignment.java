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

import beast.core.Description;
import beast.core.Input;
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
public class SimulatedNetworkAlignment extends Alignment {

    public Input<RecombinationNetwork> networkInput = new Input<>(
            "recombinationNetwork",
            "Tree down which to simulate sequence evolution.",
            Input.Validate.REQUIRED);

    public Input<SiteModel> siteModelInput = new Input<>(
            "siteModel",
            "Site model to use in simulation.",
            Input.Validate.REQUIRED);

    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName",
            "Name of file (if any) simulated alignment should be saved to.");

    private RecombinationNetwork network;
    private SiteModel siteModel;
    private int seqLength;
    private DataType dataType;
    private int[][] alignment;

    private String ancestralSeqStr;

    public SimulatedNetworkAlignment() {
        sequenceInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

    	network = networkInput.get();
        siteModel = siteModelInput.get();
        seqLength = network.totalLength;

        sequences.clear();

        grabDataType();

        simulate();

        super.initAndValidate();

        // Write simulated alignment to disk if required
        if (outputFileNameInput.get() != null) {
            try (PrintStream pstream = new PrintStream(outputFileNameInput.get())) {
                NexusBuilder nb = new NexusBuilder();
                nb.append(new TaxaBlock(new TaxonSet(this)));
                nb.append(new CharactersBlock(this));
                nb.write(pstream);
            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error writing to file "
                        + outputFileNameInput.get() + ".");
            }
        }
    }

    /**
     * Perform actual sequence simulation.
     */
    private void simulate() {
        List<RecombinationNetworkNode> networkNodes = new ArrayList<>(network.getNodes());
        // enumerate nodes
        int c=0;
        for (RecombinationNetworkNode n : networkNodes) {
        	n.setTaxonIndex(c);
        	c++;
        }

        int nNodes = network.getNodes().size();
        Node dummyNode = new Node();

        double[] categoryProbs = siteModel.getCategoryProportions(dummyNode);

        int nCategories = siteModel.getCategoryCount();
        int nStates = dataType.getStateCount();
        double[][] transitionProbs = new double[nCategories][nStates*nStates];

        alignment = new int[nNodes][seqLength];

        int[] categories = new int[seqLength];
        for (int i=0; i<seqLength; i++)
            categories[i] = Randomizer.randomChoicePDF(categoryProbs);

        RecombinationNetworkEdge edge = network.getRootEdge();

        int[] parentSequence = new int[seqLength];
        double[] frequencies = siteModel.getSubstitutionModel().getFrequencies();
        for (int i=0; i<parentSequence.length; i++)
        	alignment[network.getRootEdge().childNode.getTaxonIndex()][i] = Randomizer.randomChoicePDF(frequencies);
        
        
        ancestralSeqStr = dataType.encodingToString(parentSequence);
        traverse(edge,
                categories, transitionProbs,
                edge.breakPoints);
        
    	List<RecombinationNetworkNode> leafs = networkNodes.stream()
                .filter(e -> e.isLeaf())
                .collect(Collectors.toList());


        for (RecombinationNetworkNode l : leafs) {
            String seqString = dataType.encodingToString(alignment[l.getTaxonIndex()]);

            String taxonName;
            if (l.getTaxonLabel() != null)
                taxonName = l.getTaxonLabel();
            else
                taxonName = "t" + l.getTaxonIndex();

            sequenceInput.setValue(new Sequence(taxonName, seqString), this);
        }
    }

    /**
     * Traverse a tree, simulating a sequence alignment down it.
     *
     * @param node Node of the tree
     * @param parentSequence Sequence at the parent node in the tree
     * @param categories Mapping from sites to categories
     * @param transitionProbs transition probabilities
     * @param regionAlignment alignment for particular region
     */
    private void traverse(RecombinationNetworkEdge node,
            int[] categories, double[][] transitionProbs,
            BreakPoints computeFor) {

        Node dummyNode = new Node();

        for (RecombinationNetworkEdge child : node.childNode.getChildEdges()) {
            BreakPoints childBP = computeFor.copy();
            childBP.and(child.breakPoints);
            
            if (!childBP.isEmpty()) {            	

	            // Calculate transition probabilities
	            for (int i=0; i<siteModel.getCategoryCount(); i++) {
	                siteModel.getSubstitutionModel().getTransitionProbabilities(
	                		dummyNode, node.childNode.getHeight(), child.childNode.getHeight(),
	                        siteModel.getRateForCategory(i, dummyNode),
	                        transitionProbs[i]);
	            }
	            
	            
	            
	
	
	            // Draw characters on child sequence
	            int nStates = dataType.getStateCount();
	            double[] charProb = new double[nStates];
	            for (int i=0; i<network.totalLength; i++) {
	            	if (childBP.contains(i)) {
		                int category = categories[i];
		                System.arraycopy(transitionProbs[category],
		                		alignment[child.parentNode.getTaxonIndex()][i]*nStates, charProb, 0, nStates);
		                alignment[child.childNode.getTaxonIndex()][i] = Randomizer.randomChoicePDF(charProb);
	            	}
	            }
	
	            if (!child.childNode.isLeaf()) {
	                traverse(child,
	                        categories, transitionProbs,
	                        childBP);
	            }
            }
        }
    }

    /**
     * HORRIBLE function to identify data type from given description.
     */
    private void grabDataType() {
        if (userDataTypeInput.get() != null) {
            dataType = userDataTypeInput.get();
        } else {

            List<String> dataTypeDescList = new ArrayList<>();
            List<String> classNames = PackageManager.find(beast.evolution.datatype.DataType.class, "beast.evolution.datatype");
            for (String className : classNames) {
                try {
                    DataType thisDataType = (DataType) BEASTClassLoader.forName(className).newInstance();
                    if (dataTypeInput.get().equals(thisDataType.getTypeDescription())) {
                        dataType = thisDataType;
                        break;
                    }
                    dataTypeDescList.add(thisDataType.getTypeDescription());
                } catch (ClassNotFoundException
                    | InstantiationException
                    | IllegalAccessException e) {
                }
            }
            if (dataType == null) {
                throw new IllegalArgumentException("Data type + '"
                        + dataTypeInput.get()
                        + "' cannot be found.  Choose one of "
                        + Arrays.toString(dataTypeDescList.toArray(new String[0])));
            }
        }
    }

}