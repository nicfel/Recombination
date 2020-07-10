package recombination.operator;

import beast.util.Randomizer;
import coalre.CoalReTestClass;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkNode;
import recombination.operators.DivertLociOperator;

import org.junit.Assert;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

public class DivertLociTest extends CoalReTestClass {

    String networkString = "((((#H1[&split={7357-9999},loci={7357-9999},length=2643.0]:0.029028801150026373,(((t7[&loci={0-9999},length=10000.0]:0.3408331510142315,(t8[&loci={0-9999},length=10000.0]:0.14470426889934196,((t4[&loci={0-9999},length=10000.0]:0.008921363131720263,t3[&loci={0-9999},length=10000.0]:0.008921363131720263)[&loci={0-9999},length=10000.0]:0.07368958619017191,(t5[&loci={0-9999},length=10000.0]:0.05651985106014486,((t0[&loci={0-9999},length=10000.0]:0.01371370813629106,t2[&loci={0-9999},length=10000.0]:0.01371370813629106)[&loci={0-9999},length=10000.0]:0.01770276986295235,t1[&loci={0-9999},length=10000.0]:0.03141647799924341)[&loci={0-9999},length=10000.0]:0.02510337306090145)[&loci={0-9999},length=10000.0]:0.02609109826174731)[&loci={0-9999},length=10000.0]:0.06209331957744979)[&loci={0-9999},length=10000.0]:0.19612888211488955)[&loci={0-9999},length=10000.0]:0.2532640786850209)#H2[&split={499-9999},loci={499-9999},length=9501.0]:0.2994016137375858)#H1[&split={0-7356},loci={499-7356},length=6858.0]:0.029028801150026373)[&loci={499-9999},length=9501.0]:0.34436738009841983)#H0[&split={3613-9999},loci={3613-9999},length=6387.0]:0.16198962047105492,#H0[&split={0-3612},loci={499-3612},length=3114.0]:0.16198962047105492)[&loci={499-9999},length=9501.0]:0.2105449461319744,(#H3[&split={0-397},loci={0-397},length=398.0]:0.3596040657302466,(#H2[&split={0-498},loci={0-498},length=499.0]:0.36308767696719646,((t6[&loci={0-9999},length=10000.0]:0.08009733219216959,t9[&loci={0-9999},length=10000.0]:0.08009733219216959)[&loci={0-9999},length=10000.0]:0.703681773516321)#H3[&split={398-9999},loci={398-9999},length=9602.0]:0.1734058009579582)[&loci={0-9999},length=10000.0]:0.1861982647722884)[&loci={0-9999},length=10000.0]:0.4960464198495764)[&loci={0-9999},length=10000.0]:0.0;";

    @Test
    public void testAddRemoveSegment() {

        RecombinationNetwork network = new RecombinationNetwork(networkString);
        RecombinationNetworkNode leafNode = new ArrayList<>(network.getLeafNodes()).get(0);

        DivertLociOperator operator = new DivertLociOperator();
        
        System.out.println(network);

        List<Integer> listToAdd = new ArrayList<Integer>();
        BreakPoints lociToAdd = new BreakPoints();
        
        listToAdd.add(500);
        listToAdd.add(800);
        lociToAdd.init(listToAdd);
        
        
        double logPremove = operator.removeLociFromAncestors(
              leafNode.getParentEdges().get(0), lociToAdd);
        
        System.out.println(network);
        double logPadd = operator.addLociToAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);
        
        System.out.println(network);

//        BitSet allSegments = new BitSet();
//        allSegments.set(0, 15);
//        Assert.assertEquals(allSegments, network.getRootEdge().hasSegments);
//
//        BitSet segmentsToRemove = new BitSet();
//        segmentsToRemove.set(8, 15);
//
//        double logPremove = operator.removeSegmentsFromAncestors(
//                leafNode.getParentEdges().get(0), segmentsToRemove);
//
//        Assert.assertEquals(networkString, network.toString());
//        Assert.assertEquals(logPadd, logPremove, 1e-10);
    }

//    @Test
//    public void testAddRemoveReassortmentEdge() {
//        Network network = getContempNetwork(2, 8, 0.0);
//
//        AddRemoveReassortment operator = new AddRemoveReassortment();
//        operator.initByName("alpha", 1.0,
//                "network", network,
//                "weight", 1.0);
//
//        NetworkNode origRoot = network.getRootEdge().childNode;
//
//        NetworkEdge sourceEdge = network.getRootEdge().childNode.getChildEdges().get(0);
//        double sourceTime = sourceEdge.getLength()/2.0;
//        NetworkEdge destEdge = network.getRootEdge();
//        double destTime = destEdge.childNode.getHeight() + 1.0;
//
//        double logP1 = operator.addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);
//
//        NetworkEdge edgeToRemove = sourceEdge.parentNode.getParentEdges().get(0).parentNode == origRoot
//                ? sourceEdge.parentNode.getParentEdges().get(1)
//                : sourceEdge.parentNode.getParentEdges().get(0);
//
//        double logP2 = operator.removeReassortmentEdge(edgeToRemove);
//
//        Assert.assertEquals(logP1, -logP2, 1e-10);
//
//        sourceEdge = network.getRootEdge().childNode.getChildEdges().get(1);
//        sourceTime = sourceEdge.getLength()/4.0;
//        destEdge = sourceEdge;
//        destTime = sourceEdge.getLength()*3.0/4.0;
//
//        logP1 = operator.addReassortmentEdge(sourceEdge, sourceTime, destEdge, destTime);
//
//        edgeToRemove = sourceEdge.parentNode.getParentEdges().get(0);
//
//        logP2 = operator.removeReassortmentEdge(edgeToRemove);
//
//        Assert.assertEquals(logP1, -logP2, 1e-10);
//    }
//
//    @Test
//    public void testRemoveReassortment() {
//        // TODO Flesh out this test
//
//        Network network = new Network(networkString);
//
//        AddRemoveReassortment operator = new AddRemoveReassortment();
//        operator.initByName("network", network,
//                "alpha", 1.0,
//                "weight", 1.0);
//
//        System.out.println(network.getExtendedNewickVerbose());
//
//        operator.removeReassortment();
//
//        System.out.println(network.getExtendedNewickVerbose());
//    }
//
//    @Test
//    public void testAddReassortment() {
//        // TODO Flesh out this test
//
//        Network network = new Network(networkString);
//
//        AddRemoveReassortment operator = new AddRemoveReassortment();
//        operator.initByName("network", network,
//                "alpha", 1.0,
//                "weight", 1.0);
//
//        System.out.println(network.getExtendedNewickVerbose());
//
//        double logHR = operator.addReassortment();
//
//        System.out.println(network.getExtendedNewickVerbose());
//
//        System.out.println(logHR);
//    }
}
