package recombination.simulator;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.util.Randomizer;
import recombination.CoalReTestClass;
import recombination.statistics.RecombinationNetworkStatsLogger;

import org.junit.Assert;
import org.junit.Test;
import test.beast.beast2vs1.trace.DiscreteStatistics;


public class SimulatedCoalescentRecombinationNetworkTest extends CoalReTestClass {

    @Test
    public void testSimulator() {
        Randomizer.setSeed(1);

        TraitSet dateTrait = getContempDateTraitSet(getTaxonSet(10));

        ConstantPopulation populationFunction = new ConstantPopulation();
        populationFunction.initByName("popSize", new RealParameter("1.0"));

        int N = 10000;
        double[] reassortmentNodeCounts = new double[N];
        double[] networkHeights = new double[N];
        double[] networkLengths = new double[N];

        for (int i = 0; i < N; i++) {

            SimulatedCoalescentRecombinationNetwork network = new SimulatedCoalescentRecombinationNetwork();
            network.initByName(
                    "recombinationRate", new RealParameter("0.0001"),
                    "populationModel", populationFunction,
                    "traitSet", dateTrait,
                    "totalLength", 10000);

            reassortmentNodeCounts[i] = RecombinationNetworkStatsLogger.getRecombinationCount(network);
            networkHeights[i] = RecombinationNetworkStatsLogger.getTotalHeight(network);
            networkLengths[i] = RecombinationNetworkStatsLogger.getTotalEdgeLength(network);            
        }

        double meanCount = DiscreteStatistics.mean(reassortmentNodeCounts);
        double meanHeight = DiscreteStatistics.mean(networkHeights);
        double meanLength = DiscreteStatistics.mean(networkLengths);

        System.out.println(meanCount);
        System.out.println(meanHeight);
        System.out.println(meanLength);

        Assert.assertEquals(6.84, meanCount, 0.1);
        Assert.assertEquals(2.78, meanHeight, 0.1);
        Assert.assertEquals(9.56, meanLength, 0.5);
    }
}
