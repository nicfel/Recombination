package recombination.distribution;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.coalescent.ConstantPopulation;
import coalre.CoalReTestClass;
import junit.framework.Assert;
import recombination.network.RecombinationNetwork;

import org.junit.Test;

public class CoalescentWithRecombinationTest extends CoalReTestClass {

    @Test
    public void testDensity() {
        RecombinationNetwork network = new RecombinationNetwork(
                "((((t0[&loci={0-9999},length=10000]:0.5045351708916129)#H0[&split={4867-9999},loci={4867-9999},length=5133]"+
                		":0.2890189030749166,#H0[&split={0-4866},loci={0-4866},length=4867]:0.2890189030749166)[&loci={0-9999},length=10000]"+
                		":0.6475481073244761,((t2[&loci={0-9999},length=10000]:0.903539606702267,t1[&loci={0-9999},length=10000]:0.903539606702267)"+
                		"[&loci={0-9999},length=10000]:0.29229242390613575)#H1[&split={0-5630},loci={0-5630},length=5631]:0.24527015068260294)"+
                		"[&loci={0-9999},length=10000]:1.2362149206344832,((#H1[&split={5631-9999},loci={5631-9999},length=4369]"+
                		":0.05103233536942797)#H2[&split={7546-9999},loci={7546-9999},length=2454]:0.3230648558859137"+
                		",#H2[&split={0-7545},loci={5631-7545},length=1915]:0.3230648558859137)[&loci={5631-9999},length=4369]:"+
                		"1.1073878800617445)[&loci={0-9999},length=10000]:0.0;");

        RecombinationNetworkIntervals networkIntervals = new RecombinationNetworkIntervals();
        networkIntervals.initByName("recombinationNetwork", network);

        ConstantPopulation populationFunction = new ConstantPopulation();
        populationFunction.initByName("popSize", new RealParameter("1.0"));

        CoalescentWithRecombination coalWR = new CoalescentWithRecombination();
        coalWR.initByName("networkIntervals", networkIntervals,
                "recombinationRate", new RealParameter("0.0001"),
                "populationModel", populationFunction);
        
        
        Assert.assertEquals(-13.072190353529038, coalWR.calculateLogP(), 1e-10);
    }
}
