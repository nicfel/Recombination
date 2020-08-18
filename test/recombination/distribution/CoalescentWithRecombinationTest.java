package recombination.distribution;

import beast.core.parameter.RealParameter;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.coalescent.ConstantPopulation;
import coalre.CoalReTestClass;
import coalre.distribution.CoalescentWithReassortment;
import coalre.distribution.NetworkIntervals;
import coalre.network.Network;
import junit.framework.Assert;
import recombination.network.RecombinationNetwork;

import org.junit.Test;

public class CoalescentWithRecombinationTest extends CoalReTestClass {

    @Test
    public void testDensity() {
        RecombinationNetwork network = new RecombinationNetwork(
        		"((#H0[&split={0-929},loci={0-929},length=930]:0.4860702162314561,((#H2[&split={7633-9999},loci={7633-9999},length=2367]:0.036380107508342974,(t1[&loci={0-9999},length=10000]:0.29041085418573037,t0[&loci={0-9999},length=10000]:0.29041085418573037)[&loci={0-9999},length=10000]:0.1528435079144485)[&loci={0-9999},length=10000]:0.47005038739308824)#H1[&split={0-5172},loci={0-5172},length=5173]:1.5817575814045295)[&loci={0-5172},length=5173]:0.35688825595463936,(((t2[&loci={0-9999},length=10000]:0.4068742545918359)#H2[&split={0-7632},loci={0-7632},length=7633]:0.604218532359179,#H1[&split={5173-9999},loci={5173-9999},length=4827]:0.09778803745774778)[&loci={0-9999},length=10000]:0.9978993277153256)#H0[&split={930-9999},loci={930-9999},length=9070]:0.8429584721860954)[&loci={0-9999},length=10000]:0.0;");



        RecombinationNetworkIntervals networkIntervals = new RecombinationNetworkIntervals();
        networkIntervals.initByName("recombinationNetwork", network);

        ConstantPopulation populationFunction = new ConstantPopulation();
        populationFunction.initByName("popSize", new RealParameter("1.0"));
        

        CoalescentWithRecombination coalWR = new CoalescentWithRecombination();
        coalWR.initByName("networkIntervals", networkIntervals,
                "recombinationRate", new RealParameter("0.0001"),
                "populationModel", populationFunction);

        Assert.assertEquals(-9.611019582493, coalWR.calculateLogP(), 1e-10);
        
        //matlab code for control
//        intervals = [0.0,0.0,0.2904108541857302,0.11646340040610559,0.036380107508342974,0.47005038739308835,0.09778803745774778,0.9978993277153256,0.4860702162314561,0.35688825595463936];
//        lins = [1,2,3,2,3,2,3,2,3,2];
//        rec = [9999.0,19998.0,29997.0,19998.0,19997.0,17631.0,17630.0,15171.0,15170.0,14241.0];
//
//        prob = 5*log(1);
//        prob = prob + 3*log(0.0001*9999);
//
//        for i  = 1 : length(lins)
//            prob = prob - 0.0001 *rec(i)*intervals(i);
//            prob = prob - lins(i)*(lins(i)-1)/(2)*intervals(i)
//        end
        
        
        
        RecombinationNetwork network2 = new RecombinationNetwork(
                "(t1[&loci={0-2},length=3]:2.940457397943351,(t4[&loci={0-2},length=3]:0.4837103008790285,(t2[&loci={0-2},length=3]:0.46866995123588273,(t5[&loci={0-2},length=3]:0.09369026202003683,t3[&loci={0-2},length=3]:0.29369026202003684)[&loci={0-2},length=3]:0.07497968921584586)[&loci={0-2},length=3]:0.21504034964314578)[&loci={0-2},length=3]:2.1567470970643225)[&loci={0-2},length=3]:0.0;");
        
        
        RecombinationNetwork network3 = new RecombinationNetwork(
                "(t1[&loci={0-1},length=2]:2.940457397943351,(t4[&loci={0-1},length=2]:0.4837103008790285,(t2[&loci={0-1},length=2]:0.46866995123588273,(t5[&loci={0-1},length=2]:0.09369026202003683,t3[&loci={0-1},length=2]:0.29369026202003684)[&loci={0-1},length=2]:0.07497968921584586)[&loci={0-1},length=2]:0.21504034964314578)[&loci={0-1},length=2]:2.1567470970643225)[&loci={0-1},length=2]:0.0;");

        
        RecombinationNetworkIntervals networkIntervals2 = new RecombinationNetworkIntervals();
        networkIntervals2.initByName("recombinationNetwork", network2);

        RecombinationNetworkIntervals networkIntervals3 = new RecombinationNetworkIntervals();
        networkIntervals3.initByName("recombinationNetwork", network3);
        
        CoalescentWithRecombination coalWR2 = new CoalescentWithRecombination();
        coalWR2.initByName("networkIntervals", networkIntervals2,
                "recombinationRate", new RealParameter("0.5"),
                "populationModel", populationFunction);
        
        CoalescentWithRecombination coalWR3 = new CoalescentWithRecombination();
        coalWR3.initByName("networkIntervals", networkIntervals3,
                "recombinationRate", new RealParameter("1.0"),
                "populationModel", populationFunction);
        
        Assert.assertEquals(coalWR2.calculateLogP(), coalWR3.calculateLogP(), 1e-10);



    }
}
