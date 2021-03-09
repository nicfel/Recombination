package recombination.operator;

import beast.util.Randomizer;
import recombination.CoalReTestClass;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;
import recombination.operators.AddRemoveRecombination;
import recombination.operators.DivertLociOperator;

import org.junit.Assert;
import org.junit.Test;

import java.util.*;

public class DivertLociTest extends CoalReTestClass {

    String networkString = "((#H0[&split={0-929},loci={0-929},length=930]:0.4860702162314561,((#H2[&split={7633-9999},loci={7633-9999},length=2367]:0.036380107508342974,(t1[&loci={0-9999},length=10000]:0.29041085418573037,t0[&loci={0-9999},length=10000]:0.29041085418573037)[&loci={0-9999},length=10000]:0.1528435079144485)[&loci={0-9999},length=10000]:0.47005038739308824)#H1[&split={0-5172},loci={0-5172},length=5173]:1.5817575814045295)[&loci={0-5172},length=5173]:0.35688825595463936,(((t2[&loci={0-9999},length=10000]:0.4068742545918359)#H2[&split={0-7632},loci={0-7632},length=7633]:0.604218532359179,#H1[&split={5173-9999},loci={5173-9999},length=4827]:0.09778803745774778)[&loci={0-9999},length=10000]:0.9978993277153256)#H0[&split={930-9999},loci={930-9999},length=9070]:0.8429584721860954)[&loci={0-9999},length=10000]:0.0;";

    

    @Test
    public void testAddRemoveSegment() {

        RecombinationNetwork network = new RecombinationNetwork(networkString);
        RecombinationNetworkNode leafNode = new ArrayList<>(network.getLeafNodes()).get(0);

        DivertLociOperator operator = new DivertLociOperator();
        operator.totalLength = 10000;
        
        List<Integer> listToAdd = new ArrayList<Integer>();
        BreakPoints lociToAdd = new BreakPoints();
        
        listToAdd.add(7000);
        listToAdd.add(7999);
        lociToAdd.init(listToAdd);
        
        
        Randomizer.setSeed(10);

        
        String net1 = network.toString();
        
        double logPremove = operator.removeLociFromAncestors(
              leafNode.getParentEdges().get(0), lociToAdd);
        
        Assert.assertEquals(logPremove,-6.908754779315, 1e-10); 
        
        
        double logPadd = operator.addLociToAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);
        
        Assert.assertEquals(logPadd,-6.908754779315, 1e-10); 

        /*
         * next one
         */
        listToAdd = new ArrayList<Integer>();
        lociToAdd = new BreakPoints();
        
        listToAdd.add(8000);
        listToAdd.add(9000);
        lociToAdd.init(listToAdd);
        
        logPremove = operator.removeLociFromAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);
        logPadd = operator.addLociToAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);
        
        Assert.assertEquals(logPadd,0.0, 1e-10); 
        Assert.assertEquals(logPremove,0.0, 1e-10); 
        
        
        /*
         * next one
         */
        listToAdd = new ArrayList<Integer>();
        lociToAdd = new BreakPoints();
        
        listToAdd.add(5000);
        listToAdd.add(9999);
        lociToAdd.init(listToAdd);
        
        
        logPremove = operator.removeLociFromAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);

        logPadd = operator.addLociToAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);
        Assert.assertEquals(logPadd,-8.517193191416, 1e-10); 
        Assert.assertEquals(logPremove,-8.517193191416, 1e-10); 


        

        /*
         * next one
         */
        listToAdd = new ArrayList<Integer>();
        lociToAdd = new BreakPoints();
        
        listToAdd.add(0);
        listToAdd.add(500);

        listToAdd.add(5000);
        listToAdd.add(9999);
        lociToAdd.init(listToAdd);
        

        logPremove = operator.removeLociFromAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);

        logPadd = operator.addLociToAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);
   
      
        
        Assert.assertEquals(logPadd,-9.210340371976, 1e-10); 
        Assert.assertEquals(logPremove,-9.210340371976, 1e-10); 
        
        
        /*
         * next one
         */
        listToAdd = new ArrayList<Integer>();
        lociToAdd = new BreakPoints();
        
        listToAdd.add(0);
        listToAdd.add(9999);
        lociToAdd.init(listToAdd); 
        
        System.out.println("..............");
        System.out.println("..............");
        System.out.println("..............");
        System.out.println("..............");
        System.out.println("..............");
        logPremove = operator.removeLociFromAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);
        System.out.println(network);
        Assert.assertEquals(logPremove,Math.log(0.5)+Math.log(1/10000.0)+Math.log(1/5173.0), 1e-10);     
        
        Randomizer.setSeed(10);

        System.out.println(lociToAdd);
        logPadd = operator.addLociToAncestors(
                leafNode.getParentEdges().get(0), lociToAdd);
        System.out.println(network);

        System.out.println(logPadd);
        System.out.println(Math.log(0.5)+Math.log(1/9999.0)+Math.log(1/5173.0));
        Assert.assertEquals(logPadd,Math.log(0.5)+Math.log(1/9999.0)+Math.log(1/5173.0), 1e-10); 


        
        
//        BreakPoints bp = operator.getRangeToDivert(leafNode.getParentEdges().get(0).parentNode.getParentEdges().get(0),
//        		leafNode.getParentEdges().get(0).parentNode.getParentEdges().get(1));
//       
//        BreakPoints lociToDivert = bp.copy();
//    	lociToDivert.and(leafNode.getParentEdges().get(0).parentNode.getParentEdges().get(0).breakPoints);
//
//        
//        operator.addLociToAncestors(leafNode.getParentEdges().get(0).parentNode.getParentEdges().get(1), lociToDivert);
//        operator.removeLociFromAncestors(leafNode.getParentEdges().get(0).parentNode.getParentEdges().get(0), lociToDivert);
//        
//        leafNode.getParentEdges().get(0).parentNode.getParentEdges().get(0).passingRange.andNot(bp);
//
//        String controlNet = "((((((#H2[&split={9240-9999},loci={9240-9999},length=760]:0.3579693014979668,(#H3[&split={6573-9999},loci={7302-7756},length=455]:0.4198547464044533,#H4[&split={0-7717},loci={6648-7717},length=1070]:0.25563256204649143)[&loci={6648-7756},length=1109]:1.1619129921459237)[&loci={6648-7756,9240-9999},length=1869]:0.34286205275424475)#H1[&split={7372-9999},loci={7372-7756,9240-9999},length=1145]:0.23453108442206716,((((((#H7[&split={3116-9999},loci={3116-6068},length=2953]:0.04387002757521419,((#H8[&split={0-990},loci={0-902},length=903]:0.20086820188684906,((((#H10[&split={6130-9999},loci={6130-7425},length=1296]:0.021261671248240965,(((#H12[&split={0-4308},loci={1942-4308},length=2367]:0.017015062714915707)#H11[&split={0-3255},loci={1942-3255},length=1314]:0.09607673081432422,(((((((t7[&loci={0-9999},length=10000]:0.1967853484864821,(t3[&loci={0-9999},length=10000]:0.053221593679751766,t6[&loci={0-9999},length=10000]:0.053221593679751766)[&loci={0-9999},length=10000]:0.14356375480673034)[&loci={0-9999},length=10000]:0.051513614445772404)#H16[&split={2275-9999},loci={2275-9999},length=7725]:0.07715170782013026)#H15[&split={5532-9999},loci={5532-9999},length=4468]:0.018553890999346123,(((t4[&loci={0-9999},length=10000]:0.17649068406003154)#H17[&split={0-7301},loci={0-7301},length=7302]:0.04039091402996009,(#H18[&split={0-1941},loci={0-1941},length=1942]:0.0013853304961970458,(t9[&loci={0-9999},length=10000]:0.18579011385075272)#H19[&split={0-7244},loci={0-7244},length=7245]:0.01577528181297172)[&loci={0-7244},length=7245]:0.01531620242626719)[&loci={0-7301},length=7302]:0.10731442948560099,#H20[&split={8739-9999},loci={8739-9999},length=1261]:0.16478852729506777)[&loci={0-7301,8739-9999},length=8563]:0.019808534176138703)[&loci={0-9999},length=10000]:0.0666631589084492,#H21[&split={8920-9999},loci={8920-9999},length=1080]:0.017738595522331035)[&loci={0-9999},length=10000]:0.014016778432702548)#H14[&split={0-7425},loci={0-7425},length=7426]:0.028992203555865625)#H13[&split={2098-9999},loci={2098-7425},length=5328]:0.07487331027241328)[&loci={1942-7425},length=5484]:0.32310917034218045)#H10[&split={0-6129},loci={1942-6129},length=4188]:0.021261671248240965)[&loci={1942-7425},length=5484]:0.033310006810186366,#H22[&split={0-6338},loci={903-2097,3256-4308},length=2248]:0.04657187850738387)[&loci={903-7425},length=6523]:0.13453814563512445,((((((#H19[&split={7245-9999},loci={7245-9999},length=2755]:0.014671122587377905,(#H23[&split={6330-9999},loci={6330-9999},length=3670]:0.09734477025221899)#H20[&split={0-8738},loci={6330-8738},length=2409]:0.04105373615760577)[&loci={6330-9999},length=3670]:0.2719204176664496,#H14[&split={7426-9999},loci={7426-9999},length=2574]:0.047697155011697134)[&loci={6330-9999},length=3670]:0.19505659579144785,(((((((t5[&loci={0-9999},length=10000]:0.04149497715939221,t1[&loci={0-9999},length=10000]:0.04149497715939221)[&loci={0-9999},length=10000]:0.023470009401864278,(#H26[&split={0-4446},loci={0-4446},length=4447]:0.005486880263416705)#H25[&split={1773-9999},loci={1773-4446},length=2674]:0.053501093010555945)[&loci={0-9999},length=10000]:0.11378162648069745,(t8[&loci={0-9999},length=10000]:0.06206273002830587)#H23[&split={0-6329},loci={0-6329},length=6330]:0.11668388301364852)[&loci={0-9999},length=10000]:0.18840329381588017,#H16[&split={0-2274},loci={0-2274},length=2275]:0.1188509439255796)[&loci={0-9999},length=10000]:0.02577921828001495)#H21[&split={0-8919},loci={0-8919},length=8920]:0.2317721230555394,((((t0[&loci={0-9999},length=10000]:0.2001800651675274)#H18[&split={1942-9999},loci={1942-9999},length=8058]:0.1568564789210014)#H28[&split={0-9487},loci={1942-9487},length=7546]:0.05842167530339326)#H12[&split={4309-9999},loci={4309-9487},length=5179]:0.14856286090476134)#H27[&split={5095-9999},loci={5095-9487},length=4393]:0.060680167896705495)[&loci={0-9487},length=9488]:0.024804534978506876)#H24[&split={0-6647},loci={0-6647},length=6648]:0.01793246672413229)[&loci={0-9999},length=10000]:0.037937854132429116,#H29[&split={9676-9999},loci={9676-9999},length=324]:0.012880858126575134)[&loci={0-9999},length=10000]:0.0021093878245594944,(((#H25[&split={0-1772},loci={0-1772},length=1773]:0.1665090691103175,(t2[&loci={0-9999},length=10000]:0.005977013287283839)#H26[&split={4447-9999},loci={4447-9999},length=5553]:0.1719959493737342)[&loci={0-1772,4447-9999},length=7326]:0.37566893634034537)#H31[&split={0-8695},loci={0-1772,4447-8695},length=6022]:0.15004789222410242)#H30[&split={0-8483},loci={0-1772,4447-8483},length=5810]:0.0037957006275504135)[&loci={0-9999},length=10000]:0.017935772521958437,#H30[&split={8484-9999},loci={8484-8695},length=212]:0.02173147314950885)[&loci={0-9999},length=10000]:0.3153477425819191)[&loci={0-9999},length=10000]:0.13327931242701974)#H9[&split={0-6068},loci={0-6068},length=6069]:0.5322482481945647)[&loci={0-6068},length=6069]:0.480143525875959)#H7[&split={0-3115},loci={0-3115},length=3116]:0.04387002757521419)[&loci={0-6068},length=6069]:0.379855344642841,(((((((((#H33[&split={0-902},loci={0-902},length=903]:0.30734918020033053,#H34[&split={6993-9999},loci={7302-7756},length=455]:0.1948828154439055)[&loci={0-902,7302-7756},length=1358]:0.06246296606424773)#H3[&split={0-6572},loci={0-902},length=903]:0.16521924978447844,((((#H37[&split={4696-9999},loci={9488-9999},length=512]:0.0335526741633525,((#H27[&split={0-5094},loci={4309-5094},length=786]:0.12142467648671396,#H31[&split={8696-9999},loci={8696-9999},length=1304]:0.13180385778203352)[&loci={4309-5094,8696-9999},length=2090]:0.03215337284480002,((((#H15[&split={0-5531},loci={2275-5531},length=3257]:0.050861127364038694,#H17[&split={7302-9999},loci={7302-9999},length=2698]:0.19982111405639236)[&loci={2275-5531,7302-9999},length=5955]:0.18444917909635272)#H38[&split={0-7314},loci={2275-5531,7302-7314},length=3270]:0.11078295782418879,#H38[&split={7315-9999},loci={7315-9999},length=2685]:0.11078295782418879)[&loci={2275-5531,7302-9999},length=5955]:0.02095131086491664)#H29[&split={0-9675},loci={2275-5531,7302-9675},length=5631]:0.025103883726315335)[&loci={2275-5531,7302-9999},length=5955]:0.039060252127992356)[&loci={2275-5531,7302-9999},length=5955]:0.16350514165220442)#H36[&split={0-7756},loci={2275-5531,7302-7756},length=3712]:0.014311106406990959)#H34[&split={0-6992},loci={2275-5531},length=3257]:0.0074120040152894084)#H35[&split={3294-9999},loci={3294-5531},length=2238]:0.41515302727734227)[&loci={0-902,3294-5531},length=3141]:0.04680788491294807,((#H36[&split={7757-9999},loci={7757-9999},length=2243]:0.1231585867741507,((((((#H28[&split={9488-9999},loci={9488-9999},length=512]:0.11833012494094675,#H11[&split={3256-9999},loci={3256-4308},length=1053]:0.04289338692263778)[&loci={3256-4308,9488-9999},length=1565]:0.046333283917160184,#H13[&split={0-2097},loci={0-2097},length=2098]:0.06802325029788703)[&loci={0-2097,3256-4308,9488-9999},length=3663]:0.2014067546462015)#H37[&split={0-4695},loci={0-2097,3256-4308},length=3151]:0.0651053177459966,#H24[&split={6648-9999},loci={6648-9487},length=2840]:0.13870624216693805)[&loci={0-2097,3256-4308,6648-9487},length=5991]:0.03379723972012627)#H33[&split={903-9999},loci={903-2097,3256-4308,6648-9487},length=5088]:0.037649717755425804)#H22[&split={6339-9999},loci={6648-9487},length=2840]:0.18366412736815896)[&loci={6648-9999},length=3352]:0.31272048549895537)#H4[&split={7718-9999},loci={7718-9999},length=2282]:0.047804950339464636)[&loci={0-902,3294-5531,7718-9999},length=5423]:0.10157981967066476)#H8[&split={991-9999},loci={3294-5531,7718-9999},length=4520]:0.11373471539674496,#H9[&split={6069-9999},loci={6069-9999},length=3931]:0.4451147617044606)[&loci={3294-5531,6069-9999},length=6169]:0.012728616236560297)#H32[&split={6233-9999},loci={6233-9999},length=3767]:0.5280998221021458,(#H32[&split={0-6232},loci={3294-5531,6069-6232},length=2402]:0.14509037171222605,#H35[&split={0-3293},loci={2275-3293},length=1019]:0.8350944352064864)[&loci={2275-5531,6069-6232},length=3421]:0.38300945038991974)[&loci={2275-5531,6069-9999},length=7188]:0.2556283289488679)#H2[&split={0-9239},loci={2275-5531,6069-9239},length=6428]:0.1945456172965443)[&loci={0-9239},length=9240]:0.3456122897294751)#H6[&split={0-8604},loci={0-8604},length=8605]:0.06276462209584555)#H5[&split={0-7661},loci={0-7661},length=7662]:0.03841238216566989,#H6[&split={8605-9999},loci={8605-9239},length=635]:0.10117700426151544)[&loci={0-7661,8605-9239},length=8297]:0.1841877682963844,#H1[&split={0-7371},loci={6648-7371},length=724]:0.12469132533170768)[&loci={0-7661,8605-9239},length=8297]:0.10983975909035948)[&loci={0-7756,8605-9999},length=9152]:0.015448252837242826)#H0[&split={0-6768},loci={0-6768},length=6769]:0.1191130035471124,#H0[&split={6769-9999},loci={6769-7756,8605-9999},length=2383]:0.1191130035471124)[&loci={0-7756,8605-9999},length=9152]:0.13688107726724486,#H5[&split={7662-9999},loci={7662-8604},length=943]:0.6038822432040138)[&loci={0-9999},length=10000]:0.0;";
//        Assert.assertEquals(network.toString(), controlNet);
//		
//

    }

    @Test
    public void testAddRemoveReassortmentEdge() {
        RecombinationNetwork network = new RecombinationNetwork(networkString);

        AddRemoveRecombination operator = new AddRemoveRecombination();
        operator.totalLength=10000;
        operator.network = network;
        network.totalLength=10000;

        Randomizer.setSeed(1);
        
        RecombinationNetworkNode origRoot = network.getRootEdge().childNode;

        RecombinationNetworkEdge sourceEdge = network.getRootEdge().childNode.getChildEdges().get(0);
        double sourceTime = sourceEdge.getLength()/2.0 + sourceEdge.childNode.getHeight();
        RecombinationNetworkEdge destEdge = network.getRootEdge();
        double destTime = destEdge.childNode.getHeight() + 1.0;

        double logP1 = operator.addRecombinationEdge(sourceEdge, sourceTime, destEdge, destTime);
        Assert.assertEquals(logP1,-9.244161920452, 1e-5); 
        
        double logP2 = operator.removeRecombinationEdge(network.getRootEdge().childNode.getChildEdges().get(1));
        Assert.assertEquals(logP2,9.244161920452, 1e-5); 

        
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
    }
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
