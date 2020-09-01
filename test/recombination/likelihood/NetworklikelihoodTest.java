package recombination.likelihood;

import org.junit.Test;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.HKY;
import beast.evolution.substitutionmodel.JukesCantor;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.TestCase;
import recombination.alignment.RecombinationAlignment;
import recombination.network.RecombinationNetwork;
import test.beast.BEASTTestCase;

public class NetworklikelihoodTest extends TestCase {

    public NetworklikelihoodTest() {
        super();
    }

    protected NetworkLikelihood newNetworkLikelihood() {
    	System.setProperty("java.only","true");
        return new NetworkLikelihood();
    }
    
    protected TreeLikelihood newTreeLikelihood() {
    	System.setProperty("java.only","true");
        return new TreeLikelihood();
    }


    @Test
    public void testJC69Likelihood() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
    	RecombinationAlignment data = getAlignmentShort();
        RecombinationNetwork network = getNetworkShort();
                       
        
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        NetworkLikelihood likelihood = newNetworkLikelihood();
        likelihood.initByName("recombinationData", data, "recombinationNetwork", network, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        
        
        // compute the tree likelihoods for each positions individually
        double treeLog = 0.0;
        for (int i = 0; i < 4; i++) {
            Alignment data_pos = getAlignmentPosition(i);
            Node root = network.getLocusChildren(network.getRootEdge().childNode, i);
            TreeParser t = new TreeParser();

            t.initByName("taxa", data_pos,
                    "newick", root.toNewick(false),
                    "IsLabelledNewick", true);            
           
            TreeLikelihood treelikelihood = newTreeLikelihood();
            treelikelihood.initByName("data", data_pos, "tree", t, "siteModel", siteModel);
            treeLog += treelikelihood.calculateLogP();
        }
        
        
        
        assertEquals(logP, treeLog, BEASTTestCase.PRECISION);
        assertEquals(logP, likelihood.calculateLogP(), BEASTTestCase.PRECISION);
    }
    
    @Test
    public void testHKYLikelihood() throws Exception {
        // Set up JC69 model: uniform freqs, kappa = 1, 0 gamma categories
    	RecombinationAlignment data = getAlignmentShort();
        RecombinationNetwork network = getNetworkShort();
                       
    	Alignment data_freqs = getNormalAlignment();
       
        Frequencies freqs = new Frequencies();
        freqs.initByName("data", data_freqs);

        
        HKY hky = new HKY();
        hky.initByName("kappa", "29.739445", "frequencies", freqs);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", hky);

        NetworkLikelihood likelihood = newNetworkLikelihood();
        likelihood.initByName("recombinationData", data, "recombinationNetwork", network, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        
        
        // compute the tree likelihoods for each positions individually
        double treeLog = 0.0;
        for (int i = 0; i < 4; i++) {
            Alignment data_pos = getAlignmentPosition(i);
            Node root = network.getLocusChildren(network.getRootEdge().childNode, i);
            TreeParser t = new TreeParser();

            t.initByName("taxa", data_pos,
                    "newick", root.toNewick(false),
                    "IsLabelledNewick", true);            
           
            TreeLikelihood treelikelihood = newTreeLikelihood();
            treelikelihood.initByName("data", data_pos, "tree", t, "siteModel", siteModel);
            treeLog += treelikelihood.calculateLogP();
        }    
        
        
        assertEquals(logP, treeLog, BEASTTestCase.PRECISION);
    }

    
    
    static public RecombinationAlignment getAlignment() throws Exception {
        Sequence human = new Sequence("human", "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCCCCTTATTTCTACTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTG-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCGGAGCTTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCT-AAATCCTATATATCTTA-CACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGA");
        Sequence chimp = new Sequence("chimp", "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTCAAATCCCCTTATTTCTACTAGGACTATAAGAATCGAACTCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTAAGCACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCTCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AACCCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCGACACCTCTTTACAGTGA");
        Sequence bonobo = new Sequence("bonobo", "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTTAAATCCCCTTATTTCTACTAGGACTATGAGAGTCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAACACAGACCAAGAGCCTTCAAAGCTCTCAGTAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAGC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTTGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCCCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCAACACCTCTTTACAGTGA");
        Sequence gorilla = new Sequence("gorilla", "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGAGGTTTAAACCCCCTTATTTCTACTAGGACTATGAGAATTGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTGTCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTCACATCCTTCCCGTACTAAGAAATTTAGGTTAAACATAGACCAAGAGCCTTCAAAGCCCTTAGTAAGTTA-CAACACTTAATTTCTGTAAGGACTGCAAAACCCTACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCAATGGGACTCAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAT-TCACCTCGGAGCTTGGTAAAAAGAGGCCCAGCCTCTGTCTTTAGATTTACAGTCCAATGCCTTA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGACCTTCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAGGTTAAATTACGGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCGTTAACCTTTTAAGTTAAAGATTAAGAGTATCGGCACCTCTTTGCAGTGA");
        Sequence orangutan = new Sequence("orangutan", "AGAAATATGTCTGACAAAAGAGTTACTTTGATAGAGTAAAAAATAGAGGTCTAAATCCCCTTATTTCTACTAGGACTATGGGAATTGAACCCACCCCTGAGAATCCAAAATTCTCCGTGCCACCCATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTA--CACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTCA-CAGCACTTAATTTCTGTAAGGACTGCAAAACCCCACTTTGCATCAACTGAGCGCAAATCAGCCACTTTAATTAAGCTAAGCCCTCCTAGACCGATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAT-TGGCTTCAGTCCAAAGCCCCGGCAGGCCTTAAAGCTGCTCCTTCGAATTTGCAATTCAACATGACAA-TCACCTCAGGGCTTGGTAAAAAGAGGTCTGACCCCTGTTCTTAGATTTACAGCCTAATGCCTTAACTCGGCCATTTTACCGCAAAAAAGGAAGGAATCGAACCTCCTAAAGCTGGTTTCAAGCCAACCCCATAACCCCCATGACTTTTTCAAAAGGTACTAGAAAAACCATTTCGTAACTTTGTCAAAGTTAAATTACAGGTC-AGACCCTGTGTATCTTA-CATTGCAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGACTAAGAGAACCAGCCTCTCTTTGCAATGA");
        Sequence siamang = new Sequence("siamang", "AGAAATACGTCTGACGAAAGAGTTACTTTGATAGAGTAAATAACAGGGGTTTAAATCCCCTTATTTCTACTAGAACCATAGGAGTCGAACCCATCCTTGAGAATCCAAAACTCTCCGTGCCACCCGTCGCACCCTGTTCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCATACTAAGAAATTTAGGTTAAACACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTAACAAAACTTAATTTCTGCAAGGGCTGCAAAACCCTACTTTGCATCAACCGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCGATGGGACTTAAACCCATAAAAATTTAGTTAACAGCTAAACACCCTAAACAACCTGGCTTCAATCTAAAGCCCCGGCAGA-GTTGAAGCTGCTTCTTTGAACTTGCAATTCAACGTGAAAAATCACTTCGGAGCTTGGCAAAAAGAGGTTTCACCTCTGTCCTTAGATTTACAGTCTAATGCTTTA-CTCAGCCACTTTACCACAAAAAAGGAAGGAATCGAACCCTCTAAAACCGGTTTCAAGCCAGCCCCATAACCTTTATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATCACAGGTCCAAACCCCGTATATCTTATCACTGTAGAGCTAGACCAGCATTAACCTTTTAAGTTAAAGACTAAGAGAACTACCGCCTCTTTACAGTGA");

        RecombinationAlignment data = new RecombinationAlignment();
        data.initByName("sequence", human, "sequence", chimp, "sequence", bonobo, "sequence", gorilla, "sequence", orangutan, "sequence", siamang,
                "dataType", "nucleotide"
        );
        return data;
    }
    
    static public RecombinationAlignment getAlignmentShort() throws Exception {
        Sequence t0 = new Sequence("t0", "AGAA");
        Sequence t1 = new Sequence("t1", "AGAT");
        Sequence t2 = new Sequence("t2", "AGAA");

        RecombinationAlignment data = new RecombinationAlignment();
        data.initByName("sequence", t0, "sequence", t1, "sequence", t2,
                "dataType", "nucleotide"
        );
        return data;
    }
    
    static public Alignment getAlignmentPosition(int i) throws Exception {
    	String t0_st = "AGAA";
    	String t1_st = "AGAT";
    	String t2_st = "AGAA";
        Sequence t0 = new Sequence("t0", t0_st.substring(i, i+1));
        Sequence t1 = new Sequence("t1", t1_st.substring(i, i+1));
        Sequence t2 = new Sequence("t2", t2_st.substring(i, i+1));

        Alignment data = new Alignment();
        data.initByName("sequence", t0, "sequence", t1, "sequence", t2,
                "dataType", "nucleotide"
        );
        return data;
    }
    
    static public Alignment getNormalAlignment() throws Exception {
        Sequence t0 = new Sequence("t0", "AGAA");
        Sequence t1 = new Sequence("t1", "AGAT");
        Sequence t2 = new Sequence("t2", "AGAA");

        Alignment data = new Alignment();
        data.initByName("sequence", t0, "sequence", t1, "sequence", t2,
                "dataType", "nucleotide"
        );
        return data;
    }


    
    static public RecombinationNetwork getNetworkShort() throws Exception {
        RecombinationNetwork network = new RecombinationNetwork(
        		"((#H0[&split={0-0},loci={0-0},length=4]:0.4860702162314561,((#H2[&split={1-3},loci={1-3},length=3]:0.036380107508342974,(t1[&loci={0-3},length=4]:0.29041085418573037,t0[&loci={0-3},length=4]:0.29041085418573037)[&loci={0-3},length=4]:0.1528435079144485)[&loci={0-3},length=4]:0.47005038739308824)#H1[&split={0-0},loci={0-0},length=1]:1.5817575814045295)[&loci={0-0},length=1]:0.35688825595463936,(((t2[&loci={0-3},length=4]:0.4068742545918359)#H2[&split={0-0},loci={0-0},length=1]:0.604218532359179,#H1[&split={1-3},loci={1-3},length=3]:0.09778803745774778)[&loci={0-3},length=4]:0.9978993277153256)#H0[&split={1-3},loci={1-3},length=3]:0.8429584721860954)[&loci={0-3},length=4]:0.0;");
        return network;
    }



}
