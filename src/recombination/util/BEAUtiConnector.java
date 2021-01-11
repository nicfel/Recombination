package recombination.util;

import beast.app.beauti.BeautiDoc;
import beast.core.BEASTInterface;
import beast.core.BEASTObject;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.operators.UpDownOperator;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import coalre.distribution.CoalescentWithReassortment;
import coalre.operators.NetworkScaleOperator;
import coalre.simulator.SimulatedCoalescentNetwork;
import recombination.likelihood.NetworkLikelihood;
import recombination.simulator.SimulatedCoalescentRecombinationNetwork;

import java.util.*;

/**
 * Class containing a static method used in the BEAUti template to tidy
 * up some loose ends.
 */
public class BEAUtiConnector {

    public static boolean customConnector(BeautiDoc doc) {

        MCMC mcmc = (MCMC) doc.mcmc.get();
        
        // replace the tree likelihood with a Network likelihood
        for (Distribution distr : ((CompoundDistribution) mcmc.posteriorInput.get()).pDistributions.get()) {
	   		 boolean hasNetworkLikelihood = false;
			 for (Distribution likelihood : ((CompoundDistribution) distr).pDistributions.get()) {
				 if(likelihood instanceof NetworkLikelihood){
					 hasNetworkLikelihood = true;
				 }
			 }
			 if (!hasNetworkLikelihood) {
	        	if (distr.getID().contentEquals("likelihood")) {
	        		int i = 0;
	        		 for (Distribution likelihood : ((CompoundDistribution) distr).pDistributions.get()) {
		                SimulatedCoalescentRecombinationNetwork network = (SimulatedCoalescentRecombinationNetwork) doc.pluginmap.get("network");
		                if (!hasNetworkLikelihood) {
			                try{
				                NetworkLikelihood nlk = new NetworkLikelihood();
				                nlk.initByName("data", ((GenericTreeLikelihood) likelihood).dataInput.get(),
				                		"siteModel", ((GenericTreeLikelihood) likelihood).siteModelInput.get(), 
				                		"branchRateModel",((GenericTreeLikelihood) likelihood).branchRateModelInput.get(),
				                		"recombinationNetwork", network,
				                		"tree", ((GenericTreeLikelihood) likelihood).treeInput.get());
				                nlk.setID("NetworkLikelihood");
				                ((CompoundDistribution) distr).pDistributions.get().add(nlk);
				                ((CompoundDistribution) distr).pDistributions.get().remove(i);
				                
			                }catch (Exception e){
			                	System.out.println(e);
			                }
			                i++;
			                hasNetworkLikelihood = true;
			                break;
		                }
	        		 }
        		 }

        	}            	
        }

        List<Logger> removeLogger = new ArrayList<>();
        for (Logger log : mcmc.loggersInput.get()) {
        	if (log.getID().contains("treelog")) {
        		removeLogger.add(log);
        	}
        }
        for (Logger log : removeLogger) {
        	mcmc.loggersInput.get().remove(mcmc.loggersInput.get().indexOf(log));
        }



//        for (BEASTInterface p : doc.getPartitions("tree")) {
//            String pId = BeautiDoc.parsePartition(p.getID());
//
//            DummyTreeDistribution dummy = (DummyTreeDistribution)doc.pluginmap.get("CoalescentWithReassortmentDummy.t:" + pId);
//
//            if (dummy == null || !dummy.getOutputs().contains(doc.pluginmap.get("prior")))
//                continue;

//            GenericTreeLikelihood likelihood = (GenericTreeLikelihood) p;
////            likelihood.
//            NetworkLikelihood nlk = (NetworkLikelihood) likelihood;
//            System.out.println(likelihood);
//            System.out.println(nlk);
            
            

            // Remove segment trees from standard up/down operators.

//            for (Operator operator : mcmc.operatorsInput.get()) {
//                if (!(operator instanceof UpDownOperator))
//                    continue;
//
//                UpDownOperator upDown = (UpDownOperator) operator;
//
//                boolean segmentTreeScaler = upDown.upInput.get().contains(segmentTree)
//                        || upDown.downInput.get().contains(segmentTree);
//
//                if (segmentTreeScaler) {
//                    upDown.upInput.get().remove(segmentTree);
//                    upDown.downInput.get().remove(segmentTree);
//                }
//
//            }
//
//            // Clock rates to add to network up/down:
//            for (Operator operator : mcmc.operatorsInput.get()) {
//                if (!(operator instanceof UpDownOperator))
//                    continue;
//
//                UpDownOperator upDown = (UpDownOperator) operator;
//
//                // Note: built-in up/down operators scale trees _down_ while
//                // ours scales trees _up_, hence the up/down reversal.
//                for (BEASTObject o : upDown.upInput.get()) {
//                    if (o instanceof RealParameter) {
//                        if (o.getID().contains("clock"))
//                            parametersToScaleDown.add((RealParameter)o);
//                    }
//                }
//
//                for (BEASTObject o : upDown.downInput.get()) {
//                    if (o instanceof RealParameter) {
//                        if (o.getID().contains("clock"))
//                            parametersToScaleUp.add((RealParameter)o);
//                    }
//                }
//            }
//
//            // Extract trait set from one of the trees to use for network.
//
//            if (traitSet == null && segmentTree.hasDateTrait())
//                traitSet = segmentTree.getDateTrait();
//        }


        // Add clock rates to network up/down operator.

//        NetworkScaleOperator networkUpDown = (NetworkScaleOperator)doc.pluginmap.get("networkUpDownCwR.alltrees");
//        if (networkUpDown != null) {
//            List<RealParameter> paramsCurrent;
//
//            paramsCurrent = networkUpDown.upParametersInput.get();
//            for (RealParameter param : paramsCurrent) {
//                if (param.getID().toLowerCase().contains("clock"))
//                    networkUpDown.upParametersInput.get().remove(param);
//            }
//
//            paramsCurrent = networkUpDown.downParametersInput.get();
//            for (RealParameter param : paramsCurrent) {
//                if (param.getID().toLowerCase().contains("clock"))
//                    networkUpDown.downParametersInput.get().remove(param);
//            }
//
//            networkUpDown.upParametersInput.get().addAll(parametersToScaleUp);
//            networkUpDown.downParametersInput.get().addAll(parametersToScaleDown);
//        }

        // Update network initializer:

//        if (doc.pluginmap.containsKey("networkCwR.alltrees")) {
//            SimulatedCoalescentNetwork network = (SimulatedCoalescentNetwork) doc.pluginmap.get("networkCwR.alltrees");
//
//            // Update number of segments for initializer.
//            network.nSegmentsInput.setValue(segTreeCount, network);
//
//            // Provide trait set from first segment tree to network initializer:
//            if (traitSet != null)
//                network.traitSetInput.setValue(traitSet, network);
//        }



        return false;
    }
}
