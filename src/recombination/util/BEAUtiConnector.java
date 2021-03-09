package recombination.util;

import beast.app.beauti.BeautiDoc;
import beast.core.Distribution;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.util.CompoundDistribution;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.operators.UpDownOperator;
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

        // remove up down operator
        for (Operator operator : mcmc.operatorsInput.get()) {
            if (!(operator instanceof UpDownOperator))
                continue;

            UpDownOperator upDown = (UpDownOperator) operator;
            mcmc.operatorsInput.get().remove(mcmc.operatorsInput.get().indexOf(operator));
        }
        return false;
    }
}
