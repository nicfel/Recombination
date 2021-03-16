package recombination.operators;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

public class GaussianRecombinationNetworkScaler extends RecombinationNetworkOperator {

    final public Input<Double> sizeInput = new Input<>("size", "size of the slide, default 1.0", 1.0);
    final public Input<Boolean> gaussianInput = new Input<>("gaussian", "Gaussian (=true=default) or uniform delta", true);
    final public Input<Double> limitInput = new Input<>("limit", "limit on step size, default disable, " +
            "i.e. -1. (when positive, gets multiplied by network-height/log2(n-taxa).", -1.0);
    
    
    public Input<List<RealParameter>> upParametersInput = new Input<>(
            "upParameter",
            "Parameters to scale in the SAME direction as the network.",
            new ArrayList<>());

    public Input<List<RealParameter>> downParametersInput = new Input<>(
            "downParameter",
            "Parameters to scale in the OPPOSITE direction as the network.",
            new ArrayList<>());

    
    final public Input<Boolean> optimiseInput = new Input<>(
    		"optimise", 
    		"flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", 
    		true);

    boolean scaleRootOnly;
    List<RealParameter> upParameters, downParameters;

    // shadows size
    double size;
    private double limit;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        size = sizeInput.get();
        limit = limitInput.get();
        upParameters = upParametersInput.get();
        downParameters = downParametersInput.get();

    }

    @Override
    public double networkProposal() {

        int count = 0;

        final double delta = getDelta();

        network.startEditing(this);

        if (scaleRootOnly) {

            // Scale root

        	RecombinationNetworkNode rootNode = network.getRootEdge().childNode;

            rootNode.setHeightFilty(rootNode.getHeight() + delta);
            count += 1;

            if (delta<0) {

                for (RecombinationNetworkEdge childEdge : rootNode.getChildEdges())
                    if (rootNode.getHeight() < childEdge.childNode.getHeight())
                        return Double.NEGATIVE_INFINITY;

            }

        } else {

            // Scale network nodes

            for (RecombinationNetworkNode node : network.getInternalNodes()) {
                node.setHeightFilty(node.getHeight() + delta);
                count += 1;
            }
            if (delta<0) {
	            for (RecombinationNetworkNode leaf : network.getLeafNodes()) {
	                if (leaf.getParentEdges().get(0).parentNode.getHeight() < leaf.getHeight())
	                    return Double.NEGATIVE_INFINITY;
	            }
            }
            

        }

        // Scale parameters

        try {

            for (RealParameter param : upParameters) {
                param.startEditing(this);
                count += param.scale(Math.exp(delta));
            }

            for (RealParameter param : downParameters) {
                param.startEditing(this);
                count -= param.scale(Math.exp(-delta));
            }

        } catch (IllegalArgumentException ex) {
            return Double.NEGATIVE_INFINITY;
        }

        return 0;
    }
    

    private double getDelta() {
        if (!gaussianInput.get()) {
            return (Randomizer.nextDouble() * size) - (size / 2.0);
        } else {
            return Randomizer.nextGaussian() * size;
        }
    }

    
    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(size);
            final double f = Math.exp(delta);
            if( limit > 0 ) {
                final RecombinationNetwork network = networkInput.get();
                final double h = network.getRootEdge().childNode.getHeight();
                final double k = Math.log(network.getLeafNodes().size()) / Math.log(2);
                final double lim = (h / k) * limit;
                if( f <= lim ) {
                    size = f;
                }
            } else {
               size = f;
            }
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return size;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        size = value;
    }
    
    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;

        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        final double newDelta = size * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try decreasing size to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try increasing size to about " + formatter.format(newDelta);
        } else return "";
    }
}
