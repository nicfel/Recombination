package recombination.distribution;



import beast.base.inference.CalculationNode;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Nicola Felix Mueller
 */
public class RecombinationNetworkIntervals extends CalculationNode {
    public Input<RecombinationNetwork> recombinationNetworkInput = new Input<>("recombinationNetwork",
            "recombinationNetwork for which to calculate the intervals", Validate.REQUIRED);
    
	public Input<String> recombinationRatesChangePointsInput = new Input<>(
	        "recombinationRatesChangePoints",
            "if true, only coalescent events are allowed after the .",
            Input.Validate.OPTIONAL);


    private RecombinationNetwork recombinationNetwork;

    private List<RecombinationNetworkEvent> recombinationNetworkEventList, storedRecombinationNetworkEventList;

    public boolean eventListDirty = true;    
    public boolean hasBreakPoints = false;
    
    public BreakPoints[] recBP;


    
    @Override
    public void initAndValidate() {
        recombinationNetwork = recombinationNetworkInput.get();
        storedRecombinationNetworkEventList = new ArrayList<>();
        
        
        if (recombinationRatesChangePointsInput.get()!=null) {
        	
        	hasBreakPoints = true;
        	String[] tmp = recombinationRatesChangePointsInput.get().trim().split("\\s+");
        	recBP = new BreakPoints[tmp.length+1];
        	List<Integer> bpList = new ArrayList<>();
        	bpList.add(0);
        	for (int i = 0; i < tmp.length; i++) {            	
            	bpList.add(Integer.parseInt(tmp[i]));
            	bpList.add(Integer.parseInt(tmp[i]));
        		if (Integer.parseInt(tmp[i])>recombinationNetworkInput.get().totalLength) {
            		throw new IllegalArgumentException("recombination rate break point value is larger than the network size");
        		}
        	}   
        	bpList.add(recombinationNetworkInput.get().totalLength-1);
        	
        	for (int i = 0; i < recBP.length; i++) 
        		recBP[i] = new BreakPoints(bpList.get(2*i), bpList.get(2*i+1));
        }else {
        	recBP = new BreakPoints[1];
    		recBP[0] = new BreakPoints(0, recombinationNetworkInput.get().totalLength-1);
        }
    }

    public List<RecombinationNetworkEvent> getRecombinationNetworkEventList() {
        update();

        return recombinationNetworkEventList;
    }

    void update() {
//        if (!eventListDirty)
//            return;
        
        recombinationNetworkEventList = recombinationNetwork.getNodes().stream().map(n -> {
            RecombinationNetworkEvent event = new RecombinationNetworkEvent();
            event.time = n.getHeight();
            event.node = n;
            switch(n.getChildCount()) {
                case 0:
                    event.type = RecombinationNetworkEvent.RecombinationNetworkEventType.SAMPLE;
                    break;

                case 1:
                    event.type = RecombinationNetworkEvent.RecombinationNetworkEventType.RECOMBINATION;
                    break;

                case 2:
                    event.type = RecombinationNetworkEvent.RecombinationNetworkEventType.COALESCENCE;
                    break;

                default:
                    throw new RuntimeException("RecombinationNetwork node has illegal number of children.");
            }
            return event;
        }).sorted(Comparator.comparingDouble(e -> e.time)).collect(Collectors.toList());

        int lineages = 0;
        double[] totalReassortmentObsProb = new double[recBP.length];

        for (RecombinationNetworkEvent event : recombinationNetworkEventList) {
            switch(event.type) {
                case SAMPLE:
                    lineages += 1;
//                    System.out.println(Arrays.toString(recBP));

                    for (int i = 0; i < recBP.length;i++) {
                    	totalReassortmentObsProb[i] += event.node.getParentEdges().get(0).breakPoints.getNullLength(recBP[i]);
                    }
//                    System.out.println(Arrays.toString(totalReassortmentObsProb));
//                    System.exit(0);
                    break;

                case RECOMBINATION:
                    lineages += 1;

                    for (int i = 0; i < recBP.length;i++) {
                    	totalReassortmentObsProb[i] -= event.node.getChildEdges().get(0).breakPoints.getNullLength(recBP[i]);
                    	totalReassortmentObsProb[i] += event.node.getParentEdges().get(0).breakPoints.getNullLength(recBP[i]);
                    	totalReassortmentObsProb[i] += event.node.getParentEdges().get(1).breakPoints.getNullLength(recBP[i]);
                    }

                    event.lociToSort = event.node.getChildEdges().get(0).getRecombinationLength();
                    
                    event.breakPoint = Math.max(event.node.getParentEdges().get(0).passingRange.getMin(), 
                    		event.node.getParentEdges().get(1).passingRange.getMin());
                    
                    break;

                case COALESCENCE:
                    lineages -= 1;
                    for (int i = 0; i < recBP.length;i++) {
                    	totalReassortmentObsProb[i] -= event.node.getChildEdges().get(0).breakPoints.getNullLength(recBP[i]);
                    	totalReassortmentObsProb[i] -= event.node.getChildEdges().get(1).breakPoints.getNullLength(recBP[i]);
                    	totalReassortmentObsProb[i] += event.node.getParentEdges().get(0).breakPoints.getNullLength(recBP[i]);
                    }

                    
                    
//                    totalReassortmentObsProb -= event.node.getChildEdges().get(0).getRecombinationLength();
//                    totalReassortmentObsProb -= event.node.getChildEdges().get(1).getRecombinationLength();
//                    totalReassortmentObsProb += event.node.getParentEdges().get(0).getRecombinationLength();
                    break;
            }
            event.lineages = lineages;
            event.totalRecombinationObsProb = new double[totalReassortmentObsProb.length];
            System.arraycopy(totalReassortmentObsProb, 0, event.totalRecombinationObsProb, 0, totalReassortmentObsProb.length);
        }
        eventListDirty = false;
    }

    @Override
    protected boolean requiresRecalculation() {
        eventListDirty = true;

        return true;
    }

    @Override
    protected void restore() {
        List<RecombinationNetworkEvent> tmp = recombinationNetworkEventList;
        recombinationNetworkEventList = storedRecombinationNetworkEventList;
        storedRecombinationNetworkEventList = tmp;

        super.restore();
    }

    @Override
    protected void store() {
        storedRecombinationNetworkEventList.clear();
        storedRecombinationNetworkEventList.addAll(recombinationNetworkEventList);

        super.store();
    }
}