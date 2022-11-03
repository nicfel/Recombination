/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package recombination.annotator;

import beast.base.core.Log;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

//import bacter.Conversion;
//import bacter.ConversionGraph;
//import bacter.Locus;
//import bacter.util.BacterACGLogReader;

/**
 * A rewrite of TreeAnnotator that outputs reassortment distances 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class RecombinationDistance extends RecombinationAnnotator {


    List<Double> leaveDistance;

    
    private static class NetworkAnnotatorOptions {

        File inFile;
        File outFile = new File("recombination_distances.txt");
        double burninPercentage = 10.0;
        BreakPoints breakPoints = new BreakPoints();
        String sequence = null;

        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
            		"Remove Loci for summary: " + breakPoints + "\n";
        }
    }

    public RecombinationDistance(NetworkAnnotatorOptions options) throws IOException {

        // Display options:
        System.out.println(options + "\n");

        // Initialise reader
        RecombinationLogReader logReader = new RecombinationLogReader(options.inFile,
                options.burninPercentage);

        System.out.println(logReader.getNetworkCount() + " Networks in file.");

        System.out.println("The first " + logReader.getBurnin() +
                 " (" + options.burninPercentage + "%) ACGs will be discarded " +
                "to account for burnin.");

	      System.out.println("\nWriting output to " + options.outFile.getName()
	      + "...");
        
        // compute the pairwise reassortment distances 
	    int i =0;
        try (PrintStream ps = new PrintStream(options.outFile)) {
	        for (RecombinationNetwork network : logReader){
	        	pruneLociFromNetwork(network, options.breakPoints);
	        	
	        	if (options.sequence!=null) {
	        		//prune all breakpoints not ancestral to options.sequence
	                List<RecombinationNetworkNode> leafs = network.getNodes().stream()
	                        .filter(e -> e.isLeaf())
	                        .collect(Collectors.toList());    	

	        		boolean foundLeaf = false;
	        		for (RecombinationNetworkNode leaf : leafs) {
	        			if (leaf.getTaxonLabel().contentEquals(options.sequence)) {
	        				pruneUp(leaf.getParentEdges().get(0), leaf.getParentEdges().get(0).breakPoints);
	        				foundLeaf = true;
	        			}
	        		}
	        		
	        		if (!foundLeaf)
	        			throw new IllegalArgumentException("leaf with name " + options.sequence + " not found");
	        		
	        	}
	        	computeRecombinationDistance(network, ps, i, options.sequence!=null);
	        	i++;
	        }
	        ps.close();
        }
        System.out.println("\nDone!");
    }
    
    private void pruneUp(RecombinationNetworkEdge edge, BreakPoints bp) {
    	if (edge.isRootEdge())
    		return;
    	
    	BreakPoints bp_tmp = bp.copy(); 
    	bp_tmp.and(edge.breakPoints);
    	
    	if (bp_tmp.isEmpty())
    		return;
    	
    	edge.carryingRange = bp_tmp.copy();

    	
    	for (RecombinationNetworkEdge e : edge.parentNode.getParentEdges())
    		pruneUp(e, bp_tmp);    	
		
	}

	private void computeRecombinationDistance(RecombinationNetwork network, PrintStream ps, int samplenr, boolean hasSeq){    	
    	// get all reassortment nodes    	
        List<RecombinationNetworkNode> recombinationNodes = network.getNodes().stream()
                .filter(e -> e.isRecombination())
                .collect(Collectors.toList());    	
        
        int j=0;
        
        for (RecombinationNetworkNode node : recombinationNodes){
        	
        	if (hasSeq && node.getChildEdges().get(0).carryingRange==null)
        		continue;
        	
        	if (node.getParentEdges().get(0).breakPoints.isEmpty())
        		throw new IllegalArgumentException("Empty parent edge");
        	if (node.getParentEdges().get(1).breakPoints.isEmpty())
        		throw new IllegalArgumentException("Empty parent edge");

        	// follow the nodes uphill and always denote which
        	Map<Integer, BreakPoints> uphill = new HashMap<>();     
        	
        	
        	getNodesUphill(uphill, node.getParentEdges().get(0).breakPoints.copy(), node.getParentEdges().get(0));
        	List<Double> heights = new ArrayList<>();
        	List<BreakPoints> breaks = new ArrayList<>();
        	getDistances(uphill, heights, breaks, node.getParentEdges().get(1).breakPoints.copy(), node.getParentEdges().get(1));
        	// compute weighted height
        	double height = 0.0;
        	double weight = 0.0;
        	for (int i = 0; i < breaks.size(); i++) {
        		height += heights.get(i) * breaks.get(i).getGeneticLength();
        		if (hasSeq) {
        			BreakPoints bp_tmp = breaks.get(i).andCopy(node.getChildEdges().get(0).carryingRange);
        			weight += bp_tmp.getGeneticLength();
        		}else {
        			weight += breaks.get(i).getGeneticLength();
        		}
        	}
        	height/=weight;
        	double minedge = -1;
        	
        	if (hasSeq) {
    			BreakPoints bp_tmp1 = node.getParentEdges().get(0).breakPoints.andCopy(node.getChildEdges().get(0).carryingRange);
    			BreakPoints bp_tmp2 = node.getParentEdges().get(1).breakPoints.andCopy(node.getChildEdges().get(0).carryingRange);
    			minedge = Math.min(bp_tmp1.getGeneticLength(), bp_tmp2.getGeneticLength());    		
    		}else {
    			minedge = Math.min(node.getParentEdges().get(0).breakPoints.getGeneticLength(), node.getParentEdges().get(1).breakPoints.getGeneticLength());    		
			}
        	
        	if (weight==0 || minedge==0)
        		continue;
       
        	       	
        	
            ps.print(Math.max(node.getParentEdges().get(0).passingRange.getMin(), node.getParentEdges().get(1).passingRange.getMin()));
            ps.print("\t");
            ps.print(height);
            ps.print("\t");
            ps.print(node.getHeight());
            ps.print("\t");
            ps.print(minedge);
            ps.print("\t");

            ps.print(samplenr);
            j++;
        	ps.print("\n"); 
        }  
    }
    
    private void getNodesUphill(Map<Integer, BreakPoints> uphill, BreakPoints cp, RecombinationNetworkEdge edge) {
    	if (cp.isEmpty())
    		return;

    	RecombinationNetworkNode node = edge.parentNode;
    	if (uphill.get(node.ID)==null) {
    		uphill.put(node.ID, cp.copy());
    	}else {
    		uphill.get(node.ID).or(cp);
    	}
    	
    	for (RecombinationNetworkEdge e : node.getParentEdges()) {
    		if (e.isRootEdge())
    			return;
    		
    		BreakPoints cp_up = cp.copy();
    		cp_up.and(e.breakPoints);
        	getNodesUphill(uphill, cp_up, e);
    	}    	
    }
    
    private void getDistances(Map<Integer, BreakPoints> uphill, List<Double> heights, List<BreakPoints> breaks, BreakPoints cp, RecombinationNetworkEdge edge) {
    	if (cp.isEmpty())
    		return;
    	
    	
    	RecombinationNetworkNode node = edge.parentNode;
    	if (uphill.get(node.ID)!=null) {    		
    		heights.add(node.getHeight());
    		breaks.add(cp.copy());
    		return;
    	}
    	
    	for (RecombinationNetworkEdge e : node.getParentEdges()) {
    		BreakPoints cp_up = cp.copy();
    		cp_up.and(e.breakPoints);
    		getDistances(uphill, heights, breaks, cp_up, e);
    	}
    	
    }

    

    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(NetworkAnnotatorOptions options) {

        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("Reassortment Network Annotator");

        JLabel logFileLabel = new JLabel("Reassortment Network log file:");
        JLabel outFileLabel = new JLabel("Output file:");
        JLabel burninLabel = new JLabel("Burn-in percentage:");
        JLabel summaryMethodLabel = new JLabel("Position summary method:");
        JLabel thresholdLabel = new JLabel("Posterior conversion support threshold:");
        JCheckBox geneFlowCheckBox = new JCheckBox("Record gene flow");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");

        JTextField gfOutFilename = new JTextField(20);
        gfOutFilename.setEditable(false);
        gfOutFilename.setEnabled(false);
        JButton gfOutFileButton = new JButton("Choose File");
        gfOutFileButton.setEnabled(false);

        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.burninPercentage));
        burninSlider.setMajorTickSpacing(50);
        burninSlider.setMinorTickSpacing(10);
        burninSlider.setPaintTicks(true);
        burninSlider.setPaintLabels(true);
        burninSlider.setSnapToTicks(true);


        Container cp = dialog.getContentPane();
        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
        cp.setLayout(boxLayout);

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);
        layout.setAutoCreateGaps(true);
        layout.setAutoCreateContainerGaps(true);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(outFileLabel)
                        .addComponent(burninLabel)
                        .addComponent(summaryMethodLabel)
                        .addComponent(thresholdLabel)
                        .addComponent(geneFlowCheckBox))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename)
                        .addComponent(burninSlider)
                        .addComponent(gfOutFilename))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFileButton)
                        .addComponent(outFileButton)
                        .addComponent(gfOutFileButton)));

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(inFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(inFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(burninLabel)
                        .addComponent(burninSlider,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(geneFlowCheckBox)
                        .addComponent(gfOutFilename)
                        .addComponent(gfOutFileButton)));

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
            dialog.setVisible(false);
        });
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

        JFileChooser inFileChooser = new JFileChooser();
        inFileButton.addActionListener(e -> {
            inFileChooser.setDialogTitle("Select ACG log file to summarize");
            if (options.inFile == null)
                inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = inFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.inFile = inFileChooser.getSelectedFile();
                inFilename.setText(inFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
            outFileChooser.setDialogTitle("Select output file name.");
            if (options.inFile != null)
                outFileChooser.setCurrentDirectory(options.inFile);
            else
                outFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
            }
        });

        geneFlowCheckBox.addActionListener(e -> {
            boolean newValue = geneFlowCheckBox.isSelected();
            gfOutFilename.setEnabled(newValue);
            gfOutFileButton.setEnabled(newValue);
        });

        JFileChooser gfOutFileChooser = new JFileChooser();
        gfOutFileButton.addActionListener(e -> {
            gfOutFileChooser.setDialogTitle("Select gene flow output file name.");
            if (options.inFile != null)
                gfOutFileChooser.setCurrentDirectory(options.inFile);
            else
                gfOutFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

            int returnVal = gfOutFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                gfOutFilename.setText(gfOutFileChooser.getSelectedFile().getName());
            }
        });

        cp.add(buttonPanel);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        return !canceled[0];
    }

    /**
     * Prepare JFrame to which ACGAnnotator output streams will be
     * directed.
     */
    private static void setupGUIOutput() {

        JFrame frame = new JFrame();
        frame.setTitle("ACGAnnotator");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JTextArea textArea = new JTextArea(25, 80);
        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
        textArea.setEditable(false);
        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> System.exit(0));
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

        // Redirect streams to output window:
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                SwingUtilities.invokeLater(() -> {
                    if ((char)b == '\r') {
                        int from = textArea.getText().lastIndexOf("\n") + 1;
                        int to = textArea.getText().length();
                        textArea.replaceRange(null, from, to);
                    } else
                        textArea.append(String.valueOf((char) b));
                });
            }
        };

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));

        frame.pack();
        frame.setVisible(true);
    }

    public static String helpMessage =
            "ACGAnnotator - produces summaries of Bacter ACG log files.\n"
                    + "\n"
                    + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-positions {mean,median} Choose position summary method.\n"
                    + "                         (default mean)\n"
                    + "-burnin percentage       Choose _percentage_ of log to discard\n"
                    + "                         in order to remove burn-in period.\n"
                    + "                         (Default 10%)\n"
                    + "-threshold percentage    Choose minimum posterior probability\n"
                    + "                         for including conversion in summary.\n"
                    + "                         (Default 50%)\n"
                    + "-recordGeneFlow gfFile   Record posterior distribution of gene\n"
                    + "                         flow in given file.\n"
                    + "\n"
                    + "If no output file is specified, output is written to a file\n"
                    + "named 'summary.tree'.";

    /**
     * Print usage info and exit.
     */
    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    /**
     * Display error, print usage and exit with error.
     */
    public static void printUsageAndError(String errMsg) {
        System.err.println(errMsg);
        System.err.println(helpMessage);
        System.exit(1);
    }

    /**
     * Retrieve ACGAnnotator options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, NetworkAnnotatorOptions options) {
        int i=0;
        while (args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

                case "-burnin":
                    if (args.length<=i+1)
                        printUsageAndError("-burnin must be followed by a number (percent)");

                    try {
                        options.burninPercentage = Double.parseDouble(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
                    }

                    if (options.burninPercentage<0 || options.burninPercentage>100) {
                        printUsageAndError("Burnin percentage must be >= 0 and < 100.");
                    }

                    i += 1;
                    break;

                case "-subsetRange":
                    if (args.length<=i+1) {
                        printUsageAndError("-subsetRange must be a range in the format of 0-100.");
                    }

                    try {
                    	String[] argarray = args[i + 1].split(",");
                    	List<Integer> bp_list = new ArrayList<>();
                    	for (int j = 0; j < argarray.length; j++) {
                    		String[] tmp = argarray[j].split("-");
                    		bp_list.add(Integer.parseInt(tmp[0]));
                    		bp_list.add(Integer.parseInt(tmp[1]));
                    	}
                		options.breakPoints.init(bp_list);
                    } catch (NumberFormatException e) {
                        printUsageAndError("removeSegments must be an array of integers separated by commas if more than one");
                     }

                    i += 1;
                    break;
                    
                case "-sequence":
                    if (args.length<=i+1) {
                        printUsageAndError("-sequence must be tha name of a sequence in the network file.");
                    }

                    try {
                		options.sequence = args[i + 1];
                    } catch (NumberFormatException e) {
                        printUsageAndError("sequence must be a sequence name");
                     }

                    i += 1;
                    break;

                    
                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }

        if (i >= args.length)
            printUsageAndError("No input file specified.");
        else
            options.inFile = new File(args[i]);

        if (i+1<args.length)
            options.outFile = new File(args[i+1]);
    }

    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
    	NetworkAnnotatorOptions options = new NetworkAnnotatorOptions();

        if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
                Log.warning.println("Error setting cross-platform look and feel.");
            }

            try {
                SwingUtilities.invokeAndWait(() -> {
                    if (!getOptionsGUI(options))
                        System.exit(0);

                    setupGUIOutput();
                });
            } catch (InterruptedException | InvocationTargetException e) {
                e.printStackTrace();
            }


        } else {
            getCLIOptions(args, options);
        }

        // Run ACGAnnotator
        try {
            new RecombinationDistance(options);

        } catch (Exception e) {
            if (args.length == 0) {
                JOptionPane.showMessageDialog(null, e.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                System.err.println("Error: " + e.getMessage());
                e.printStackTrace();
                System.err.println();
                System.err.println(helpMessage);
            }

            System.exit(1);
        }
    }
}