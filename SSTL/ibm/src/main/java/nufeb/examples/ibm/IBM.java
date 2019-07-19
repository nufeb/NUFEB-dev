package nufeb.examples.ibm;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Locale;




//import eu.quanticol.jsstl.core.dsl.ScriptLoader;
import eu.quanticol.jsstl.core.formula.Signal;
import eu.quanticol.jsstl.core.formula.SignalStatistics;
import eu.quanticol.jsstl.core.formula.jSSTLScript;
import eu.quanticol.jsstl.core.signal.BooleanSignal;
import eu.quanticol.jsstl.core.space.GraphModel;
import eu.quanticol.jsstl.core.io.TraGraphModelReader;
import eu.quanticol.jsstl.core.monitor.SpatialBooleanSignal;

public class IBM {

	public static void main(String[] args) throws Exception {
		long tStart = System.currentTimeMillis();
		//// /// %%%%%%% GRAPH %%%%% /////////
		TraGraphModelReader graphread = new TraGraphModelReader();
		GraphModel graph = graphread.read("models/sstl_model.txt");
		long start = System.currentTimeMillis();

		graph.dMcomputation();
		double dmtime = (System.currentTimeMillis()-start)/1000.0;
		int runs = 1;

		//// /// %%%%%%% PROPERTY %%%%% //////////
		jSSTLScript script = new IBMScript();
		// /// %%%%%%% DATA import %%%%%%%%%%%%/////////

		double t1 = check(script, graph, runs);

		System.out.println("phi1: "+t1);
		System.out.println("DM Computation: "+dmtime);

	}


	public static Signal readSystemTraj(GraphModel graph, String filename) throws IOException{
		FileReader fr = new FileReader(filename);
		BufferedReader br = new BufferedReader(fr);
		ArrayList<String> strings = new ArrayList<String>();
		String line = null;
		while ((line = br.readLine()) != null) {
			strings.add(line);
		}
		int timeSize = strings.size();

		int locSize = 2;
		double[][][] data = new double[locSize][timeSize][1];
		double[] time = new double[strings.size()];
		for (int t = 0; t < timeSize; t++) {
			String[] splitted = strings.get(t).split(",");
			time[t] = Double.valueOf(splitted[0]);
		}
	
		for (int t = 0; t < timeSize; t++) {
			String[] splitted = strings.get(t).split(",");
			for (int i = 0; i < locSize; i++) {		
				data[i][t][0] = Double.valueOf(splitted[1]);
			}
		}

		br.close();
		fr.close();
		Signal s = new Signal(graph, time, data);
		return s;
 }
	
	public static Signal readMeshTraj(GraphModel graph, String filename) throws IOException{
		FileReader fr = new FileReader(filename);
		BufferedReader br = new BufferedReader(fr);
		ArrayList<String> strings = new ArrayList<String>();
		String line = null;
		while ((line = br.readLine()) != null) {
			strings.add(line);
		}
		String[] header = strings.get(0).split(",");

		int timeSize = Integer.valueOf(header[0]);
		int locSize = Integer.valueOf(header[1]);
		
		double[][][] data = new double[locSize][timeSize][1];
		double[] time = new double[timeSize];
		
		for (int t = 0; t < timeSize; t++) {
			String[] splitted = strings.get(t*locSize+t+1).split(",");
			time[t] = Double.valueOf(splitted[1]);
		}

		for (int t = 0; t < timeSize; t++) {
			for (int i = 0; i < locSize; i++) {	
				String[] splitted = strings.get(t*locSize+i+t+2).split(",");
				int loc = Integer.valueOf(splitted[0])-1;
				data[loc][t][0] = Double.valueOf(splitted[1]);
				//System.out.println("loc: "+loc + " value: "+  data[loc][t][0]);
			}
		}

		br.close();
		fr.close();
		Signal s = new Signal(graph, time, data);
		return s;
 }

	public static double check(jSSTLScript script, GraphModel graph, int runs) throws IOException {
		PrintWriter stprinter = new PrintWriter("data/"+"surface"+"boolen.txt");
		
		for (int i = 0; i < runs; i++) {
			
			Signal s = readMeshTraj(graph, "models/trajectories/surface/volf_all.txt");

			SpatialBooleanSignal boolSign = script.getFormula("surface").booleanCheck(null, graph, s);
			for (int j = 0; j < graph.getNumberOfLocations(); j++) {
				BooleanSignal b = boolSign.spatialBoleanSignal.get(graph.getLocation(j));
				int value = b.getValueAt(10.0) ? 1 : 0;
				stprinter.println(graph.getLocation(j)+ "," + value);
			}

		}
		stprinter.close();
		return runs;
	}
	
}
