/**
 * 
 */
package nufeb.examples.ibm;

import java.util.Map;

import eu.quanticol.jsstl.core.formula.AndFormula;
import eu.quanticol.jsstl.core.formula.AtomicFormula;
import eu.quanticol.jsstl.core.formula.EventuallyFormula;
import eu.quanticol.jsstl.core.formula.GloballyFormula;
import eu.quanticol.jsstl.core.formula.NotFormula;
import eu.quanticol.jsstl.core.formula.ParametricExpression;
import eu.quanticol.jsstl.core.formula.ParametricInterval;
import eu.quanticol.jsstl.core.formula.ReferencedFormula;
import eu.quanticol.jsstl.core.formula.SignalExpression;
import eu.quanticol.jsstl.core.formula.SomewhereFormula;
import eu.quanticol.jsstl.core.formula.jSSTLScript;

/**
 * @author Bowen Li
 *
 */
public class IBMScript extends jSSTLScript {

	public static final double bioDt_CONST_ = 4;
	public static final double maxT_CONST_ = 10;
	
	// pressure
	public static final double minP_CONST_ = 0.99;
	public static final int P_VAR_ = 0;
	
	
	// surface 
	public static final double minVF_CONST_ = 0.01;
	public static final int VF_VAR_ = 0;
	
	public IBMScript() {
		super( 
			new String[] {
				"P"
			}
		);	
		
		// pressure: G[0,maxT](!((P > minP) & !F[0,bioDt](P < minP)))
		addFormula( "pressure-s2" ,
				new EventuallyFormula( 
					new ParametricInterval( 
						new ParametricExpression() {
						
							public SignalExpression eval( final Map<String,Double> parameters ) {
					
								return new SignalExpression() {
									
									public double eval( double ... variables ) {
										return 0;
									}
									
								};					
								
							}
							
						} , 
						new ParametricExpression() {
						
							public SignalExpression eval( final Map<String,Double> parameters ) {
					
								return new SignalExpression() {
									
									public double eval( double ... variables ) {
										return bioDt_CONST_;
									}
									
								};					
								
							}
							
						} 		
					),
					new AtomicFormula( 
							new ParametricExpression( ) {
							
								public SignalExpression eval( Map<String, Double> parameters ) {
									
									return new SignalExpression() {						
												
										public double eval(double... variables) {
											return minP_CONST_ - variables[getIndex(P_VAR_)];
										}	
															
									};	
												
								}
							
							} , 
					false
					)		
				)		
				 ,
		null );
		
		addFormula( "pressure-s1" ,
				new AtomicFormula( 
						new ParametricExpression( ) {
						
							public SignalExpression eval( Map<String, Double> parameters ) {
								
								return new SignalExpression() {						
											
									public double eval(double... variables) {
										return variables[getIndex(P_VAR_)] - minP_CONST_;
									}	
														
								};	
											
							}
						
						} , 
				false
				),
		null );
		
		addFormula( "pressure" ,
				new GloballyFormula( 
						new ParametricInterval( 
							new ParametricExpression() {
							
								public SignalExpression eval( final Map<String,Double> parameters ) {
						
									return new SignalExpression() {
										
										public double eval( double ... variables ) {
											return 0;
										}
										
									};					
									
								}
								
							} , 
							new ParametricExpression() {
							
								public SignalExpression eval( final Map<String,Double> parameters ) {
						
									return new SignalExpression() {
										
										public double eval( double ... variables ) {
											return maxT_CONST_;
										}
										
									};					
									
								}
								
							} 		
						),
						new NotFormula( 
								new AndFormula(
										new ReferencedFormula( 
											this ,
											"pressure-s1"
											),
										new NotFormula( 
											new ReferencedFormula( 
												this ,
												"pressure-s2"
											)	
										)	
								)				
						)		
					),
		null );

		
		// surface
		addFormula( "surface-s1" ,
				new AtomicFormula( 
						new ParametricExpression( ) {
						
							public SignalExpression eval( Map<String, Double> parameters ) {
								
								return new SignalExpression() {						
											
									public double eval(double... variables) {
										return variables[getIndex(VF_VAR_)] - minVF_CONST_;
									}	
														
								};	
											
							}
						
						} , 
				false
				)			
				 ,
		null );
		
		addFormula( "surface-s2" ,
				new SomewhereFormula( 
					new ParametricInterval( 
						new ParametricExpression() {
						
							public SignalExpression eval( final Map<String,Double> parameters ) {
					
								return new SignalExpression() {
									
									public double eval( double ... variables ) {
										return 0;
									}
									
								};					
								
							}
							
						} , 
						new ParametricExpression() {
						
							public SignalExpression eval( final Map<String,Double> parameters ) {
					
								return new SignalExpression() {
									
									public double eval( double ... variables ) {
										return 1;
									}
									
								};					
								
							}
							
						} 		
					)		
					 ,
					new AtomicFormula( 
							new ParametricExpression( ) {
							
								public SignalExpression eval( Map<String, Double> parameters ) {
									
									return new SignalExpression() {						
												
										public double eval(double... variables) {
											return minVF_CONST_ - variables[getIndex(VF_VAR_)];
										}	
															
									};	
												
								}
							
							} , 
					false
					)	
				)		
				 ,
		null );
		
		addFormula( "surface-s3" ,
				new NotFormula( 
					new AndFormula(
							new ReferencedFormula( 
								this ,
								"surface-s1"
								),
							new NotFormula( 
								new ReferencedFormula( 
									this ,
									"surface-s2"
								)	
							)	
					)				
				),
		null );
		
		addFormula( "surface" ,
				new AndFormula(
					new ReferencedFormula( 
						this ,
						"surface-s1"
						),
					new ReferencedFormula( 
						this ,
						"surface-s3"
					)
				),
		null );
	}
}
