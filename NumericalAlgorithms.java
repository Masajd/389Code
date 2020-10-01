
/**
 * An interface of constants used in classes which require numerical integration algorithms
 * @author Ian Bailey
 * @version 1.2
 */
 public interface NumericalAlgorithms{
 	 
 	 final int nEulerMidPoint=0;
	 final int nEuler=1;
 	 final int nEulerCromer=2;
 	 final String[] naNames = {"Euler MidPoint", "Euler (FPA)", "Euler-Cromer (LPA)"};
 }
