/**
* A Class to keep track of a bunch's orbit
* @author Ian Bailey
* @author Sam Dicker
* @version 1.4
*/

import java.lang.Math;
import java.util.*;

public class OrbitTracker<T extends Particle>{

	public boolean converge=false; // flag to indicate whether the bunch is approaching or diverging from the origin
	double displacement=0.0; // distance between bunch origin and current position
	PhysicsVector origin;
	
	double displacement2 = 0.0; // distance between bunch origin and current position
	public boolean outField = false; //flag to indicate when the bunch is out of the electric field, and nearing the maxima 
	
	private double c=3e8; // speed of light
	
	public OrbitTracker()
	{
		converge=false;
		displacement=0.0;
		origin = new PhysicsVector();
	}

	public OrbitTracker(PhysicsVector originIn)
	{
		converge=false;
		displacement=0.0;
		origin = new PhysicsVector(originIn);
	}
	
	/** 
	* Tests whether an orbit has been completed by assuming that the bunch
	* will return to its approximate initial position. 
	* @return true if the bunch has completed an orbit
	*/
	public boolean hasOrbited(Bunch<T> theBunch)
	{
	 	double currentDis=(PhysicsVector.subtract(theBunch.getPosition(),origin)).magnitude();
	 	boolean anOrbit=false;
	 	
	 	if (currentDis>=displacement)
		{
	 		if (converge)
			{
	 			// passed through closest approach and is now diverging from origin
	 			anOrbit=true;
	 		}
	 		converge=false;
	 	}
	 	else
		{
	 		converge=true; // bunch is approaching origin
	 	}
	 	displacement=currentDis;
	 	return anOrbit ;	 
	}
	
	/**
	* Tests whether the proton has reached the edge of the cyclotron
	* @param cyclotronRad radius of the cyclotron
	* @param orbitRadius radius of the protons orbit
	* @return true if the proton is at a distance greater than the cyclotron radius
	*/
	public boolean hasReachedEdge(Bunch<T> theBunch, double cyclotronRad, double orbitRadius)
	{
		boolean theEdge = false;
		
		if(orbitRadius >= cyclotronRad)
		{
			theEdge = true;
		}
	
	return theEdge;
	}
	
	/**
	* Tests whether the proton is at the maxima (in terms of Y) of its orbit
	* @param Radius the distance of the proton to the centre of its orbit
	* @return true if the proton is at the maxima
	*/
	public boolean MiddleOrbit(Bunch<T> theBunch, double Radius)
	{
		boolean MidOrbit = false;
		double currentDisplacement = Math.abs(theBunch.getPosition().getX()-origin.getX() - Radius); //Distance from centre of orbit
		
		if(currentDisplacement >= displacement2)
		{
			if(outField)
			{	// passed through closest approach and is now diverging from maxima
				MidOrbit = true;		
			}
			outField = false;
		}
		else
		{
			outField = true; // bunch is approaching maxima and out of the electric field
		}
		displacement2 = currentDisplacement;
		return MidOrbit;
	}
}
