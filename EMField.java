public class EMField extends GeneralEMField{
	/**
	* Represents a uniform electric and/or magnetic field.  SI units are used throughout.
	*
	* @author Ian Bailey
	* @version 1.0
	*/
	
 	protected PhysicsVector electric; // electric field strength
 	protected PhysicsVector magnetic; // magnetic flux density
	
 	
 	/**
 	*  Default constructor. Set data members to zero.
 	*
 	*/
 	public EMField(){
 		electric = new PhysicsVector();
 		magnetic = new PhysicsVector();
 	}
 	
 	/**
 	*  Constructor with two inputs - the electric field strength and magnetic flux density
 	*
 	* @param electricIn The electric field strength
 	* @param magneticIn The magnetic flux density
 	*/
 	public EMField(PhysicsVector electricIn, PhysicsVector magneticIn){
 		electric = new PhysicsVector(electricIn);
 		magnetic = new PhysicsVector(magneticIn);
 	}
 	
 	
 	/**
 	*  Set the electric field strength
 	*
 	* @param electricIn The electric field strength
 	*/
 	public void setElectric(PhysicsVector electricIn){
 		electric = new PhysicsVector(electricIn);
 	}
 	
 	/**
 	*  Set the magnetic flux density
 	*
 	* @param magneticIn The magnetic flux density
 	*/
 	public void setMagnetic(PhysicsVector magneticIn){
 		magnetic = new PhysicsVector(magneticIn);
 	}
 	
 	/**
 	*  Return the electric field strength
 	*
 	* @return The current value of the electric field strength
 	*/
 	public PhysicsVector getElectric(){
 		return new PhysicsVector(electric);
 	}
 	
 	
 	/**
 	*  Get the magnetic flux density
 	*
 	* @return The current value of the magnetic flux density
 	*/
 	public PhysicsVector getMagnetic(){
 		return new PhysicsVector(magnetic);
 	}
 	
 	
 	/**
 	*  Return the electric field strength
 	*
 	* @return The current value of the electric field strength
 	*/
 	public PhysicsVector getElectric(PhysicsVector aPosition){
 		return getElectric();
 	}
 	
 	
 	/**
 	*  Get the magnetic flux density
 	*
 	* @return The current value of the magnetic flux density
 	*/
 	public PhysicsVector getMagnetic(PhysicsVector aPosition){
 		return getMagnetic();
 	}
 	
 	/**
 	*  Return the electric field strength
 	*
 	* @return The current value of the electric field strength
 	*/
 	public PhysicsVector getElectric(Particle aParticle){
 		return getElectric();
 	}
 	
 	
 	/**
 	*  Get the magnetic flux density
 	*
 	* @return The current value of the magnetic flux density
 	*/
 	public PhysicsVector getMagnetic(Particle aParticle){
 		return getMagnetic();
 	}
 	
 	/**
 	*  Return the potential energy of a particle due to the uniform field (assumed to be zero).
 	*
 	* @param aParticle The charged particle whose energy is being calculated
 	* @return The potential energy of the particle in J
 	*/
 	public double getPotentialE(ChargedParticle aParticle){
 		return 0.0;
 		
 		
 	} 	
}
