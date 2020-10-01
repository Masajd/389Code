public abstract class GeneralEMField{
	/**
	* Represents an arbitrary electromagnetic field.  SI units are used throughout.
	*
	* @author Ian Bailey
	* @version 1.0
	*/

 	
 
 	/**
 	*  Return the electric field strength being experienced by a particle
 	*
 	* @return The electric field strength
 	*/
 	public abstract PhysicsVector getElectric(Particle aParticle);
 	
 	/**
 	*  Get the magnetic flux density being experienced by a particle
 	*
 	* @return The magnetic flux density
 	*/
 	public abstract PhysicsVector getMagnetic(Particle aParticle);
 	
 	
 	/**
 	*  Return the electric field strength at a position
 	*
 	* @return The electric field strength 
 	*/
 	public abstract PhysicsVector getElectric(PhysicsVector aPosition);
 	
 	/**
 	*  Get the magnetic flux density being experienced by a particle
 	*
 	* @return The magnetic flux density
 	*/
 	public abstract PhysicsVector getMagnetic(PhysicsVector aPosition);
 	
 	/**
 	*  Get the potential energy of a particle due to the field
 	*
 	* @return The potential energy in J
 	*/
 	public abstract double getPotentialE(ChargedParticle aParticle);
 	
 	/**
 	* Returns the acceleration experienced by a charged particle according to the Lorentz force law (non-relativistic).
 	* @param theParticle - the charged particle moving in the field
 	* @return the acceleration calculated from (qE + vXB / m)
 	*/
 	public PhysicsVector getAcceleration(ChargedParticle theParticle)
 	{
 		PhysicsVector electric = this.getElectric(theParticle);
 		PhysicsVector magnetic = this.getMagnetic(theParticle);
 		PhysicsVector lorentz=new PhysicsVector(electric); // E
 		lorentz.increaseBy(PhysicsVector.cross(theParticle.getVelocity(),magnetic)); //+ v cross B
 		lorentz.scale(theParticle.getCharge()/theParticle.getMass()); // multiply by charge divided by mass
 		return lorentz;
 		
 	}
 	
 	
}
