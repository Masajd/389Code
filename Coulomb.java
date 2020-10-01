public class Coulomb extends GeneralEMField{
	/**
	* Represents the Coulomb field around a point charge. SI units are used throughout.
	*
	* @author Ian Bailey
	* @version 1.0
	*/
	
	protected ChargedParticle sourceParticle;
	public final static double k= 1/(4*Math.PI*8.85418782E-12);// 1/4\pi\epsilon_0
	
 	/**
 	*  Default constructor. Associate the Coulomb field with a default charged
 	* particle
 	*
 	*/
 	public Coulomb(){
 		sourceParticle=new ChargedParticle();
 		
 	}
 	
 	/**
 	*  Default constructor. Associate a Coulomb field with a charged particle
 	* @param aParticle the particle with which to associate the Coulomb field
 	*/
 	public Coulomb(ChargedParticle aParticle){
 		sourceParticle=aParticle;
 		
 	}
 	
 	
 	/**
 	*  Return the electric field strength due to the Coulomb field
 	*
 	* @param aPosition The position at which the electric field should be calculated
 	* @return The value of the electric field strength from the Coulomb field
 	*/
 	public PhysicsVector getElectric(PhysicsVector aPosition){
 		PhysicsVector sourcePos = sourceParticle.getPosition();
 		double sourceQ = sourceParticle.getCharge();
 		
 		PhysicsVector r = PhysicsVector.subtract(aPosition, sourcePos);
 		return PhysicsVector.scale(k*sourceQ/(r.magnitude()*r.magnitude()),r.getUnitVector());
 		
 		
 	}
 	
 		
 	/**
 	*  Return the potential energy of a particle in the Coulomb field
 	*
 	* @param aParticle The charged particle whose energy is being calculated
 	* @return The potential energy of the particle in J
 	*/
 	public double getPotentialE(ChargedParticle aParticle){
 		PhysicsVector sourcePos = sourceParticle.getPosition();
 		PhysicsVector particlePos = aParticle.getPosition();
 		double sourceQ = sourceParticle.getCharge();
 		double particleQ= aParticle.getCharge();
 		
 		PhysicsVector r = PhysicsVector.subtract(particlePos, sourcePos);
 		return k*particleQ*sourceQ/r.magnitude();
 		
 		
 	}
 	
 	/**
 	*  Get the magnetic flux density due to the Coulomb field(i.e. zero if relativity is neglected)
 	*
 	* @param aPosition The position at which the magnetic flux density should be calculated
 	* @return The  value of the magnetic flux density associated with the Coulomb field
 	*/
 	public PhysicsVector getMagnetic(PhysicsVector aPosition){
 		return new PhysicsVector();
 	}
 	
 		
 	/**
 	*  Return the electric field strength due to the Coulomb field
 	*
 	* @param  aParticle the particle denoting the position at which the field is calculated
 	* @return The  value of the electric field strength at the position of the particle
 	*/
 	public PhysicsVector getElectric(Particle aParticle){
 		
 		return getElectric(aParticle.getPosition());
 	}
 	
 	
 	/**
 	*  Get the magnetic flux density due to the Coulomb field (i.e. zero if relativity is neglected)
 	*
 	* @param aParticle aParticle the particle denoting the position at which the field is calculated
 	* @return  The value of the magnetic flux density at the position of the particle
 	*/
 	public PhysicsVector getMagnetic(Particle aParticle){
 		return new PhysicsVector();
 	}
 	
 	
}
