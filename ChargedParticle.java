
/**
 * A Class to represent a massive, electrically charged particle
 * @author Ian Bailey
 * @version 1.6
 */
public class ChargedParticle extends Particle implements Comparable{

	protected double charge; //the electric charge of the particle in Coulombs
        private double defaultCharge=0.0;
        
	/**
	* The Default Constructor. Sets everything to zero.
	*
	*/
	public ChargedParticle(){
		super();
		setCharge(this.defaultCharge);
	}
    
	/**
	* Constructor with two inputs - the charge and mass of the particle. Set everything else to zero.
	* @param mIn the mass
	* @param qIn the charge
	*/
	public ChargedParticle(double mIn, double qIn){
		super(mIn);
		setCharge(qIn);
		
	}
    
	/**
	*  Constructor that sets mass, position and velocity
	*  @param mIn mass of the charged particle
	*  @param positionIn initial position of charged particle
	*  @param velocityIn initial velocity of charged particle 
	*/
	public ChargedParticle(double mIn,PhysicsVector positionIn,PhysicsVector velocityIn)
	{
		super(mIn, positionIn, velocityIn);
	}
	
	/**
	*  Constructor that sets mass, charge, position and velocity
	*  @param mIn mass of the particle
	*  @param qIn charge of the particle
	*  @param positionIn initial position of particle
	*  @param velocityIn initial velocity of particle 
	*/
	public ChargedParticle(double mIn,double qIn, PhysicsVector positionIn,PhysicsVector velocityIn)
	{
		super(mIn, positionIn, velocityIn);
		setCharge(qIn);
		
	}
    
	/**
	*  Copy Constructor 
	*  @param particleIn charged particle whose properties are to be copied to the new particle
	*/
	public ChargedParticle(ChargedParticle particleIn)
	{
		super(particleIn);
		setCharge(particleIn.charge);
		
	}
	
	/**
	*  Method to set the properties of the particle equal to those of another particle
	*  @param particleIn particle whose properties are to be copied to 'this' particle
	*/
	public void setAll(ChargedParticle particleIn)
	{
		super.setAll(particleIn);
		setCharge(particleIn.charge);
	}
	
	/**
	* Return the charge
	*
	* @return charge
	*/
	public double getCharge()
	{
		return charge;
	}
    
	/**
	* Set the charge
	*
	* @param chargeIn The new charge
	*/
	public void setCharge(double qIn)
	{
		charge=qIn;
	}
    
    
	/**
	* Create a string containing the mass, charge, position, velocity, and acceleration of the particle.
	* This method is called automatically by System.out.println(someparticle)
	* @return string with the format
	* " mass "+mass+" charge "+charge+" Position: "+position+" Velocity: "+velocity+" Acceleration: "+acceleration
	*/
	@Override
	public String toString()
	{
		return " mass "+mass+ " charge " +charge+" Position: "+position.returnSimpleString()+" Velocity: "+velocity.returnSimpleString()+" Acceleration: "+acceleration.returnSimpleString();
	}
	
	/**
	* Compares a ChargedParticle with position (p2) to 'this' ChargedParticle's position (p1)
	* and determines which has the largest absolute value of the 
	* y coordinate of position.
	* @return -1  if |p2|>|p1|, 0 if |p2|==|p1|, 1 otherwise
	*/
	public int compareTo(Object other) throws ClassCastException
	{
		double myY=Math.abs((this.getPosition()).getY());
		double otherY=Math.abs((((ChargedParticle)other).getPosition()).getY());
    	
		if (myY< otherY){
			return(-1);
		}
		else if (myY==otherY){
			return(0);
		}
		else
		{
			return(1);
		}			
	}
	
}


