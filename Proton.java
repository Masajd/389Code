
/**
 * A Class to represent a massive, electrically charged particle
 * @author Ian Bailey
 * @version 1.0
 */
public final class Proton extends ChargedParticle{

	public static final double pCharge=1.60217657e-19; //the electric charge of the proton in Coulombs
	public static final double pMass=1.67262178E-27; // the mass of the proton in kilograms
	
	private double defaultCharge=pCharge;
	private double defaultMass=pMass;
	
	
	/**
	* The Default Constructor. Sets everything to zero except the mass and charge of the proton.
	*
	*/
	public Proton(){
		super();
		setDefaults();
	}
    
	/**
	* Constructor with two inputs - the charge and mass of the proton. Inputs are ignored.
	* @param mIn the mass (ignored)
	* @param qIn the charge (ignored)
	*/
	public Proton(double mIn, double qIn){
		super();
		setDefaults();
		
	}
    
	/**
	*  Constructor that sets mass, charge, position and velocity
	*  @param mIn mass of the proton (ignored)
	*  @param qIn charge of the proton (ignored)
	*  @param positionIn initial position of proton
	*  @param velocityIn initial velocity of proton 
	*/
	public Proton(double mIn,double qIn, PhysicsVector positionIn,PhysicsVector velocityIn)
	{
		super(mIn, positionIn, velocityIn);
		setDefaults();
	}
    
	/**
	*  Copy Constructor 
	*  @param particleIn proton whose properties are to be copied to the new proton
	*/
	public Proton(Proton particleIn)
	{
		super(particleIn.mass, particleIn.position, particleIn.velocity);
		setDefaults();	
		
	}
	
	/**
	* Set the default charge and mass for a proton
	*
	*/
	private void setDefaults()
	{
		super.setCharge(defaultCharge);
		super.setMass(defaultMass);
	}
	
	/**
	* Set the charge
	*
	* @param chargeIn The new charge
	*/
	@Override
	public void setCharge(double qIn)
	{
		
	}
    
    	/**
	* Set the mass
	*
	* @param massIn The new mass
	*/
	@Override
	public void setMass(double massIn)
	{
		
	}
	
	/**
	*  Method to set the properties of the proton equal to those of another proton
	*  @param protonIn proton whose properties are to be copied to 'this' proton
	*/
	public void setAll(Proton protonIn)
	{
		setPosition(protonIn.position);
		setVelocity(protonIn.velocity);
		setAcceleration(protonIn.acceleration);
		
	}
		
}


