
/**
* A Class to represent a massive particle
*  It can have position, velocity, acceleration and mass
* @author  Alex Finch
* @author Ian Bailey
* @author Sam Dicker
* @version 1.7
*/
public class Particle implements NumericalAlgorithms{
	// NB Private => only visible inside this class
	// Protected => visible inside this class AND any subclasses
	protected double mass; //the mass of the particle 
	protected PhysicsVector position, velocity, acceleration;
	private double defaultMass=0.0;
	protected int algorithm=nEulerMidPoint; // Numerical integration algorithm to use
	
	/**
	* The Default Constructor. Sets everything to zero.
	*
	*/
	public Particle(){
		mass = this.defaultMass;
		position = new PhysicsVector();
		velocity = new PhysicsVector();
		acceleration= new PhysicsVector();
	}
    
	/**
	* Constructor with one input, the mass of the particle. Set everything else to zero.
	* @param mIn mass of the particle
	*/
	public Particle(double mIn){
		this();
		mass = mIn;
	}
    
	/**
	*  Constructor that sets mass, position and velocity
	*  @param mIn mass of the particle
	*  @param positionIn initial position of particle
	*  @param velocityIn initial velocity of particle 
	*/
	public Particle(double mIn,PhysicsVector positionIn,PhysicsVector velocityIn)
	{
		this(mIn);
		position = new PhysicsVector(positionIn);
		velocity = new PhysicsVector(velocityIn);
	}
    
    
	/**
	*  Copy Constructor 
	*  @param particleIn particle whose properties are to be copied to the new particle
	*/
	public Particle(Particle particleIn)
	{
		this(particleIn.mass,particleIn.position,particleIn.velocity);
		acceleration= new PhysicsVector(particleIn.acceleration);
	}
	
	/**
	*  Method to set the properties of the particle equal to those of another particle
	*  @param particleIn particle whose properties are to be copied to 'this' particle
	*/
	public void setAll(Particle particleIn)
	{
		setMass(particleIn.mass);
		setPosition(particleIn.position);
		setVelocity(particleIn.velocity);
		setAcceleration(particleIn.acceleration);
	}
	
	/**
	* Return the kinetic energy
	*
	* @return the kinetic energy in J
	*/
	public double getKE()
	{
		double speed=velocity.magnitude();
		return 0.5*mass*speed*speed;
	}
    
	
	/**
	* Return the position
	*
	* @return position
	*/
	public PhysicsVector getPosition()
	{
		return new PhysicsVector(position);
	}
    
	/**
	* Return the velocity
	*
	* @return velocity
	*/
	public PhysicsVector getVelocity()
	{
		return new PhysicsVector(velocity);
	}
	
	/**
	* Return the acceleration
	*
	* @return acceleration
	*/
	public PhysicsVector getAcceleration()
	{
		return new PhysicsVector(acceleration);
	}
	
	
	/**
	* Return the mass
	*
	* @return mass
	*/
	public double getMass()
	{
		return mass;
	}
	
	
	/**
	* Set the numerical algorithm for integration
	*
	* @param algoIn The algorithm to use as defined in the NumericalAlgorithms interface
	*/
	public void setAlgo(int algoIn)
	{
		algorithm=algoIn;
	}
	
	/**
	* Set the mass
	*
	* @param massIn The new mass
	*/
	public void setMass(double massIn)
	{
		mass=massIn;
	}
    
	/**
	* Set the position
	*
	* @param pIn The new position
	*/
	public void setPosition(PhysicsVector pIn)
	{
		position = new PhysicsVector(pIn);
		return;
	}
    
	
	/**
	* Set the velocity
	*
	* @param vecIn The new velocity
	*/
	public void setVelocity(PhysicsVector velocityIn)
	{
		velocity = new PhysicsVector(velocityIn);
		return;
	}
	
	/**
	* Set the acceleration
	*
	* @param accIn The new acceleration
	*/
	public void setAcceleration(PhysicsVector accIn)
	{
		acceleration= new PhysicsVector(accIn);
		return;
	}
    
	
	/**
	* Update the position and velocity of the particle subject to a constant acceleration for a time.
	* After the time has passed the acceleration reverts to its previous value.
	*                                      
	* @param deltaTime  The change in time
	* @param accelIn    The applied acceleration
	*/
	public void update(double deltaTime, PhysicsVector accelIn)
	{
		PhysicsVector savedAcceleration = acceleration;
		
		// apply the new acceleration for a short time
		acceleration = new PhysicsVector(accelIn);
		update(deltaTime);
		
		// revert acceleration to previous value 
		acceleration = savedAcceleration;
		
		return;
	}
    
	/**
	* Update the position and velocity of the particle after a short time has
	* passed when the particle is experiencing the acceleration stored in the class.
	* Applies the formula s = ut + 1/2 at**2 to the position
	* Applies the formula v=u+at to the velocity
	* @param deltaTime  The change in time
	*/
	public void update(double deltaTime)
	{
		switch (algorithm){
			
		case nEuler: updateFPA(deltaTime);
			break;
		case nEulerCromer: updateEC(deltaTime);
			break;
			default: updateEulerMidPoint(deltaTime);
		}
		
		return;
		
	}
    
	/**
	* Update the position and velocity of the particle after a short time has
	* passed when the particle is experiencing the acceleration stored in the class.
	* Use the Euler FPA algorithm
	* @param deltaTime  The change in time
	*/
	private void updateFPA(double deltaTime)
	{
		position.increaseBy(PhysicsVector.scale(deltaTime,velocity)); // old position + new velocity * time
		velocity.increaseBy(PhysicsVector.scale(deltaTime,acceleration)); // change the velocity due the acceleration 
		
		return;	
	}	
	
	/**
	* Update the position and velocity of the particle after a short time has
	* passed when the particle is experiencing the acceleration stored in the class.
	* Use the Euler Cromer algorithm
	* @param deltaTime  The change in time
	*/
	private void updateEC(double deltaTime)
	{
		velocity.increaseBy(PhysicsVector.scale(deltaTime,acceleration)); // change the velocity due the acceleration 
		position.increaseBy(PhysicsVector.scale(deltaTime,velocity)); // old position + new velocity * time
		
		return;	
	}
	
	/**
	* Update the position and velocity of the particle after a short time has
	* passed when the particle is experiencing the acceleration stored in the class.
	* Applies the formula s = ut + 1/2 at**2 to the position
	* Applies the formula v=u+at to the velocity
	* @param deltaTime  The change in time
	*/
	private void updateEulerMidPoint(double deltaTime)
	{
		PhysicsVector oldvelocity = velocity;
		velocity.increaseBy(PhysicsVector.scale(deltaTime, acceleration)); // v_(n+1) = v_n + at
		
		PhysicsVector sumVelocity = PhysicsVector.add(oldvelocity, velocity);
		position.increaseBy(PhysicsVector.scale(0.5*deltaTime, sumVelocity)); // old position + 0.5(v_n + v(n+1))t
		
		return;
	}
	
	/**
	* Create a string containing the mass, position, velocity, and acceleration of the particle.
	* This method is called automatically by System.out.println(someparticle)
	* @return string with the format
	* " mass "+mass+" Position: "+position+" Velocity: "+velocity+" Acceleration: "+acceleration
	*/
	public String toString()
	{
		return " mass "+mass+" Position: "+position.returnSimpleString()+" Velocity: "+velocity.returnSimpleString()+" Acceleration: "+acceleration.returnSimpleString();
	}
}


