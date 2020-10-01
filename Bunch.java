import java.util.*;
/**
* A Class to represent a bunch of particles.
* This is a generic class of Particle objects
* 
* @author Ian Bailey
* @version 1.3
*/
public  class Bunch <T extends Particle> implements Iterable{
	
	
	private ArrayList<T> particles; // the class of particles in the bunch
	private char rDist='U'; // the random distribution to use for generating random vectors, etc
	
	/**
	* Constructor for the Bunch of Particles
	* 
	*/
	public Bunch(){
		particles= new ArrayList<T>();
		
	}
	
	/**
	* add a particle to the bunch of particles
	* @param particle the particle to add to the bunch
	* 
	*/
	public void addParticle(T particle){
		particles.add(particle);
		
	}
	
	/**
	* Returns the iterator for the collection of particles
	* @return iterator of the bunch
	*/
    @Override
    public Iterator<T> iterator() {
        return particles.iterator();
    }
    
	
    /**
	* set the numerical algorithm inside all particles in the bunch
	* @param the numerical algorithm as defined in the NumericalAlgorithms interface
	*/
	public void setAlgo(int algoIn){
		for (T aParticle : this.particles){
			aParticle.setAlgo(algoIn);
		}
	}
	
	/**
	* set the distribution from which to draw random numbers
	* @param rDist the distribution to use ('U' => uniform, 'G' => Gaussian)
	*/
	public void setDist(char rDist){
		this.rDist=rDist;
	}
	
	/**
	* get the distribution from which random numbers are drawn
	* @return the distribution being used ('U' => uniform, 'G' => Gaussian)
	*/
	public char getDist(){
		return this.rDist;
	}
	
	/**
	* setPosition
	* @param positionIn the nominal position of the bunch
	* @param spreadIn the spread around the nominal position of the bunch
	*/
	public void setPosition(PhysicsVector positionIn, PhysicsVector spreadIn){
		for (T aParticle : this.particles){
			aParticle.setPosition(randomVector3D(positionIn,spreadIn,rDist));
		}
	}
	
	
	/**
	* setVelocity
	* @param velocityIn the nominal velocity direction of the bunch
	* @param spreadIn the relative spread in the direction of the bunch
	* @param energy the nominal energy of the bunch
	* @param energySpread the absolute spread in the energy
	*/
	public void setVelocity(PhysicsVector velocityIn, PhysicsVector spreadIn, double energy, double energySpread){
		
		for (T aParticle : this.particles){
			double speed = Math.sqrt(2.0*randomScalar(energy,energySpread,rDist)/aParticle.getMass());
			PhysicsVector velocity = PhysicsVector.scale(speed,randomVector3D((velocityIn.getUnitVector()),spreadIn,rDist));
			aParticle.setVelocity(velocity);
		}
	}
	
	/** 
	* Return a random vector with a specified average and maximum spread. 
	* A vector of random numbers uniformally distributed along the line average +/- spread  
	* @param average the centre of the line segment from which the random vectors are drawn
	* @param spread half the line segment length
	* @return the random vector.
	*/
	private PhysicsVector randomVector(PhysicsVector average, PhysicsVector spread){
		return PhysicsVector.add(average,PhysicsVector.scale(randomScalar(0.0,1.0), spread));
	}
	
	
	/** 
	* Return a random vector with a specified average and  spread. 
	* @param average mean of the distribution from which the vector is to be drawn
	* @param spread the spread of the distribution from which the vector is to be drawn
	* @param dist 'U' => uniform distribution, 'G' => Gaussian distribution, return 0 otherwise
	* @return the random vector.
	*/
	private PhysicsVector randomVector3D(PhysicsVector average, PhysicsVector spread, char dist){
		if (dist=='U'){	
			return randomVector3D(average, spread);
		}
		else if(dist=='G'){
			return randomGVector3D(average, spread);
		}
		else{
			return new PhysicsVector();
		}
	}
	
	
	/** 
	* Return a random vector with a specified average and maximum spread. 
	* A vector of random numbers distributed between average +/- spread  
	* @param average the centre of the cube from which the random vectors are drawn
	* @param spread half the cube edge length
	* @return the random vector.
	*/
	private PhysicsVector randomVector3D(PhysicsVector average, PhysicsVector spread){
		return new PhysicsVector(randomScalar(average.getX(),spread.getX()),
			randomScalar(average.getY(),spread.getY()),
			randomScalar(average.getZ(),spread.getZ()));
	}
	
	/** 
	* Return a random vector with a specified average and s.d. with a Gaussian distribution  
	* @param average the mean of the 3d Gaussian distribution
	* @param sd the standard deviation of the 3d Gaussian distribution
	* @return the random vector.
	*/
	private PhysicsVector randomGVector3D(PhysicsVector average, PhysicsVector spread){
		return new PhysicsVector(randomScalar(average.getX(),spread.getX(),'G'),
			randomScalar(average.getY(),spread.getY(),'G'),
			randomScalar(average.getZ(),spread.getZ(),'G'));
	}
	
	/** 
	* Return a random scalar with a specified mean and range.
	* I.e. the number is drawn from a uniform distribution with minimum and maximum values(av+/- spread)
	* @param average the centre of the uniform distribution from which the random number is drawn 
	* @param spread the half spread of the uniform distribution from which the random number is drawn
	* @return the random scalar.
	*/
	private double randomScalar(double average, double spread){
		return average+(2*Math.random()-1.0)*spread;
	}
	
	/** 
	* Return a random scalar with a specified mean and standard deviation
	* The number is drawn from a Gaussian distribution.
	* @param average the centre of the Gaussian distribution from which the random number is drawn 
	* @param sd the standard deviation of the Gaussian distribution
	* @return the random scalar.
	*/
	private double randomGScalar(double average, double sd){
		Random aRandom = new Random();
		return average+aRandom.nextGaussian()*sd;
		
	}
	
	/** 
	* Return a  random scalar with a specified average and spread.
	* @param average the centre of the distribution from which the random number is drawn 
	* @param spread the half spread of the distribution from which the random number is drawn
	* @param dist the distribution to use 'U' => uniform, and 'G' => Gaussian (return 0 otherwise)
	* @return the random scalar.
	*/
	private double randomScalar(double average, double spread, char dist){
		if (dist=='U'){
			return randomScalar(average, spread);
		}
		else if (dist=='G'){
			return randomGScalar(average, spread);
		}
		else{
			return 0;
		}	
	}
	
	
	/** 
	* Find the average position of the particles in the bunch.
	* @return the average position of particles in the bunch.
	*/
	public PhysicsVector getPosition(){
		PhysicsVector position = new PhysicsVector();
		for (T aParticle : this.particles){
			position.increaseBy(aParticle.getPosition());
		}
		return PhysicsVector.scale(1.0/this.particles.size(),position);
	}
	
	
	/** 
	* Find the average velocity of the particles in the bunch.
	* @return the average velocity of particles in the bunch.
	*/
	public PhysicsVector getVelocity(){
		PhysicsVector velocity = new PhysicsVector();
		for (T aParticle : this.particles){
			velocity.increaseBy(aParticle.getVelocity());
		}
		return PhysicsVector.scale(1.0/this.particles.size(),velocity);
	}
	
	/** 
	* Find the total kinetic energy of the bunch
	* @return the total kinetic energy of the bunch in J.
	*/
	public double getTotalKE(){
		double kinE=0.0;
		for (T aParticle : this.particles){
			kinE+=aParticle.getKE();
		}
		return kinE;
	}
	
	/** 
	* Find the standard deviation of the bunch distribution in x, y and z.
	* @return the standard deviation in x, y and z.
	*/
	public PhysicsVector getSpreadSD(){
		PhysicsVector average = new PhysicsVector(); 
	 	PhysicsVector sqAverage = new PhysicsVector();
	 	for (T aParticle: this.particles){
	 		PhysicsVector position=aParticle.getPosition();
	 		average.increaseBy(position);
	 		sqAverage.increaseBy(PhysicsVector.multiply(position,position));
	 	}
	 	average.scale(1.0/this.particles.size());
	 	sqAverage.scale(1.0/this.particles.size());
	 	PhysicsVector variance=PhysicsVector.subtract(sqAverage,PhysicsVector.multiply(average,average));
	 	return new PhysicsVector(Math.sqrt(variance.getX()), Math.sqrt(variance.getY()), Math.sqrt(variance.getZ()));
	}
	
	/** 
	* Find the full spread of the bunch in x, y and z
	* I.e. find the particles which is furthest away from the average position of the bunch in the +ve q direction
	* and the particle which is furthest away from the average position of the bunch in the -ve q direction and take
	* the difference.
	* Return this distance where q=x,y,z
	* @return the maximum spread in x, y and z.
	*/
	public PhysicsVector getFullSpread(){
		PhysicsVector av = this.getPosition();
		PhysicsVector max = new PhysicsVector(av);
		PhysicsVector min = new PhysicsVector(av);
		for (T aParticle: this.particles){
			PhysicsVector pos = aParticle.getPosition();
			if (pos.getX()>max.getX()) max.setX(pos.getX());
			if (pos.getY()>max.getY()) max.setY(pos.getY());
			if (pos.getZ()>max.getZ()) max.setZ(pos.getZ());
			if (pos.getX()<min.getX()) min.setX(pos.getX());
			if (pos.getY()<min.getY()) min.setY(pos.getY());
			if (pos.getZ()<min.getZ()) min.setZ(pos.getZ());
		}
		
		return(PhysicsVector.subtract(max,min));
	}
	
	/** 
	* Find the maximum spread of the bunch in x, y and z (furthest deviations from the average position of the bunch)
	* I.e. find the particle which is furthest away from the average position of the bunch in the q direction 
	* and return this distance where q=x,y,z
	* @return the maximum spread in x, y and z.
	*/
	public PhysicsVector getSpreadMax(){
		PhysicsVector average = this.getPosition();
		PhysicsVector max = new PhysicsVector();
		for (T aParticle: this.particles){
			PhysicsVector diff=PhysicsVector.subtract(aParticle.getPosition(),average);
			if (Math.abs(diff.getX())>Math.abs(max.getX())) max.setX(Math.abs(diff.getX()));
			if (Math.abs(diff.getY())>Math.abs(max.getY())) max.setY(Math.abs(diff.getY()));
			if (Math.abs(diff.getZ())>Math.abs(max.getZ())) max.setZ(Math.abs(diff.getZ()));
		}
		
		return(max);
	}
	
	
	/** 
	* Find the maximum spread of the bunch in v_x, v_y and v_z.
	* I.e. find the particle whose speed is the furthest from the average speed of the bunch in the q direction 
	* and return this speed where q=x,y,z
	* @return the maximum spread in v_x, v_y and v_z.
	*/
	public PhysicsVector getVSpreadMax(){
		PhysicsVector average = this.getVelocity();
		PhysicsVector max = new PhysicsVector();
		for (T aParticle: this.particles){
			PhysicsVector diff=PhysicsVector.subtract(aParticle.getVelocity(),average);
			if (Math.abs(diff.getX())>Math.abs(max.getX())) max.setX(Math.abs(diff.getX()));
			if (Math.abs(diff.getY())>Math.abs(max.getY())) max.setY(Math.abs(diff.getY()));
			if (Math.abs(diff.getZ())>Math.abs(max.getZ())) max.setZ(Math.abs(diff.getZ()));
		}
		
		return(max);
	}
	
	/** 
	* Find the particle in the bunch which is furthest from the y axis.
	* This is not a very useful method, but it does illustrate a use for ordered sets!
	* @return the position in y of the particle furthest from the y axis.
	*/
	public double getSpreadY(){
		TreeSet<T> myTree = new TreeSet<T>(this.particles);
		
		return( ((myTree.last()).getPosition()).getY());
	}
	
	/**
	* Create a string containing the mass, charge, position, velocity, and acceleration of the particle.
	* This method is called automatically by System.out.println(someparticle)
	* @return string with the format
	*/
	@Override
	public String toString()
	{
		return String.format(
			"Number of particles: %d \n"+ 
			"Bunch position (m):                 %s \n"+ 
			"Bunch velocity (m/s):               %s \n"+ 
			"Bunch spread (sd in m):             %s\n"+
			"Bunch spread (max in m):            %s\n"+
			"Bunch spread (full in m):           %s\n"+
			"Bunch spread (y from TreeSet in m): %10.6f\n"+
			"Bunch v spread (max in m/s):        %s\n",
			this.particles.size(),
			(getPosition().formatString()),
			(getVelocity().formatString()),
			(getSpreadSD().formatString()),
			(getSpreadMax().formatString()),
			(getFullSpread().formatString()),
			(getSpreadY()),
			(getVSpreadMax().formatString()));
		
	}
	
	
 	
}


