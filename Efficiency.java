/**
* A class to represent the efficiency of the electric field
* SI units used throughout.
*
* @author Sam Dicker
  @version 1.0
**/
public class Efficiency
{
	private double phi;
	private int rev;
	private double simPeriod;
	private double finalRad;
	private double simKE;
	/**
	*Default efficiency constructor
	**/
	public Efficiency()
	{
		phi=0;
		rev=0;
		simPeriod=0;
		finalRad=0;
		simKE=0;
	}
	 /**
	 * Constructor with 2 inputs - the Phi value and the revolution number
	 * @param phiIn the new Phi value
	 * @param revIn the new Revolution value
	 **/
	public Efficiency(double phiIn, int revIn, double simPeriodIn, double finalRadIn, double simKEIn)
	{
		phi = phiIn;
		rev = revIn;
		simPeriod = simPeriodIn;
		finalRad = finalRadIn;
		simKE = simKEIn;
	}
	
	/**
	* Return the phi value
	* @return phi
	*/
	public double getPhi()
	{
		return phi;
	}

	/**
	* Return Revolution Number
	* @return rev
	*/
	public int getRev()
	{
		return rev;
	}
	
	/**
	* Return simulation period
	* @return simPeriod
	*/
	public double getSimPeriod()
	{
		return simPeriod;
	}
	
	/**
	* Return the final radius
	* @return finalRad
	*/
	public double getFinalRad()
	{
		return finalRad;
	}
	
	/**
	* Return the simulated KE
	* @return simKE
	*/
	public double getSimKE()
	{
		return simKE;
	}
	/**
	* Set the Phi value of the object
	* @param phiIn new Phi value
	*/
	public void setPhi(double phiIn)
	{
		phi = phiIn;
	}
	
	/**
	* Set the Revolution value of the object
	* @param revIn new Revolution value
	*/
	public void setRev(int revIn)
	{
		rev= revIn;
	}
	
	/**
	* Set the simulation period of the object
	* @param simPeriodIn new simulation period
	*/
	public void setSimPeriod(double simPeriodIn)
	{
		simPeriod = simPeriodIn;
	}
	
	/**
	* Set the final radius of the object
	* @param finalRadIn new final radius
	*/
	public void setFinalRad(double finalRadIn)
	{
		finalRad= finalRadIn;
	}
	
	/**
	* Set the simulated KE of the object
	* @param simKEIn new kinetic energy
	*/
	public void setSimKE(double simKEIn)
	{
		simKE=simKEIn;
	}
	
	/**
	*Print values of the Efficiency object to the screen.
	*/
	public void printHere()
	{
		System.out.println( this.getPhi() + ", Revolution of " + this.getRev());
		System.out.println("Simulated period of " + this.getSimPeriod() + "s");
		System.out.println("Final radius of " + this.getFinalRad() + "m");
		System.out.println("Final simulated KE of orbit " + this.getSimKE() + "J"); 
	}
}
