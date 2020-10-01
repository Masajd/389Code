/**
* A Class to simulate a bunch of protons in a (mostly) uniform magnetic field.
* @author Ian Bailey
* @author Sam Dicker
* @version 2.1
*/

import java.lang.Math;
import java.util.*;
import java.io.*;

public class BunchSim implements NumericalAlgorithms{
	
	public final static double k= 1/(4*Math.PI*8.85418782E-12);// 1/4\pi\epsilon_0
	public final static double c = 3e8; //speed of light
	
	/**
	* Main method to simulate the motion of a single bunch of protons in a uniform B field 
	* or an equivalent Coulomb field.
	*
	*/
	public static void main (String[] args) 
	{
   	 	//Method selector
		Scanner s1= new Scanner(System.in);
		int Method=0;
		int iAlgo=0; //the numerical algorithm to use
		do
		{
			try
			{
				System.out.println("Select Method");
				System.out.println("1. Midpoint Euler" + "\t" + "2.Euler" + "\t" + "3.Euler-Cromer");  //Do loop ensures one of methods is selected
				Method=s1.nextInt();
			}
			catch(InputMismatchException e)
			{
				System.out.println("Invalid Input!");
				System.exit(4);	
			}
		}
		while(Method!=1 && Method!=2 && Method !=3);
		if(Method == 1)
		{
			iAlgo = nEulerMidPoint;
		}
		if(Method == 2)
		{
			iAlgo = nEuler;
		}
		if(Method == 3)
		{
			iAlgo = nEulerCromer;
		}
   	 	
		// Simulation controls
		double timeStep=0;
		do
		{
			try
			{
				System.out.println("Timestep (in seconds i.e 1e-5)");  // time step in seconds
				timeStep = s1.nextDouble();
			}
			catch(InputMismatchException e)
			{
				System.out.println("Invalid Input!");
				System.exit(4);	
			}
		}
		while(timeStep <= 0);
   	    
   	 	long maxTime=3000; // maximum simulation time in seconds
   	 	int maxRev=10000; // maximum number of orbital revolutions for proton bunch
   	 	char rDist='U'; // the random distribution to use (U for uniform, G for Gaussian)
   	 	int pStep=1000; // print output on every pStep time steps
   	 	boolean checkSpread=false; // whether to check for spreadX == spreadY
   	 	double spreadTol=0.1; // the tolerance to use in checking for the spreads to be equal
   	 	boolean useEField=true; // flag to switch between using the uniform B field and an equivalent (in 2d) Coulomb E field.
		
		double eFieldRadius = 0.05;
		boolean doCyclotron =  true; // flag to switch accelerating E field on and off
		boolean CyclotronHasEdge = true; //flag to see if Cyclotron has a finite radius
		double cyclotronRad = 0.1; //Radius of Cyclotr
		
		boolean files = true;     //Whether to print to files or not
		
		boolean WARPSPEEEEEED= false; //check for whether certain speed (as a percentage of c) has been achieved
		double speedLimit = 0.1; // Percentage of c
   	 	
   	 	// E and B parameters
   	 	double bMag = 1.0e-7; // magnetic flux density in Tesla
		double bError = 1.0e-6; // error in b field in Tesla
		boolean bErrorTrue = true;
		
		double eMag = 1.0e-7; // Electric field strength in NC^-1
		double phi =  0.5* Math.PI;     //alteration of the electric field
		double phiCounter = 0;  //What phi increases by
		double phiLimit = 2*Math.PI;  //Limit of phi
		ArrayList<Efficiency> Eff = new ArrayList<Efficiency>(); //Array list of Efficiency objects
   	 	
		// Bunch parameters 
		int nProtons=0;	// Number of protons
		do
		{
			try
			{
				System.out.println("Number of Protons?");
				nProtons = s1.nextInt();
			}
			catch(InputMismatchException e)
			{
				System.out.println("Invalid Input!");
				System.exit(4);	
			}
		}
		while(nProtons<0); 
		
   	 	double pSpeed=0.2; // initial average speed in ms^-1
   	 	double pEnergy = 0.5*Proton.pMass*(pSpeed)*(pSpeed); // initial average energy of particles in bunch (non-relativistic)
   	 	double pESpread= 0.00*pEnergy; // absolute spread in energy of particles in the bunch
   	 	double frequency=(Proton.pCharge*bMag/(2*Math.PI*Proton.pMass)); // expected frequency of bunch orbit
		double radius=Proton.pMass*pSpeed/(Proton.pCharge*bMag); // expected radius of bunch orbit
   	 	
   	 	PhysicsVector pDirn = new PhysicsVector(0,1,0); // direction of bunch at start of simulation
   	 	PhysicsVector pDirnSpread=new PhysicsVector(0,0,0); // relative spread in the direction of the bunch at the start of the simulation
   	 	PhysicsVector pOrigin= new PhysicsVector(); // nominal average position of bunch at start of simulation
   	 	PhysicsVector pSpread = new PhysicsVector(0*radius/100, 0*radius/100, 0); // absolute spread in position of particles in bunch
		
		//Final Values
		double calcKE = Math.pow(Proton.pCharge*bMag*cyclotronRad,2)/(2*Proton.pMass);
		double calcPeriod=1/frequency;
   	 	
		
		String outFileName="Position, Method " + Method + " with timestep " + timeStep + ", efield " + useEField + ", phi value " + phi + ".txt"; //output file name for bunch trajectory data
		String EnergyoutFileName = "Energy, Method " + Method + " with timestep " + timeStep + ", efield " + useEField + ", phi value " + phi + ".txt";  //output file name for the energy of the bunch
		String RandomoutFileName = "Random final values, Method " + Method + " with timestep " + timeStep + ", efield " + useEField + ", phi value " + phi + ".txt"; //file name for random final values and periods
   	 	
		// Open a file to save the bunch positions in
		final PrintWriter outFile;
		final PrintWriter EnergyOutFile;
		final PrintWriter RandomOutFile;
		PrintWriter tryFile=null;
		PrintWriter EnergytryFile = null;
		PrintWriter RandomtryFile = null;			
	
		try
		{
			tryFile = new PrintWriter(outFileName);
			EnergytryFile = new PrintWriter(EnergyoutFileName);
			RandomtryFile = new PrintWriter(RandomoutFileName);
		}
		catch (IOException e)
		{
			System.out.println("Exception opening file: " + e.getMessage());
			System.exit(4);	
		}
		finally
		{
			outFile=tryFile;
			EnergyOutFile = EnergytryFile;
			RandomOutFile = RandomtryFile;
		}
		
   	 	
   	 	//Make a bunch of protons 
   	 	Proton aProton = new Proton();
   	 	Bunch<Particle> pBunch = new Bunch<Particle>();
   	 	
   	 	for (int i=1; i<=nProtons; i++)
		{
   	 		pBunch.addParticle(new Proton(aProton));
   	 	}
   	 	   	 	
   	 	// Set up the simulation
		double time=0.0; 
		int nRev=0; 			// number of orbits the proton bunch completes
		double lastTime=0.0; // The time at which the last 'turn' ended
		double spreadX=0.0;
		double spreadY=0.0;
		OrbitTracker<Particle> bunchOrbit = new OrbitTracker<Particle>(pOrigin); // Start tracking the orbit of the bunch
		
		while(phi<phiLimit)
		{
			time=0;
			nRev=0;
			pBunch.setAlgo(iAlgo);
			pBunch.setDist(rDist);
			pBunch.setPosition(pOrigin,pSpread);
			pBunch.setVelocity(pDirn,pDirnSpread,pEnergy,pESpread);
			System.out.println(pBunch);
			
			while(nRev<maxRev && time < maxTime) // Loop over time
			{
				time+=timeStep;
				double eTimeDependent = eMag* Math.sin((2*Math.PI*frequency*time) + phi); // calculates time dependent Electric field
				
				//Move all particles in the bunch 
				Iterator<Particle> bunchIt = pBunch.iterator();
				while(bunchIt.hasNext())
				{
					Particle aParticle=bunchIt.next();
					PhysicsVector acceleration = new PhysicsVector();
					// Loop over all fields
				
					// Set up the cyclotron fields, and magnetic field
					ArrayList<GeneralEMField> cyclotron = new ArrayList<GeneralEMField>();	
					cyclotron.add(addField(bMag,eTimeDependent,eMag, aParticle, useEField, doCyclotron,radius, eFieldRadius, bError, bErrorTrue));
					
					for (GeneralEMField field : cyclotron)
					{					
						if (aParticle instanceof ChargedParticle)
						{
							acceleration.increaseBy(field.getAcceleration((ChargedParticle)aParticle));
						}
					
					}
					aParticle.update(timeStep, acceleration);
				}
			
				if (((int)(time/timeStep))%pStep==1 && files) // Write out the bunch position at intervals
				{	 
					outFile.println(time + " " + (pBunch.getPosition()).returnSimpleString()); 
					outFile.flush();
				}
			
				if(bunchOrbit.MiddleOrbit(pBunch, radius) && files) //writes to energy file whenever protons have moved outside the electric field
				{
					EnergyOutFile.println(pBunch.getTotalKE() + "\t" + (pBunch.getPosition()).returnSimpleString());
					EnergyOutFile.flush();
				}
			
				if (bunchOrbit.hasOrbited(pBunch))
				{
					// bunch has completed an orbit
					nRev+=1;
					
						System.out.printf("Revolution number %3d at time %10.6f s\n", nRev, time);
						System.out.printf("Period of this revolution is  %10.6f s\n",(time-lastTime)); 
					if(files)
					{
						RandomOutFile.println((time-lastTime));
						RandomOutFile.flush();
					}
						System.out.println(pBunch);
						System.out.println();
					
					lastTime=time;
					PhysicsVector spread= pBunch.getSpreadMax();
					spreadX=spread.getX();
					spreadY=spread.getY();
	
					if (checkSpread && spreadTest(spreadX, spreadY,spreadTol))// Check the spreads and quit if they are equal within tolerance.
					{	 
						break;
					}
				}	
			
				double radius2=Proton.pMass*(pBunch.getVelocity()).magnitude()/(Proton.pCharge*bMag);
				//double radius3=(PhysicsVector.subtract(pBunch.getPosition(),centre)).magnitude();
				if(bunchOrbit.hasReachedEdge(pBunch, cyclotronRad, radius2) && CyclotronHasEdge) //Check whether the protons have reached the last orbit in the cyclotron
				{
					break;				
				}
				
				if(WARPSPEEEEEED && pBunch.getVelocity().magnitude() > speedLimit*c)
				{
					break;
				}
			}                          
		
			// After simulation, calculate the average radius, KE and periodic orbit 
			double simPeriod=time/nRev;
			double radius2=Proton.pMass*(pBunch.getVelocity()).magnitude()/(Proton.pCharge*bMag);
			double simKE = Math.pow(Proton.pCharge*bMag*radius2,2)/(2*Proton.pMass);
			if(files)
			{
				RandomOutFile.println(simPeriod + "\t" + radius2);
				RandomOutFile.println(pBunch.getTotalKE()/nProtons + "\t" + simKE);
			}

			Eff.add(new Efficiency(phi, nRev, simPeriod, radius2, simKE));
			if(files || phiCounter==0)
			{
				break;
			}
			
			phi+=phiCounter;
		}
		
		System.out.println("\n\n");
		System.out.println("Calculated period of orbit " + calcPeriod + " s");
		System.out.println("Initial radius of orbit (calculated) " + radius + " m");
		System.out.println("Final kinetic energy of orbit (calculated) " + calcKE + " J");
		
		getSmallestRev(Eff, maxRev);
		
		outFile.close(); 
		RandomOutFile.close();
		EnergyOutFile.close();	
	}
	
	/**
	* Method to test whether spreads in X and Y are equal within tolerance
	*
	* @param spreadX 	Bunch spread in X
	* @param spreadY	Bunch spread in Y
	* @param spreadTol	relative difference in spreadX and spreadY
	* @return true if spreads are equal within tolerance. False otherwise.
	*/
	public static boolean spreadTest(double spreadX, double spreadY, double spreadTol)
	{
		return ((Math.abs(spreadX-spreadY)<=(spreadTol*spreadX))&&
			(Math.abs(spreadX-spreadY)<=(spreadTol*spreadY)));
	}
	
	
	
	/**
	* Method to calculate the total potential energy of a bunch of particles
	* in the fields of an accelerator
	*
	* @param acceleratorFields the e/m fields of the accelerator
	* @param aBunch the bunch of particles
	* @return the potential energy in J
	*/
	public static double totalPE(ArrayList<GeneralEMField> acceleratorFields, Bunch<Particle> aBunch)
	{
		Iterator<Particle> bunchIt = aBunch.iterator();
		double potE=0.0;
		while(bunchIt.hasNext())
		{
			Particle aParticle=bunchIt.next();
			// Loop over all fields
			for (GeneralEMField field : acceleratorFields){
				if (aParticle instanceof ChargedParticle){
					potE+=(field.getPotentialE((ChargedParticle)aParticle));
				}
				
			}
		}
		
		return potE;
	}
	
	
	/**
	* Method creates an electric and magnetic field depending on whether the electric field is being used, or if it is a cyclotron.
	*
	* @param bField the magnetic flux density
	* @param eField1 the time dependent electric field strength
	* @param eField2 the static electric field strength
	* @param aParticle a particle who's poistion will be used to turn the cyclotron electric field on
	* @param useEField boolean to show whether an electric field will be turned on
	* @param doCyclotron boolean to show whether simulation is of a cyclotron
	* @param radius radius of a charged particle moving just in the magnetic field
	* @return The EMfield present at the particles location
	**/
	public static GeneralEMField addField(double bField, double eField1, double eField2, Particle aParticle, boolean useEField, boolean doCyclotron, double radius, double eFieldRadius, double bError, boolean bErrorTrue) 
	{
		GeneralEMField bendingField;
		
		bendingField = new EMField(); 
		if(bErrorTrue)
		{
			((EMField)bendingField).setMagnetic(new PhysicsVector(0,0,(bField-bError+(Math.random()*2*bError))));
		}
		else
		{
			((EMField)bendingField).setMagnetic(new PhysicsVector(0,0,bField));
		}
		
		if(useEField)
		{					
			if(doCyclotron)
			{
				if(aParticle.getPosition().getY() >= (-eFieldRadius*radius) && aParticle.getPosition().getY() <= (eFieldRadius*radius))
				{
					((EMField)bendingField).setElectric(new PhysicsVector(0, eField1, 0));
				}
				else
				{
					((EMField)bendingField).setElectric(new PhysicsVector());
				}
			}
			else
			{
				((EMField)bendingField).setElectric(new PhysicsVector(0, eField2, 0));
			}
		}
		else
		{
			((EMField)bendingField).setElectric(new PhysicsVector()); 
		}
		
		return bendingField;
	}
	
	/**
	* Finds the smallest Revolution when applying phi condition, 
	* then prints the smallest rev and its associated values.
	*
	* @param List An arraylist of efficiency objects
	* @param nRevLimit the max number of revs
	*/
	public static void getSmallestRev(ArrayList<Efficiency> List, int nRevLimit)
	{
		int maxRev=nRevLimit;
		double maxPhi=0;
		double maxPeriod=0;
		double maxRad=0;
		double maxKE=0;
		
		for(int i=0; i< List.size(); i++)
		{
			int S = List.get(i).getRev();
			
			if( S < maxRev)
			{
				maxRev = S;
				maxPhi = List.get(i).getPhi();
				maxPeriod = List.get(i).getSimPeriod();
				maxRad = List.get(i).getFinalRad();
				maxKE = List.get(i).getSimKE();
				
			}		
		}
		
		Efficiency maxEff = new Efficiency(maxPhi, maxRev, maxPeriod, maxRad, maxKE);
		maxEff.printHere();
	}
}
