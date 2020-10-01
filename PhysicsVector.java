/**
 * A 3-dimensional vector class to be used in the development of
 * computer simulations of various physical systems.
 *
 * 
 * The 3-dimensional vector represents a standard cartesian vector
 * given by <b> xi + yj +zk</b> where <b>i</b>, <b>j</b>, and <b>k</b> are unit
 * vectors in the x, y, and z direction respectiveley.
 *
 * @author Iain A. Bertram
 * @version 1.2
 **/
public class PhysicsVector {
	// Fix the dimension of the array representing the vectors 
	private static final int  vectorSize=3;
	
	// In this case we have a three dimensional vector, the x component is [0] with y[1] and z [2]
	private double[]  vectorComponents = new double[vectorSize];
	
    /**
     * Default contructor that creates a PhysicsVector with zero
     * magnitude
     **/
    public PhysicsVector(){
	for (int i=0; i<vectorComponents.length; i++) {
		vectorComponents[i] =0.;
	}
    }
	

    /**
     * Contructor that takes two arguments to create a PhysicsVector
     * with the first argument x giving the component in the <b>i</b>
     * direction, the second argument y giving the component in the
     * <b>j</b> direction, the z component is set to zero.
     **/
    public PhysicsVector(double x, double y){
	vectorComponents[0] = x;
	vectorComponents[1] = y;
	vectorComponents[2] = 0.;
    }


    public PhysicsVector(double x, double y, double z) {
        vectorComponents[0] = x;
        vectorComponents[1] = y;
        vectorComponents[2] = z;
    }
	
	
    /**
     * Contructor that takes an array of length 3 and sets the components of the vector 
     * to the values stored in the array. If array length 2 it sets the x and y components.
     * Returns a default vector if the array is not 
     * length 2 or 3 and prints a warning message.      
     **/
    public PhysicsVector(double [] x){
    	    if (x.length == vectorSize ) {
		for (int i=0; i<vectorComponents.length; i++) {
			vectorComponents[i] = x[i];			
		}
	    }
	    else if (x.length == vectorSize-1 ) {
	    	for (int i=0; i<x.length; i++) {
	    		vectorComponents[i] = x[i];
		}
		vectorComponents[vectorComponents.length-1] = 0.;
	    }
	    else {
	    	this.setVector(new PhysicsVector());
		System.out.println(" WARNING: PhysicsVector(double [] x) requires an array of length " + vectorSize);
	    }
    }
	
	
    /**
     * Constructor that creates a vector with same value as the
     * PhysicsVector passed as an argument.
     **/
    public PhysicsVector(PhysicsVector v){
    // Note: need to copy values not the array itself. 
    	for (int i=0; i<vectorComponents.length; i++) {
    	    vectorComponents[i] = v.vectorComponents[i];
    	}
		
    }
    
    
       /**
    * return the x-component of the vector
    *
    * @return x-component
    **/
    public double getX(){
    	    return vectorComponents[0];
    }

    /**
    * set the x-component of the vector
    *
    * @param x  x-component
    **/
    public void setX(double x) {
    	    vectorComponents[0] = x;
    }
	
	
    /**
     * return the y-component of the vector
     *
     * @return y-component
     **/
    public double getY(){
    	    return vectorComponents[1];
    }

    /**
    * set the y-component of the vector
    *
    * @param y  y-component
    **/ 
    public void setY(double y) {
        vectorComponents[1] = y;
    }
	
    /**
     * return the z-component of the vector
     *
     * @return z-component
     **/
    public double getZ(){
	return vectorComponents[2];
    }

    /**
    * set the z-component of the vector
    *
    * @param z  z-component
    **/
    public void setZ(double z) {
        vectorComponents[2] = z;
    }
	
    /**
     * Method to return a String which represents the vector and has
     * the form <b> xi + yj</b>.
     * @return String <b> xi + yj</b>
     **/
    public String return2DString(){
	String text = ""+vectorComponents[0]+"i ";
	if (vectorComponents[1] < 0){
		text += " -"+(-1*vectorComponents[1]);
	}else{
		text += " +"+(vectorComponents[1]);
	}
	text+="j ";
	return text;
    }
	
    /**
     * Method to return a String which represents the vector and has
     * the form <b> x y</b> (I.e. separated with spaces).
     * @return String <b>x y</b>
     **/
    public String returnSimple2DString(){
	String text = "";
	text += vectorComponents[0] +" ";
	text += vectorComponents[1] +" ";
	return text;
    }
	
	
    /**
     * Method to return a String which represents the vector and has
     * the form <b> xi + yj +zk</b>.
     * @return String <b> xi + yj + zk</b>
     **/
    public String returnString(){
	String text = ""+vectorComponents[0]+"i ";
	if (vectorComponents[1] < 0){
		text += " -"+(-1*vectorComponents[1]);
	}else{
		text += " +"+(vectorComponents[1]);
	}
	text+="j ";
	if (vectorComponents[2] < 0){
		text += " -"+(-1*vectorComponents[2]);
	}else{
		text += " +"+(vectorComponents[2]);
	}
		
	text+="k ";
	return text;
    }

    /**
     * Method to return a String which represents the vector and has
     * the form <b>x y z</b> (I.e. only seperated by spaces).
     * @return String <b>"x y z"</b>
     **/
    public String returnSimpleString(){
	String text = "";
	for (int i=0; i<vectorComponents.length; i++) {
		text+=vectorComponents[i]+" ";
	}
	return text;
    }
	
	
    /**
     * prints the vector to standard out using the
     * <code>System.out.println</code> command. Note: it will include
     * a carriage return. Prints: <b> xi + yj</b>
     *
     **/
    public void print2D(){
	String text = this.return2DString();
	System.out.println(text);
    }
	
	
    /**
     * prints the vector to standard out using the
     * <code>System.out.println</code> command. Note: it will include
     * a carriage return. Prints: <b> xi + yj + zk</b>
     *
     **/
    public void print(){
	String text = this.returnString();
	System.out.println(text);
    }
	
      /**
     * Method to return a formatted String which represents the vector and has
     * the form <b>x y z</b> (I.e. only seperated by spaces).
     * @return String <b>"x y z"</b>
     **/
    public String formatString(){
        String text = "";
        for (int i=0; i<vectorComponents.length; i++) {
            text+=String.format("%10.6f ", vectorComponents[i]);
        }
        return new String(text);
    }
	
    /**
     * Change the value of the vector to <b> xi + yj + zk</b>
     *
     * @param x <code>x-component</code> of the vector
     * @param y <code>y-component</code> of the vector
     * @param z <code>z-component</code> of the vector	
     **/
    public void setVector(double x, double y, double z){
	vectorComponents[0] = x;
	vectorComponents[1] = y;
	vectorComponents[2] = z;
		
    }
	
    /**
     * Change the value of the vector to <b> xi + yj</b>. The z-component is set to zero. 
     *
     * @param x <code>x-component</code> of the vector
     * @param y <code>y-component</code> of the vector
     **/
    public void setVector(double x, double y){
	vectorComponents[0] = x;
	vectorComponents[1] = y;
	vectorComponents[2] = 0.;
		
    }
	
	
    /**
     * Modify the vector so that it is equal to vector <b>v</b>
     *
     * @param v the vector that the used to set the values of the
     * modified vector.
     **/
    public void setVector(PhysicsVector v){
	// Note: need to copy values not the array itself. 
	for (int i=0; i<vectorComponents.length; i++) {
		vectorComponents[i] = v.vectorComponents[i];
	}
		
    }
	
    /**
     * Add a vector <b>v</b> to the original vector. Normal vector
     * addition is carried out. I.e. the x-components are added and
     * the y components are added.
     *
     * @param v vector to be added to original vector. 
     **/
    public void increaseBy(PhysicsVector v) {
	for (int i=0; i<vectorComponents.length; i++) {
		vectorComponents[i] += v.vectorComponents[i];
	}
    }
    
   /**
     * Subtract a vector <b>v</b> from the original vector. 
     * @param v vector to be subtracted from original vector. 
     **/
    public void decreaseBy(PhysicsVector v) {
        for (int i = 0; i < vectorComponents.length; i++) {
            vectorComponents[i] -= v.vectorComponents[i];
        }
    }
    
    /**
     * returns the magnitude of the vector. If the vector is given by
     * <b>xi +yj +zk</b> Then the magnitude is Math.sqrt(x*x + y*y + zz)
     *
     * @return a scalar with the magnitude
     **/
    public double magnitude(){
	double mag = dot(this,this);
	mag = Math.sqrt(mag);
	return mag;
    }
			
    /**
     * returns a vector with magnitude of one and the original
     * direction of the original vector. I.e. scale the vector by
     * 1/magnitude.
     *
     * If vector is a zero vector then returns zero length vector
     *
     * @return unit vector with direction of original vector
     **/
    public PhysicsVector getUnitVector(){
	PhysicsVector unit = new PhysicsVector(this);
	double magnitude = this.magnitude();
	if (Math.abs(magnitude) > 1.e-34){
		unit.scale(1/magnitude);}
	else {
		unit = new PhysicsVector();
	}
	return unit;
    }
	
    /**
     * scale the components of the vector by the scalar. I.e. If the
     * scalar is s then the vector transfoms to <b>s*xComponent i +
     * s*yComponent j</b>.
     *
     * @param x the factor to scale the vector by.  
     **/
    public void scale(double x){
	for (int i=0; i<vectorComponents.length; i++) {
		vectorComponents[i] *= x;
	}	
    }
	
	
    /**
     * static scale the components of the vector by the scalar. I.e. If the
     * scalar is s then the vector transfoms to <b>s*xComponent i +
     * s*yComponent j</b>.
     *
     * @param x the factor to scale the vector by.  
     * @param v the vector to scale
     **/
    public static PhysicsVector scale(double x, PhysicsVector v){
	PhysicsVector scaled = new PhysicsVector(v);
	for (int i=0; i<scaled.vectorComponents.length; i++) {
		scaled.vectorComponents[i] *=x;
	}
	return scaled;
    }
	
	
    /**
     * dot is a standard vector scalar product. E.g. If <b> v = xi + yj</b>
     * and <b>u = wi + zy</b>, then the dot product returns <b>xw +
     * yz</b>
     *
     * @param v first vector in product
     * @param u second vector in product
     * @return return scalar product
     **/
    public static double dot(PhysicsVector v, PhysicsVector u){
	double product = 0d;
        for (int i = 0; i < v.vectorComponents.length; i++) {
            product += v.vectorComponents[i] * u.vectorComponents[i];
        }
        return product;
    }
   /**
     * cross performs a standard vector cross product.
     *
     * @param v first vector in product
     * @param u second vector in product
     * @return return cross product
     **/
    
    public static PhysicsVector cross(PhysicsVector v, PhysicsVector u) {
        return new PhysicsVector(
            v.vectorComponents[1] * u.vectorComponents[2] - v.vectorComponents[2] * u.vectorComponents[1],
            v.vectorComponents[2] * u.vectorComponents[0] - v.vectorComponents[0] * u.vectorComponents[2],
            v.vectorComponents[0] * u.vectorComponents[1] - v.vectorComponents[1] * u.vectorComponents[0]);
  
    }
	
    /**
     * standard vector addition. If <b> v = xi + yj</b>
     * and <b>u = wi + zy</b>, then the addition returns a vector
     * <b>(x+w)i + (y+z)j</b>
     *
     * @param v first vector in sum
     * @param u second vector in sum
     * @return return summed vector
     **/
    public static PhysicsVector add(PhysicsVector v, PhysicsVector u){
	PhysicsVector sum = new PhysicsVector(v);
	sum.increaseBy(u);
	return sum;
    }
	
    /**
     * vector subtraction. If <b> v = xi + yj</b> and <b>u = wi +
     * zy</b>, then the subtraction method returns a vector <b>(x-w)i +
     * (y-z)j</b>
     *
     * @param v first vector
     * @param u vector to subtract from first vector
     * @return return v-u vector
     **/
    public static PhysicsVector subtract(PhysicsVector v, PhysicsVector u){
	PhysicsVector diff = new PhysicsVector(v);
	diff.decreaseBy(u);
	return diff;
    }
	
     /**
     * component by component multiplication of two vectors. If <b> v = xi + yj</b> and <b>u = wi +
     * zy</b>, then the multiplication method returns a vector <b>(x*w)i +
     * (y*z)j</b>
     *
     * @param v first vector
     * @param u vector to 'multiply' by first vector
     * @return return v_i*u_i vector
     **/
    public static PhysicsVector multiply(PhysicsVector v, PhysicsVector u){
	PhysicsVector mult = new PhysicsVector(u.getX()*v.getX(), u.getY()*v.getY(), u.getZ()*v.getZ());
	return mult;
    }
    
    /** 
     * returns true if the vector components have equal length to
     * within 1<sup>-34</sup> (chosen since exact equivalence for
     * doubles is difficult to obtain exact equivalence. <par>
     *
     * If <b> v = xi + yj</b> and <b>u = wi + zy</b>, then return true
     * if <b>abs(x-w) &lt; 1<sup>-34</sup></b> and <b>abs(y-z) &lt;
     * 1<sup>-34</sup></b>
     *
     * @param v first vector 
     * @param u second vector 
     * @return return true if vectors are sufficiently equal
     * 
     **/
    public static boolean equals(PhysicsVector v, PhysicsVector u){
	boolean equal = true;
	for (int i=0; i<v.vectorComponents.length; i++) {
		if (Math.abs(v.vectorComponents[i]-u.vectorComponents[i]) >=1.e-34) equal = false;
	}
	return equal;
    }
    
    /**
     * Method to return a String which represents the vector and has
     * the form <b> xi + yj +zk</b>.
     * @return String <b> xi + yj + zk</b>
     **/
    public String toString(){
	
	return returnString();
    }
	
}


