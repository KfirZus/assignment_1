package assignments.Ex1;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {

            if (lx == 2){
               double x1 = xx[0], y1 = yy[0];

               double x2 = xx[1], y2 = yy[1];

               double m = (y2 - y1) / (x2 - x1);
               double b = y1 - m * x1;
               ans = new double[]{m, b};
            }
            else{
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                double x3 = xx[2], y3 = yy[2];

                double denom = (x1 - x2)*(x1 - x3)*(x2 - x3);

                double a = (x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3 - y2)) / denom;
                double b = (x3*x3*(y1 - y2) + x1*x1*(y2 - y3) + x2*x2*(y3 - y1)) / denom;
                double c = (x2*x3*(x2 - x3)*y1 + x3*x1*(x3 - x1)*y2 + x1*x2*(x1 - x2)*y3) / denom;

                ans = new double[]{a,b, c};

            }
		}

		return ans;
	}

	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
        if (p1 == null || p2 == null) return false;

        int n = Math.max(p1.length, p2.length) - 1;

        for (int x = 0; x <= n; x++) {
            double y1 = f(p1, x);
            double y2 = f(p2, x);
            if (Math.abs(y1 - y2) > EPS) {
                return false;
            }
        }
        return true;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
    public static String poly(double[] poly) {
        if (poly == null || poly.length == 0) {
            return "0";
        }

        StringBuilder ans = new StringBuilder();
        boolean first = true; // if nothing is yet printed

        // from highest to lowest power
        for (int power = poly.length - 1; power >= 0; power--) {
            double coef = poly[power];
            if (coef == 0) {
                continue; // skipping '0' coefs
            }

            // sign before next coef ( if not first)
            if (first) {
                first = false;
            } else {
                if (coef > 0) {
                    ans.append(" +");
                } else {
                    ans.append(" ");
                }
            }

            // adding the coef
            if (power == 0) {
                // last coef
                ans.append(coef);
            } else {
                // there's an x
                ans.append(coef).append("x");
                if (power > 1) {
                    ans.append("^").append(power);
                }
            }
        }

        // if skips all (arr values are 0)
        if (first) {
            return "0";
        }

        return ans.toString();
    }
    /**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */

	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {

        double v1 = g(p1, p2, x1);
        double xm = (x1 + x2) / 2;
        double vm = g(p1, p2, xm);

        //polynoms almost equal on point xm - stops recursion
        if (Math.abs(vm) < eps){
            return xm;
        }
        //root is on the left side of xm
        if((v1 * vm) <= 0) {
            return sameValue(p1, p2, x1, xm, eps);
        }
        //root is on the right side of xm
        return sameValue(p1, p2, xm, x2, eps);
    }


	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = 0;
        double dx = (x2 - x1) / numberOfSegments;
        for(int i=0; i<numberOfSegments; i++){
            double xl = x1 + i*dx;
            double xr = xl + dx;
            double yl = f(p, xl);
            double yr = f(p, xr);
            double dy = yr - yl;
            double len = Math.sqrt(dx*dx + dy*dy);
            ans+=len;

        }
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
        if (p1 == null || p2 == null || x1 == x2) return 0;

        // make sure x1 bigger than x2
        if (x2 < x1) {
            double tmp = x1;
            x1 = x2;
            x2 = tmp;
        }

        // for better accuracy
        int n = Math.max(numberOfTrapezoid, 1000);

        double ans = 0;
        double dx = (x2 - x1) / n;

        for (int i = 0; i < n; i++) {
            double xi = x1 + i * dx;
            double nextx = x1 + (i + 1) * dx;

            double y1p1 = f(p1, xi);
            double y1p2 = f(p2, xi);
            double y2p1 = f(p1, nextx);
            double y2p2 = f(p2, nextx);

            double dLeft = Math.abs(y1p1 - y1p2);
            double dRight = Math.abs(y2p1 - y2p2);

            double A = (dLeft + dRight) * 0.5 * dx;
            ans += A;
        }
		return ans;
	}
	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return polynom array
	 */
	public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0

        if (p != null && !p.trim().isEmpty()){

            String s = p.replace(" ",""); // "-1.0x^2+3.0x+2.0"
            s = s.replace("-", "+-");// "+-1.0x^2+3.0x+2.0"

            if (s.charAt(0) == '+') {
                s = s.substring(1); // "-1.0x^2+3.0x+2.0"
            }

            String[] terms = s.split("\\+");// ["-1.0x^2","3.0x","2.0"]
            int maxP = 0;
            for (String term : terms) {
                if (term.isEmpty()) continue;

                int power;
                if (term.contains("x^")) {
                    power = Integer.parseInt(term.substring(term.indexOf('^') + 1));
                } else if (term.contains("x")) {
                    power = 1;
                } else {
                    power = 0;
                }
                if (power > maxP) {
                    maxP = power;
                }
            }
            double[] coeffs = new double[maxP + 1];

            for (String term : terms) {
                if (term.isEmpty()) continue;

                int power;
                double coef;

                if (term.contains("x")) {
                    int xIndex = term.indexOf('x');
                    String coefStr = term.substring(0, xIndex);  // what's before x

                    if (coefStr.isEmpty() || coefStr.equals("+")) {
                        coef = 1.0;
                    } else if (coefStr.equals("-")) {
                        coef = -1.0;
                    } else {
                        coef = Double.parseDouble(coefStr);
                    }

                    if (term.contains("^")) {
                        power = Integer.parseInt(term.substring(term.indexOf('^') + 1));
                    } else {
                        power = 1;
                    }
                } else {
                    // no x
                    coef = Double.parseDouble(term);
                    power = 0;
                }

                // in case of multiples with same power sum
                coeffs[power] += coef;
            }
            ans = coeffs;
        }

		return ans;
	}
	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1 - first polynom
	 * @param p2 - second polynom
	 * @return sum of p1 and p2
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//

        if(p1 != null && p2 !=null){

            int p1Len = p1.length;
            int p2Len = p2.length;
            int minLen = Math.min(p1Len,p2Len);
            ans = new double[Math.max(p1Len,p2Len)];
            for (int i = 0; i < minLen; i++) {
                ans[i] = p1[i] + p2[i];
            }
            if(p1Len > p2Len){
                for(int i = minLen; i < p1Len; i ++){
                    ans[i] = p1[i];
                }
            }else if(p2Len > p1Len){
                for(int i = minLen; i < p2Len; i ++){
                    ans[i] = p2[i];
                }
            }
        }
		return ans;
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1 - first polynom
	 * @param p2 - second polynom
	 * @return ans = multiplication of p1 and p2
	 */
	public static double[] mul(double[] p1, double[] p2) {

        if (p1 == null || p2 == null) return ZERO;

        int n1 = p1.length;
        int n2 = p2.length;


        double[] ans = new double[n1 + n2 - 1];

        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                ans[i + j] += p1[i] * p2[j];
            }
        }

        return ans;
	}
	/**
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
        if (po == null || po.length <= 1) return ZERO;

        double[] ans = new double[po.length - 1];

        for (int i = 1; i < po.length; i++) {
            ans[i - 1] = po[i] * i;
        }

        return ans;
	}

    private static int effectiveDegree(double[] p) {
        int deg = p.length - 1;
        while (deg > 0 && Math.abs(p[deg]) <= EPS) {
            deg--;
        }
        return deg;
    }
    private static double g(double[] p1,double[] p2 , double x){
        return f(p1, x) - f(p2, x);
    }
}
