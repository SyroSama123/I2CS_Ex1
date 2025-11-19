import java.util.Arrays;

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
            /** add you code below

             /////////////////// */
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
        boolean ans = true;
        /** add you code below

         /////////////////// */
        return ans;
    }

    /**
     * Computes a String representing the polynomial function.
     * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
     * @param poly the polynomial function represented as an array of doubles
     * @return String representing the polynomial function:
     */
    public static String poly(double[] poly) {
        String ans = "";
        int l = poly.length;
        if(l==0) {
            ans="0";
        }

        else {
            for(int i=1; i <= l; i++) {
                if(i<l) {
                    if(poly[l - i - 1] > 0) {
                        ans += poly[l - i] + "x^" + (l-i) + " +";
                    }
                    else if (poly[l - i - 1] < 0) {
                        ans += poly[l - i] + "x^" + (l-i) + " ";
                    }
                }
                else {
                    ans += poly[0];
                }

            }
        }
        return ans;
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
        double ans;
        double x_1 = Math.min(x1, x2), x_2 = Math.max(x1, x2);
        double deltaY1 = Math.abs( f(p1, x_1) - f(p2, x_1) ), deltaY2 = Math.abs( f(p1, x_2) - f(p2, x_2) );

        if(deltaY1 > deltaY2) {
            ans = x_1;
            for(int i=0; ans <= x_2; i++) {
                ans = x_1 + (i * eps);
                double absDelta = Math.abs(f(p1, ans) - f(p2, ans));
                if(absDelta < eps) {
                    return ans;
                }
            }
        }
        else {
            ans = x_2;
            for(int i=0; ans >= x_1; i--) {
                ans = x_1 + (i * eps);
                double absDelta = Math.abs(f(p1, ans) - f(p2, ans));
                if(absDelta < eps) {
                    return ans;
                }
            }
        }


        return ans;
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
        double ans = x1;
        /** add you code below

         /////////////////// */
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
        double ans = 0;
        /** add you code below

         /////////////////// */
        return ans;
    }
    /**
     * This function computes the array representation of a polynomial function from a String
     * representation. Note:given a polynomial function represented as a double array,
     * getPolynomFromString(poly(p)) should return an array equals to p.
     *
     * @param p - a String representing polynomial function.
     * @return
     */
    public static double[] getPolynomFromString(String p) {
        double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        /** add you code below

         /////////////////// */
        return ans;
    }
    /**
     * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
     * @param p1
     * @param p2
     * @return
     */
    public static double[] add(double[] p1, double[] p2) {
        double [] ans = new double[Math.max(p1.length, p2.length)];//
        for(int i=0; i<ans.length; i++) {
            if(i<Math.min(p1.length, p2.length)) {
                ans[i] = p1[i] + p2[i];
            }
            else if(p1.length > p2.length) {
                ans[i] = p1[i];
            }
            else {
                ans[i] = p2[i];
            }

        }
        return ans;
    }
    /**
     * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
     *
     * {1, 2, 3, 4, 5} * {1,2, 3, 4, 5} =
     *
     * @param p1
     * @param p2
     * @return
     */
    public static double[] mul(double[] p1, double[] p2) {
        int l1 = p1.length, l2 = p2.length;
        double [] ans = new double[l1+l2+1];//
        double [][] pArr = new double[l1][l1+l2+1];

        for(int i=0; i<l1; i++) {
            for(int j=0; j<l2; j++) {
                pArr[i][j+i] = p1[i] * p2[j];
            }
        }

        for(int i = 0; i<l1; i++) {
            ans = add(ans, pArr[i]);
        }

        return ans;
    }
    /**
     * This function computes the derivative of the p0 polynomial function.
     *
     * 3+2x+5x^2+6x^3 = {3, 2, 5, 6}
     *
     * 2+10x+18x^2 = {2 * 1, 5 * 2, 6 * 3} = {2, 10, 18}
     *
     * @param po
     * @return
     */
    public static double[] derivative (double[] po) {
        double [] ans = new double[po.length];//
        for(int i = 1; i < po.length; i++) {
            ans[i-1] = po[i] * (i);
        }

        return ans;
    }
}
