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
     * Computes the polynomial function that passes through the given points.
     *
     * Given arrays of x-values and y-values of equal length, this method returns
     * the coefficients of the unique polynomial of degree 1 (for 2 points) or
     * degree 2 (for 3 points) that interpolates all the points.
     *
     * The returned array is in the format:
     *   poly[i] = coefficient of x^i
     *
     * Restrictions:
     * - xx and yy must be non-null
     * - both arrays must have the same length
     * - the number of points must be 2 or 3
     *
     * For 2 points: returns a linear polynomial (ax + b).
     * For 3 points: returns a quadratic polynomial (ax^2 + bx + c).
     *
     * @param xx - array of x-coordinates
     * @param yy - array of y-coordinates
     * @return a coefficient array of the interpolating polynomial,
     *         or null if the input is invalid.
     */
    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        double [] ans = null;
        int lx = xx.length;
        int ly = yy.length;
        if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
            if (lx == 2) {
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];

                double a = (y2 - y1) / (x2 - x1);
                double b = y1 - a * x1;

                ans = new double[2];
                ans[0] = b;
                ans[1] = a;
            }

            else if (lx == 3) {
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                double x3 = xx[2], y3 = yy[2];

                double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);

                double a = ( x3*(y2 - y1) + x2*(y1 - y3) + x1*(y3 - y2) ) / denom;
                double b = ( x3*x3*(y1 - y2) + x2*x2*(y3 - y1) + x1*x1*(y2 - y3) ) / denom;
                double c = ( x2*x3*(x2 - x3)*y1 + x3*x1*(x3 - x1)*y2 + x1*x2*(x1 - x2)*y3 ) / denom;

                ans = new double[3];
                ans[0] = c;
                ans[1] = b;
                ans[2] = a;
            }
        }
        return ans;
    }
    /**
     * Checks whether two polynomial functions are equal.
     *
     * The method compares the values of the two polynomials p1 and p2 at all
     * integer points from 0 up to the maximum degree of the two polynomials.
     * If, for every such point, the absolute difference between the function
     * values is smaller than EPS, the polynomials are considered equal.
     *
     * Notes:
     * - The comparison is based on function values, not on coefficient arrays.
     * - The degree of each polynomial is extracted from its string representation.
     * - Equality tolerance is determined by the constant EPS.
     *
     * @param p1 - first polynomial (coefficients)
     * @param p2 - second polynomial (coefficients)
     * @return true if the polynomials are equal within EPS tolerance,
     *         false otherwise.
     */
    public static boolean equals(double[] p1, double[] p2) {
        boolean ans = true;
        String p1Str = poly(p1), p2Str = poly(p2);
        int maxDeg = Math.max(polSizeFromString(p1Str) , polSizeFromString(p2Str));
        double absDelta;

        for(int i=0; i<maxDeg+1; i++) {
            absDelta = Math.abs( f(p1, i) - f(p2, i) );
            if(absDelta >= EPS) {
                ans = false;
            }
        }
        return ans;
    }

    /**
     * Converts a polynomial (given as a coefficient array) into a string.
     *
     * The input array poly[] is interpreted such that:
     *     poly[i] = coefficient of x^i
     *
     * The method constructs a string representation of the polynomial in
     * descending powers of x, including terms of the form:
     *     ax^n, ax, or a
     *
     *
     * @param poly - polynomial coefficients where poly[i] corresponds to x^i
     * @return a string representation of the polynomial
     */
    public static String poly(double[] poly) {
        String ans = "";
        int l = poly.length;
        if(l==0) {
            ans="0";
        }

        else {
            for(int i=1; i <= l; i++) {
                if(i<l-1) {
                    if(poly[l - i - 1] >= 0) {
                        ans += poly[l - i] + "x^" + (l-i) + " +";
                    }
                    else if (poly[l - i - 1] < 0) {
                        ans += poly[l - i] + "x^" + (l-i) + " ";
                    }
                }
                else if(i==l-1) {
                    if(poly[l - i - 1] > 0) {
                        ans += poly[l - i] + "x" + " +";
                    }
                    else if (poly[l - i - 1] < 0) {
                        ans += poly[l - i] + "x" + " ";
                    }
                    else {
                        ans += poly[l - i] + "x";
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
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps.
     * This function computes an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps,
     * assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
     * This function is implemented recursively.
     *
     * @param p1 - the first polynomial function
     * @param p2 - the second polynomial function
     * @param x1 - minimal value of the range
     * @param x2 - maximal value of the range
     * @param eps - epsilon (positive small value, often 10^-3 or 10^-6)
     * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
     */

    public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        double deltaF1 = f(p1, x1) - f(p2, x1), x12 = (x1 + x2) / 2;
        double deltaF12 = f(p1, x12) - f(p2, x12);

        if (Math.abs(deltaF12) < eps) {
            return x12;
        }

        if (deltaF1 * deltaF12 <= 0) {
            return sameValue(p1, p2, x1, x12, eps);
        } else {
            return sameValue(p1, p2, x12, x2, eps);
        }
    }

    /**
     * Computes the approximate arc length of a polynomial function over an interval.
     *
     * The method divides the interval [x1, x2] into a given number of equal segments
     * and approximates the total curve length by summing the lengths of straight-line
     * segments connecting consecutive points on the graph of the polynomial.
     *
     * For each segment, the method computes:
     *   dx = horizontal step size
     *   dy = difference in function values between the segment endpoints
     *   length += sqrt(dx^2 + dy^2)
     *
     *
     * @param p                 - polynomial coefficients
     * @param x1                - start of the interval
     * @param x2                - end of the interval
     * @param numberOfSegments  - number of straight segments used in the approximation
     * @return the approximate length of the polynomial curve between x1 and x2
     */

    public static double length(double[] p, double x1, double x2, int numberOfSegments) {
        double length = 0;

        double dx = (x2-x1) / numberOfSegments, dy, l;

        for(int i=0; i<numberOfSegments; i++) {
            dy = f(p, x1+(i*dx)) - f(p, x1+((i+1)*dx));
            l = Math.sqrt(dy*dy + dx*dx);
            length += l;
        }

        return length;
    }

    /**
     * Computes the approximate area between two polynomial functions over an interval.
     *
     * The method divides the interval [x1, x2] into a given number of equal trapezoids
     * and approximates the total area by summing the areas of trapezoids formed between
     * the two polynomial curves. If the polynomials cross within the interval, the method
     * recursively splits the interval at the intersection point to ensure correct calculation.
     *
     * For each trapezoid, the method computes:
     *   dy1 = vertical distance between polynomials at the left endpoint
     *   dy2 = vertical distance between polynomials at the right endpoint
     *   A   = 0.5 * dx * (dy1 + dy2)   // area of the trapezoid
     *   area += |A|
     *
     * If the polynomials cross, the interval is split at the crossing point and the areas
     * of the subintervals are computed recursively to avoid negative contributions.
     *
     * @param p1                  - coefficients of the first polynomial
     * @param p2                  - coefficients of the second polynomial
     * @param x1                  - start of the interval
     * @param x2                  - end of the interval
     * @param numberOfTrapezoid   - number of trapezoids used in the approximation
     * @return the approximate area between the two polynomials over [x1, x2]
     */


    public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
        double area = 0;
        double dx = (x2-x1) / numberOfTrapezoid, dy1, dy2, A;

        double deltaY1 = f(p1, x1) - f(p2, x1), deltaY2 = f(p1, x2) - f(p2, x2);

        if(deltaY1*deltaY2 >= 0) {
            for(int i=0; i<numberOfTrapezoid; i++) {
                dy1 = Math.abs( f(p1, x1+(i*dx)) - f(p2, x1+(i*dx)) );
                dy2 = Math.abs( f(p1, x1+((i+1)*dx)) - f(p2, x1+((i+1)*dx)) );

                A = 0.5 * dx * (dy1+dy2);
                area += Math.abs(A);
            }
        }

        else {
            double x12 = sameValue(p1, p2, x1, x2, EPS);
            area += area(p1, p2, x1, x12, numberOfTrapezoid) + area(p1, p2, x12, x2, numberOfTrapezoid);
        }
        return area;
    }


    /**
     * This function computes the array representation of a polynomial function from a String
     * representation.
     *
     * The input string p should represent a polynomial in standard format, for example:
     *     "-1.0x^2 + 3.0x + 2.0"
     *
     * The method parses each term, determines its coefficient and power, and fills
     * the resulting array such that:
     *    ans[i] = coefficient of x^i
     *
     * Notes:
     * -  given a polynomial function represented as a double array,
     * -  getPolynomFromString(poly(p)) should return an array equals to p.
     *
     * @param p - a String representing polynomial function.
     * @return
     */
    public static double[] getPolynomFromString(String p) {
        int polSize = polSizeFromString(p), pow;
        double [] ans = new double[polSize];//  -1.0x^2 +3.0x +2.0
        String temp = "";
        boolean atNumber = false;

        for(int i=0; i<p.length(); i++) {
            if( (isNumber(p.charAt(i)) || p.charAt(i) == '-') && !atNumber) {
                temp += p.charAt(i);
                atNumber = true;
            }

            else if(atNumber && (p.charAt(i) == '.' || isNumber(p.charAt(i)))) {
                temp += p.charAt(i);
            }

            if(i < p.length() - 2 && atNumber && (p.charAt(i) == 'x' && p.charAt(i+1) == '^')) {
                pow  = getPowFromString(p , i+2);
                ans[pow] = Double.parseDouble(temp);
                temp = "";
                atNumber = false;
                i+=3;
            }

            if(i<p.length() && atNumber && (p.charAt(i) == 'x' && p.charAt(i+1) != '^')) {
                ans[1] = Double.parseDouble(temp);
                temp = "";
                atNumber = false;
                i+=1;
            }
            if(i == p.length()-1 && atNumber) {
                ans[0] = Double.parseDouble(temp);
            }
            /* Finish this later */

        }
        return ans;
    }
    /**
     * Computes the sum of two polynomials.
     *
     * Given two polynomials represented by coefficient arrays p1 and p2, this method
     * returns a new array representing their sum. The coefficients are aligned by
     * power such that ans[i] = p1[i] + p2[i]. If the polynomials have different
     * degrees, the extra coefficients from the longer polynomial are copied as-is.
     *
     * Notes:
     * - The length of the returned array equals the maximum degree of p1 and p2.
     *
     * @param p1 - first polynomial
     * @param p2 - second polynomial
     * @return a new coefficient array representing p1 + p2
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
     * Computes the product of two polynomials.
     *
     * Given two polynomials represented by coefficient arrays p1 and p2, this method
     * returns a new array representing their product. The coefficients are combined
     * according to the distributive property:
     *     (p1[0] + p1[1]x + ... ) * (p2[0] + p2[1]x + ... )
     *
     * Notes:
     * - The length of the resulting array is at most p1.length + p2.length + 1.
     * - The method uses an intermediate array for partial products, then sums them
     *   using the add() method.
     *
     * @param p1 - first polynomial
     * @param p2 - second polynomial
     * @return a new coefficient array representing the product p1 * p2
     */
    public static double[] mul(double[] p1, double[] p2) {
        int l1 = p1.length, l2 = p2.length;

        if(l1 == 0 || l2 ==0) {
            return new double[] {};
        }

        double [] ans = new double[l1+l2-1];//
        double [][] pArr = new double[l1][l1+l2-1];



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
     * Computes the derivative of a polynomial.
     *
     * Given a polynomial represented by a coefficient array po, this method returns
     * a new array representing its derivative. For each term a*x^i, the derivative
     * is computed as i*a*x^(i-1) and placed in the corresponding index of the result.
     *
     * Special case:
     * - If the input polynomial is constant (length 1), the derivative is defined
     *   as 0, and the result is an array containing a single element {0}.
     *
     * Notes:
     * - For polynomials of degree n ≥ 1, the resulting array has length po.length - 1.
     *
     * @param po - polynomial coefficients
     * @return a new coefficient array representing the derivative of po
     */

    public static double[] derivative (double[] po) {
        double [] ans = new double[po.length-1];

        for(int i = 1; i < po.length; i++) {
            ans[i-1] = po[i] * (i);
        }

        if(po.length == 1) {
            double [] ansForConstant = {0};
            ans = ansForConstant;
        }

        return ans;
    }


    /**
     *
     * My helper functions.
     *
     * */


    /**
     * Checks whether a character represents a numeric digit.
     *
     * The method determines if the input character x is one of the digits '0' through '9'.
     *
     * Notes:
     * - Only decimal digits are considered numeric.
     * - Returns false for all other characters, including letters, symbols, or whitespace.
     *
     * @param x - the character to check
     * @return true if x is a digit ('0' to '9'), false otherwise
     */
    public static boolean isNumber(char x) {
        char[] nums = {'0','1','2','3','4','5','6','7','8','9'};
        boolean ans = false;

        for(int i=0; i<nums.length; i++) {
            if(nums[i] == x) {
                ans = true;
            }
        }

        return ans;
    }

    /**
     * Determines the size of the coefficient array needed for a polynomial string.
     *
     * Given a string representation of a polynomial, this method returns the number
     * of coefficients required to store the polynomial in an array. The size is
     * determined by the highest power of x in the string.
     *
     * Rules:
     * - If the highest term is x^n, the returned size is n + 1.
     * - If the highest term is x without an explicit power, the returned size is 2.
     * - If the string contains no x terms, the returned size is 1 (for a constant).
     *
     *
     * @param p - string representation of the polynomial
     * @return the required size of the coefficient array to represent the polynomial
     */

    public static int polSizeFromString(String p) {
        int size = 0;
        String temp = "";

        for(int i=0; i<p.length(); i++) {
            if(p.charAt(i) == 'x') {
                if(p.charAt(i+1) == '^') {
                    for(int j=i+2; j< p.length(); j++) {
                        if(isNumber(p.charAt(j))) {
                            temp += p.charAt(j);
                        }
                        else {
                            break;
                        }
                    }
                    size = Integer.parseInt(temp) + 1;

                    return size;
                }
                else {
                    size = 2;

                    return size;
                }
            }
            if(i == p.length() - 1) {
                size = 1;
            }
        }

        return size;
    }

    /**
     * Extracts the exponent of a term from a polynomial string starting at a given index.
     *
     * The method reads consecutive numeric characters starting from the specified index
     * and converts them into an integer representing the power of x for that term.
     *
     *
     * @param p     - the polynomial string
     * @param index - starting index in the string where the exponent begins
     * @return the integer exponent parsed from the string
     */

    public static int getPowFromString(String p, int index) {
        int power = 0;
        String temp = "";

        for(int i = index; i<p.length(); i++) {
            if(isNumber(p.charAt(i))) {
                temp += p.charAt(i);
            }
            else {
                break;
            }
        }

        power += Integer.parseInt(temp);

        return power;
    }


}



