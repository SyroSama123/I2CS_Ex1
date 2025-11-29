import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2026, Ariel University,
 *  * Ex1: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as define in Ex1.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex1Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};


    /**
     * JUnit test for isSameArray().
     *
     * Checks whether two arrays of doubles are equal within a given tolerance.
     *
     * The method compares each corresponding element of arr1 and arr2. If the absolute
     * difference between any pair of elements is greater than or equal to Ex1.EPS, or
     * if the arrays have different lengths, the arrays are considered not equal.
     *
     *
     * @param arr1 - first array of doubles
     * @param arr2 - second array of doubles
     * @return true if both arrays have the same length and all elements are equal
     *         within the EPS tolerance, false otherwise
     */

    boolean isSameArray(double[] arr1, double[] arr2) {
        boolean ans = true;

        if(arr1.length == arr2.length) {
            for(int i=0; i<arr1.length; i++) {
                if( Math.abs( arr1[i] - arr2[i] ) >= Ex1.EPS) {
                    ans = false;
                }
            }
        }
        else {
            ans = false;
        }

        return ans;
    }


    @Test

    /**
     * JUnit test for isSameArray().
     *
     * Tests that the function returns correct boolean value in each of the following cases:
     *  1. Two identical arrays should be considered equal - should return true.
     *  2. Arrays that differ only slightly (within EPS) should be considered equal - should return true.
     *  3. Arrays of different lengths should be considered not equal - should return false.
     */

    void testIsSameArray() {
        double[] arr1 = {-3.234, 0, 5.2, 3.456, 2, 1}, arr2 = {-3.234, 0, 5.2, 3.456, 2, 1}, arr3 = {-3.234, 0, 5.2, 3.456, 2, 1, 0}, arr4 = {-3.233, 0, 5.2, 3.455, 2, 1};

        if(!isSameArray(arr1, arr2)) {
            fail("Arrays are equal but function thinks they're not.");
        }
        if(!isSameArray(arr1, arr4)) {
            fail("Arrays are equal (to EPS) but function thinks they're not.");
        }
        if(isSameArray(arr1, arr3)) {
            fail("Arrays are not equal but function thinks they are.");
        }
    }


    @Test

    /**
     * JUnit test for Ex1.poly().
     *
     * Tests that the function correctly converts an array of polynomial coefficients into a string in the following cases:
     *  1. Empty array - should return "0".
     *  2. Single coefficient - should return the coefficient as a string.
     *  3. Linear polynomial - should return the correct "ax + b" format.
     *  4. Quadratic polynomial  - should return the correct "ax^2 + bx + c" format.
     *  5. Polynomial with negative coefficients - ensures negative signs are correctly formatted.
     *  6. Polynomial with zero coefficients in the middle - ensures zero coefficients are included in the string.
     *
     * Each assertion includes a message describing the specific case being tested.
     *
     */

    void testPoly() {
        double[] emptyPoly = {}, singleCoefPoly = {5}, linearPoly = {1, 2}, quadPoly = {1, 2, 3}, negPoly = {1, -2, 3}, contZeroPoly = {1, 0, 3};

        assertEquals("0", Ex1.poly(emptyPoly), "Empty array should return 0");
        assertEquals("5.0", Ex1.poly(singleCoefPoly), "Single coefficient should return just the number");
        assertEquals("2.0x +1.0", Ex1.poly(linearPoly), "Linear polynomial failed");
        assertEquals("3.0x^2 +2.0x +1.0", Ex1.poly(quadPoly), "Quadratic polynomial failed");
        assertEquals("3.0x^2 -2.0x +1.0", Ex1.poly(negPoly), "Polynomial with negative coefficient failed");
        assertEquals("3.0x^2 +0.0x +1.0", Ex1.poly(contZeroPoly), "Polynomial with zero coefficient failed");

    }


    @Test

    /**
     * JUnit test for Ex1.length().
     *
     * Tests that the function correctly approximates the arc length of a polynomial over an interval using a segmented approach in the following cases:
     *  1. Linear polynomial with varying starting points (i from 0 to 10000 in steps of 50) - ensures that the computed length is within Ex1.EPS of the expected value.
     *
     * The test fails if the calculated length deviates from the expected value by more than Ex1.EPS at any iteration.
     *
     */

    void testLength() {
        double[] p1 = {-9, 83.5};
        int n = 19;
        double realLength = 167.0119;

        for(int i=0; i<10000; i+=50) {
            double length = Ex1.length(p1, i, i+2, n*(i+1));
            boolean isCorrect = Math.abs(length - realLength) < Ex1.EPS;
            if(!isCorrect) {
                fail("function failed at i = " + i);
            }
        }

    }

    @Test

    /**
     * JUnit test for Ex1.polSizeFromString().
     *
     * Tests that the function correctly determines the size of the polynomial (maximum exponent + 1) in the following case:
     *  1. Polynomial string "ax^n + .... + z" - should return n+1 as the size of the polynomial array.
     *
     * The test fails if the returned size does not match the expected value.
     *
     */

    void testPolSize() {
        String poly = "29x^78 + 3";

        int polSize = Ex1.polSizeFromString(poly);

        if(polSize != 79) {
            fail("polSize is:" + polSize + " and not 79");
        }

    }

    @Test

    /**
     * JUnit test for Ex1.equals().
     *
     * Tests that the function correctly determines whether two polynomials are equal in the following cases:
     *  1. Polynomials with different first coefficients but identical remaining coefficients - should return false.
     *  2. Polynomials that are exactly equal - should return true.
     *
     * The test fails if Ex1.equals incorrectly reports equality or inequality.
     *
     */

    void testEquals_() {
        double[] p1 = {2, 33, -15, 52, -30};
        double[] p2 = {1, 33, -15, 52, -30};
        double[] p3 = {2, 33, -15, 52, -30};

        if(Ex1.equals(p1 , p2)) {
           fail("Functions are not equal but func thinks they are");
        }
        if(!Ex1.equals(p1, p3)) {
            fail("Functions are equal but func thinks they aren't");
        }

    }

    @Test

    /**
     * JUnit test for Ex1.getPowFromString().
     *
     * Tests that the function correctly extracts the exponent from a polynomial string at given indices.
     * Each assertion verifies that the returned exponent matches the expected value within Ex1.EPS tolerance.
     *
     */

    void testGetPowFromString() {
        String poly = "29x^78 +12x^32 -2x^23 +3x^2";

        int pow1 = Ex1.getPowFromString(poly, 4);
        int pow2 = Ex1.getPowFromString(poly, 12);
        int pow3 = Ex1.getPowFromString(poly, 19);
        int pow4 = Ex1.getPowFromString(poly, 26);

        assertEquals(pow1 , 78 , Ex1.EPS);
        assertEquals(pow2 , 32 , Ex1.EPS);
        assertEquals(pow3 , 23 , Ex1.EPS);
        assertEquals(pow4 , 2 , Ex1.EPS);

    }

    @Test

    /**
     * JUnit test for Ex1.root_rec().
     *
     * Tests that the function correctly computes a root of a polynomial using a recursive method.
     * The test fails if the computed root deviates from the expected root by more than Ex1.EPS.
     *
     */

    void testRoot_Rec() {
        double[] p1 = {2.0, 3.5, -34.2, 398.345};
        double actualRoot = -0.132774, rootX = Ex1.root_rec(p1, -0.2, 0, Ex1.EPS);

        double absDelta = Math.abs( rootX - actualRoot );

        if(absDelta >= Ex1.EPS){
            fail("Function returned an incorrect root.");
        }
    }

    @Test

    /**
     * JUnit test for Ex1.derivative().
     *
     * Tests that the function correctly computes the derivative of a polynomial represented as an array of coefficients.
     * The test fails if the computed derivative does not match the expected derivative.
     *
     */

    void testDerivative() {
        double[] p1 = {2.0, 3.5, -34.2, 398.345}, actualDer = {3.5, -68.4, 1195.035}, der;

        der = Ex1.derivative(p1);

        if (!isSameArray(der, actualDer)) {
            fail("Derivative is incorrect");
        }
    }


    @Test

    /**
     *  JUnit test for Ex1.PolynomfromPoints().
     *
     *  Tests that the function returns the correct polynomial for each of the following cases:
     *   1. 2 points given - returns a linear polynomial {a,b}.
     *   2. 3 points given - returns a quadratic polynomial {a,b,c}.
     *   3. Array xx and yy are not of equal size - returns null.
     *   4. More than 3 points given - returns null.
     *
     * */

    void testPolynomFromPoints() {
        double[] x3x = {-1.34403, 0.74403, -0.3}, y3y = {0, 0, -2.725}, actualPoly3P = {-2.5, 1.5, 2.5};
        double[] polyFrom3Points = Ex1.PolynomFromPoints(x3x, y3y);

        double[] x2x = {2,6}, y2y = {2, 6}, actualPoly2P = {0,1};
        double[] polyFrom2Points = Ex1.PolynomFromPoints(x2x, y2y);
        double[] polyFrom2_3Points = Ex1.PolynomFromPoints(x2x, y3y);

        if(!isSameArray(actualPoly3P, polyFrom3Points)) {
            fail("Calculated polynome (quadratic) is incorrect.");
        }
        if(!isSameArray(actualPoly2P, polyFrom2Points)) {
            fail("Calculated polynome (linear) is incorrect.");
        }
        if(polyFrom2_3Points != null) {
            fail("Function should return null.");
        }

    }

    @Test

    /**
     * JUnit test for Ex1.getPolynomFromString().
     *
     * Tests that the function correctly parses a polynomial string into an array of coefficients in the following cases:
     *  1. Linear polynomial (bx +a) - returns {a, b}.
     *  2. Quadratic polynomial (cx^2 +bx +c) - returns {a, b, c}.
     *  3. Constant polynomial (a) - returns {a}.
     *  4. Cubic polynomial with negative coefficients (dx^3 + cx^2 -bx -c) - returns {-a, -b, c, d}.
     *  5. Polynomial with floating point coefficients (c.fffx^2 +d.FFx +a) - returns {a, d.FF, c.fff}.
     *
     */

    void testGetPolynomFromString() {

        double[] expectedLinear = {3.0, 2.0};
        assertArrayEquals(expectedLinear, Ex1.getPolynomFromString("2x +3"), Ex1.EPS);

        double[] expectedQuadratic = {2.0, 3.0, -1.0};
        assertArrayEquals(expectedQuadratic, Ex1.getPolynomFromString("-1x^2 +3x +2"), Ex1.EPS);

        double[] expectedConstant = {5.0};
        assertArrayEquals(expectedConstant, Ex1.getPolynomFromString("5"), Ex1.EPS);

        double[] expectedCubic = {-1.0, 4.0, 0.0, -2.0};
        assertArrayEquals(expectedCubic, Ex1.getPolynomFromString("-2x^3 +0x^2 +4x -1"), Ex1.EPS);

        double[] expectedFloat = {3.0, -1.25, 0.5};
        assertArrayEquals(expectedFloat, Ex1.getPolynomFromString("0.5x^2 -1.25x +3.0"), Ex1.EPS);

    }



    @Test

    /**
     * JUnit test for Ex1.add().
     *
     * Tests that the function correctly adds two polynomials represented as arrays of coefficients in the following cases:
     *  1. Polynomials of equal length - sums each corresponding coefficient.
     *  2. First polynomial longer than second - sums overlapping coefficients and appends remaining coefficients from the first polynomial.
     *  3. Second polynomial longer than first - sums overlapping coefficients and appends remaining coefficients from the second polynomial.
     *  4. One or both polynomials are empty arrays - returns the non-empty polynomial or an empty array if both are empty.
     *
     */

    void testAdd_() {

        double[] p1 = {1, 2, 3};
        double[] p2 = {4, 5, 6};
        double[] expectedEqual = {5, 7, 9};
        assertArrayEquals(expectedEqual, Ex1.add(p1, p2), Ex1.EPS);

        double[] p3 = {1, 2, 3, 4};
        double[] p4 = {5, 6};
        double[] expectedFirstLonger = {6, 8, 3, 4};
        assertArrayEquals(expectedFirstLonger, Ex1.add(p3, p4), Ex1.EPS);

        double[] p5 = {1, 2};
        double[] p6 = {3, 4, 5};
        double[] expectedSecondLonger = {4, 6, 5};
        assertArrayEquals(expectedSecondLonger, Ex1.add(p5, p6), Ex1.EPS);

        double[] empty = {};
        double[] nonEmpty = {7, 8, 9};
        assertArrayEquals(nonEmpty, Ex1.add(empty, nonEmpty), Ex1.EPS);
        assertArrayEquals(nonEmpty, Ex1.add(nonEmpty, empty), Ex1.EPS);
        assertArrayEquals(empty, Ex1.add(empty, empty), Ex1.EPS);
    }


    @Test

    /**
     * JUnit test for Ex1.mul().
     *
     * Tests that the function correctly multiplies two polynomials represented as arrays of coefficients in the following cases:
     *  1. Polynomials of equal length (e.g., {1,2} * {3,4}) - returns {3,10,8}.
     *  2. First polynomial longer than second (e.g., {1,2,3} * {4,5}) - returns {4,13,22,15}.
     *  3. Second polynomial longer than first (e.g., {2,1} * {1,2,3}) - returns {2,5,8,3}.
     *  4. One or both polynomials are empty arrays - returns an empty array.
     *  5. Polynomials with negative or floating-point coefficients.
     *
     */


    void testMul() {
        double[] p1 = {1, 2};
        double[] p2 = {3, 4};
        double[] expectedEqual = {3, 10, 8}; // 1*3 + (1*4 + 2*3)x + 2*4 x^2
        assertArrayEquals(expectedEqual, Ex1.mul(p1, p2), Ex1.EPS);

        double[] p3 = {1, 2, 3};
        double[] p4 = {4, 5};
        double[] expectedFirstLonger = {4, 13, 22, 15}; // (manual calculation)
        assertArrayEquals(expectedFirstLonger, Ex1.mul(p3, p4), Ex1.EPS);

        double[] p5 = {2, 1};
        double[] p6 = {1, 2, 3};
        double[] expectedSecondLonger = {2, 5, 8, 3};
        assertArrayEquals(expectedSecondLonger, Ex1.mul(p5, p6), Ex1.EPS);

        double[] empty = {};
        double[] nonEmpty = {1, 2};
        assertArrayEquals(empty, Ex1.mul(empty, nonEmpty), Ex1.EPS);
        assertArrayEquals(empty, Ex1.mul(nonEmpty, empty), Ex1.EPS);
        assertArrayEquals(empty, Ex1.mul(empty, empty), Ex1.EPS);

        double[] p7 = {0.5, -1};
        double[] p8 = {2, -0.5};
        double[] expectedFloat = {1.0, -2.25, 0.5};
        assertArrayEquals(expectedFloat, Ex1.mul(p7, p8), Ex1.EPS);
    }

    @Test

    /**
     * JUnit test for Ex1.area().
     *
     * Tests that the function correctly computes the area between two polynomials over a given range using the trapezoidal method in the following cases:
     *  1. Polynomials do not cross in the interval - area is sum of trapezoids between p1 and p2.
     *  2. Polynomials cross in the interval - splits at intersection point and sums areas of subintervals.
     *  3. Polynomials are equal - area is 0.
     *  4. Polynomials are constant - area equals the absolute difference times the interval length.
     *
     */


    void testArea_() {
        double[] f1 = {0,2};
        double[] f2 = {0,1};
        double computedArea1 = Ex1.area(f1, f2, 0, 3, 1000);
        double expectedArea1 = 4.5; // integral of |2x+1 - (x+2)| dx = integral |x-1| dx from 0 to 3
        assertEquals(expectedArea1, computedArea1, 0.01);

        double[] g1 = {0,1};
        double[] g2 = {0,0.5};
        double computedArea2 = Ex1.area(g1, g2, 0, 4, 1000);
        double expectedArea2 = 4.0; // integral of |x-0.5x| dx = integral of 0.5x dx from 0 to 4
        assertEquals(expectedArea2, computedArea2, 0.01);

        double[] h1 = {0,0,1};
        double[] h2 = {0,0,1};
        assertEquals(0.0, Ex1.area(h1, h2, -2, 2, 1000), 1e-9);

        double[] c1 = {5};
        double[] c2 = {2};
        double expectedArea4 = 9.0; // |5-2|*3
        assertEquals(expectedArea4, Ex1.area(c1, c2, 0, 3, 1000), 1e-9);

    }


    @Test

    /**
     * JUnit test for Ex1.isNumber().
     *
     * Tests that the function correctly identifies numeric characters ('0' to '9') and distinguishes them from non-numeric characters in the following cases:
     *  1. Characters '0' through '9' - should return true.
     *  2. Alphabetic characters 'a', 'Z' - should return false.
     *  3. Special characters '@', '-' - should return false.
     *  4. Space character ' ' - should return false.
     *
     */

    void testIsNumber() {
        assertTrue(Ex1.isNumber('0'), "Character '" + 0 + "' should be recognized as a number");
        assertTrue(Ex1.isNumber('1'), "Character '" + 1 + "' should be recognized as a number");
        assertTrue(Ex1.isNumber('2'), "Character '" + 2 + "' should be recognized as a number");
        assertTrue(Ex1.isNumber('3'), "Character '" + 3 + "' should be recognized as a number");
        assertTrue(Ex1.isNumber('4'), "Character '" + 4 + "' should be recognized as a number");
        assertTrue(Ex1.isNumber('5'), "Character '" + 5 + "' should be recognized as a number");
        assertTrue(Ex1.isNumber('6'), "Character '" + 6 + "' should be recognized as a number");
        assertTrue(Ex1.isNumber('7'), "Character '" + 7 + "' should be recognized as a number");
        assertTrue(Ex1.isNumber('8'), "Character '" + 8 + "' should be recognized as a number");
        assertTrue(Ex1.isNumber('9'), "Character '" + 9 + "' should be recognized as a number");

        assertFalse(Ex1.isNumber('a'));
        assertFalse(Ex1.isNumber('Z'));
        assertFalse(Ex1.isNumber('@'));
        assertFalse(Ex1.isNumber('-'));
        assertFalse(Ex1.isNumber(' '));
    }


    @Test

    /**
     * JUnit test for Ex1.sameValue().
     *
     * Tests that the function correctly finds an x-value within a given interval where two polynomials have approximately the same value.
     *
     * The test fails if the absolute difference between the polynomials at the returned x exceeds Ex1.EPS.
     *
     */

    void testSameValue() {
        double[] p1 = {0, 1};
        double[] p2 = {0, 0.5};
        double x1 = Ex1.sameValue(p1, p2, 0, 4, Ex1.EPS);
        assertTrue(Math.abs(Ex1.f(p1, x1) - Ex1.f(p2, x1)) < Ex1.EPS, "Returned x does not satisfy f(p1,x1) ≈ f(p2,x1)");

        double[] p3 = {-1, 1};
        double[] p4 = {-0.5, 0.5};
        double x2 = Ex1.sameValue(p3, p4, 0, 4, Ex1.EPS);
        assertTrue(Math.abs(Ex1.f(p3, x2) - Ex1.f(p4, x2)) < Ex1.EPS, "Returned x does not satisfy f(p3,x2) ≈ f(p4,x2)");
    }












    /**
     *
     * ___________________________________________________________________________________________
     *
     * Original Unchanged tests:
     *
     *
     * */




 	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex1.f(po1, 0);
		double fx1 = Ex1.f(po1, 1);
		double fx2 = Ex1.f(po1, 2);
		assertEquals(fx0, 2, Ex1.EPS);
		assertEquals(fx1, 4, Ex1.EPS);
		assertEquals(fx2, 6, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex1.add(po1, po2);
		double f1x = Ex1.f(po1, x);
		double f2x = Ex1.f(po2, x);
		double f12x = Ex1.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex1.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex1.mul(po2, minus1);
		double[] p1 = Ex1.add(p12, pp2);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex1.add(po1, po2);
		double[] p21 = Ex1.add(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex1.add(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex1.mul(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, Ex1.ZERO));
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex1.mul(po1, po2);
		double[] p21 = Ex1.mul(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex1.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex1.f(po1, x);
			double f2x = Ex1.f(po2, x);
			double f12x = Ex1.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex1.EPS);
		}
	}
	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex1.derivative(p); // 2x + 6
		double[] dp2 = Ex1.derivative(dp1); // 2
		double[] dp3 = Ex1.derivative(dp2); // 0
		double[] dp4 = Ex1.derivative(dp3); // 0
		assertTrue(Ex1.equals(dp1, pt));
		assertTrue(Ex1.equals(Ex1.ZERO, dp3));
		assertTrue(Ex1.equals(dp4, dp3));
	}
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "3.1x^2 +2.3x -1.1";
		String sp = Ex1.poly(p);
		double[] p1 = Ex1.getPolynomFromString(sp);
		double[] p2 = Ex1.getPolynomFromString(sp2);
		boolean isSame1 = Ex1.equals(p1, p);
		boolean isSame2 = Ex1.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex1.poly(p1));
	}
	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex1.ZERO, {1+ Ex1.EPS/2}, {1,2}};
		double[][] xx = {{-2* Ex1.EPS}, {1+ Ex1.EPS*1.2}, {1,2, Ex1.EPS/2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex1.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex1.equals(d1[i], xx[i]));
		}
	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1=-4, x2=0;
		double rs1 = Ex1.sameValue(po1,po2, x1, x2, Ex1.EPS);
		double rs2 = Ex1.sameValue(po2,po1, x1, x2, Ex1.EPS);
		assertEquals(rs1,rs2, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=-4, x2=0;
		double a1 = Ex1.area(po1, po2, x1, x2, 100);
		double a2 = Ex1.area(po2, po1, x1, x2, 100);
		assertEquals(a1,a2, Ex1.EPS);
}
	@Test
	/**
	 * Test the area f1(x)=0, f2(x)=x;
	 */
	public void testArea2() {
		double[] po_a = Ex1.ZERO;
		double[] po_b = {0,1};
		double x1 = -1;
		double x2 = 2;
		double a1 = Ex1.area(po_a,po_b, x1, x2, 1);
		double a2 = Ex1.area(po_a,po_b, x1, x2, 2);
		double a3 = Ex1.area(po_a,po_b, x1, x2, 3);
		double a100 = Ex1.area(po_a,po_b, x1, x2, 100);
		double area =2.5;
		assertEquals(a1,area, Ex1.EPS);
		assertEquals(a2,area, Ex1.EPS);
		assertEquals(a3,area, Ex1.EPS);
		assertEquals(a100,area, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function.
	 */
	public void testArea3() {
		double[] po_a = {2,1,-0.7, -0.02,0.02};
		double[] po_b = {6, 0.1, -0.2};
		double x1 = Ex1.sameValue(po_a,po_b, -10,-5, Ex1.EPS);
		double a1 = Ex1.area(po_a,po_b, x1, 6, 8);
		double area = 58.5658;
		assertEquals(a1,area, Ex1.EPS);
	}






}
