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

    /**
     *
     * ___________________________________________________________________________________________
     *
     * My tests:
     *
     *
     * */


    @Test

    /**
     *
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
     * My test -- remove this line before submitting
     *
     * Tests that Poly() returns correct string.
     */
    void testPoly() {
        double[] p1 = {9.0, 7.3, -2.0, 4.0, -13.25};
        String polyStr = "-13.25x^4 +4.0x^3 -2.0x^2 +7.3x +9.0";

        String p1Str = Ex1.poly(p1);

        if(!p1Str.equals(polyStr)) {
            fail("Strings not equal");
        }
    }

    @Test
    /**
     * My test -- remove this line before submitting
     *
     * Tests that length() returns correct string.
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
     * My test -- remove this line before submitting
     *
     * Tests that polSizeFromString() returns correct polSize.
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
     * My test -- remove this line before submitting
     *
     * Tests that equals() returns correct boolean value.
     */
    void testEquals_1() {
        double[] p1 = {2, 33, -15, 52, -30};
        double[] p2 = {1, 33, -15, 52, -30};
        double[] p3 = {2, 33, -15, 52, -30};

        if(Ex1.equals(p1, p2)) {
           fail("Functions are not equal but func thinks they are");
        }
        if(!Ex1.equals(p1, p3)) {
            fail("Functions are equal but func thinks they aren't");
        }

    }

    @Test
    /**
     * My test -- remove this line before submitting
     *
     * Tests that getPowFromString() returns correct power.
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
     * My test -- remove this line before submitting
     *
     * Tests that root_rec() return an x value such that | x - actualRoot | < EPS  .
     *
     */
    void testRoot_Rec() {
        double[] p1 = {2.0, 3.5, -34.2, 398.345};
        double actualRoot = -0.132774, rootX = Ex1.root_rec(p1, -0.2, 0, Ex1.EPS);

        double absDelta = Math.abs( rootX - actualRoot );

        if(absDelta >= Ex1.EPS){
            fail();
        }
    }

    @Test

    /**
     * My test -- remove this line before submitting
     *
     * Tests that root_rec() return an x value such that | x - actualRoot | < EPS  .
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
     *
     *
     *
     * */

    void testPolynomFromPoints() {
        double[] xx = {-1.34403, 0.74403, -0.3}, yy = {0, 0, -2.725}, actualPoly = {-2.5, 1.5, 2.5};
        double[] polyFromPoints = Ex1.PolynomFromPoints(xx, yy);

        if(!isSameArray(actualPoly, polyFromPoints)) {
            fail("Calculated polynome is incorrect.");
        }
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
