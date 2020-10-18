"""
find_roots.py


OVERVIEW

This module contains functions that implement four algorithms for finding roots
of 1-D functions.

AUTHOR

Dr. Phillip M. Feldman
"""

def find_root_bisection(f, a, b, ftol=1.e-6, xtol=1.e-6, both=True,
  max_steps=3000, verbose=False):
   """
   OVERVIEW

   This function finds a root (zero) of a function `f` via bisection search.
   `f` is required to be everywhere continuous, but is not required to have a
   derivative anywhere.  The calling program provides the function, a pair of
   starting values `a` and `b` for the independent variable, and optionally
   other inputs as well (see below).  This root-finding algorithm requires that
   f(a) and f(b) have opposite signs.


   INPUTS

   `f` is a function that takes a single real-valued argument and returns a
   single real values result.

   `a` and `b` are the starting values of the independent variable `x`.

   `ftol` and `xtol` are convergence thresholds; both default to 1.e-6.

   `both` is a bool value that specifies whether both convergence thresholds
   must be simultaneously satisfied; the default is `True`.  If `both` equals
   `False`, the algorithm stops as soon as either convergence threshold is
   satisfied.

   `max_steps` is the maximum allowed number of iterations.  If convergence
   does not occur within this number of steps, an `AlgorithFailure` exception is
   raised.

   `verbose`: Setting this input to `True` causes the function to display the
   numbers of function calls and iterations required for convergence.
   """

   if xtol <= 0.0:
      raise ValueError("If specified, `xtol` must be positive.")
   if ftol <= 0.0:
      raise ValueError("If specified, `ftol` must be positive.")

   f_a= f(a)
   f_b= f(b)

   if f_a * f_b > 0.0:
      raise AlgorithmFailure("This root-finding algorithm requires that f(a) "
        "and f(b) have opposite signs.  You may try using `find_root`, which "
        "does not have this restriction.")

   calls= 2
   steps= 1

   if both:
      if abs(b - a) <= xtol:
         if abs(f_a) <= ftol:
            if verbose:
               print("Convergence achieved after %d steps and %d calls."
                 % (steps, calls))
            return a

         if abs(f_b) <= ftol:
            if verbose:
               print("Convergence achieved after %d steps and %d calls."
                 % (steps, calls))
            return b

   else:
      if abs(f_a) <= ftol:
         if verbose:
            print("Convergence achieved after %d steps and %d calls."
              % (steps, calls))
         return a

      if abs(f_b) <= ftol:
         if verbose:
            print("Convergence achieved after %d steps and %d calls."
              % (steps, calls))
         return b


   while True:
      c= 0.5 * (a+b)
      f_c= f(c)
      calls+= 1

      if both:
         if abs(c - b) <= xtol and abs(f_c) <= ftol:
            if verbose:
               print("Convergence achieved after %d steps and %d calls."
                 % (steps, calls))
            return c

      else:
         if abs(c - b) <= xtol or abs(f_c) <= ftol:
            if verbose:
               print("Convergence achieved after %d steps and %d calls."
                 % (steps, calls))
            return c

      # One might expect the following test to be `if f_a * f_c < 0`, but this
      # fails to produce the desired result when f_c equals zero.
      if f_a * f_c < f_b * f_c:
         b, f_b= c, f_c
      else:
         a, f_a, b, f_b= b, f_b, c, f_c

      steps+= 1
      if steps > max_steps:
         raise AlgorithmFailure("Limit of %d iterations has been reached."
           % max_steps)

   # end while True


def find_root_secant(f, a, b, ftol=1.e-6, xtol=1.e-6, both=True,
  max_steps=3000, min_slope=1.e-60, verbose=False):
   """
   OVERVIEW

   This function uses the secant method to find a root (zero) of a function
   starting from two values of x that do not necessarily enclose a root.

   Like Newton's method, the closely-related secant and modified Regula Falsi
   methods use successive linear approximations to the function, but unlike
   Newton's method, they do not require knowledge of the derivative.


   INPUTS

   `f` is a function that takes a single real-valued argument and returns a
   single real values result.

   `a` and `b` are the starting values of the independent variable `x`.

   `ftol` and `xtol` are convergence thresholds; both default to 1.e-6.

   `both` is a bool value that specifies whether both convergence thresholds
   must be simultaneously satisfied; the default is `True`.  If `both` equals
   `False`, the algorithm stops as soon as either convergence threshold is
   satisfied.

   `max_steps` is the maximum allowed number of iterations.  If convergence
   does not occur within this number of steps, an `AlgorithFailure` exception is
   raised.

   `min_slope` is the minimum absolute value of the slope that is considered
   adequately different from zero for the algorithm to proceed.  The default is
   1.e-60.

   `verbose`: Setting this input to `True` causes the function to display the
   number of steps required for convergence.


   REFERENCE

   http://en.wikipedia.org/wiki/Secant_method
   """

   if xtol <= 0.0:
      raise ValueError("If specified, `xtol` must be positive.")
   if ftol <= 0.0:
      raise ValueError("If specified, `ftol` must be positive.")

   calls= 0
   steps= 1

   if both:
      f_a= f(a)
      f_b= f(b)
      calls+= 2

      if abs(b - a) <= xtol:
         if abs(f_a) <= ftol:
            if verbose:
               print("Convergence achieved after %d steps and %d calls."
                 % (steps, calls))
            return a

         if abs(f_b) <= ftol:
            if verbose:
               print("Convergence achieved after %d steps and %d calls."
                 % (steps, calls))
            return b

   else:
      if abs(b - a) <= xtol:
         raise ValueError("When `both` equals `False`, the initial values of "
           "x (a and b) must differ by more than xtol= %f." % xtol)

      f_a= f(a)
      calls+= 1
      if abs(f_a) <= ftol:
         if verbose:
            print("Convergence achieved after %d steps and %d calls."
              % (steps, calls))
         return a

      f_b= f(b)
      calls+= 1
      if abs(f_b) <= ftol:
         if verbose:
            print("Convergence achieved after %d steps and %d calls."
              % (steps, calls))
         return b


   while True:
      slope= (f_b - f_a)/float(b - a)
      if abs(slope) < min_slope:

         # The slope is essentially zero.  Attempt to fix this by shifting `b`
         # towards `a` (the new `b` is a weighted average of `a` and the old
         # `b`):
         b= 0.3679*a + 0.6321*b
         f_b= f(b)
         calls+= 1

         slope= (f_b - f_a)/float(b - a)
         if abs(slope) < min_slope:
            if abs(b - a) <= xtol:
               break
            raise AlgorithmFailure("At step %d, the slope is too close to zero!"
              % steps)

      # Calculate c= updated estimate of root location and f_c= corresponding
      # function value:
      c= b - f_b/slope
      f_c= f(c)
      calls+= 1

      if both:
         if abs(c - b) <= xtol and abs(f_c) <= ftol:
            break
      else:
         if abs(c - b) <= xtol or  abs(f_c) <= ftol:
            break

      steps+= 1
      if steps > max_steps:
         raise AlgorithmFailure(
           "Limit of %d iterations has been reached." % max_steps)

      a, f_a= b, f_b
      b, f_b= c, f_c

   # end while True

   if verbose:
      print("Convergence achieved after %d steps and %d calls."
        % (steps, calls))

   return c


def find_root_Regula_Falsi(f, a, b, ftol=1.e-6, xtol=1.e-6, both=True,
  max_steps=3000, min_slope=1.e-60, verbose=False):
   """
   OVERVIEW

   This function finds a root (zero) of a function `f` using the modified Regula
   Falsi method (Ref. 2), which we henceforth call 'Regula Falsi'.  `f` is
   required to be everywhere continuous, but is not required to have a
   derivative anywhere.

   The calling program provides the function, a pair of starting values `a` and
   `b` for the independent variable, and optionally other inputs as well (see
   below).  If `a` and `b` bracket a root--more precisely, f(a) and f(b) have
   opposite signs--or if at any point the algorithm discovers a pair of values
   that bracket at root, subsequent intervals will continue to bracket that
   root.  Note that this behavior is different from that of the secant method,
   which given an interval [a, b] that contains a root, can converge to a root
   outside that interval or fail to converge.


   INPUTS

   `f` is a function that takes a single real-valued argument and returns a
   single real values result.

   `a` and `b` are the starting values of the independent variable `x`.

   `ftol` and `xtol` are convergence thresholds; both default to 1.e-6.

   `both` is a bool value that specifies whether both convergence thresholds
   must be simultaneously satisfied; the default is `True`.  If `both` equals
   `False`, the algorithm stops as soon as either convergence threshold is
   satisfied.

   `max_steps` is the maximum allowed number of iterations.  If convergence
   does not occur within this number of steps, an `AlgorithFailure` exception is
   raised.

   `min_slope` is the minimum absolute value of the slope that is considered
   adequately different from zero for the algorithm to proceed.  The default is
   1.e-60.

   `verbose`: Setting this input to `True` causes the function to display the
   numbers of function calls and iterations required for convergence.
   """

   if xtol <= 0.0:
      raise ValueError("If specified, `xtol` must be positive.")
   if ftol <= 0.0:
      raise ValueError("If specified, `ftol` must be positive.")
   if not 0.5 <= contraction_factor <= 1.0:
      raise ValueError("If specified, `contraction` factor must be in the "
        "interval [0, 1].")

   if not isinstance(use_bisection, int) or use_bisection < 0:
      raise ValueError("`use_bisection` must be a non-negative integer.")

   calls= steps= 0

   if both:
      f_a= f(a)
      f_b= f(b)
      calls+= 2

      if abs(b - a) <= xtol:
         if abs(f_a) <= ftol:
            if verbose:
               print("Convergence achieved after %d steps and %d calls."
                 % (steps, calls))
            return a

         if abs(f_b) <= ftol:
            if verbose:
               print("Convergence achieved after %d steps and %d calls."
                 % (steps, calls))
            return b

   else:
      if abs(b - a) <= xtol:
         raise ValueError("When `both` equals `False`, the starting values of "
           "x (a and b) must differ by more than xtol= %f." % xtol)

      f_a= f(a)
      calls+= 1

      if abs(f_a) <= ftol:
         if verbose:
            print("Convergence achieved after %d steps and %d calls."
              % (steps, calls))
         return a

      f_b= f(b)
      calls+= 1

      if abs(f_b) <= ftol:
         if verbose:
            print("Convergence achieved after %d steps and %d calls."
              % (steps, calls))
         return b


   while True:
      steps+= 1
      if steps > max_steps:
         raise AlgorithmFailure("Limit of %d iterations has been reached."
           % max_steps)

      slope= (f_b - f_a) / float(b - a)
      if abs(slope) < min_slope:

         # The slope is essentially zero.  Attempt to fix this by shifting `b`
         # towards `a` (the new `b` is a weighted average of `a` and the old
         # `b`):
         b= 0.3679*a + 0.6321*b
         f_b= f(b)
         calls+= 1

         slope= (f_b - f_a)/float(b - a)
         if abs(slope) < min_slope:
            if abs(b - a) <= xtol:
               break
            raise AlgorithmFailure("At step %d, the slope is too close to zero!"
              % steps)

      # Calculate c= updated estimate of root location and f_c= corresponding
      # function value:
      c= b - f_b/slope
      f_c= f(c)
      calls+= 1

      if both:
         if abs(c - b) <= xtol and abs(f_c) <= ftol:
            break
      else:
         if abs(c - b) <= xtol or  abs(f_c) <= ftol:
            break

      if f_a * f_b < 0.0:

         # a and b bracket a root of the function.  (More accurately, `a` and
         # `b` bracket n roots, where n is an odd number).

         x_min= min(a, b)
         x_max= max(a, b)

         if not x_min < c < x_max:
            if verbose:
               print("At step %d, following %d calls, f(a) and f(b) have "
                 "opposite signs but c is not between a and b.  This should "
                 "not happen.  a=%f, b=%f, c=%f" % (steps, calls, a, b, c))
            continue

         x_rng= x_max - x_min

         if f_b * f_c >= 0.0:

            # a and b bracket a root.  c is between a and b.  b and c do
            # not bracket a root.  ==> a and c must bracket a root.  (All that
            # we know for certain is that a and c bracket an odd number of
            # roots).
            b, f_b= c, f_c

         else:

            # b and c bracket a root.
            a, f_a, b, f_b= b, f_b, c, f_c

      else:

         # It appears that a and b do not bracket a root of the function.
         a, f_a, b, f_b= b, f_b, c, f_c

   # end while True

   if verbose:
      print("Convergence achieved after %d steps and %d calls."
        % (steps, calls))

   return c


def find_root(f, a, b, ftol=1.e-6, xtol=1.e-6, both=True,
  contraction_factor=0.7071, max_steps=3000, min_slope=1.e-60,
  use_bisection=0, verbose=False):
   """
   OVERVIEW

   This function finds a root (zero) of a function `f` using an algorithm that
   is a hybrid of the modified Regula Falsi method (Ref. 2), which we henceforth
   call 'Regula Falsi', and bisection search.  Each iteration of the algorithm
   uses either bisection or Regula Falsi.  `f` is required to be everywhere
   continuous, but is not required to have a derivative anywhere.  (Compare with
   Newton's method, which requires that `f` have a derivative).

   The calling program provides the function, a pair of starting values `a` and
   `b` for the independent variable, and optionally other inputs as well (see
   below).  If `a` and `b` bracket a root--more precisely, f(a) and f(b) have
   opposite signs--or if at any point the algorithm discovers a pair of values
   that bracket at root, subsequent intervals will continue to bracket that
   root.  Note that this behavior is different from that of the secant method,
   which given an interval [a, b] that contains a root, can converge to a root
   outside that interval or fail to converge.

   In the initial `use_bisection` steps--this argument equals zero by default--
   the algorithm uses the bisection method.  (The default behavior is that the
   first step uses the Regula Falsi method).

   After any step in which bisection is used, the value of `use_bisection` is
   decremented by 1.  When `use_bisection` reaches zero, the algorithm switches
   to the Regula Falsi method.

   After any step in which the Regula Falsi method is used, and in which the
   starting values of x bracket a root, the lengths of the new and old intervals
   are compared.  If the ratio exceeds a threshold (see below), the algorithm
   will set `use_bisection` to 1, causing the bisection method to be used at the
   next step.


   INPUTS

   `f` is a function that takes a single real-valued argument and returns a
   single real values result.

   `a` and `b` are the starting values of the independent variable `x`.

   `ftol` and `xtol` are convergence thresholds; both default to 1.e-6.

   `both` is a bool value that specifies whether both convergence thresholds
   must be simultaneously satisfied; the default is `True`.  If `both` equals
   `False`, the algorithm stops as soon as either convergence threshold is
   satisfied.

   `contraction_factor` controls the application of the bisection method.  When
   a pair of consecutive x's brackets a root, and the ratio of the new
   interval's length to that of the old exceeds this factor, the algorithm
   applies the bisection method at the next step and then reverts to its normal
   behavior.  Allowed values for the `contraction_factor` are 0.5 to 1.  The
   default value is 0.7071= sqrt(2), chosen so that the algorithm will take no
   more than twice as many steps as bisection in the worst case.  This value
   appears to typically work well, but may be tuned as needed.

   `max_steps` is the maximum allowed number of iterations.  If convergence
   does not occur within this number of steps, an `AlgorithFailure` exception is
   raised.

   `min_slope` is the minimum absolute value of the slope that is considered
   adequately different from zero for the algorithm to proceed.  The default is
   1.e-60.

   `use_bisection` specifies the number of bisection steps to be performed
   before the algorithm tries the modified Regula Falsi method.  The default
   value is zero.

   `verbose`: Setting this input to `True` causes the function to display the
   numbers of function calls and iterations required for convergence.


   NOTE

   My goals in designing this function were the following:

   (1) to provide more reliable convergence than the secant method, and in
   particular, to provide guaranteed convergence to a root when f(a) and f(b)
   have opposite signs (which means that `a` and `b` bracket a root).

   (2) to provide faster convergence than bisection search when f() is
   close to a straight line, or becomes close to a straight line after a small
   number of bisection steps.

   (3) to provide faster convergence than the Modified Regula Falsi method when
   f() is irregular on a large scale but smooth at smaller scales.

   Unless the calling program specifies that the hybrid algorithm should use
   bisection search initially, f(a) and f(b) are not required to have opposite
   signs.

   If f(a) and f(b) have opposite signs, or after any step that produces a
   pair of x's (x_{i}, x_{i+1}) such that f(x_{i}) and f(x_{i+1}) have opposite
   signs, the algorithm tends to converge more rapidly than bisection search,
   but less rapidly than secant search.  As per Ref. 2, 'The safeguarding
   modification of the secant method that leads to the Modified Regula Falsi
   method promotes convergence, but destroys its fast convergence'.


   AUTHOR

   This algorithm was designed by Phillip M. Feldman.


   REFERENCES

   1. http://en.wikipedia.org/wiki/Secant_method

   2. Joanna M. Papakonstantinou and Richard A. Tapia, 'Origin and Evolution of
   the Secant Method in One Dimension', The American Mathematical Monthly, Vol.
   120, No. 6, June-July 2013, pp. 500-518.


   EXAMPLE #1

   Consider the following function, which has a simple root at x= 0.7:

      f(x)= 2(x-0.7) + 0.03(x-0.7)^3


   Code:

   from mymath import *

   def f(x):
      return 2.*(x-0.7) + 0.03*(x-0.7)**3

   print('Testing `find_root_bisection` ...')
   x= find_root_bisection(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
   print('%6f' % x)

   print('Testing `find_root_secant` ...')
   x= find_root_secant(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
   print('%6f' % x)

   print('Testing `find_root` with default value of `contraction_factor` ...')
   x= find_root(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
   print('%6f' % x)

   print('Testing `find_root` with `contraction_factor= 1.0` ...')
   x= find_root(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, contraction_factor=1.0,
     verbose=True)
   print('%6f' % x)


   Output:

   Testing `find_root_bisection` ...
   Convergence achieved after 23 steps and 25 calls.
   0.700000

   Testing `find_root_secant` ...
   Convergence achieved after 4 steps and 6 calls.
   0.700000

   Testing `find_root` with default value of `contraction_factor` ...
   Convergence achieved after 14 steps and 16 calls.
   0.700000

   Testing `find_root` with `contraction_factor=1.0` ...
   Convergence achieved after 11 steps and 13 calls.
   0.700000


   Discussion:

   `find_root_secant` achieves the most rapid convergence, followed by
   `find_root`, with `find_root_bisection` in last place.  The performance of
   `find_root` is slightly improved by using a non-default value of the
   `contraction_factor` parameter.


   EXAMPLE #2

   Consider the following function, which has a root of multiplicity 4 at
   x= 0.7:

      f(x)= (x-0.7)^4

   The function touches but does not cross the x-axis.


   Code:

   from mymath import *

   def f(x):
      return (x-0.7)**4

   print('\nTesting `find_root_bisection` ...')
   try:
      x= find_root_bisection(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
      print('%6f' % x)
   except Exception as ex:
      print(ex)

   print('\nTesting `find_root_secant` ...')
   x= find_root_secant(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
   print('%6f' % x)

   print('\nTesting `find_root` with default value of `contraction_factor` ...')
   x= find_root(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
   print('%6f' % x)


   Output:

   Testing `find_root_bisection` ...
   This function requires that f(a) and f(b) have opposite signs.  You may try
   using `find_root`, which does not have this restriction.

   Testing `find_root_secant` ...
   Convergence achieved after 52 steps and 54 calls.
   0.699996

   Testing `find_root` with default value of `contraction_factor` ...
   Convergence achieved after 52 steps and 54 calls.
   0.699996


   Discussion:

   `find_root_secant` and `find_root` converge in exactly the same number of
   steps, while `find_root_bisection` fails because f(a) and f(b) are both
   positive.  So far, `find_root` appears to have no advantage over
   `find_root_secant`.


   EXAMPLE #3

   Consider the following function, which is piecewise-linear with three
   sections:

         -1, x < 1
   f(x)=  0, -1 <= x <= 1
          1, x > -1


   Code:

   from mymath import *

   def f(x):
      return numpy.clip(x, -1, 1)

   print('\nTesting `find_root_bisection` ...')
   try:
      x= find_root_bisection(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
      print('%6f' % x)
   except Exception as ex:
      print(ex)

   print('\nTesting `find_root_secant` ...')
   try:
      x= find_root_secant(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
      print('%6f' % x)
   except Exception as ex:
      print(ex)

   print('\nTesting `find_root` with default value of `contraction_factor` ...')
   x= find_root(f, 0.6, 6.0, xtol=1e-6, ftol=1e-6, verbose=True)
   print('%6f' % x)


   Output:

   Testing `find_root_bisection` ...
   This function requires that f(a) and f(b) have opposite signs.  You may try
   using `find_root`, which does not have this restriction.

   Testing `find_root_secant` ...
   At step 5, the slope is too close to zero!

   Testing `find_root` with default value of `contraction_factor` ...
   Convergence achieved after 6 steps and 8 calls.
   0.000000


   Discussion:

   `find_root_bisection` and `find_root_secant` both fail.  `find_root` achieves
   rapid convergence.


   EXAMPLE #4:

   Consider the following function:

      f(x)= x e^(-|x|)

   f() is everywhere continuous and differentiable.  It has peaks on either side
   of the single root, which occurs at x= 0.


   Code:

   from mymath import *

   def f(x):
      return x * exp(-abs(x))

   print('\nTesting `find_root_bisection` ...')
   x= find_root_bisection(f, -0.5, 10., xtol=1e-6, ftol=1e-6, verbose=True)
   print('%6f' % x)

   print('\nTesting `find_root_secant` ...')
   try:
      x= find_root_secant(f, -0.5, 10., xtol=1e-6, ftol=1e-6, verbose=True)
      print('%6f' % x)
   except Exception as ex:
      print(ex)

   print('\nTesting `find_root` with default value of `contraction_factor` ...')
   x= find_root(f, -0.5, 10., xtol=1e-6, ftol=1e-6, verbose=True)
   print('%6f' % x)


   Output:

   Testing `find_root_bisection` ...
   Convergence achieved after 24 steps and 26 calls.
   -0.000000

   Testing `find_root_secant` ...
   At step 191, the slope is too close to zero!

   Testing `find_root` with default value of `contraction_factor` ...
   After step 1.  Converging too slowly.
     ==> Temporarily switching to bisection.
   At step 2, following 3 calls, using bisection.
   After step 3.  Converging too slowly.
     ==> Temporarily switching to bisection.
   At step 4, following 5 calls, using bisection.
   Convergence achieved after 11 steps and 13 calls.
   -0.000000


   Discussion:

   `find_root_secant` fails.  `find_root` and `find_root_bisection` both
   converge, but `find_root` requires half as many function calls.


   Conclusions:

   1. `find_root` and `find_root_secant` are more flexible and convenient to
   use than `find_root_bisection` because they do not require that the initial
   values of x bracket a root.

   2. `find_root` is more reliable than `find_root_secant`.

   3. In situations where `find_root` and `find_root_bisection` are both
   applicable, `find_root` often converges more rapidly.
   """

   if xtol <= 0.0:
      raise ValueError("If specified, `xtol` must be positive.")
   if ftol <= 0.0:
      raise ValueError("If specified, `ftol` must be positive.")
   if not 0.5 <= contraction_factor <= 1.0:
      raise ValueError("If specified, `contraction` factor must be in the "
        "interval [0, 1].")

   if not isinstance(use_bisection, int) or use_bisection < 0:
      raise ValueError("`use_bisection` must be a non-negative integer.")

   calls= steps= 0

   if both:
      f_a= f(a)
      f_b= f(b)
      calls+= 2

      if abs(b - a) <= xtol:
         if abs(f_a) <= ftol:
            if verbose:
               print("Convergence achieved after %d steps and %d calls."
                 % (steps, calls))
            return a

         if abs(f_b) <= ftol:
            if verbose:
               print("Convergence achieved after %d steps and %d calls."
                 % (steps, calls))
            return b

   else:
      if abs(b - a) <= xtol:
         raise ValueError("When `both` equals `False`, the starting values of "
           "x (a and b) must differ by more than xtol= %f." % xtol)

      f_a= f(a)
      calls+= 1

      if abs(f_a) <= ftol:
         if verbose:
            print("Convergence achieved after %d steps and %d calls."
              % (steps, calls))
         return a

      f_b= f(b)
      calls+= 1

      if abs(f_b) <= ftol:
         if verbose:
            print("Convergence achieved after %d steps and %d calls."
              % (steps, calls))
         return b


   if use_bisection and f_a * f_b > 0.0:
      raise AlgorithmFailure("When `use_bisection` is positive, f(a) and "
        "f(b) must have opposite signs.")


   while True:
      steps+= 1
      if steps > max_steps:
         raise AlgorithmFailure("Limit of %d iterations has been reached."
           % max_steps)

      if use_bisection:

         # The current step uses the bisection method.

         if verbose:
            print("At step %d, following %d calls, using bisection."
              % (steps, calls))

         c= 0.5 * (a+b)
         f_c= f(c)
         calls+= 1

         if f_a * f_c < 0.0:
            b, f_b= c, f_c

         else:

            # In the following statement, it might seem that we are doing
            # some unnecessary copying, but checking the convergence
            # conditions requires that we keep track of the order in which
            # results were generated.
            a, f_a, b, f_b= b, f_b, c, f_c

         use_bisection-= 1
         continue

      # end if use_bisection


      # The current step uses the modified Regula Falsi method.

      slope= (f_b - f_a) / float(b - a)
      if abs(slope) < min_slope:

         # The slope is essentially zero.  Attempt to fix this by shifting `b`
         # towards `a` (the new `b` is a weighted average of `a` and the old
         # `b`):
         b= 0.3679*a + 0.6321*b
         f_b= f(b)
         calls+= 1

         slope= (f_b - f_a)/float(b - a)
         if abs(slope) < min_slope:
            if abs(b - a) <= xtol:
               break
            raise AlgorithmFailure("At step %d, the slope is too close to zero!"
              % steps)

      # Calculate c= updated estimate of root location and f_c= corresponding
      # function value:
      c= b - f_b/slope
      f_c= f(c)
      calls+= 1

      if both:
         if abs(c - b) <= xtol and abs(f_c) <= ftol:
            break
      else:
         if abs(c - b) <= xtol or  abs(f_c) <= ftol:
            break

      if f_a * f_b < 0.0:

         # a and b bracket a root of the function.  (More accurately, `a` and
         # `b` bracket n roots, where n is an odd number).

         x_min= min(a, b)
         x_max= max(a, b)

         if not x_min < c < x_max:
            if verbose:
               print("At step %d, following %d calls, f(a) and f(b) have "
                 "opposite signs but c is not between a and b.  This should "
                 "not happen.  a=%f, b=%f, c=%f" % (steps, calls, a, b, c))
            use_bisection= 1
            continue

         x_rng= x_max - x_min

         if f_b * f_c >= 0.0:

            # a and b bracket a root.  c is between a and b.  b and c do
            # not bracket a root.  ==> a and c must bracket a root.  (All that
            # we know for certain is that a and c bracket an odd number of
            # roots).
            b, f_b= c, f_c

         else:

            # b and c bracket a root.
            a, f_a, b, f_b= b, f_b, c, f_c

         if max(a, b) - min(a, b) > contraction_factor * x_rng:

            # We are converging too slowly.  ==> Perform a bisection step.
            if verbose:
               print("After step %d.  Converging too slowly.  ==> Temporarily "
                 "switching to bisection." % steps)
            use_bisection= 1

      else:

         # It appears that a and b do not bracket a root of the function.
         a, f_a, b, f_b= b, f_b, c, f_c

   # end while True

   if verbose:
      print("Convergence achieved after %d steps and %d calls."
        % (steps, calls))

   return c

# end def find_root