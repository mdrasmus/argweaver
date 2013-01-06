# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

# include "r8lib.h"

extern "C" {

/******************************************************************************/

int i4_log_10 ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_LOG_10 returns the integer part of the logarithm base 10 of an I4.

  Example:

        I  I4_LOG_10
    -----  --------
        0    0
        1    0
        2    0
        9    0
       10    1
       11    1
       99    1
      100    2
      101    2
      999    2
     1000    3
     1001    3
     9999    3
    10000    4

  Discussion:

    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number whose logarithm base 10 is desired.

    Output, int I4_LOG_10, the integer part of the logarithm base 10 of
    the absolute value of X.
*/
{
  int i_abs;
  int ten_pow;
  int value;

  if ( i == 0 )
  {
    value = 0;
  }
  else
  {
    value = 0;
    ten_pow = 10;

    i_abs = abs ( i );

    while ( ten_pow <= i_abs )
    {
      value = value + 1;
      ten_pow = ten_pow * 10;
    }
  }
  return value;
}
/******************************************************************************/

int i4_max ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MAX returns the maximum of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, are two integers to be compared.

    Output, int I4_MAX, the larger of I1 and I2.
*/
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:

    I4_MIN returns the smaller of two I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 August 2006

  Author:

    John Burkardt

  Parameters:

    Input, int I1, I2, two integers to be compared.

    Output, int I4_MIN, the smaller of I1 and I2.
*/
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
/******************************************************************************/

int i4_modp ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_MODP returns the nonnegative remainder of I4 division.

  Discussion:

    If
      NREM = I4_MODP ( I, J )
      NMULT = ( I - NREM ) / J
    then
      I = J * NMULT + NREM
    where NREM is always nonnegative.

    The MOD function computes a result with the same sign as the
    quantity being divided.  Thus, suppose you had an angle A,
    and you wanted to ensure that it was between 0 and 360.
    Then mod(A,360) would do, if A was positive, but if A
    was negative, your result would be between -360 and 0.

    On the other hand, I4_MODP(A,360) is between 0 and 360, always.

  Example:

        I         J     MOD  I4_MODP   I4_MODP Factorization

      107        50       7       7    107 =  2 *  50 + 7
      107       -50       7       7    107 = -2 * -50 + 7
     -107        50      -7      43   -107 = -3 *  50 + 43
     -107       -50      -7      43   -107 =  3 * -50 + 43

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 January 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number to be divided.

    Input, int J, the number that divides I.

    Output, int I4_MODP, the nonnegative remainder when I is
    divided by J.
*/
{
  int value;

  if ( j == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4_MODP - Fatal error!\n" );
    fprintf ( stderr, "  I4_MODP ( I, J ) called with J = %d\n", j );
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
/******************************************************************************/

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    I4_POWER returns the value of I^J.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the base and the power.  J should be nonnegative.

    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
/******************************************************************************/

int i4_sign ( int i )

/******************************************************************************/
/*
  Purpose:

    I4_SIGN returns the sign of an I4.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2007

  Author:

    John Burkardt

  Parameters:

    Input, int I, the integer whose sign is desired.

    Output, int I4_SIGN, the sign of I.
*/
{
  int value;

  if ( i < 0 )
  {
    value = -1;
  }
  else
  {
    value = 1;
  }
  return value;
}
/******************************************************************************/

int i4_uniform_ab ( int a, int b, int *seed )

/******************************************************************************/
/*
  Purpose:

    I4_UNIFORM_AB returns a pseudorandom I4 scaled to [A,B].

  Discussion:

    The pseudorandom number should be uniformly distributed
    between A and B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2006

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, int I4_UNIFORM_AB, a number between A and B.
*/
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
/*
  Scale R to lie between A-0.5 and B+0.5.
*/
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
/*
  Use rounding to convert R to an integer between A and B.
*/
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
/******************************************************************************/

int i4_wrap ( int ival, int ilo, int ihi )

/******************************************************************************/
/*
  Purpose:

    I4_WRAP forces an I4 to lie between given limits by wrapping.

  Example:

    ILO = 4, IHI = 8

    I   Value

    -2     8
    -1     4
     0     5
     1     6
     2     7
     3     8
     4     4
     5     5
     6     6
     7     7
     8     8
     9     4
    10     5
    11     6
    12     7
    13     8
    14     4

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int IVAL, an integer value.

    Input, int ILO, IHI, the desired bounds for the integer value.

    Output, int I4_WRAP, a "wrapped" version of IVAL.
*/
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
/******************************************************************************/

double i4int_to_r8int ( int imin, int imax, int i, double rmin, double rmax )

/******************************************************************************/
/*
  Purpose:

    I4INT_TO_R8INT maps an I4 interval to an R8 interval.

  Discussion:

    The formula is

      R := RMIN + ( RMAX - RMIN ) * ( I - IMIN ) / ( IMAX - IMIN )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2005

  Author:

    John Burkardt

  Parameters:

    Input, int IMIN, IMAX, the range.

    Input, int I, the integer to be converted.

    Input, double RMIN, RMAX, the range.

    Output, double R, the corresponding value in [RMIN,RMAX].
*/
{
  double r;

  if ( imax == imin )
  {
    r = 0.5 * ( rmin + rmax );
  }
  else
  {
    r = ( ( double ) ( imax - i        ) * rmin
        + ( double ) (        i - imin ) * rmax )
        / ( double ) ( imax     - imin );
  }

  return r;
}
/******************************************************************************/

void i4vec_copy ( int n, int a1[], int a2[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_COPY copies an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 April 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, int A1[N], the vector to be copied.

    Input, int A2[N], the copy of A1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
/******************************************************************************/

int *i4vec_indicator_new ( int n )

/******************************************************************************/
/*
  Purpose:

    I4VEC_INDICATOR_NEW sets an I4VEC to the indicator vector.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Output, int I4VEC_INDICATOR_NEW[N], the array.
*/
{
  int *a;
  int i;

  a = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }
  return a;
}
/******************************************************************************/

void i4vec_permute ( int n, int p[], int base, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PERMUTE permutes an I4VEC in place.

  Discussion:

    An I4VEC is a vector of I4's.

    This routine permutes an array of integer "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.

  Example:

    Input:

      N = 5
      P = (   1,   3,   4,   0,   2 )
      A = (   1,   2,   3,   4,   5 )

    Output:

      A    = (   2,   4,   5,   1,   3 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of objects.

    Input, int P[N], the permutation.  P(I) = J means
    that the I-th element of the output array should be the J-th
    element of the input array.

    Input, int BASE, is 0 for a 0-based permutation and 1 for
    a 1-based permutation.

    Input/output, int A[N], the array to be permuted.
*/
{
  int a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm_check ( n, p, base ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "I4VEC_PERMUTE - Fatal error!\n" );
    fprintf ( stderr, "  PERM_CHECK rejects this permutation.\n" );
    exit ( 1 );
  }
/*
  In order for the sign negation trick to work, we need to assume that the
  entries of P are strictly positive.  Presumably, the lowest number is BASE.
  So temporarily add 1-BASE to each entry to force positivity.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }
/*
  Search for the next element of the permutation that has not been used.
*/
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;
/*
  Copy the new value into the vacated entry.
*/
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          fprintf ( stderr, "\n" );
          fprintf ( stderr, "I4VEC_PERMUTE - Fatal error!\n" );
          fprintf ( stderr, "  Entry IPUT = %d of the permutation has\n", iput );
          fprintf ( stderr, "  an illegal value IGET = %d.\n", iget );
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }
/*
  Restore the signs of the entries.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
/*
  Restore the base of the entries.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 + base;
  }

  return;
}
/******************************************************************************/

void i4vec_print ( int n, int a[], char *title )

/******************************************************************************/
/*
  Purpose:

    I4VEC_PRINT prints an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 November 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, int A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %6d: %8d\n", i, a[i] );
  }
  return;
}
/******************************************************************************/

void i4vec_zero ( int n, int a[] )

/******************************************************************************/
/*
  Purpose:

    I4VEC_ZERO zeroes an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Output, int A[N], a vector of zeroes.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return;
}
/******************************************************************************/

int *i4vec_zero_new ( int n )

/******************************************************************************/
/*
  Purpose:

    I4VEC_ZERO_NEW creates and zeroes an I4VEC.

  Discussion:

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
*/
{
  int *a;
  int i;

  a = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return a;
}
/******************************************************************************/

double *legendre_zeros ( int order )

/******************************************************************************/
/*
  Purpose:

    LEGENDRE_ZEROS returns the zeros of the Legendre polynomial of degree N.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 June 2011

  Author:

    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt.

  Reference:

    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.

  Parameters:

    Input, int ORDER, the order.
    ORDER must be greater than 0.

    Output, double LEGENDRE_ZEROS[ORDER], the zeros.
*/
{
  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pi = 3.141592653589793;
  double pk;
  double pkm1;
  double pkp1;
  double t;
  double u;
  double v;
  double x0;
  double *xtab;
  double xtemp;

  xtab = ( double * ) malloc ( order * sizeof ( double ) );

  e1 = ( double ) ( order * ( order + 1 ) );

  m = ( order + 1 ) / 2;

  for ( i = 1; i <= m; i++ )
  {
    mp1mi = m + 1 - i;

    t = ( double ) ( 4 * i - 1 ) * pi / ( double ) ( 4 * order + 2 );

    x0 = cos ( t ) * ( 1.0 - ( 1.0 - 1.0 / ( double ) ( order ) ) 
      / ( double ) ( 8 * order * order ) );

    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= order; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }

    d1 = ( double ) ( order ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;
/*
  Initial approximation H:
*/
    h = -u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn / ( 3.0 * dpn ) ) ) );
/*
  Refine H using one step of Newton's method:
*/
    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0 
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;

    xtemp = x0 + h;

    xtab[mp1mi-1] = xtemp;

    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );
  }

  if ( ( order % 2 ) == 1 )
  {
    xtab[0] = 0.0;
  }
/*
  Shift the data up.
*/
  nmove = ( order + 1 ) / 2;
  ncopy = order - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = order + 1 - i;
    xtab[iback-1] = xtab[iback-ncopy-1];
  }
/*
  Reflect values for the negative abscissas.
*/
  for ( i = 1; i <= order - nmove; i++ )
  {
    xtab[i-1] = - xtab[order-i];
  }

  return xtab;
}
/******************************************************************************/

int perm_check ( int n, int p[], int base )

/******************************************************************************/
/*
  Purpose:

    PERM_CHECK checks that a vector represents a permutation.

  Discussion:

    The routine verifies that each of the integers from BASE to
    to BASE+N-1 occurs among the N entries of the permutation.

    Set the input quantity BASE to 0, if P is a 0-based permutation,
    or to 1 if P is a 1-based permutation.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries.

    Input, int P[N], the array to check.

    Input, int BASE, the index base.

    Output, int PERM_CHECK, is TRUE if the permutation is OK.
*/
{
  int found;
  int i;
  int seek;

  for ( seek = base; seek < base + n; seek++ )
  {
    found = 0;

    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == seek )
      {
        found = 1;
        break;
      }
    }

    if (( found = 0 ))
    {
      return 0;
    }

  }

  return 1;
}
/******************************************************************************/

int *perm_uniform_new ( int n, int base, int *seed )

/******************************************************************************/
/*
  Purpose:

    PERM_UNIFORM_NEW selects a random permutation of N objects.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 June 2010

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the number of objects to be permuted.

    Input, int BASE, is 0 for a 0-based permutation and 1 for
    a 1-based permutation.

    Input/output, int *SEED, a seed for the random number generator.

    Output, int PERM_UNIFORM_NEW[N], a permutation of
    (BASE, BASE+1, ..., BASE+N-1).
*/
{
  int i;
  int j;
  int k;
  int *p;

  p = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    p[i] = i + base;
  }

  for ( i = 0; i < n; i++ )
  {
    j = i4_uniform_ab ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }

  return p;
}
/******************************************************************************/

int r4_nint ( float x )

/******************************************************************************/
/*
  Purpose:

    R4_NINT returns the nearest integer to an R4.

  Example:

        X         R4_NINT

      1.3         1
      1.4         1
      1.5         1 or 2
      1.6         2
      0.0         0
     -0.7        -1
     -1.1        -1
     -1.6        -2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, float X, the value.

    Output, int R4_NINT, the nearest integer to X.
*/
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = - 1;
  }
  else
  {
    s = + 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ABS returns the absolute value of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the quantity whose absolute value is desired.

    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
/******************************************************************************/

double r8_acos ( double c )

/******************************************************************************/
/*
  Purpose:

    R8_ACOS computes the arc cosine function, with argument truncation.

  Discussion:

    If you call your system ACOS routine with an input argument that is
    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
    This routine truncates arguments outside the range.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    13 June 2002

  Author:

    John Burkardt

  Parameters:

    Input, double C, the argument, the cosine of an angle.

    Output, double R8_ACOS, an angle whose cosine is C.
*/
{
  double angle;
  double pi = 3.141592653589793;

  if ( c <= -1.0 )
  {
    angle = pi;
  } 
  else if ( 1.0 <= c )
  {
    angle = 0.0;
  }
  else
  {
    angle = acos ( c );
  }
  return angle;
}
/******************************************************************************/

double r8_add ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_ADD returns the sum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the numbers to be added.

    Output, double R8_ADD, the sum of X and Y.
*/
{
  double value;

  value = x + y;

  return value;
}
/******************************************************************************/

double r8_aint ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_AINT truncates an R8 argument to an integer.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    1 September 2011

  Author:

    John Burkardt.

  Parameters:

    Input, double X, the argument.

    Output, double R8_AINT, the truncated version of X.
*/
{
  double value;

  if ( x < 0.0E+00 )
  {
    value = - ( double ) ( ( int ) ( r8_abs ( x ) ) );
  }
  else
  {
    value =   ( double ) ( ( int ) ( r8_abs ( x ) ) );
  }

  return value;
}
/******************************************************************************/

double r8_atan ( double y, double x )

/******************************************************************************/
/*
  Purpose:

    R8_ATAN computes the inverse tangent of the ratio Y / X.

  Discussion:

    R8_ATAN returns an angle whose tangent is ( Y / X ), a job which
    the built in functions ATAN and ATAN2 already do.

    However:

    * R8_ATAN always returns a positive angle, between 0 and 2 PI,
      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
      and [-PI,+PI] respectively;

    * R8_ATAN accounts for the signs of X and Y, (as does ATAN2).  The ATAN
     function by contrast always returns an angle in the first or fourth
     quadrants.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 May 2003

  Author:

    John Burkardt

  Parameters:

    Input, double Y, X, two quantities which represent the tangent of
    an angle.  If Y is not zero, then the tangent is (Y/X).

    Output, double R8_ATAN, an angle between 0 and 2 * PI, whose tangent is
    (Y/X), and which lies in the appropriate quadrant so that the signs
    of its cosine and sine match those of X and Y.
*/
{
  double abs_x;
  double abs_y;
  double pi = 3.141592653589793;
  double theta;
  double theta_0;
/*
  Special cases:
*/
  if ( x == 0.0 )
  {
    if ( 0.0 < y )
    {
      theta = pi / 2.0;
    }
    else if ( y < 0.0 )
    {
      theta = 3.0 * pi / 2.0;
    }
    else if ( y == 0.0 )
    {
      theta = 0.0;
    }
  }
  else if ( y == 0.0 )
  {
    if ( 0.0 < x )
    {
      theta = 0.0;
    }
    else if ( x < 0.0 )
    {
      theta = pi;
    }
  }
/*
  We assume that ATAN2 is correct when both arguments are positive.
*/
  else
  {
    abs_y = fabs ( y );
    abs_x = fabs ( x );

    theta_0 = atan2 ( abs_y, abs_x );

    if ( 0.0 < x && 0.0 < y )
    {
      theta = theta_0;
    }
    else if ( x < 0.0 && 0.0 < y )
    {
      theta = pi - theta_0;
    }
    else if ( x < 0.0 && y < 0.0 )
    {
      theta = pi + theta_0;
    }
    else if ( 0.0 < x && y < 0.0 )
    {
      theta = 2.0 * pi - theta_0;
    }
  }

  return theta;
}
/******************************************************************************/

double r8_cas ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_CAS returns the "casine" of an R8.

  Discussion:

    The "casine", used in the discrete Hartley transform, is abbreviated
    CAS(X), and defined by:

      CAS(X) = cos ( X ) + sin( X )
             = sqrt ( 2 ) * sin ( X + pi/4 )
             = sqrt ( 2 ) * cos ( X - pi/4 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose casine is desired.

    Output, double R8_CAS, the casine of X, which will be between
    plus or minus the square root of 2.
*/
{
  double value;

  value = cos ( x ) + sin ( x );

  return value;
}
/******************************************************************************/

double r8_ceiling ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_CEILING rounds an R8 up to the nearest integral R8.

  Discussion:

    The "ceiling" of X is the value of X rounded towards plus infinity.

  Example:

    X        R8_CEILING(X)

   -1.1      -1.0
   -1.0      -1.0
   -0.9       0.0
   -0.1       0.0
    0.0       0.0
    0.1       1.0
    0.9       1.0
    1.0       1.0
    1.1       2.0
    2.9       3.0
    3.0       3.0
    3.14159   4.0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose ceiling is desired.

    Output, double R8_CEILING, the ceiling of X.
*/
{
  double value;

  value = ( double ) ( ( int ) x );

  if ( value < x )
  {
    value = value + 1.0;
  }

  return value;
}
/******************************************************************************/

double r8_choose ( int n, int k )

/******************************************************************************/
/*
  Purpose:

    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.

  Discussion:

    The value is calculated in such a way as to avoid overflow and
    roundoff.  The calculation is done in R8 arithmetic.

    The formula used is:

      C(N,K) = N! / ( K! * (N-K)! )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Reference:

    ML Wolfson, HV Wright,
    Algorithm 160:
    Combinatorial of M Things Taken N at a Time,
    Communications of the ACM,
    Volume 6, Number 4, April 1963, page 161.

  Parameters:

    Input, int N, K, the values of N and K.

    Output, double R8_CHOOSE, the number of combinations of N
    things taken K at a time.
*/
{
  int i;
  int mn;
  int mx;
  double value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }

  return value;
}
/******************************************************************************/

double r8_chop ( int place, double x )

/******************************************************************************/
/*
  Purpose:

    R8_CHOP chops an R8 to a given number of binary places.

  Example:

    3.875 = 2 + 1 + 1/2 + 1/4 + 1/8.

    The following values would be returned for the 'chopped' value of
    3.875:

    PLACE  Value

       1      2
       2      3     = 2 + 1
       3      3.5   = 2 + 1 + 1/2
       4      3.75  = 2 + 1 + 1/2 + 1/4
       5+     3.875 = 2 + 1 + 1/2 + 1/4 + 1/8

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int PLACE, the number of binary places to preserve.
    PLACE = 0 means return the integer part of X.
    PLACE = 1 means return the value of X, correct to 1/2.
    PLACE = 2 means return the value of X, correct to 1/4.
    PLACE = -1 means return the value of X, correct to 2.

    Input, double X, the number to be chopped.

    Output, double R8_CHOP, the chopped number.
*/
{
  double fac;
  int temp;
  static double two = 2.0;
  double value;

  temp = ( int ) ( r8_log_2 ( x ) );
  fac = pow ( two, temp - place + 1 );
  value = ( double ) ( ( int ) ( x / fac ) ) * fac;

  return value;
}
/******************************************************************************/

double r8_csc ( double theta )

/******************************************************************************/
/*
  Purpose:

    R8_CSC returns the cosecant of X.

  Discussion:

    R8_CSC ( THETA ) = 1.0 / SIN ( THETA )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, double THETA, the angle, in radians, whose cosecant is desired.
    It must be the case that SIN ( THETA ) is not zero.

    Output, double R8_CSC, the cosecant of THETA.
*/
{
  double value;

  value = sin ( theta );

  if ( value == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_CSC - Fatal error!\n" );
    fprintf ( stderr, "  Cosecant undefined for THETA = %g\n", theta );
    exit ( 1 );
  }

  value = 1.0 / value;

  return value;
}
/******************************************************************************/

double r8_cube_root ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_CUBE_ROOT returns the cube root of an R8.

  Discussion:

    This routine is designed to avoid the possible problems that can occur
    when formulas like 0.0**(1/3) or (-1.0)**(1/3) are to be evaluated.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose cube root is desired.

    Output, double R8_CUBE_ROOT, the cube root of X.
*/
{
  double e;
  double value;

  e = 1.0 / 3.0;

  if ( 0.0 < x )
  {
    value = pow ( x, e );
  }
  else if ( x == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = - pow ( - x , e );
  }

  return value;
}
/******************************************************************************/

double r8_diff ( double x, double y, int n )

/******************************************************************************/
/*
  Purpose:

    R8_DIFF computes (X-Y) to a specified accuracy.

  Discussion:

    The user controls how many binary digits of accuracy
    are to be used.

    N determines the accuracy of the value.  If N = 10,
    for example, only 11 binary places will be used in the arithmetic.
    In general, only N+1 binary places will be used.

    N may be zero.  However, a negative value of N should
    not be used, since this will cause both X and Y to look like 0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the two values whose difference is desired.

    Input, int N, the number of binary digits to use.

    Output, double R8_DIFF, the value of X-Y.
*/
{
  double cx;
  double cy;
  double pow2;
  double size;
  static double two = 2.0;
  double value;

  if ( x == y )
  {
    value = 0.0;
    return value;
  }

  pow2 = pow ( two, n );
/*
  Compute the magnitude of X and Y, and take the larger of the
  two.  At least one of the two values is not zero!
*/
  size = r8_max ( fabs ( x ), fabs ( y ) );
/*
  Make normalized copies of X and Y.  One of the two values will
  actually be equal to 1.
*/
  cx = x / size;
  cy = y / size;
/*
  Here's where rounding comes in.  We know that the larger of the
  the two values equals 1.  We multiply both values by 2**N,
  where N+1 is the number of binary digits of accuracy we want
  to use, truncate the values, and divide back by 2**N.
*/
  cx = ( double ) ( ( int ) ( cx * pow2 + 0.5 * r8_sign ( cx ) ) ) / pow2;
  cy = ( double ) ( ( int ) ( cy * pow2 + 0.5 * r8_sign ( cy ) ) ) / pow2;
/*
  Take the difference now.
*/
  value = cx - cy;
/*
  Undo the scaling.
*/
  value = value * size;

  return value;
}
/******************************************************************************/

int r8_digit ( double x, int idigit )

/******************************************************************************/
/*
  Purpose:

    R8_DIGIT returns a particular decimal digit of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose IDIGIT-th decimal digit is desired.
    Note that if X is zero, all digits will be returned as 0.

    Input, int IDIGIT, the position of the desired decimal digit.
    A value of 1 means the leading digit, a value of 2 the second digit
    and so on.

    Output, int R8_DIGIT, the value of the IDIGIT-th decimal digit of X.
*/
{
  int digit;
  int i;
  int ival;

  if ( x == 0.0 )
  {
    digit = 0;
    return digit;
  }

  if ( idigit <= 0 )
  {
    digit = 0;
    return digit;
  }
/*
  Force X to lie between 1 and 10.
*/
  x = fabs ( x );

  while ( x < 1.0 )
  {
    x = x * 10.0;
  }

  while ( 10.0 <= x )
  {
    x = x / 10.0;
  }

  for ( i = 1; i <= idigit; i++ )
  {
    ival = ( int ) ( x );
    x = ( x - ( double ) ival ) * 10.0;
  }

  digit = ival;

  return digit;
}
/******************************************************************************/

double r8_divide_i4 ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    R8_DIVIDE_I4 returns an I4 fraction as an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the numerator and denominator.

    Output, double R8_DIVIDE_I4, the value of (I/J).
*/
{
  double value;

  value = ( double ) ( i ) / ( double ) ( j );

  return value;
}
/******************************************************************************/

double r8_epsilon ( void )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
  static double value = 2.220446049250313E-016;

  return value;
}
/******************************************************************************/

double r8_epsilon_compute ( void )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON_COMPUTE computes the R8 round off unit.

  Discussion:

    The round off unit is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON_COMPUTE, the R8 round-off unit.
*/
{
  double one;
  double temp;
  double test;
  static double value = 0.0;

  if ( value == 0.0 )
  {
    one = ( double ) ( 1 );

    value = one;
    temp = value / 2.0;
    test = r8_add ( one, temp );

    while ( one < test )
    {
      value = temp;
      temp = value / 2.0;
      test = r8_add ( one, temp );
    }
  }

  return value;
}
/******************************************************************************/

double r8_exp ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_EXP computes the exponential function, avoiding overflow and underflow.

  Discussion:

    My experience with the G95 compiler has included many unpleasant
    floating point exceptions when very small arguments are given to
    the exponential function.

    This routine is designed to avoid such problems.

    Ideally, the rule would be:

                    X <= log ( TINY ) => R8_EXP ( X ) = 0
    log ( HUGE ) <= X                 => R8_EXP ( X ) = HUGE

    However, the G95 math library seems to produce infinity for
    EXP ( LOG ( HUGE ( X ) ), rather than HUGE ( X ), so we've
    included a fudge factor.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the exponential function.

    Output, double R8_EXP, the value of exp ( X ).
*/
{
  double log_max = 709.711;
  double log_min = -708.467;
  double value;

  if ( x <= log_min )
  {
    value = 0.0;
  }
  else if ( x < log_max )
  {
    value = exp ( x );
  }
  else
  {
    value = r8_huge ( );
  }

  return value;
}
/******************************************************************************/

double r8_factorial ( int n )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL computes the factorial of N.

  Discussion:

    factorial ( N ) = product ( 1 <= I <= N ) I

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the argument of the factorial function.
    If N is less than 1, the function value is returned as 1.

    Output, double R8_FACTORIAL, the factorial of N.
*/
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
/******************************************************************************/

double r8_factorial2 ( int n )

/******************************************************************************/
/*
  Purpose:

    R8_FACTORIAL2 computes the double factorial function.

  Discussion:

    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)

  Example:

     N    Factorial2(N)

     0     1
     1     1
     2     2
     3     3
     4     8
     5    15
     6    48
     7   105
     8   384
     9   945
    10  3840

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the argument of the double factorial
    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.

    Output, double R8_FACTORIAL2, the value of Factorial2(N).
*/
{
  int n_copy;
  double value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( double ) n_copy;
    n_copy = n_copy - 2;
  }

  return value;
}
/******************************************************************************/

double r8_floor ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_FLOOR rounds an R8 down to the nearest integral R8.

  Example:

    X        R8_FLOOR(X)

   -1.1      -2.0
   -1.0      -1.0
   -0.9      -1.0
   -0.1      -1.0
    0.0       0.0
    0.1       0.0
    0.9       0.0
    1.0       1.0
    1.1       1.0
    2.9       2.0
    3.0       3.0
    3.14159   3.0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose floor is desired.

    Output, double R8_FLOOR, the floor of X.
*/
{
  double value;

  value = ( double ) ( ( int ) x );

  if ( x < value )
  {
    value = value - 1.0;
  }

  return value;
}
/******************************************************************************/

double r8_fraction ( int i, int j )

/******************************************************************************/
/*
  Purpose:

    R8_FRACTION uses real arithmetic on an integer ratio.

  Discussion:

    Given integer variables I and J, both FORTRAN and C will evaluate
    an expression such as "I/J" using what is called "integer division",
    with the result being an integer.  It is often convenient to express
    the parts of a fraction as integers but expect the result to be computed
    using real arithmetic.  This function carries out that operation.

  Example:

       I     J   I/J  R8_FRACTION

       1     2     0  0.5
       7     4     1  1.75
       8     4     2  2.00
       9     4     2  2.25

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int I, J, the arguments.

    Output, double R8_FRACTION, the value of the ratio.
*/
{
  double value;

  value = ( double ) ( i ) / ( double ) ( j );

  return value;
}
/*****************************************************************************/

double r8_fractional ( double x )

/*****************************************************************************/
/*
  Purpose:

    R8_FRACTIONAL returns the fractional part of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument.

    Output, double R8_FRACTIONAL, the fractional part of X.
*/
{
  double value;

  value = r8_abs ( x ) - ( double ) ( ( int ) r8_abs ( x ) );

  return value;
}
/******************************************************************************/

double r8_gamma ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_GAMMA evaluates Gamma(X) for a real argument.

  Discussion:

    The C math library includes the GAMMA ( X ) function which should generally
    be used instead of this function.

    This routine calculates the gamma function for a real argument X.

    Computation is based on an algorithm outlined in reference 1.
    The program uses rational functions that approximate the gamma
    function to at least 20 significant decimal digits.  Coefficients
    for the approximation over the interval (1,2) are unpublished.
    Those for the approximation for 12 <= X are from reference 2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by William Cody, Laura Stoltz.
    C version by John Burkardt.

  Reference:

    William Cody,
    An Overview of Software Development for Special Functions,
    in Numerical Analysis Dundee, 1975,
    edited by GA Watson,
    Lecture Notes in Mathematics 506,
    Springer, 1976.

    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    Charles Mesztenyi, John Rice, Henry Thatcher,
    Christoph Witzgall,
    Computer Approximations,
    Wiley, 1968,
    LC: QA297.C64.

  Parameters:

    Input, double X, the argument of the function.

    Output, double R8_GAMMA, the value of the function.
*/
{
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  int i;
  int n;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  int parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = 0;
  fact = 1.0;
  n = 0;
  y = x;
/*
  Argument is negative.
*/
  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = 1;
      }

      fact = - pi / sin ( pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Argument is positive.
*/
  if ( y < eps )
  {
/*
  Argument < EPS.
*/
    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
/*
  0.0 < argument < 1.0.
*/
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
/*
  1.0 < argument < 12.0.
  Reduce argument if necessary.
*/
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }
/*
  Evaluate approximation for 1.0 < argument < 2.0.
*/
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;
/*
  Adjust result for case  0.0 < argument < 1.0.
*/
    if ( y1 < y )
    {
      res = res / y1;
    }
/*
  Adjust result for case 2.0 < argument < 12.0.
*/
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {
/*
  Evaluate for 12.0 <= argument.
*/
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Final adjustments and return.
*/
  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
/******************************************************************************/

double r8_huge ( void )

/******************************************************************************/
/*
  Purpose:

    R8_HUGE returns a "huge" R8.

  Discussion:

    The value returned by this function is NOT required to be the
    maximum representable R8.  This value varies from machine to machine,
    from compiler to compiler, and may cause problems when being printed.
    We simply want a "very large" but non-infinite number.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2007

  Author:

    John Burkardt

  Parameters:

    Output, double R8_HUGE, a "huge" R8 value.
*/
{
  double value;

  value = 1.0E+30;

  return value;
}
/******************************************************************************/

double r8_hypot ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_HYPOT returns the value of sqrt ( X^2 + Y^2 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the arguments.

    Output, double R8_HYPOT, the value of sqrt ( X^2 + Y^2 ).
*/
{
  double a;
  double b;
  double value;

  if ( r8_abs ( x ) < r8_abs ( y ) )
  {
    a = r8_abs ( y );
    b = r8_abs ( x );
  }
  else
  {
    a = r8_abs ( x );
    b = r8_abs ( y );
  }
/*
  A contains the larger value.
*/
  if ( a == 0.0 )
  {
    value = 0.0;
  }
  else
  {
    value = a * sqrt ( 1.0 + ( b / a ) * ( b / a ) );
  }

  return value;
}
/******************************************************************************/

int r8_in_01 ( double a )

/******************************************************************************/
/*
  Purpose:

    R8_IN_01 is 1 if an R8 is in the range [0,1].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double A, the value.

    Output, int R8_IN_01, is 1 if A is between 0 and 1.
*/
{
  if ( a < 0.0 || 1.0 < a )
  {
    return 0;
  }

  return 1;
}
/******************************************************************************/

int r8_insignificant ( double r, double s )

/******************************************************************************/
/*
  Purpose:

    R8_INSIGNIFICANT determines if an R8 is insignificant.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, double R, the number to be compared against.

    Input, double S, the number to be compared.

    Output, int R8_INSIGNIFICANT, is TRUE if S is insignificant
    compared to R.
*/
{
  double t;
  double tol;
  int value;

  value = 1;

  t = r + s;
  tol = r8_epsilon ( ) * r8_abs ( r );

  if ( tol < r8_abs ( r - t ) )
  {
    value = 0;
  }
  
  return value;
}
/******************************************************************************/

int r8_is_int ( double r )

/******************************************************************************/
/*
  Purpose:

    R8_IS_INT is 1 if an R8 represents an integer value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double R, the number to be checked.

    Output, int R8_IS_INT, is 1 if R is an integer value.
*/
{
  int i4_huge = 2147483647;
  int value;

  if ( ( double ) ( i4_huge ) < r )
  {
    value = 0;
  }
  else if ( r < - ( double ) ( i4_huge ) )
  {
    value = 0;
  }
  else if ( r == ( double ) ( ( int ) ( r ) ) )
  {
    value = 1;
  }
  else
  {
    value = 0;
  }
  return value;
}
/******************************************************************************/

double r8_log_10 ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_LOG_10 returns the logarithm base 10 of an R8.

  Discussion:

    R8_LOG_10 ( X ) = Log10 ( |X| )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose base 2 logarithm is desired.
    X should not be 0.

    Output, double R8_LOG_10, the logarithm base 10 of the absolute
    value of X.  It should be true that |X| = 10**R_LOG_10.
*/
{
  double value;

  if ( x == 0.0 )
  {
    value = - r8_huge ( );
  }
  else
  {
    value = log10 ( fabs ( x ) );
  }

  return value;
}
/******************************************************************************/

double r8_log_2 ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_LOG_2 returns the logarithm base 2 of an R8.

  Discussion:

    R8_LOG_2 ( X ) = Log ( |X| ) / Log ( 2.0 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose base 2 logarithm is desired.
    X should not be 0.

    Output, double R8_LOG_2, the logarithm base 2 of the absolute
    value of X.  It should be true that |X| = 2**R_LOG_2.
*/
{
  double value;

  if ( x == 0.0 )
  {
    value = - r8_huge ( );
  }
  else
  {
    value = log ( fabs ( x ) ) / log ( 2.0 );
  }

  return value;
}
/******************************************************************************/

double r8_log_b ( double x, double b )

/******************************************************************************/
/*
  Purpose:

    R8_LOG_B returns the logarithm base B of an R8.

  Discussion:

    R8_LOG_B ( X ) = log ( |X| ) / log ( |B| )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 May 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose base B logarithm is desired.
    X should not be 0.

    Input, double B, the base, which should not be 0, 1 or -1.

    Output, double R8_LOG_B, the logarithm base B of the absolute
    value of X.  It should be true that |X| = |B|**R_LOG_B.
*/
{
  double value;

  if ( b == 0.0 || b == 1.0 || b == -1.0 )
  {
    value = - r8_huge ( );
  }
  else if ( fabs ( x ) == 0.0 )
  {
    value = - r8_huge ( );
  }
  else
  {
    value = log ( fabs ( x ) ) / log ( fabs ( b ) );
  }

  return value;
}
/******************************************************************************/

void r8_mant ( double x, int *s, double *r, int *l )

/******************************************************************************/
/*
  Purpose:

    R8_MANT computes the "mantissa" or "fraction part" of an R8.

  Formula:

    X = S * R * 2^L

    S is +1 or -1,
    R is a real between 1.0 and 2.0,
    L is an integer.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 April 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, the real number to be decomposed.

    Output, int *S, the "sign" of the number.
    S will be -1 if X is less than 0, and +1 if X is greater
    than or equal to zero.

    Output, double *R, the mantissa of X.  R will be greater
    than or equal to 1, and strictly less than 2.  The one
    exception occurs if X is zero, in which case R will also
    be zero.

    Output, int *L, the integer part of the logarithm (base 2) of X.
*/
{
/*
  Determine the sign.
*/
  if ( x < 0.0 )
  {
    *s = -1;
  }
  else
  {
    *s = 1;
  }
/*
  Set R to the absolute value of X, and L to zero.
  Then force R to lie between 1 and 2.
*/
  if ( x < 0.0 )
  {
    *r = -x;
  }
  else
  {
    *r = x;
  }

  *l = 0;
/*
  Time to bail out if X is zero.
*/
  if ( x == 0.0 )
  {
    return;
  }

  while ( 2.0 <= *r )
  {
    *r = *r / 2.0;
    *l = *l + 1;
  }

  while ( *r < 1.0 )
  {
    *r = *r * 2.0;
    *l = *l - 1;
  }

  return;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MAX returns the maximum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

double r8_min ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MIN returns the minimum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to compare.

    Output, double R8_MIN, the minimum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
/******************************************************************************/

double r8_mod ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MOD returns the remainder of R8 division.

  Discussion:

    If
      REM = R8_MOD ( X, Y )
      RMULT = ( X - REM ) / Y
    then
      X = Y * RMULT + REM
    where REM has the same sign as X, and abs ( REM ) < Y.

  Example:

        X         Y     R8_MOD   R8_MOD  Factorization

      107        50       7     107 =  2 *  50 + 7
      107       -50       7     107 = -2 * -50 + 7
     -107        50      -7    -107 = -2 *  50 - 7
     -107       -50      -7    -107 =  2 * -50 - 7

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 June 2007

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number to be divided.

    Input, double Y, the number that divides X.

    Output, double R8_MOD, the remainder when X is divided by Y.
*/
{
  double value;

  if ( y == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_MOD - Fatal error!\n" );
    fprintf ( stderr, "  R8_MOD ( X, Y ) called with Y = %g\n", y );
    exit ( 1 );
  }

  value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

  if ( x < 0.0 && 0.0 < value )
  {
    value = value - r8_abs ( y );
  }
  else if ( 0.0 < x && value < 0.0 )
  {
    value = value + r8_abs ( y );
  }

  return value;
}
/******************************************************************************/

double r8_modp ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_MODP returns the nonnegative remainder of R8 division.

  Formula:

    If
      REM = R8_MODP ( X, Y )
      RMULT = ( X - REM ) / Y
    then
      X = Y * RMULT + REM
    where REM is always nonnegative.

  Discussion:

    The MOD function computes a result with the same sign as the
    quantity being divided.  Thus, suppose you had an angle A,
    and you wanted to ensure that it was between 0 and 360.
    Then mod(A,360.0) would do, if A was positive, but if A
    was negative, your result would be between -360 and 0.

    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.

  Example:

        I         J     MOD  R8_MODP   R8_MODP Factorization

      107        50       7       7    107 =  2 *  50 + 7
      107       -50       7       7    107 = -2 * -50 + 7
     -107        50      -7      43   -107 = -3 *  50 + 43
     -107       -50      -7      43   -107 =  3 * -50 + 43

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number to be divided.

    Input, double Y, the number that divides X.

    Output, double R8_MODP, the nonnegative remainder when X is divided by Y.
*/
{
  double value;

  if ( y == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_MODP - Fatal error!\n" );
    fprintf ( stderr, "  R8_MODP ( X, Y ) called with Y = %g\n", y );
    exit ( 1 );
  }

  value = x - ( ( double ) ( ( int ) ( x / y ) ) ) * y;

  if ( value < 0.0 )
  {
    value = value + fabs ( y );
  }

  return value;
}
/******************************************************************************/

double r8_mop ( int i )

/******************************************************************************/
/*
  Purpose:

    R8_MOP returns the I-th power of -1 as an R8 value.

  Discussion:

    An R8 is an double value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int I, the power of -1.

    Output, double R8_MOP, the I-th power of -1.
*/
{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = + 1.0;
  }
  else
  {
    value = - 1.0;
  }

  return value;
}
/******************************************************************************/

int r8_nint ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_NINT returns the nearest integer to an R8.

  Example:

        X         R8_NINT

      1.3         1
      1.4         1
      1.5         1 or 2
      1.6         2
      0.0         0
     -0.7        -1
     -1.1        -1
     -1.6        -2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the value.

    Output, int R8_NINT, the nearest integer to X.
*/
{
  int s;
  int value;

  if ( x < 0.0 )
  {
    s = - 1;
  }
  else
  {
    s = + 1;
  }
  value = s * ( int ) ( fabs ( x ) + 0.5 );

  return value;
}
/******************************************************************************/

double r8_normal ( double a, double b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_NORMAL returns a scaled pseudonormal R8.

  Discussion:

    The normal probability distribution function (PDF) is sampled,
    with mean A and standard deviation B.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 November 2004

  Author:

    John Burkardt

  Parameters:

    Input, double A, the mean of the PDF.

    Input, double B, the standard deviation of the PDF.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8_NORMAL, a sample of the normal PDF.
*/
{
  double value;

  value = a + b * r8_normal_01 ( seed );

  return value;
}
/******************************************************************************/

double r8_normal_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_NORMAL_01 samples the standard normal probability distribution.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

    The Box-Muller method is used, which is efficient, but
    generates two values at a time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8_NORMAL_01, a normally distributed random value.
*/
{
  double pi = 3.141592653589793;
  double r1;
  double r2;
  static int used = -1;
  double x;
  static double y = 0.0;

  if ( used == -1 )
  {
    used = 0;
  }
/*
  If we've used an even number of values so far, generate two more, return one,
  and save one.
*/
  if ( ( used % 2 )== 0 )
  {
    for ( ; ; )
    {
      r1 = r8_uniform_01 ( seed );
      if ( r1 != 0.0 )
      {
        break;
      }
    }

    r2 = r8_uniform_01 ( seed );

    x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );
    y = sqrt ( -2.0 * log ( r1 ) ) * sin ( 2.0 * pi * r2 );
  }
  else
  {

    x = y;

  }

  used = used + 1;

  return x;
}
/******************************************************************************/

double r8_pi ( void )

/******************************************************************************/
/*
  Purpose:

    R8_PI returns the value of PI as an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Output, double R8_PI, the value of PI.
*/
{
  double value = 3.141592653589793;

  return value;
}
/******************************************************************************/

double r8_pi_sqrt ( void )

/******************************************************************************/
/*
  Purpose:

    R8_PI_SQRT returns the square root of PI as an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_PI_SQRT, the square root of PI.
*/
{
  double value = 1.7724538509055160273;

  return value;
}
/******************************************************************************/

double r8_power ( double r, int p )

/******************************************************************************/
/*
  Purpose:

    R8_POWER computes an integer power of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, double R, the base.

    Input, int P, the power, which may be negative.

    Output, double R8_POWER, the value of R^P.
*/
{
  double value;
/*
  Special case.  R^0 = 1.
*/
  if ( p == 0 )
  {
    value = 1.0;
  }
/*
  Special case.  Positive powers of 0 are 0.
  We go ahead and compute negative powers, relying on the software to complain.
*/
  else if ( r == 0.0 )
  {
    if ( 0 < p )
    {
      value = 0.0;
    }
    else
    {
      value = pow ( r, p );
    }
  }
  else if ( 1 <= p )
  {
    value = pow ( r, p );
  }
  else
  {
    value = pow ( r, p );
  }

  return value;
}
/******************************************************************************/

double r8_power_fast ( double r, int p, int *mults )

/******************************************************************************/
/*
  Purpose:

    R8_POWER_FAST computes the P-th power of R, for real R and integer P.

  Discussion:

    Obviously, R^P can be computed using P-1 multiplications.

    However, R^P can also be computed using at most 2*LOG2(P) multiplications.
    To do the calculation this way, let N = LOG2(P).
    Compute A, A^2, A^4, ..., A^N by N-1 successive squarings.

    Start the value of R^P at A, and each time that there is a 1 in
    the binary expansion of P, multiply by the current result of the squarings.

    This algorithm is not optimal.  For small exponents, and for special
    cases, the result can be computed even more quickly.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, double R, the base.

    Input, int P, the power, which may be negative.

    Output, int *MULTS, the number of multiplications and divisions.

    Output, double R8_POWER_FAST, the value of R^P.
*/
{
  int p_mag;
  int p_sign;
  double r2;
  double value;

  *mults = 0;
/*
  Special bases.
*/
  if ( r == 1.0 )
  {
    value = 1.0;
    return value;
  }

  if ( r == -1.0 )
  {
    if ( ( p % 2 ) == 1 )
    {
      value = -1.0;
    }
    else
    {
      value = 1.0;
    }
    return value;
  }

  if ( r == 0.0 )
  {
    if ( p <= 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8_POWER_FAST - Fatal error!\n" );
      fprintf ( stderr, "  Base is zero, and exponent is negative.\n" );
      exit ( 1 );
    }

    value = 0.0;
    return value;
  }
/*
  Special powers.
*/
  if ( p == -1 )
  {
    value = 1.0 / r;
    *mults = *mults + 1;
    return value;
  }
  else if ( p == 0 )
  {
    value = 1.0;
    return value;
  }
  else if ( p == 1 )
  {
    value = r;
    return value;
  }
/*
  Some work to do.
*/
  p_mag = abs ( p );
  p_sign = i4_sign ( p );

  value = 1.0;
  r2 = r;

  while ( 0 < p_mag )
  {
    if ( ( p_mag % 2 ) == 1 )
    {
      value = value * r2;
      *mults = *mults + 1;
    }

    p_mag = p_mag / 2;
    r2 = r2 * r2;
    *mults = *mults + 1;
  }

  if ( p_sign == -1 )
  {
    value = 1.0 / value;
    *mults = *mults + 1;
  }

  return value;
}
/******************************************************************************/

double r8_pythag ( double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8_PYTHAG computes sqrt ( A^2 + B^2 ), avoiding overflow and underflow.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 January 2002

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the values for which sqrt ( A^2 + B^2 ) is desired.

    Output, double R8_PYTHAG, the value of sqrt ( A^2 + B^2 ).
*/
{
  double a_abs;
  double b_abs;
  double value;

  a_abs = r8_abs ( a );
  b_abs = r8_abs ( b );

  if ( b_abs < a_abs )
  {
    value = a_abs * sqrt ( 1.0 + pow ( b_abs / a_abs, 2 ) );
  }
  else if ( b_abs == 0.0 )
  {
    value = 0.0;
  }
  else if ( a_abs <= b_abs )
  {
    value = b_abs * sqrt ( 1.0 + pow ( a_abs / b_abs, 2 ) );
  }

  return value;
}
/******************************************************************************/

double r8_reverse_bytes ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_REVERSE_BYTES reverses the bytes in an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, a value whose bytes are to be reversed.

    Output, R8_REVERSE_BYTES, a value with bytes in reverse order;
*/
{
  char c;
  union
  {
    double ydouble;
    char ychar[8];
  } y;

  y.ydouble = x;

  c = y.ychar[0];
  y.ychar[0] = y.ychar[7];
  y.ychar[7] = c;

  c = y.ychar[1];
  y.ychar[1] = y.ychar[6];
  y.ychar[6] = c;

  c = y.ychar[2];
  y.ychar[2] = y.ychar[5];
  y.ychar[5] = c;

  c = y.ychar[3];
  y.ychar[3] = y.ychar[4];
  y.ychar[4] = c;

  return ( y.ydouble );
}
/******************************************************************************/

double r8_round ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_ROUND rounds an R8 to the nearest integral value.

  Example:

        X         Value

      1.3         1.0
      1.4         1.0
      1.5         1.0 or 2.0
      1.6         2.0
      0.0         0.0
     -0.7        -1.0
     -1.1        -1.0
     -1.6        -2.0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, the value.

    Output, double R8_ROUND, the rounded value.
*/
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( double ) floor ( - x + 0.5 );
  }
  else
  {
    value =   ( double ) floor (   x + 0.5 );
  }

  return value;
}
/******************************************************************************/

double r8_round2 ( int nplace, double x )

/******************************************************************************/
/*
  Purpose:

    R8_ROUND2 rounds an R8 in base 2.

  Discussion:

    Assume that the input quantity X has the form

      X = S * J * 2^L

    where S is plus or minus 1, L is an integer, and J is a binary
    mantissa which is either exactly zero, or greater than or equal
    to 0.5 and less than 1.0.

    Then on return, XROUND = R8_ROUND2 ( NPLACE, X ) will satisfy

      XROUND = S * K * 2^L

    where S and L are unchanged, and K is a binary mantissa which
    agrees with J in the first NPLACE binary digits and is zero
    thereafter.

    If NPLACE is 0, XROUND will always be zero.

    If NPLACE is 1, the mantissa of XROUND will be 0 or 0.5.

    If NPLACE is 2, the mantissa of XROUND will be 0, 0.25, 0.50,
    or 0.75.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int NPLACE, the number of binary digits to
    preserve.  NPLACE should be 0 or positive.

    Input, double X, the real number to be decomposed.

    Output, double R8_ROUND2, the rounded value of X.
*/
{
  int iplace;
  int l;
  int s;
  static double two = 2.0;
  double xmant;
  double xtemp;
  double value;

  value = 0.0;
/*
  1: Handle the special case of 0.
*/
  if ( x == 0.0 )
  {
    return value;
  }

  if ( nplace <= 0 )
  {
    return value;
  }
/*
  2: Determine the sign S.
*/
  if ( 0.0 < x )
  {
    s = 1;
    xtemp = x;
  }
  else
  {
    s = -1;
    xtemp = -x;
  }
/*
  3: Force XTEMP to lie between 1 and 2, and compute the
  logarithm L.
*/
  l = 0;

  while ( 2.0 <= xtemp )
  {
    xtemp = xtemp / 2.0;
    l = l + 1;
  }

  while ( xtemp < 1.0 )
  {
    xtemp = xtemp * 2.0;
    l = l - 1;
  }
/*
  4: Strip out the digits of the mantissa as XMANT, and decrease L.
*/
  xmant = 0.0;
  iplace = 0;

  for ( ; ; )
  {
    xmant = 2.0 * xmant;

    if ( 1.0 <= xtemp )
    {
      xmant = xmant + 1.0;
      xtemp = xtemp - 1.0;
    }

    iplace = iplace + 1;

    if ( xtemp == 0.0 || nplace <= iplace )
    {
      value = s * xmant * pow ( two, l );
      break;
    }

    l = l - 1;
    xtemp = xtemp * 2.0;
  }

  return value;
}
/******************************************************************************/

double r8_roundb ( int base, int nplace, double x )

/******************************************************************************/
/*
  Purpose:

    R8_ROUNDB rounds an R8 in a given base.

  Discussion:

    The code does not seem to do a good job of rounding when
    the base is negative

    Assume that the input quantity X has the form

      X = S * J * BASE^L

    where S is plus or minus 1, L is an integer, and J is a
    mantissa base BASE which is either exactly zero, or greater
    than or equal to (1/BASE) and less than 1.0.

    Then on return, XROUND will satisfy

      XROUND = S * K * BASE^L

    where S and L are unchanged, and K is a mantissa base BASE
    which agrees with J in the first NPLACE digits and is zero
    thereafter.

    Note that because of rounding, for most bases, most numbers
    with a fractional quantities cannot be stored exactly in the
    computer, and hence will have trailing "bogus" digits.

    If NPLACE is 0, XROUND will always be zero.

    If NPLACE is 1, the mantissa of XROUND will be 0,
    1/BASE, 2/BASE, ..., (BASE-1)/BASE.

    If NPLACE is 2, the mantissa of XROUND will be 0,
    BASE/BASE^2, (BASE+1)/BASE^2, ...,
    BASE^2-2/BASE^2, BASE**2-1/BASE^2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int BASE, the base of the arithmetic.
    BASE must not be zero.  Theoretically, BASE may be negative.

    Input, int NPLACE, the number of digits base BASE to
    preserve.  NPLACE should be 0 or positive.

    Input, double X, the number to be decomposed.

    Output, double R8_ROUNDB, the rounded value of X.
*/
{
  int iplace;
  int is;
  int js;
  int l;
  double r8_base;
  double value;
  double xmant;
  double xtemp;

  value = 0.0;
  r8_base = ( double ) base;
/*
  0: Error checks.
*/
  if ( base == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_ROUNDB - Fatal error!\n" );
    fprintf ( stderr, "  The base BASE cannot be zero.\n" );
    exit ( 1 );
  }
/*
  1: Handle the special case of 0.
*/
  if ( x == 0.0 )
  {
    return value;
  }

  if ( nplace <= 0 )
  {
    return value;
  }
/*
  2: Determine the sign IS.
*/
  if ( 0.0 < x )
  {
    is = 1;
    xtemp = x;
  }
  else
  {
    is = -1;
    xtemp = -x;
  }
/*
  3: Force XTEMP to lie between 1 and ABS(BASE), and compute the
  logarithm L.
*/
  l = 0;

  while ( r8_abs ( r8_base ) <= r8_abs ( xtemp ) )
  {
    xtemp = xtemp / r8_base;

    if ( xtemp < 0.0 )
    {
      is = -is;
      xtemp = -xtemp;
    }
    l = l + 1;
  }

  while ( r8_abs ( xtemp ) < 1.0 )
  {
    xtemp = xtemp * r8_base;

    if ( xtemp < 0.0 )
    {
      is = -is;
      xtemp = -xtemp;
    }

    l = l - 1;
  }
/*
  4: Now strip out the digits of the mantissa as XMANT, and
  decrease L.
*/
  xmant = 0.0;
  iplace = 0;
  js = is;

  for ( ; ; )
  {
    xmant = r8_base * xmant;

    if ( xmant < 0.0 )
    {
      js = -js;
      xmant = -xmant;
    }

    if ( 1.0 <= xtemp )
    {
      xmant = xmant + ( int ) ( xtemp );
      xtemp = xtemp - ( int ) ( xtemp );
    }

    iplace = iplace + 1;

    if ( xtemp == 0.0 || nplace <= iplace )
    {
      value = ( double ) js * xmant * pow ( r8_base, l );
      break;
    }

    l = l - 1;
    xtemp = xtemp * r8_base;

    if ( xtemp < 0.0 )
    {
      is = -is;
      xtemp = -xtemp;
    }
  }

  return value;
}
/******************************************************************************/

double r8_roundx ( int nplace, double x )

/******************************************************************************/
/*
  Purpose:

    R8_ROUNDX rounds an R8 in base 10.

  Discussion:

    Assume that the input quantity X has the form

      X = S * J * 10^L

    where S is plus or minus 1, L is an integer, and J is a decimal
    mantissa which is either exactly zero, or greater than or equal
    to 0.1 and less than 1.0.

    Then on return, XROUND will satisfy

      XROUND = S * K * 10^L

    where S and L are unchanged, and K is a decimal mantissa which
    agrees with J in the first NPLACE decimal digits and is zero
    thereafter.

    Note that because of rounding, most decimal fraction quantities
    cannot be stored exactly in the computer, and hence will have
    trailing "bogus" digits.

    If NPLACE is 0, XROUND will always be zero.

    If NPLACE is 1, the mantissa of XROUND will be 0, 0.1,
    0.2, ..., or 0.9.

    If NPLACE is 2, the mantissa of XROUND will be 0, 0.01, 0.02,
    0.03, ..., 0.98, 0.99.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int NPLACE, the number of decimal digits to
    preserve.  NPLACE should be 0 or positive.

    Input, double X, the number to be decomposed.

    Output, double R8_ROUNDX, the rounded value of X.
*/
{
  int iplace;
  int is;
  int l;
  static double ten = 10.0;
  double xmant;
  double xround;
  double xtemp;

  xround = 0.0;
/*
  1: Handle the special case of 0.
*/
  if ( x == 0.0 )
  {
    return xround;
  }

  if ( nplace <= 0 )
  {
    return xround;
  }
/*
  2: Determine the sign IS.
*/
  if ( 0.0 < x )
  {
    is = 1;
    xtemp = x;
  }
  else
  {
    is = -1;
    xtemp = -x;
  }
/*
  3: Force XTEMP to lie between 1 and 10, and compute the
  logarithm L.
*/
  l = 0;

  while ( 10.0 <= x )
  {
    xtemp = xtemp / 10.0;
    l = l + 1;
  }

  while ( xtemp < 1.0 )
  {
    xtemp = xtemp * 10.0;
    l = l - 1;
  }
/*
  4: Now strip out the digits of the mantissa as XMANT, and
  decrease L.
*/
  xmant = 0.0;
  iplace = 0;

  for ( ; ; )
  {
    xmant = 10.0 * xmant;

    if ( 1.0 <= xtemp )
    {
      xmant = xmant + ( int ) xtemp;
      xtemp = xtemp - ( int ) xtemp;
    }

    iplace = iplace + 1;

    if ( xtemp == 0.0 || nplace <= iplace )
    {
      xround = is * xmant * pow ( ten, l );
      break;
    }

    l = l - 1;
    xtemp = xtemp * 10.0;
  }

  return xround;
}
/******************************************************************************/

double r8_sech ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SECH evaluates the hyperbolic secant, while avoiding COSH overflow.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    21 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the function.

    Output, double R8_SECH, the value of the function.
*/
{
  double log_huge = 80.0;
  double value;

  if ( log_huge < r8_abs ( x ) )
  {
    value = 0.0;
  }
  else
  {
    value = 1.0 / cosh ( x );
  }
  return value;
}
/******************************************************************************/

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, double R8_SIGN, the sign of X.
*/
{
  double value;

  if ( x < 0.0 )
  {
    value = - 1.0;
  }
  else
  {
    value = + 1.0;
  }
  return value;
}
/******************************************************************************/

char r8_sign_char ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN_CHAR returns a character indicating the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 April 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, char R8_SIGN_CHAR, the sign of X, '-', '0' or '+'.
*/
{
  char value;

  if ( x < 0.0 )
  {
    value = '-';
  }
  else if ( x == 0.0 )
  {
    value = '0';
  }
  else
  {
    value = '+';
  }
  return value;
}
/******************************************************************************/

int r8_sign_match ( double r1, double r2 )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN_MATCH is TRUE if two R8's are of the same sign.

  Discussion:

    This test could be coded numerically as

      if ( 0 <= r1 * r2 ) then ...

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 April 2012

  Author:

    John Burkardt

  Parameters:

    Input, double R1, R2, the values to check.

    Output, int R8_SIGN_MATCH, is TRUE if ( R1 <= 0 and R2 <= 0 )
    or ( 0 <= R1 and 0 <= R2 ).
*/
{
  int value;

  value = ( r1 <= 0.0 && r2 <= 0.0 ) ||
          ( 0.0 <= r1 && 0.0 <= r2 );

  return value;
}
/******************************************************************************/

int r8_sign_match_strict ( double r1, double r2 )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN_MATCH_STRICT is TRUE if two R8's are of the same strict sign.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 April 2012

  Author:

    John Burkardt

  Parameters:

    Input, double R1, R2, the values to check.

    Output, int R8_SIGN_MATCH_STRICT, is TRUE if the signs match.
*/
{
  int value;

  value = ( r1 <  0.0 && r2 <  0.0 ) ||
          ( r1 == 0.0 && r2 == 0.0 ) ||
          ( 0.0 < r1 && 0.0 < r2 );

  return value;
}
/******************************************************************************/

int r8_sign_opposite ( double r1, double r2 )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN_OPPOSITE is TRUE if two R8's are not of the same sign.

  Discussion:

    This test could be coded numerically as

      if ( r1 * r2 <= 0.0 ) ...

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, double R1, R2, the values to check.

    Output, int R8_SIGN_OPPOSITE, is TRUE if ( R1 <= 0 and 0 <= R2 )
    or ( R2 <= 0 and 0 <= R1 ).
*/
{
  int value;

  value = ( r1 <= 0.0 && 0.0 <= r2 ) || ( r2 <= 0.0 && 0.0 <= r1 );

  return value;
}
/******************************************************************************/

int r8_sign_opposite_strict ( double r1, double r2 )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN_OPPOSITE_STRICT is TRUE if two R8's are strictly of opposite sign.

  Discussion:

    This test could be coded numerically as

      if ( r1 * r2 < 0.0 ) ...

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, double R1, R2, the values to check.

    Output, int R8_SIGN_OPPOSITE_STRICT, is TRUE if ( R1 < 0 and 0 < R2 )
    or ( R2 < 0 and 0 < R1 ).
*/
{
  int value;

  value = ( r1 < 0.0 && 0.0 < r2 ) || ( r2 < 0.0 && 0.0 < r1 );

  return value;
}
/******************************************************************************/

double r8_sign2 ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN2 returns the first argument with the sign of the second.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the input arguments.

    Output, double R8_SIGN2, is equal in absolute value to the absolute value 
    of X, and has the sign of Y.
*/
{
  double value;

  if ( 0.0 <= y )
  {
    value = r8_abs ( x );
  } 
  else
  {
    value = - r8_abs ( x );
  }
  return value;
}
/******************************************************************************/

double r8_sqrt_i4 ( int i )

/******************************************************************************/
/*
  Purpose:

    R8_SQRT_I4 returns the square root of an I4 as an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    04 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int I, the number whose square root is desired.

    Output, double R8_SQRT_I4, the value of sqrt(I).
*/
{
  double value;

  value = sqrt ( ( double ) ( i ) );

  return value;
}
/******************************************************************************/

double r8_sum ( double x, double y )

/******************************************************************************/
/*
  Purpose:

    R8_SUM returns the sum of two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, the quantities to add.

    Output, double R8_SUM, the sum of X and Y.
*/
{
  double value;

  value = x + y;

  return value;
}
/******************************************************************************/

void r8_swap ( double *x, double *y )

/******************************************************************************/
/*
  Purpose:

    R8_SWAP switches two R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input/output, double *X, *Y.  On output, the values of X and
    Y have been interchanged.
*/
{
  double z;

  z = *x;
  *x = *y;
  *y = z;

  return;
}
/******************************************************************************/

void r8_swap3 ( double *x, double *y, double *z )

/******************************************************************************/
/*
  Purpose:

    R8_SWAP3 swaps three R8's.

  Example:

    Input:

      X = 1, Y = 2, Z = 3

    Output:

      X = 2, Y = 3, Z = 1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input/output, double *X, *Y, *Z, three values to be swapped.
*/
{
  double w;

   w = *x;
  *x = *y;
  *y = *z;
  *z =  w;

  return;
}
/******************************************************************************/

double r8_tiny ( void )

/******************************************************************************/
/*
  Purpose:

    R8_TINY returns a "tiny" R8.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Output, double R8_TINY, a "tiny" R8 value.
*/
{
  float value;

  value = 0.4450147717014E-307;

  return value;
}
/******************************************************************************/

void r8_to_dhms ( double r, int *d, int *h, int *m, int *s )

/******************************************************************************/
/*
  Purpose:

    R8_TO_DHMS converts an R8 day value into days, hours, minutes, seconds.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, double R, a real number representing a time period measured in days.

    Output, int D, H, M, S, the equivalent number of days, hours,
    minutes and seconds.
*/
{
  int sign;

  if ( 0.0 <= r )
  {
    sign = 1;
  }
  else if ( r < 0.0 )
  {
    sign = -1;
    r = -r;
  }

  *d = ( int ) r;

  r = r - ( double ) *d;
  r = 24.0 * r;
  *h = ( int ) r;

  r = r - ( double ) *h;
  r = 60.0 * r;
  *m = ( int ) r;

  r = r - ( double ) *m;
  r = 60.0 * r;
  *s = ( int ) r;

  if ( sign == -1 )
  {
    *d = -(*d);
    *h = -(*h);
    *m = -(*m);
    *s = -(*s);
  }

  return;
}
/******************************************************************************/

int r8_to_i4 ( double x, double xmin, double xmax, int ixmin, int ixmax )

/******************************************************************************/
/*
  Purpose:

    R8_TO_I4 maps real X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].

  Discussion:

    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
    IX := min ( IX, max ( IXMIN, IXMAX ) )
    IX := max ( IX, min ( IXMIN, IXMAX ) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, the real number to be converted.

    Input, double XMIN, XMAX, the real range.  XMAX and XMIN must not be
    equal.  It is not necessary that XMIN be less than XMAX.

    Input, int IXMIN, IXMAX, the allowed range of the output
    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
    It is not necessary that IXMIN be less than IXMAX.

    Output, int R8_TO_I4, the value in the range [IXMIN,IXMAX] that
    corresponds to X.
*/
{
  int ix;
  double temp;

  if ( xmax == xmin )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_TO_I4 - Fatal error!\n" );
    fprintf ( stderr, "  XMAX = XMIN, making a zero divisor.\n" );
    fprintf ( stderr, "  XMAX = %g\n", xmax );
    fprintf ( stderr, "  XMIN = %g\n", xmin );
    exit ( 1 );
  }

  temp =
      ( ( xmax - x        ) * ( double ) ixmin
      + (        x - xmin ) * ( double ) ixmax )
      / ( xmax     - xmin );

  if ( 0.0 <= temp )
  {
    temp = temp + 0.5;
  }
  else
  {
    temp = temp - 0.5;
  }

  ix = ( int ) temp;

  return ix;
}
/******************************************************************************/

double r8_to_r8_discrete ( double r, double rmin, double rmax, int nr )

/******************************************************************************/
/*
  Purpose:

    R8_TO_R8_DISCRETE maps R to RD in [RMIN, RMAX] with NR possible values.

  Discussion:

    if ( R < RMIN ) then
      RD = RMIN
    else if ( RMAX < R ) then
      RD = RMAX
    else
      T = nint ( ( NR - 1 ) * ( R - RMIN ) / ( RMAX - RMIN ) )
      RD = RMIN + T * ( RMAX - RMIN ) / real ( NR - 1 )

    In the special case where NR = 1, when

      XD = 0.5 * ( RMAX + RMIN )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, double R, the number to be converted.

    Input, double RMAX, RMIN, the maximum and minimum
    values for RD.

    Input, int NR, the number of allowed values for XD.
    NR should be at least 1.

    Output, double RD, the corresponding discrete value.
*/
{
  int f;
  double rd;
/*
  Check for errors.
*/
  if ( nr < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_TO_R8_DISCRETE - Fatal error!\n" );
    fprintf ( stderr, "  NR = %d\n", nr );
    fprintf ( stderr, "  but NR must be at least 1.\n" );
    exit ( 1 );
  }

  if ( nr == 1 )
  {
    rd = 0.5 * ( rmin + rmax );
    return rd;
  }

  if ( rmax == rmin )
  {
    rd = rmax;
    return rd;
  }

  f = r8_nint ( ( double ) ( nr ) * ( rmax - r ) / ( rmax - rmin ) );
  f = i4_max ( f, 0 );
  f = i4_min ( f, nr );

  rd = ( ( double ) (      f ) * rmin
       + ( double ) ( nr - f ) * rmax )
       / ( double ) ( nr     );

  return rd;
}
/******************************************************************************/

double r8_uniform_ab ( double a, double b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM returns a pseudorandom R8 scaled to [A,B].

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 November 2004

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the limits of the interval.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double R8_UNIFORM_AB, a number strictly between A and B.
*/
{
  double value;

  value = a + ( b - a ) * r8_uniform_01 ( seed );

  return value;
}
/******************************************************************************/

double r8_uniform_01 ( int *seed )

/******************************************************************************/
/*
  Purpose:

    R8_UNIFORM_01 returns a pseudorandom R8 scaled to [0,1].

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      r8_uniform_01 = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    If the initial seed is 12345, then the first three computations are

      Input     Output      R8_UNIFORM_01
      SEED      SEED

         12345   207482415  0.096616
     207482415  1790989824  0.833995
    1790989824  2035175616  0.947702

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation
    edited by Jerry Banks,
    Wiley Interscience, page 95, 1998.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    P A Lewis, A S Goodman, J M Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
    0 and 1.
*/
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( ( double ) ( *seed ) ) * 4.656612875E-10;

  return r;
}
/******************************************************************************/

void r8_unswap3 ( double *x, double *y, double *z )

/******************************************************************************/
/*
  Purpose:

    R8_UNSWAP3 unswaps three real items.

  Example:

    Input:

      X = 2, Y = 3, Z = 1

    Output:

      X = 1, Y = 2, Z = 3

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input/output, double *X, *Y, *Z, three values to be swapped.
*/
{
  double w;

   w = *z;
  *z = *y;
  *y = *x;
  *x =  w;

  return;
}
/******************************************************************************/

double r8_walsh_1d ( double x, int digit )

/******************************************************************************/
/*
  Purpose:

    R8_WALSH_1D evaluates the Walsh function of a real scalar argument.

  Discussion:

    Consider the binary representation of X, and number the digits
    in descending order, from leading to lowest, with the units digit
    being numbered 0.

    The Walsh function W(J)(X) is equal to the J-th binary digit of X.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the Walsh function.

    Input, int DIGIT, the index of the Walsh function.

    Output, double R8_WALSH_1D, the value of the Walsh function.
*/
{
  int n;
  static double two = 2.0;
  double value;
/*
  Hide the effect of the sign of X.
*/
  x = r8_abs ( x );
/*
  If DIGIT is positive, divide by 2 DIGIT times.
  If DIGIT is negative, multiply by 2 (-DIGIT) times.
*/
  x = x / pow ( two, digit );
/*
  Make it an integer.
  Because it's positive, and we're using INT, we don't change the
  units digit.
*/
  n = ( int ) x;
/*
  Is the units digit odd or even?
*/
  if ( ( n % 2 ) == 0 )
  {
    value = 0.0;
  }
  else
  {
    value = 1.0;
  }

  return value;
}
/******************************************************************************/

double r8_wrap ( double r, double rlo, double rhi )

/******************************************************************************/
/*
  Purpose:

    R8_WRAP forces an R8 to lie between given limits by wrapping.

  Discussion:

    An R8 is a double value.

  Example:

    RLO = 4.0, RHI = 8.0

     R  Value

    -2     8
    -1     4
     0     5
     1     6
     2     7
     3     8
     4     4
     5     5
     6     6
     7     7
     8     8
     9     4
    10     5
    11     6
    12     7
    13     8
    14     4

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, double R, a value.

    Input, double RLO, RHI, the desired bounds.

    Output, double R8_WRAP, a "wrapped" version of the value.
*/
{
  int n;
  double rhi2;
  double rlo2;
  double rwide;
  double value;
/*
  Guarantee RLO2 < RHI2.
*/
  rlo2 = r8_min ( rlo, rhi );
  rhi2 = r8_max ( rlo, rhi );
/*
  Find the width.
*/
  rwide = rhi2 - rlo2;
/*
  Add enough copies of (RHI2-RLO2) to R so that the
  result ends up in the interval RLO2 - RHI2.
*/
  if ( rwide == 0.0 )
  {
    value = rlo;
  }
  else if ( r < rlo2 )
  {
    n = ( int ) ( ( rlo2 - r ) / rwide ) + 1;
    value = r + n * rwide;
    if ( value == rhi )
    {
      value = rlo;
    }
  }
  else
  {
    n = ( int ) ( ( r - rlo2 ) / rwide );
    value = r - n * rwide;
    if ( value == rlo )
    {
      value = rhi;
    }
  }
  return value;
}
/******************************************************************************/

double r82_dist_l2 ( double a1[2], double a2[2] )

/******************************************************************************/
/*
  Purpose:

    R82_DIST_L2 returns the L2 distance between a pair of R82's.

  Discussion:

    An R82 is a vector of type R8, with two entries.

    The vector L2 norm is defined as:

      sqrt ( sum ( 1 <= I <= N ) A(I) * A(I) ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, double A1[2], A2[2], the vectors.

    Output, double R82_DIST_L2, the L2 norm of A1 - A2.
*/
{
  double value;

  value = sqrt ( pow ( a1[0] - a2[0], 2 )
               + pow ( a1[1] - a2[1], 2 ) );

  return value;
}
/******************************************************************************/

void r82_print ( double a[2], char *title )

/******************************************************************************/
/*
  Purpose:

    R82_PRINT prints an R82.

  Discussion:

    An R82 is an R8VEC with two entries.

    A format is used which suggests a coordinate pair:

  Example:

    Center : ( 1.23, 7.45 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, double A[2], the coordinates of the vector.

    Input, char *TITLE, a title.
*/
{
  printf ( "  %s : ( %12g, %12g )\n", title, a[0], a[1] );

  return;
}
/******************************************************************************/

void r82_uniform_ab ( double a, double b, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R82_UNIFORM_AB returns a pseudorandom R82 scaled to [A,B].

  Discussion:

    An R82 is an R8VEC with two entries.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, the minimum and maximum values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R[2], the randomly chosen value.
*/
{
  int i;

  for ( i = 0; i < 2; i++ )
  {
    r[i] = r8_uniform_ab ( a, b, seed );
  }

  return;
}
/******************************************************************************/

void r82poly2_print ( double a, double b, double c, double d, double e,
  double f )

/******************************************************************************/
/*
  Purpose:

    R82POLY2_PRINT prints a second order polynomial in two variables.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, C, D, E, F, the coefficients.
*/
{
  printf ( "  %g * x^2 + %g * y^2 + %g * xy +\n", a, b, c );
  printf ( "  %g * x + %g * y + %g\n", d, e, f );

  return;
}
/******************************************************************************/

int r82poly2_type ( double a, double b, double c, double d, double e, double f )

/******************************************************************************/
/*
  Purpose:

    R82POLY2_TYPE analyzes a second order polynomial in two variables.

  Discussion:

    The polynomial has the form

      A x^2 + B y^2 + C xy + Dx + Ey + F = 0

    The possible types of the solution set are:

     1: a hyperbola;
        9x^2 -  4y^2       -36x - 24y -  36 = 0
     2: a parabola;
        4x^2 +  1y^2 - 4xy + 3x -  4y +   1 = 0;
     3: an ellipse;
        9x^2 + 16y^2       +36x - 32y -  92 = 0;
     4: an imaginary ellipse (no real solutions);
         x^2 +   y^2       - 6x - 10y + 115 = 0;
     5: a pair of intersecting lines;
                        xy + 3x -   y -   3 = 0
     6: one point;
         x^2 +  2y^2       - 2x + 16y +  33 = 0;
     7: a pair of distinct parallel lines;
                 y^2            -  6y +   8 = 0
     8: a pair of imaginary parallel lines (no real solutions);
                 y^2            -  6y +  10 = 0
     9: a pair of coincident lines.
                 y^2            -  2y +   1 = 0
    10: a single line;
                             2x -   y +   1 = 0;
    11; all space;
                                          0 = 0;
    12; no solutions;
                                          1 = 0;

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2011

  Author:

    John Burkardt

  Reference:

    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    CRC Press, 30th Edition, 1996, pages 282-284.

  Parameters:

    Input, double A, B, C, D, E, F, the coefficients.

    Output, int TYPE, indicates the type of the solution set.
*/
{
  double delta;
  double j;
  double k;
  int type;
/*
  Handle the degenerate case.
*/
  if ( a == 0.0 && b == 0.0 && c == 0.0 )
  {
    if ( d == 0.0 && e == 0.0 )
    {
      if ( f == 0.0 )
      {
        type = 11;
      }
      else
      {
        type = 12;
      }
    }
    else
    {
      type = 10;
    }
    return type;
  }

  delta =
      8.0 * a * b * f
    + 2.0 * c * e * d
    - 2.0 * a * e * e
    - 2.0 * b * d * d
    - 2.0 * f * c * c;

  j = 4.0 * a * b - c * c;

  if ( delta != 0.0 )
  {
    if ( j < 0.0 )
    {
      type = 1;
    }
    else if ( j == 0.0 )
    {
      type = 2;
    }
    else if ( 0.0 < j )
    {
      if ( r8_sign ( delta ) != r8_sign ( a + b ) )
      {
        type = 3;
      }
      else if ( r8_sign ( delta ) == r8_sign ( a + b ) )
      {
        type = 4;
      }
    }
  }
  else if ( delta == 0.0 )
  {
    if ( j < 0.0 )
    {
      type = 5;
    }
    else if ( 0.0 < j )
    {
      type = 6;
    }
    else if ( j == 0.0 )
    {
      k = 4.0 * ( a + b ) * f - d * d - e * e;

      if ( k < 0.0 )
      {
        type = 7;
      }
      else if ( 0.0 < k )
      {
        type = 8;
      }
      else if ( k == 0.0 )
      {
        type = 9;
      }
    }
  }

  return type;
}
/******************************************************************************/

void r82poly2_type_print ( int type )

/******************************************************************************/
/*
  Purpose:

    R82POLY2_TYPE_PRINT prints the meaning of the output from R82POLY2_TYPE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int TYPE, the type index returned by R82POLY2_TYPE.
*/
{
  if ( type == 1 )
  {
    printf ( "  The set of solutions forms a hyperbola.\n" );
  }
  else if ( type == 2 )
  {
    printf ( "  The set of solutions forms a parabola.\n" );
  }
  else if ( type == 3 )
  {
    printf ( "  The set of solutions forms an ellipse.\n" );
  }
  else if ( type == 4 )
  {
    printf ( "  The set of solutions forms an imaginary ellipse.\n" );
    printf ( "  (There are no real solutions).\n" );
  }
  else if ( type == 5 )
  {
    printf ( "  The set of solutions forms a pair of intersecting lines.\n" );
  }
  else if ( type == 6 )
  {
    printf ( "  The set of solutions is a single point.\n" );
  }
  else if ( type == 7 )
  {
    printf ( "  The set of solutions form a pair of distinct parallel lines.\n" );
  }
  else if ( type == 8 )
  {
    printf ( "  The set of solutions forms a pair of imaginary parallel lines.\n" );
    printf ( "  (There are no real solutions).\n" );
  }
  else if ( type == 9 )
  {
    printf ( "  The set of solutions forms a pair of coincident lines.\n" );
  }
  else if ( type == 10 )
  {
    printf ( "  The set of solutions forms a single line.\n" );
  }
  else if ( type == 11 )
  {
    printf ( "  The set of solutions is all space.\n" );
  }
  else if ( type == 12 )
  {
    printf ( "  The set of solutions is empty.\n" );
  }
  else
  {
    printf ( "  This type index is unknown.\n" );
  }
  return;
}
/******************************************************************************/

double *r82vec_max ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R82VEC_MAX returns the maximum value in an R82VEC.

  Discussion:

    An R82VEC is an array of pairs of double precision real values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[2*N], the array.

    Output, double R82VEC_MAX[2]; the largest entries in each row.
*/
{
# define DIM_NUM 2

  double *amax = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amax = ( double * ) malloc ( DIM_NUM * sizeof ( double ) );

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amax[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[0+j*DIM_NUM] )
      {
        amax[i] = a[0+j*DIM_NUM];
      }
    }
  }
  return amax;
# undef DIM_NUM
}
/******************************************************************************/

double *r82vec_min ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R82VEC_MIN returns the minimum value in an R82VEC.

  Discussion:

    An R82VEC is an array of pairs of double precision real values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[2*N], the array.

    Output, double R82VEC_MIN[2]; the smallest entries in each row.
*/
{
# define DIM_NUM 2

  double *amin = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amin = ( double * ) malloc ( DIM_NUM * sizeof ( double ) );

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amin[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( a[0+j*DIM_NUM] < amin[i] )
      {
        amin[i] = a[0+j*DIM_NUM];
      }
    }
  }
  return amin;
# undef DIM_NUM
}
/******************************************************************************/

int r82vec_order_type ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R82VEC_ORDER_TYPE finds if an R82VEC is (non)strictly ascending/descending.

  Discussion:

    An R82VEC is a vector whose entries are R82's.
    An R82 is a vector of type double precision with two entries.
    An R82VEC may be stored as a 2 by N array.

    The dictionary or lexicographic ordering is used.

    (X1,Y1) < (X2,Y2)  <=>  X1 < X2 or ( X1 = X2 and Y1 < Y2).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the array.

    Input, double A[2*N], the array to be checked.

    Output, int R82VEC_ORDER_TYPE, order indicator:
    -1, no discernable order;
    0, all entries are equal;
    1, ascending order;
    2, strictly ascending order;
    3, descending order;
    4, strictly descending order.
*/
{
  int i;
  int order;
/*
  Search for the first value not equal to A(1,1).
*/
  i = 0;

  for ( ; ; )
  {
    i = i + 1;

    if ( n <= i )
    {
      order = 0;
      return order;
    }

    if ( a[0+0*2] < a[0+i*2] || ( a[0+0*2] == a[0+i*2] && a[1+0*2] < a[1+i*2] ) )
    {
      if ( i == 2 )
      {
        order = 2;
      }
      else
      {
        order = 1;
      }
      break;
    }
    else if ( a[0+i*2] < a[0+0*2] || ( a[0+i*2] == a[0+0*2] && a[1+i*2] < a[1+0*2] ) )
    {
      if ( i == 2 )
      {
        order = 4;
      }
      else
      {
        order = 3;
      }
      break;
    }
  }
/*
  Now we have a "direction".  Examine subsequent entries.
*/
  for ( ; ; )
  {
    i = i + 1;
    if ( n <= i )
    {
      break;
    }

    if ( order == 1 )
    {
      if ( a[0+i*2] < a[0+(i-1)*2] ||
        ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] < a[1+(i-1)*2] ) )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 2 )
    {
      if ( a[0+i*2] < a[0+(i-1)*2] ||
        ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] < a[1+(i-1)*2] ) )
      {
        order = -1;
        break;
      }
      else if ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] == a[1+(i-1)*2] )
      {
        order = 1;
      }
    }
    else if ( order == 3 )
    {
      if ( a[0+(i-1)*2] < a[0+i*2] ||
        ( a[0+(i-1)*2] == a[0+i*2] && a[1+(i-1)*2] < a[1+i*2] ) )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 4 )
    {
      if ( a[0+(i-1)*2] < a[0+i*2] ||
        ( a[0+(i-1)*2] == a[0+i*2] && a[1+(i-1)*2] < a[1+i*2] ) )
      {
        order = -1;
        break;
      }
      else if ( a[0+i*2] == a[0+(i-1)*2] && a[1+i*2] == a[1+(i-1)*2] )
      {
        order = 3;
      }
    }
  }
  return order;
}
/******************************************************************************/

void r82vec_part_quick_a ( int n, double a[], int *l, int *r )

/******************************************************************************/
/*
  Purpose:

    R82VEC_PART_QUICK_A reorders an R82VEC as part of a quick sort.

  Discussion:

    The routine reorders the entries of A.  Using A(1:2,1) as a
    key, all entries of A that are less than or equal to the key will
    precede the key, which precedes all entries that are greater than the key.

  Example:

    Input:

      N = 8

      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )

    Output:

      L = 2, R = 4

      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
             -----------          ----------------------------------
             LEFT          KEY    RIGHT

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of A.

    Input/output, double A[N*2].  On input, the array to be checked.
    On output, A has been reordered as described above.

    Output, int *L, *R, the indices of A that define the three segments.
    Let KEY = the input value of A(1:2,1).  Then
    I <= L                 A(1:2,I) < KEY;
         L < I < R         A(1:2,I) = KEY;
                 R <= I    A(1:2,I) > KEY.
*/
{
  int i;
  int j;
  double key[2];
  int ll;
  int m;
  int rr;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R82VEC_PART_QUICK_A - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }

  if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key[0] = a[0+0*2];
  key[1] = a[1+0*2];
  m = 1;
/*
  The elements of unknown size have indices between L+1 and R-1.
*/
  ll = 1;
  rr = n + 1;

  for ( i = 2; i <= n; i++ )
  {
    if ( r8vec_gt ( 2, a+2*ll, key ) )
    {
      rr = rr - 1;
      r8vec_swap ( 2, a+2*(rr-1), a+2*ll );
    }
    else if ( r8vec_eq ( 2, a+2*ll, key ) )
    {
      m = m + 1;
      r8vec_swap ( 2, a+2*(m-1), a+2*ll );
      ll = ll + 1;
    }
    else if ( r8vec_lt ( 2, a+2*ll, key ) )
    {
      ll = ll + 1;
    }

  }
/*
  Now shift small elements to the left, and KEY elements to center.
*/
  for ( i = 0; i < ll - m; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = a[2*(i+m)+j];
    }
  }

  ll = ll - m;

  for ( i = ll; i < ll+m; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = key[j];
    }
  }

  *l = ll;
  *r = rr;

  return;
}
/******************************************************************************/

void r82vec_permute ( int n, int p[], int base, double a[] )

/******************************************************************************/
/*
  Purpose:

    R82VEC_PERMUTE permutes an R82VEC in place.

  Discussion:

    An R82VEC is a vector whose entries are R82's.
    An R82 is a vector of R8's with two entries.
    An R82VEC may be stored as a 2 by N array.

    This routine permutes an array of real "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.

  Example:

    Input:

      N = 5
      P = (   2,    4,    5,    1,    3 )
      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
          (11.0, 22.0, 33.0, 44.0, 55.0 )

    Output:

      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of objects.

    Input, int P[N], the permutation.  P(I) = J means
    that the I-th element of the output array should be the J-th
    element of the input array.

    Input, int BASE, is 0 for a 0-based permutation and 1 for a 1-based permutation.

    Input/output, double A[2*N], the array to be permuted.
*/
{
  double a_temp[2];
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm_check ( n, p, base ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R82VEC_PERMUTE - Fatal error!\n" );
    fprintf ( stderr, "  PERM_CHECK rejects this permutation.\n" );
    exit ( 1 );
  }
/*
  In order for the sign negation trick to work, we need to assume that the
  entries of P are strictly positive.  Presumably, the lowest number is BASE.
  So temporarily add 1-BASE to each entry to force positivity.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }
/*
  Search for the next element of the permutation that has not been used.
*/
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp[0] = a[0+(istart-1)*2];
      a_temp[1] = a[1+(istart-1)*2];
      iget = istart;
/*
  Copy the new value into the vacated entry.
*/
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          fprintf ( stderr, "\n" );
          fprintf ( stderr, "R82VEC_PERMUTE - Fatal error!\n" );
          fprintf ( stderr, "  Entry IPUT = %d of the permutation has\n", iput );
          fprintf ( stderr, "  an illegal value IGET = %d.\n", iget );
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[0+(iput-1)*2] = a_temp[0];
          a[1+(iput-1)*2] = a_temp[1];
          break;
        }
        a[0+(iput-1)*2] = a[0+(iget-1)*2];
        a[1+(iput-1)*2] = a[1+(iget-1)*2];
      }
    }
  }
/*
  Restore the signs of the entries.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
/*
  Restore the base of the entries.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 + base;
  }
  return;
}
/******************************************************************************/

void r82vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R82VEC_PRINT prints an R82VEC.

  Discussion:

    An R82VEC is a vector whose entries are R82's.
    An R82 is a vector of R8's with two entries.
    An R82VEC may be stored as a 2 by N array.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[2*N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int j;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( j = 0; j < n; j++ )
  {
    fprintf ( stdout, "  %8d: %14f  %14f\n", j, a[0+j*2], a[1+j*2] );
  }

  return;
}
/******************************************************************************/

void r82vec_print_part ( int n, double a[], int max_print, char *title )

/******************************************************************************/
/*
  Purpose:

    R82VEC_PRINT_PART prints "part" of an R82VEC.

  Discussion:

    The user specifies MAX_PRINT, the maximum number of lines to print.

    If N, the size of the vector, is no more than MAX_PRINT, then
    the entire vector is printed, one entry per line.

    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    followed by a line of periods suggesting an omission,
    and the last entry.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double A[2*N], the vector to be printed.

    Input, int MAX_PRINT, the maximum number of lines
    to print.

    Input, char *TITLE, a title.
*/
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g\n", i, a[0+i*2], a[1+i*2] );
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g\n", i, a[0+i*2], a[1+i*2] );
    }
    fprintf ( stdout, "  ........  ..............  ..............\n" );
    i = n - 1;
    fprintf ( stdout, "  %8d: %14g  %14g\n", i, a[0+i*2], a[1+i*2] );
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g\n", i, a[0+i*2], a[1+i*2] );
    }
    i = max_print - 1;
    fprintf ( stdout, "  %8d: %14g  %14g  ...more entries...\n", 
      i, a[0+i*2], a[1+i*2] );
  }

  return;
}
/******************************************************************************/

int *r82vec_sort_heap_index_a ( int n, int base, double a[] )

/******************************************************************************/
/*
  Purpose:

    R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.

  Discussion:

    An R82VEC is a vector whose entries are R82's.
    An R82 is a vector of R8's with two entries.
    An R82VEC may be stored as a 2 by N array.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    Once the index array is computed, the sorting can be carried out
    "implicitly:

      a(*,indx(*))

    or explicitly, by the call

      r82vec_permute ( n, indx, base, a )

    after which a(*,*) is sorted.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, int BASE, the desired indexing for the sort index:
    0 for 0-based indexing,
    1 for 1-based indexing.

    Input, double A[2*N], an array to be index-sorted.

    Output, int R82VEC_SORT_HEAP_INDEX_A[N], the sort index.  The
    I-th element of the sorted array is A(0:1,R82VEC_SORT_HEAP_INDEX_A(I)).
*/
{
  double aval[2];
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0] + base;
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval[0] = a[0+indxt*2];
      aval[1] = a[1+indxt*2];
    }
    else
    {
      indxt = indx[ir-1];
      aval[0] = a[0+indxt*2];
      aval[1] = a[1+indxt*2];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }
    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if (   a[0+indx[j-1]*2] <  a[0+indx[j]*2] ||
             ( a[0+indx[j-1]*2] == a[0+indx[j]*2] &&
               a[1+indx[j-1]*2] <  a[1+indx[j]*2] ) )
        {
          j = j + 1;
        }
      }

      if (   aval[0] <  a[0+indx[j-1]*2] ||
           ( aval[0] == a[0+indx[j-1]*2] &&
             aval[1] <  a[1+indx[j-1]*2] ) )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
/*
  Take care of the base.
*/
  for ( i = 0; i < n; i++ )
  {
    indx[i] = indx[i] + base;
  }

  return indx;
}
/******************************************************************************/

void r82vec_sort_quick_a ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R82VEC_SORT_QUICK_A ascending sorts an R82VEC using quick sort.

  Discussion:

    A is a two dimensional array of order N by 2, stored as a vector
    of rows: A(0,0), A(0,1),  A(1,0), A(1,1)  ...

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input/output, double A[N*2].
    On input, the array to be sorted.
    On output, the array has been sorted.
*/
{
# define LEVEL_MAX 25

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R82VEC_SORT_QUICK_A - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {
/*
  Partition the segment.
*/
    r82vec_part_quick_a ( n_segment, a+2*(base-1)+0, &l_segment, &r_segment );
/*
  If the left segment has more than one element, we need to partition it.
*/
    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "R82VEC_SORT_QUICK_A - Fatal error!\n" );
        fprintf ( stderr, "  Exceeding recursion maximum of %d\n", LEVEL_MAX );
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
/*
  The left segment and the middle segment are sorted.
  Must the right segment be partitioned?
*/
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
/*
  Otherwise, we back up a level if there is an earlier one.
*/
    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          n_segment = 0;
          break;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }

      }

    }

  }
  return;
# undef LEVEL_MAX
}
/******************************************************************************/

double r83_norm ( double x, double y, double z )

/******************************************************************************/
/*
  Purpose:

    R83_NORM returns the Euclidean norm of an R83.

  Discussion:

    An R83 is a vector of 3 R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, double X, Y, Z, the vector.

    Output, double R83_NORM, the norm of the vector.
*/
{
  double value;

  value = sqrt ( x * x + y * y + z * z );

  return value;
}
/******************************************************************************/

double *r83vec_max ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R83VEC_MAX returns the maximum value in an R83VEC.

  Discussion:

    An R83VEC is an array of triples of double precision real values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 January 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[3*N], the array.

    Output, double R83VEC_MAX[3]; the largest entries in each row.
*/
{
# define DIM_NUM 3

  double *amax = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amax = ( double * ) malloc ( DIM_NUM * sizeof ( double ) );

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amax[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[i+j*DIM_NUM] )
      {
        amax[i] = a[i+j*DIM_NUM];
      }
    }
  }
  return amax;
# undef DIM_NUM
}
/******************************************************************************/

double *r83vec_min ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R83VEC_MIN returns the minimum value in an R83VEC.

  Discussion:

    An R83VEC is an array of triples of double precision real values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[3*N], the array.

    Output, double R83VEC_MIN[3]; the smallest entries in each row.
*/
{
# define DIM_NUM 3

  double *amin = NULL;
  int i;
  int j;

  if ( n <= 0 )
  {
    return NULL;
  }

  amin = ( double * ) malloc ( DIM_NUM * sizeof ( double ) );

  for ( i = 0; i < DIM_NUM; i++ )
  {
    amin[i] = a[i+0*DIM_NUM];
    for ( j = 1; j < n; j++ )
    {
      if ( a[i+j*DIM_NUM] < amin[i] )
      {
        amin[i] = a[i+j*DIM_NUM];
      }
    }
  }
  return amin;
# undef DIM_NUM
}
/******************************************************************************/

void r83vec_part_quick_a ( int n, double a[], int *l, int *r )

/******************************************************************************/
/*
  Purpose:

    R83VEC_PART_QUICK_A reorders an R83VEC as part of a quick sort.

  Discussion:

    An R83VEC is a vector whose entries are R83's.
    An R83 is a vector of type double precision with three entries.
    An R83VEC may be stored as a 3 by N array.

    The routine reorders the entries of A.  Using A(1:3,1) as a
    key, all entries of A that are less than or equal to the key will
    precede the key, which precedes all entries that are greater than the key.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 September 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of A.

    Input/output, double A[3*N].  On input, the array to be checked.
    On output, A has been reordered as described above.

    Output, int *L, *R, the indices of A that define the three segments.
    Let KEY = the input value of A(1:3,1).  Then
    I <= L                 A(1:3,I) < KEY;
         L < I < R         A(1:3,I) = KEY;
                 R <= I    A(1:3,I) > KEY.
*/
{
  int i;
  int j;
  double key[3];
  int ll;
  int m;
  int rr;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R83VEC_PART_QUICK_A - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }

  if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key[0] = a[3*0+0];
  key[1] = a[3*0+1];
  key[2] = a[3*0+2];
  m = 1;
/*
  The elements of unknown size have indices between L+1 and R-1.
*/
  ll = 1;
  rr = n + 1;

  for ( i = 2; i <= n; i++ )
  {
    if ( r8vec_gt ( 3, a+3*ll, key ) )
    {
      rr = rr - 1;
      r8vec_swap ( 3, a+3*(rr-1), a+3*ll );
    }
    else if ( r8vec_eq ( 3, a+3*ll, key ) )
    {
      m = m + 1;
      r8vec_swap ( 3, a+3*(m-1), a+3*ll );
      ll = ll + 1;
    }
    else if ( r8vec_lt ( 3, a+3*ll, key ) )
    {
      ll = ll + 1;
    }
  }
/*
  Now shift small elements to the left, and KEY elements to center.
*/
  for ( i = 0; i < ll - m; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      a[3*i+j] = a[3*(i+m)+j];
    }
  }

  ll = ll - m;

  for ( i = ll; i < ll+m; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      a[3*i+j] = key[j];
    }
  }

  *l = ll;
  *r = rr;

  return;
}
/******************************************************************************/

void r83vec_print_part ( int n, double a[], int max_print, char *title )

/******************************************************************************/
/*
  Purpose:

    R83VEC_PRINT_PART prints "part" of an R83VEC.

  Discussion:

    The user specifies MAX_PRINT, the maximum number of lines to print.

    If N, the size of the vector, is no more than MAX_PRINT, then
    the entire vector is printed, one entry per line.

    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    followed by a line of periods suggesting an omission,
    and the last entry.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double A[3*N], the vector to be printed.

    Input, int MAX_PRINT, the maximum number of lines
    to print.

    Input, char *TITLE, a title.
*/
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g  %14g\n", i, a[0+i*3], a[1+i*3], a[2+i*3] );
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g  %14g\n", i, a[0+i*3], a[1+i*3], a[2+i*3] );
    }
    fprintf ( stdout, "  ........  ..............  ..............  ..............\n" );
    i = n - 1;
    fprintf ( stdout, "  %8d: %14g  %14g  %14g\n", i, a[0+i*3], a[1+i*3], a[2+i*3] );
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      fprintf ( stdout, "  %8d: %14g  %14g  %14g\n", i, a[0+i*3], a[1+i*3], a[2+i*3] );
    }
    i = max_print - 1;
    fprintf ( stdout, "  %8d: %14g  %14g  %14g  ...more entries...\n", 
      i, a[0+i*3], a[1+i*3], a[2+i*3] );
  }

  return;
}
/******************************************************************************/

void r83vec_sort_quick_a ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R83VEC_SORT_QUICK_A ascending sorts an R83VEC using quick sort.

  Discussion:

    An R83VEC is a vector whose entries are R83's.
    An R83 is a vector of type double precision with three entries.
    An R83VEC may be stored as a 3 by N array.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 September 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input/output, double A[N*3].
    On input, the array to be sorted.
    On output, the array has been sorted.
*/
{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R83VEC_SORT_QUICK_A - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {
/*
  Partition the segment.
*/
    r83vec_part_quick_a ( n_segment, a+3*(base-1)+0, &l_segment, &r_segment );
/*
  If the left segment has more than one element, we need to partition it.
*/
    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "R83VEC_SORT_QUICK_A - Fatal error!\n" );
        fprintf ( stderr, "  Exceeding recursion maximum of %d\n", LEVEL_MAX );
        exit ( 1 );
      }
      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
/*
  The left segment and the middle segment are sorted.
  Must the right segment be partitioned?
*/
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
/*
  Otherwise, we back up a level if there is an earlier one.
*/
    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          n_segment = 0;
          break;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }
      }
    }
  }
  return;
# undef LEVEL_MAX
}
/******************************************************************************/

void r8block_delete ( double ***a, int l, int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8BLOCK_DELETE frees memory associated with an R8BLOCK.

  Discussion:

    This function releases the memory associated with an array that was 
    created by a command like
      double ***a;
      a = r8block_new ( l, m, n );

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, double ***A, the pointer to the data.

    Input, int L, M, N, the number of rows, columns, and layers in the array.
*/
{
  int i;
  int j;

  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      free ( a[i][j] );
    }
  }

  for ( i = 0; i < l; i++ )
  {
    free ( a[i] );
  }

  free ( a );

  return;
}
/******************************************************************************/

double *r8block_expand_linear ( int l, int m, int n, double x[], int lfat,
  int mfat, int nfat )

/******************************************************************************/
/*
  Purpose:

    R8BLOCK_EXPAND_LINEAR linearly interpolates new data into a 3D block.

  Discussion:

    In this routine, the expansion is specified by giving the number
    of intermediate values to generate between each pair of original
    data rows and columns.

    The interpolation is not actually linear.  It uses the functions

      1, x, y, z, xy, xz, yz, xyz.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int L, M, N, the dimensions of the input data.

    Input, double X[L*M*N], the original data.

    Input, int LFAT, MFAT, NFAT, the number of data values to interpolate
    original data values in the first, second and third dimensions.

    Output, double XFAT[L2*M2*N2], the fattened data, where
    L2 = (L-1)*(LFAT+1)+1,
    M2 = (M-1)*(MFAT+1)+1,
    N2 = (N-1)*(NFAT+1)+1.
*/
{
  int i;
  int ihi;
  int ii;
  int iii;
  int ip1;
  int j;
  int jhi;
  int jj;
  int jjj;
  int jp1;
  int k;
  int khi;
  int kk;
  int kkk;
  int kp1;
  int l2;
  int m2;
  int n2;
  double r;
  double s;
  double t;
  double x000;
  double x001;
  double x010;
  double x011;
  double x100;
  double x101;
  double x110;
  double x111;
  double *xfat;

  l2 = ( l - 1 ) * ( lfat + 1 ) + 1;
  m2 = ( m - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( n - 1 ) * ( nfat + 1 ) + 1;

  xfat = ( double * ) malloc ( l2 * m2 * n2 * sizeof ( double ) );

  for ( i = 1; i <= l; i++ )
  {
    if ( i < l )
    {
      ihi = lfat;
    }
    else
    {
      ihi = 0;
    }

    for ( j = 1; j <= m; j++ )
    {
      if ( j < m )
      {
        jhi = mfat;
      }
      else
      {
        jhi = 0;
      }

      for ( k = 1; k <= n; k++ )
      {
        if ( k < n )
        {
          khi = nfat;
        }
        else
        {
          khi = 0;
        }

        if ( i < l )
        {
          ip1 = i + 1;
        }
        else
        {
          ip1 = i;
        }

        if ( j < m )
        {
          jp1 = j + 1;
        }
        else
        {
          jp1 = j;
        }

        if ( k < n )
        {
          kp1 = k + 1;
        }
        else
        {
          kp1 = k;
        }

        x000 = x[i-1+(j-1)*l+(k-1)*l*m];
        x001 = x[i-1+(j-1)*l+(kp1-1)*l*m];
        x100 = x[ip1-1+(j-1)*l+(k-1)*l*m];
        x101 = x[ip1-1+(j-1)*l+(kp1-1)*l*m];
        x010 = x[i-1+(jp1-1)*l+(k-1)*l*m];
        x011 = x[i-1+(jp1-1)*l+(kp1-1)*l*m];
        x110 = x[ip1-1+(jp1-1)*l+(k-1)*l*m];
        x111 = x[ip1-1+(jp1-1)*l+(kp1-1)*l*m];

        for ( ii = 0; ii <= ihi; ii++ )
        {
          r = ( double ) ( ii ) / ( double ) ( ihi + 1 );

          for ( jj = 0; jj <= jhi; jj++ )
          {
            s = ( double ) ( jj ) / ( double ) ( jhi + 1 );

            for ( kk = 0; kk <= khi; kk++ )
            {
              t = ( double ) ( kk ) / ( double ) ( khi + 1 );

              iii = 1 + ( i - 1 ) * ( lfat + 1 ) + ii;
              jjj = 1 + ( j - 1 ) * ( mfat + 1 ) + jj;
              kkk = 1 + ( k - 1 ) * ( nfat + 1 ) + kk;

              xfat[iii-1+(jjj-1)*l2+(kkk-1)*l2*m2] =
                  x000 * ( 1.0 - r ) * ( 1.0 - s ) * ( 1.0 - t )
                + x001 * ( 1.0 - r ) * ( 1.0 - s ) * (       t )
                + x010 * ( 1.0 - r ) * (       s ) * ( 1.0 - t )
                + x011 * ( 1.0 - r ) * (       s ) * (       t )
                + x100 * (       r ) * ( 1.0 - s ) * ( 1.0 - t )
                + x101 * (       r ) * ( 1.0 - s ) * (       t )
                + x110 * (       r ) * (       s ) * ( 1.0 - t )
                + x111 * (       r ) * (       s ) * (       t );
            }
          }
        }
      }
    }
  }

  return xfat;
}
/******************************************************************************/

double ***r8block_new ( int l, int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8BLOCK_NEW allocates a new R8BLOCK.

  Discussion:

    A declaration of the form
      double ***a;
    is necesary.  Then an assignment of the form:
      a = r8block_new ( l, m, n );
    allows the user to assign entries to the matrix using typical
    3D array notation:
      a[2][3][4] = 17.0;
      y = a[1][0][3];
    and so on.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 March 2012

  Author:

    John Burkardt

  Parameters:

    Input, int L, M, N, the number of rows, columns and layers.

    Output, double R8BLOCK_NEW[L][M][N], a new block.
*/
{
  double ***a;
  int i;
  int j;

  a = ( double *** ) malloc ( l * sizeof ( double ** ) );

  if ( a == NULL )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8BLOCK_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Unable to allocate row pointer array.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < l; i++ )
  {
    a[i] = ( double ** ) malloc ( m * sizeof ( double * ) );
    if ( a[i] == NULL )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8BLOCK_NEW - Fatal error!\n" );
      fprintf ( stderr, "  Unable to allocate column pointer array.\n" );
      exit ( 1 );
    }
  }

  for ( i = 0; i < l; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      a[i][j] = ( double * ) malloc ( n * sizeof ( double ) );
      if ( a[i][j] == NULL )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "R8BLOCK_NEW - Fatal error!\n" );
        fprintf ( stderr, "  Unable to allocate layer array.\n" );
        exit ( 1 );
      }
    }
  }
  return a;
}
/******************************************************************************/

void r8block_print ( int l, int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8BLOCK_PRINT prints an R8BLOCK block (a 3D matrix).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int L, M, N, the dimensions of the block.

    Input, double A[L*M*N], the matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int jhi;
  int jlo;
  int k;

  printf ( "\n" );
  printf ( "%s\n", title );

  for ( k = 1; k <= n; k++ )
  {
    printf ( "\n" );
    printf ( "  K = %d\n", k );
    printf ( "\n" );
    for ( jlo = 1; jlo <= m; jlo = jlo + 5 )
    {
      jhi = i4_min ( jlo + 4, m );
      printf ( "\n" );
      printf ( "      \n" );
      for ( j = jlo; j <= jhi; j++ )
      {
        printf ( "%7d       ", j );
      }
      printf ( "\n" );
      printf ( "\n" );
      for ( i = 1; i <= l; i++ )
      {
        printf ( "%5d:", i );
        for ( j = jlo; j <= jhi; j++ )
        {
          printf ( "  %12g", a[i-1+(j-1)*l+(k-1)*l*m] );
        }
        printf ( "\n" );
      }
    }
  }

  return;
}
/******************************************************************************/

int r8col_compare ( int m, int n, double a[], int i, int j )

/******************************************************************************/
/*
  Purpose:

    R8COL_COMPARE compares two columns in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Example:

    Input:

      M = 3, N = 4, I = 2, J = 4

      A = (
        1.  2.  3.  4.
        5.  6.  7.  8.
        9. 10. 11. 12. )

    Output:

      R8COL_COMPARE = -1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the M by N array.

    Input, int I, J, the columns to be compared.
    I and J must be between 1 and N.

    Output, int R8COL_COMPARE, the results of the comparison:
    -1, column I < column J,
     0, column I = column J,
    +1, column J < column I.
*/
{
  int k;
  int value;
/*
  Check.
*/
  if ( i < 1 || n < i )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8COL_COMPARE - Fatal error!\n" );
    fprintf ( stderr, "  Column index I is out of bounds.\n" );
    fprintf ( stderr, "  I = %d\n", i );
    exit ( 1 );
  }

  if ( j < 1 || n < j )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8COL_COMPARE - Fatal error!\n" );
    fprintf ( stderr, "  Column index J is out of bounds.\n" );
    fprintf ( stderr, "  J = %d\n", j );
    exit ( 1 );
  }

  value = 0;

  if ( i == j )
  {
    return value;
  }

  k = 0;

  while ( k < m )
  {
    if ( a[k+(i-1)*m] < a[k+(j-1)*m] )
    {
      value = -1;
      return value;
    }
    else if ( a[k+(j-1)*m] < a[k+(i-1)*m] )
    {
      value = +1;
      return value;
    }
    k = k + 1;
  }

  return value;
}
/******************************************************************************/

double *r8col_duplicates ( int m, int n, int n_unique, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8COL_DUPLICATES generates an R8COL with some duplicate columns.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    This routine generates a random R8COL with a specified number of
    duplicate columns.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in each column of A.

    Input, int N, the number of columns in A.

    Input, int N_UNIQUE, the number of unique columns in A.
    1 <= N_UNIQUE <= N.

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, double R8COL_DUPLICATES[M*N], the array.
*/
{
  double *a;
  int i;
  int j1;
  int j2;
  double temp;

  if ( n_unique < 1 || n < n_unique )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8COL_DUPLICATES - Fatal error!\n" );
    fprintf ( stderr, "  1 <= N_UNIQUE <= N is required.\n" );
    exit ( 1 );
  }

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  r8mat_uniform_01 ( m, n_unique, seed, a );
/*
  Randomly copy unique columns.
*/
  for ( j1 = n_unique; j1 < n; j1++ )
  {
    j2 = i4_uniform_ab ( 0, n_unique - 1, seed );
    for ( i = 0; i < m; i++ )
    {
      a[i+j1*m] = a[i+j2*m];
    }
  }
/*
  Permute the columns.
*/
  for ( j1 = 0; j1 < n; j1++ )
  {
    j2 = i4_uniform_ab ( j1, n - 1, seed );
    for ( i = 0; i < m; i++ )
    {
      temp      = a[i+j1*m];
      a[i+j1*m] = a[i+j2*m];
      a[i+j2*m] = temp;
    }
  }
  return a;
}
/******************************************************************************/

int r8col_find ( int m, int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_FIND seeks a column value in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Example:

    Input:

      M = 3,
      N = 4,

      A = (
        1.  2.  3.  4.
        5.  6.  7.  8.
        9. 10. 11. 12. )

      x = ( 3.,
            7.,
           11. )

    Output:

      R8COL_FIND = 3

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], a table of numbers, regarded as
    N columns of vectors of length M.

    Input, double X[M], a vector to be matched with a column of A.

    Output, int R8COL_FIND, the (one-based) index of the first column of A
    which exactly matches every entry of X, or -1 if no match
    could be found.
*/
{
  int col;
  int i;
  int j;

  col = -1;

  for ( j = 1; j <= n; j++ )
  {
    col = j;

    for ( i = 1; i <= m; i++ )
    {
      if ( x[i-1] != a[i-1+(j-1)*m] )
      {
        col = -1;
        break;
      }
    }
    if ( col != -1 )
    {
      return col;
    }
  }
  return col;
}
/******************************************************************************/

int *r8col_first_index ( int m, int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8COL_FIRST_INDEX indexes the first occurrence of values in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    For element A(1:M,J) of the matrix, FIRST_INDEX(J) is the index in A of
    the first column whose entries are equal to A(1:M,J).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of A.
    The length of an "element" of A, and the number of "elements".

    Input, double A[M*N], the array.

    Input, double TOL, a tolerance for equality.

    Output, int R8COL_FIRST_INDEX[N], the first occurrence index.
*/
{
  double diff;
  int *first_index;
  int i;
  int j1;
  int j2;

  first_index = ( int * ) malloc ( n * sizeof ( int ) );

  for ( j1 = 0; j1 < n; j1++ )
  {
    first_index[j1] = -1;
  }
  for ( j1 = 0; j1 < n; j1++ )
  {
    if ( first_index[j1] == -1 )
    {
      first_index[j1] = j1;

      for ( j2 = j1 + 1; j2 < n; j2++ )
      {
        diff = 0.0;
        for ( i = 0; i < m; i++ )
        {
          diff = r8_max ( diff, r8_abs ( a[i+j1*m] - a[i+j2*m] ) );
        }
        if ( diff <= tol )
        {
          first_index[j2] = j1;
        }
      }
    }
  }
  return first_index;
}
/******************************************************************************/

int r8col_insert ( int n_max, int m, int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_INSERT inserts a column into an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Example:

    Input:

      N_MAX = 10,
      M = 3,
      N = 4,

      A = (
        1.  2.  3.  4.
        5.  6.  7.  8.
        9. 10. 11. 12. )

      X = ( 3., 4., 18. )

    Output:

      N = 5,

      A = (
        1.  2.  3.  3.  4.
        5.  6.  4.  7.  8.
        9. 10. 18. 11. 12. )

      R8COL_INSERT = 3

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N_MAX, the maximum number of columns in A.

    Input, int M, the number of rows.

    Input/output, int N, the number of columns.
    If the new column is inserted into the table, then the output
    value of N will be increased by 1.

    Input/output, double A[M*N_MAX], a table of numbers, regarded
    as an array of columns.  The columns must have been sorted
    lexicographically.

    Input, double X[M], a vector of data which will be inserted
    into the table if it does not already occur.

    Output, int R8COL_INSERT.
    I, X was inserted into column I.
    -I, column I was already equal to X.
    0, N = N_MAX.
*/
{
  int col;
  int high;
  int i;
  int isgn;
  int j;
  int low;
  int mid;
/*
  Refuse to work if N_MAX <= N.
*/
  if ( n_max <= n )
  {
    col = 0;
    return col;
  }
/*
  Stick X temporarily in column N+1, just so it's easy to use R8COL_COMPARE.
*/
  for ( i = 0; i < m; i++ )
  {
    a[i+n*m] = x[i];
  }
/*
  Do a binary search.
*/
  low = 1;
  high = n;

  for ( ; ; )
  {
    if ( high < low )
    {
      col = low;
      break;
    }

    mid = ( low + high ) / 2;

    isgn = r8col_compare ( m, n+1, a, mid, n+1 );

    if ( isgn == 0 )
    {
      col = -mid;
      return col;
    }
    else if ( isgn == -1 )
    {
      low = mid + 1;
    }
    else if ( isgn == +1 )
    {
      high = mid - 1;
    }
  }
/*
  Shift part of the table up to make room.
*/
  for ( j = n-1; col-1 <= j; j-- )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+(j+1)*m] = a[i+j*m];
    }
  }
/*
  Insert the new column.
*/
  for ( i = 0; i < m; i++ )
  {
    a[i+(col-1)*m] = x[i];
  }

  n = n + 1;

  return col;
}
/******************************************************************************/

double *r8col_max ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_MAX returns the column maximums of an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the array to be examined.

    Output, double R8COL_MAX[N], the maximums of the columns.
*/
{
  double *amax;
  int i;
  int j;

  amax = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    amax[j] = a[0+j*m];
    for ( i = 0; i < m; i++ )
    {
      amax[j] = r8_max ( amax[j], a[i+j*m] );
    }
  }

  return amax;
}
/******************************************************************************/

int *r8col_max_index ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_MAX_INDEX returns the indices of column maximums in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the array to be examined.

    Output, int R8COL_MAX_INDEX[N]; entry I is the row of A in which
    the maximum for column I occurs.
*/
{
  double amax;
  int i;
  int *imax;
  int j;

  imax = ( int * ) malloc ( n * sizeof ( int ) );

  for ( j = 0; j < n; j++ )
  {
    imax[j] = 1;
    amax = a[0+j*m];

    for ( i = 1; i < m; i++ )
    {
      if ( amax < a[i+j*m] )
      {
        imax[j] = i+1;
        amax = a[i+j*m];
      }
    }
  }

  return imax;
}
/******************************************************************************/

void r8col_max_one ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_MAX_ONE rescales an R8COL so each column maximum is 1.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N], the array to be rescaled.
*/
{
  int i;
  int i_big;
  int j;
  double temp;

  for ( j = 0; j < n; j++ )
  {
    i_big = 0;
    for ( i = 1; i < m; i++ )
    {
      if ( r8_abs ( a[i_big+j*m] ) < r8_abs ( a[i+j*m] ) )
      {
        i_big = i;
      }
    }
    temp = a[i_big+j*m];

    if ( temp != 0.0 )
    {
      for ( i = 0; i < m; i++ )
      {
        a[i+j*m] = a[i+j*m] / temp;
      }
    }
  }
  return;
}
/******************************************************************************/

double *r8col_mean ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_MEAN returns the column means of an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded
    as an array of N columns of length M.

  Example:

    A =
      1  2  3
      2  6  7

    R8COL_MEAN =
      1.5  4.0  5.0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the array to be examined.

    Output, double R8COL_MEAN[N], the means, or averages, of the columns.
*/
{
  int i;
  int j;
  double *mean;

  mean = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    mean[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      mean[j] = mean[j] + a[i+j*m];
    }
    mean[j] = mean[j] / ( double ) ( m );
  }

  return mean;
}
/******************************************************************************/

double *r8col_min ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_MIN returns the column minimums of an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the array to be examined.

    Output, double R8COL_MIN[N], the minimums of the columns.
*/
{
  double *amin;
  int i;
  int j;

  amin = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    amin[j] = a[0+j*m];
    for ( i = 0; i < m; i++ )
    {
      amin[j] = r8_min ( amin[j], a[i+j*m] );
    }
  }

  return amin;
}
/******************************************************************************/

int *r8col_min_index ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_MIN_INDEX returns the indices of column minimums in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the array to be examined.

    Output, int R8COL_MIN_INDEX[N]; entry I is the row of A in which
    the minimum for column I occurs.
*/
{
  double amin;
  int i;
  int *imin;
  int j;

  imin = ( int * ) malloc ( n * sizeof ( int ) );

  for ( j = 0; j < n; j++ )
  {
    imin[j] = 1;
    amin = a[0+j*m];

    for ( i = 1; i < m; i++ )
    {
      if ( a[i+j*m] < amin )
      {
        imin[j] = i+1;
        amin = a[i+j*m];
      }
    }
  }
  return imin;
}
/******************************************************************************/

void r8col_normalize_li ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_NORMALIZE_LI normalizes an R8COL with the column infinity norm.

  Discussion:

    Each column is scaled so that the entry of maximum norm has the value 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N], the array to be normalized.
*/
{
  double c;
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    c = a[0+j*m];
    for ( i = 1; i < m; i++ )
    {
      if ( r8_abs ( c ) < r8_abs ( a[i+j*m] ) )
      {
        c = a[i+j*m];
      }
    }

    if ( c != 0.0 )
    {
      for ( i = 0; i < m; i++ )
      {
        a[i+m*j] = a[i+m*j] / c;
      }
    }
  }

  return;
}
/******************************************************************************/

void r8col_part_quick_a ( int m, int n, double a[], int *l, int *r )

/******************************************************************************/
/*
  Purpose:

    R8COL_PART_QUICK_A reorders the columns of an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The routine reorders the columns of A.  Using A(1:M,1) as a
    key, all entries of A that are less than or equal to the key will
    precede the key, which precedes all entries that are greater than the key.

  Example:

    Input:

      M = 2, N = 8
      A = ( 2  8  6  0 10 10  0  5
            4  8  2  2  6  0  6  8 )

    Output:

      L = 2, R = 4

      A = (  0  0  2  8  6 10 10  4
             2  6  4  8  2  6  0  8 )
             ----     -------------
             LEFT KEY     RIGHT

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, the row dimension of A, and the length of a column.

    Input, int N, the column dimension of A.

    Input/output, double A[M*N].  On input, the array to be checked.
    On output, A has been reordered as described above.

    Output, int *L, *R, the indices of A that define the three segments.
    Let KEY = the input value of A(1:M,1).  Then
    I <= L                 A(1:M,I) < KEY;
         L < I < R         A(1:M,I) = KEY;
                 R <= I    KEY < A(1:M,I).
*/
{
  int i;
  int j;
  int k;
  double *key;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8COL_PART_QUICK_A - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }

  if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    key[i] = a[i+0*m];
  }
  k = 1;
/*
  The elements of unknown size have indices between L+1 and R-1.
*/
  *l = 1;
  *r = n + 1;

  for ( j = 1; j < n; j++ )
  {
    if ( r8vec_gt ( m, a+(*l)*m, key ) )
    {
      *r = *r - 1;
      r8vec_swap ( m, a+(*r-1)*m, a+(*l)*m );
    }
    else if ( r8vec_eq ( m, a+(*l)*m, key ) )
    {
      k = k + 1;
      r8vec_swap ( m, a+(k-1)*m, a+(*l)*m );
      *l = *l + 1;
    }
    else if ( r8vec_lt ( m, a+(*l)*m, key ) )
    {
      *l = *l + 1;
    }
  }
/*
  Shift small elements to the left.
*/
  for ( j = 0; j < *l - k; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+(j+k)*m];
    }
  }
/*
  Shift KEY elements to center.
*/
  for ( j = *l-k; j < *l; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = key[i];
    }
  }
/*
  Update L.
*/
  *l = *l - k;

  free ( key );

  return;
}
/******************************************************************************/

void r8col_permute ( int m, int n, int p[], int base, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_PERMUTE permutes an R8COL in place.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    This routine permutes an array of real "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.

  Example:

    Input:

      M = 2
      N = 5
      P = (   2,    4,    5,    1,    3 )
      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
          (11.0, 22.0, 33.0, 44.0, 55.0 )
      BASE = 1

    Output:

      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, the length of objects.

    Input, int N, the number of objects.

    Input, int P[N], the permutation.  P(I) = J means
    that the I-th element of the output array should be the J-th
    element of the input array.

    Input, int BASE, is 0 for a 0-based permutation and 1 for a
    1-based permutation.

    Input/output, double A[M*N], the array to be permuted.
*/
{
  double *a_temp;
  int i;
  int iget;
  int iput;
  int istart;
  int j;

  if ( !perm_check ( n, p, base ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8COL_PERMUTE - Fatal error!\n" );
    fprintf ( stderr, "  PERM_CHECK rejects this permutation.\n" );
    exit ( 1 );
  }
/*
  In order for the sign negation trick to work, we need to assume that the
  entries of P are strictly positive.  Presumably, the lowest number is BASE.
  So temporarily add 1-BASE to each entry to force positivity.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }

  a_temp = ( double * ) malloc ( m * sizeof ( double ) );
/*
  Search for the next element of the permutation that has not been used.
*/
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      for ( i = 0; i < m; i++ )
      {
        a_temp[i] = a[i+(istart-1)*m];
      }
      iget = istart;
/*
  Copy the new value into the vacated entry.
*/
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          fprintf ( stderr, "\n" );
          fprintf ( stderr, "R8COL_PERMUTE - Fatal error!\n" );
          fprintf ( stderr, "  Entry IPUT = %d of the permutation has\n", iput );
          fprintf ( stderr, "  an illegal value IGET = %d.\n", iget );
          exit ( 1 );
        }

        if ( iget == istart )
        {
          for ( i = 0; i < m; i++ )
          {
            a[i+(iput-1)*m] = a_temp[i];
          }
          break;
        }
        for ( i = 0; i < m; i++ )
        {
          a[i+(iput-1)*m] = a[i+(iget-1)*m];
        }
      }
    }
  }
/*
  Restore the signs of the entries.
*/
  for ( j = 0; j < n; j++ )
  {
    p[j] = - p[j];
  }
/*
  Restore the base of the entries.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 +  base;
  }

  free ( a_temp );

  return;
}
/******************************************************************************/

void r8col_sort_heap_a ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORT_HEAP_A ascending heapsorts an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    In lexicographic order, the statement "X < Y", applied to two real
    vectors X and Y of length M, means that there is some index I, with
    1 <= I <= M, with the property that

      X(J) = Y(J) for J < I,
    and
      X(I) < Y(I).

    In other words, the first time they differ, X is smaller.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N].
    On input, the array of N columns of M-vectors.
    On output, the columns of A have been sorted in lexicographic order.
*/
{
  int i;
  int indx;
  int isgn;
  int j;

  if ( m <= 0 )
  {
    return;
  }

  if ( n <= 1 )
  {
    return;
  }
/*
  Initialize.
*/
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
/*
  Call the external heap sorter.
*/
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
/*
  Interchange the I and J objects.
*/
    if ( 0 < indx )
    {
      r8col_swap ( m, n, a, i, j );
    }
/*
  Compare the I and J objects.
*/
    else if ( indx < 0 )
    {
      isgn = r8col_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

int *r8col_sort_heap_index_a ( int m, int n, int base, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    A(*,J1) < A(*,J2) if the first nonzero entry of A(*,J1)-A(*,J2) is negative.

    Once the index array is computed, the sorting can be carried out
    "implicitly:

      A(*,INDX(*)) is sorted,

    Note that the index vector is 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in each column of A.

    Input, int N, the number of columns in A.

    Input, int BASE, the desired indexing for the sort index:
    0 for 0-based indexing,
    1 for 1-based indexing.

    Input, double A[M*N], the array.

    Output, int R8COL_SORT_HEAP_INDEX_A[N], contains the sort index.  The
    I-th column of the sorted array is A(*,INDX(I)).
*/
{
  double *column;
  int i;
  int *indx;
  int indxt;
  int ir;
  int isgn;
  int j;
  int k;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0] + base;
    return indx;
  }

  column = ( double * ) malloc ( m * sizeof ( double ) );

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      for ( k = 0; k < m; k++ )
      {
        column[k] = a[k+indxt*m];
      }
    }
    else
    {
      indxt = indx[ir-1];
      for ( k = 0; k < m; k++ )
      {
        column[k] = a[k+indxt*m];
      }
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        isgn = r8vec_compare ( m, a+indx[j-1]*m, a+indx[j]*m );

        if ( isgn < 0 )
        {
          j = j + 1;
        }
      }

      isgn = r8vec_compare ( m, column, a+indx[j-1]*m );

      if ( isgn < 0 )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
  free ( column );
/*
  Take care of the base.
*/
  for ( i = 0; i < n; i++ )
  {
    indx[i] = indx[i] + base;
  }

  return indx;
}
/******************************************************************************/

void r8col_sort_quick_a ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORT_QUICK_A ascending quick sorts an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, the row order of A, and the length of a column.

    Input, int N, the number of columns of A.

    Input/output, double A[M*N].
    On input, the array to be sorted.
    On output, the array has been sorted.
*/
{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( m <= 0 )
  {
    return;
  }

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8COL_SORT_QUICK_A - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    fprintf ( stderr, "  N = %d\n", n );
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  for ( ; ; )
  {
/*
  Partition the segment.
*/
    r8col_part_quick_a ( m, n_segment, a+(base-1)*m, &l_segment, &r_segment );
/*
  If the left segment has more than one element, we need to partition it.
*/
    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "R8COL_SORT_QUICK_A - Fatal error!\n" );
        fprintf ( stderr, "  Exceeding recursion maximum of %d\n", LEVEL_MAX );
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
/*
  The left segment and the middle segment are sorted.
  Must the right segment be partitioned?
*/
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
/*
  Otherwise, we back up a level if there is an earlier one.
*/
    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          return;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }
      }
    }
  }
  return;
# undef LEVEL_MAX
}
/******************************************************************************/

void r8col_sorted_tol_undex ( int m, int n, double a[], int unique_num,
  double tol, int undx[], int xdnu[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTED_TOL_UNDEX returns tolerably unique indexes for a sorted R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The goal of this routine is to determine a vector UNDX,
    which points, to the unique elements of A, in sorted order,
    and a vector XDNU, which identifies, for each entry of A, the index of
    the unique sorted element of A.

    This is all done with index vectors, so that the elements of
    A are never moved.

    Assuming A is already sorted, we examine the entries of A in order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.

    Once this process has been completed, the vector A could be
    replaced by a compressed vector XU, containing the unique entries
    of A in sorted order, using the formula

      XU(*) = A(UNDX(*)).

    We could then, if we wished, reconstruct the entire vector A, or
    any element of it, by index, as follows:

      A(I) = XU(XDNU(I)).

    We could then replace A by the combination of XU and XDNU.

    Later, when we need the I-th entry of A, we can locate it as
    the XDNU(I)-th entry of XU.

    Here is an example of a vector A, the unique sort and inverse unique
    sort vectors and the compressed unique sorted vector.

      I      A      XU  Undx  Xdnu
    ----+------+------+-----+-----+
      0 | 11.0 |  11.0    0     0
      1 | 11.0 |  22.0    4     0
      2 | 11.0 |  33.0    7     0
      3 | 11.0 |  55.0    8     0
      4 | 22.0 |                1
      5 | 22.0 |                1
      6 | 22.0 |                1
      7 | 33.0 |                2
      8 | 55.0 |                3

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the dimension of the data values.

    Input, int N, the number of data values,

    Input, double A[M*N], the data values.

    Input, int UNIQUE_NUM, the number of unique values in A.
    This value is only required for languages in which the size of
    UNDX must be known in advance.

    Input, double TOL, a tolerance for equality.

    Output, int UNDX[UNIQUE_NUM], the UNDX vector.

    Output, int XDNU[N], the XDNU vector.
*/
{
  double diff;
  int i;
  int i2;
  int i3;
  int j;
  int k;
  int unique;
/*
  Consider entry I = 0.
  It is unique, so set the number of unique items to K.
  Set the K-th unique item to I.
  Set the representative of item I to the K-th unique item.
*/
  i = 0;
  k = 0;
  undx[k] = i;
  xdnu[i] = k;
/*
  Consider entry I.

  If it is unique, increase the unique count K, set the
  K-th unique item to I, and set the representative of I to K.

  If it is not unique, set the representative of item I to a
  previously determined unique item that is close to it.
*/
  for ( i = 1; i < n; i++ )
  {
    unique = 1;

    for ( j = 0; j <= k; j++ )
    {
      i2 = undx[j];
      diff = 0.0;
      for ( i3 = 0; i3 < m; i3++ )
      {
        diff = r8_max ( diff, r8_abs ( a[i3+i*m] - a[i3+i2*m] ) );
      }
      if ( diff <= tol )
      {
        unique = 0;
        xdnu[i] = j;
        break;
      }
    }
    if ( unique )
    {
      k = k + 1;
      undx[k] = i;
      xdnu[i] = k;
    }
  }
  return;
}
/******************************************************************************/

int r8col_sorted_tol_unique ( int m, int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTED_TOL_UNIQUE keeps tolerably unique elements in a sorted R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The columns of the array can be ascending or descending sorted.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A(M,N).
    On input, the sorted array of N columns of M-vectors.
    On output, a sorted array of columns of M-vectors.

    Input, double TOL, a tolerance for equality.

    Output, int R8COL_SORTED_TOL_UNIQUE, the number of unique columns.
*/
{
  double diff;
  int i;
  int j;
  int k;
  int unique;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    unique = 1;
    for ( j = 0; j < unique_num; j++ )
    {
      diff = 0.0;
      for ( k = 0; k < m; k++ )
      {
        diff = r8_max ( diff, r8_abs ( a[k+i*m] - a[k+j*m] ) );
      }
      if ( diff < tol )
      {
        unique = 0;
        break;
      }
    }
    if ( unique )
    {
      for ( k = 0; k < m; k++ )
      {
        a[k+unique_num*m] = a[k+i*m];
      }
      unique_num = unique_num + 1;
    }
  }
  return unique_num;
}
/******************************************************************************/

int r8col_sorted_tol_unique_count ( int m, int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTED_TOL_UNIQUE_COUNT counts tolerably unique elements in a sorted R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The columns of the array may be ascending or descending sorted.

    If the tolerance is large enough, then the concept of uniqueness
    can become ambiguous.  If we have a tolerance of 1.5, then in the
    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
    one unique entry?  That would be because 1 may be regarded as unique,
    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
    be unique and so on.

    This seems wrongheaded.  So I prefer the idea that an item is not
    unique under a tolerance only if it is close to something that IS unique.
    Thus, the unique items are guaranteed to cover the space if we include
    a disk of radius TOL around each one.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], a sorted array, containing
    N columns of data.

    Input, double TOL, a tolerance for equality.

    Output, int R8COL_SORTED_TOL_UNIQUE_COUNT, the number of unique columns.
*/
{
  double diff;
  int i;
  int i2;
  int i3;
  int j;
  int k;
  int *undx;
  int unique;

  undx = ( int * ) malloc ( n * sizeof ( int ) );
/*
  Consider entry I = 0.
  It is unique, so set the number of unique items to K.
  Set the K-th unique item to I.
  Set the representative of item I to the K-th unique item.
*/
  i = 0;
  k = 0;
  undx[k] = i;
/*
  Consider entry I.

  If it is unique, increase the unique count K, set the
  K-th unique item to I, and set the representative of I to K.

  If it is not unique, set the representative of item I to a
  previously determined unique item that is close to it.
*/
  for ( i = 1; i < n; i++ )
  {
    unique = 1;

    for ( j = 0; j <= k; j++ )
    {
      i2 = undx[j];
      diff = 0.0;
      for ( i3 = 0; i3 < m; i3++ )
      {
        diff = r8_max ( diff, r8_abs ( a[i3+i*m] - a[i3+i2*m] ) );
      }
      if ( diff <= tol )
      {
        unique = 0;
        break;
      }
    }
    if ( unique )
    {
      k = k + 1;
      undx[k] = i;
    }
  }
  free ( undx );

  k = k + 1;

  return k;
}
/******************************************************************************/

void r8col_sorted_undex ( int m, int n, double a[], int unique_num,
  int undx[], int xdnu[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTED_UNDEX returns unique indexes for a sorted R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The goal of this routine is to determine a vector UNDX,
    which points, to the unique elements of A, in sorted order,
    and a vector XDNU, which identifies, for each entry of A, the index of
    the unique sorted element of A.

    This is all done with index vectors, so that the elements of
    A are never moved.

    Assuming A is already sorted, we examine the entries of A in order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.

    Once this process has been completed, the vector A could be
    replaced by a compressed vector XU, containing the unique entries
    of A in sorted order, using the formula

      XU(*) = A(UNDX(*)).

    We could then, if we wished, reconstruct the entire vector A, or
    any element of it, by index, as follows:

      A(I) = XU(XDNU(I)).

    We could then replace A by the combination of XU and XDNU.

    Later, when we need the I-th entry of A, we can locate it as
    the XDNU(I)-th entry of XU.

    Here is an example of a vector A, the unique sort and inverse unique
    sort vectors and the compressed unique sorted vector.

      I      A      XU  Undx  Xdnu
    ----+------+------+-----+-----+
      0 | 11.0 |  11.0    0     0
      1 | 11.0 |  22.0    4     0
      2 | 11.0 |  33.0    7     0
      3 | 11.0 |  55.0    8     0
      4 | 22.0 |                1
      5 | 22.0 |                1
      6 | 22.0 |                1
      7 | 33.0 |                2
      8 | 55.0 |                3

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the dimension of the data values.

    Input, int N, the number of data values,

    Input, double A[M*N], the data values.

    Input, int UNIQUE_NUM, the number of unique values in A.
    This value is only required for languages in which the size of
    UNDX must be known in advance.

    Output, int UNDX[UNIQUE_NUM], the UNDX vector.

    Output, int XDNU[N], the XDNU vector.
*/
{
  double diff;
  int i;
  int j;
  int k;
/*
  Walk through the sorted array.
*/
  i = 0;

  j = 0;
  undx[j] = i;

  xdnu[i] = j;

  for ( i = 1; i < n; i++ )
  {
    diff = 0.0;
    for ( k = 0; k < m; k++ )
    {
      diff = r8_max ( diff, r8_abs ( a[k+i*m] - a[k+undx[j]*m] ) );
    }
    if ( 0.0 < diff )
    {
      j = j + 1;
      undx[j] = i;
    }
    xdnu[i] = j;
  }

  return;
}
/******************************************************************************/

int r8col_sorted_unique ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTED_UNIQUE keeps unique elements in a sorted R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The columns of the array can be ascending or descending sorted.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A(M,N).
    On input, the sorted array of N columns of M-vectors.
    On output, a sorted array of columns of M-vectors.

    Output, int UNIQUE_NUM, the number of unique columns.
*/
{
  //double diff;
  int equal;
  int i;
  int j1;
  int j2;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    equal = 1;
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j1*m] != a[i+j2*m] )
      {
        equal = 0;
        break;
      }
    }
    if ( !equal )
    {
      j1 = j1 + 1;
      for ( i = 0; i < m; i++ )
      {
        a[i+j1*m] = a[i+j2*m];
      }
    }
  }

  unique_num = j1 + 1;

  return unique_num;
}
/******************************************************************************/

int r8col_sorted_unique_count ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTED_UNIQUE_COUNT counts unique elements in a sorted R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The columns of the array may be ascending or descending sorted.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], a sorted array, containing
    N columns of data.

    Output, int R8COL_SORTED_UNIQUE_COUNT, the number of unique columns.
*/
{
  //double diff;
  int equal;
  int i;
  int j1;
  int j2;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  unique_num = 1;
  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    equal = 1;
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j1*m] != a[i+j2*m] )
      {
        equal = 0;
        break;
      }
    }
    if ( !equal )
    {
      unique_num = unique_num + 1;
      j1 = j2;
    }
  }

  return unique_num;
}
/******************************************************************************/

void r8col_sortr_a ( int m, int n, double a[], int key )

/******************************************************************************/
/*
  Purpose:

    R8COL_SORTR_A ascending sorts one column of an R8COL, adjusting all entries.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N].
    On input, an unsorted M by N array.
    On output, rows of the array have been shifted in such
    a way that column KEY of the array is in nondecreasing order.

    Input, int KEY, the column in which the "key" value
    is stored.  On output, column KEY of the array will be
    in nondecreasing order.
*/
{
  int i;
  int indx;
  int isgn;
  int j;

  if ( m <= 0 )
  {
    return;
  }

  if ( key < 1 || n < key )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8COL_SORTR_A - Fatal error!\n" );
    fprintf ( stderr, "  The value of KEY is not a legal column index.\n" );
    fprintf ( stderr, "  KEY = %d\n", key );
    fprintf ( stderr, "  N = %d\n", n );
    exit ( 1 );
  }
/*
  Initialize.
*/
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
/*
  Call the external heap sorter.
*/
  for ( ; ; )
  {
    sort_heap_external ( m, &indx, &i, &j, isgn );
/*
  Interchange the I and J objects.
*/
    if ( 0 < indx )
    {
      r8row_swap ( m, n, a, i, j );
    }
/*
  Compare the I and J objects.
*/
    else if ( indx < 0 )
    {
      if ( a[i-1+(key-1)*m] < a[j-1+(key-1)*m] )
      {
        isgn = -1;
      }
      else
      {
        isgn = +1;
      }
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

double *r8col_sum ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_SUM sums the columns of an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the array to be examined.

    Output, double R8COL_SUM[N], the sums of the columns.
*/
{
  double *colsum;
  int i;
  int j;

  colsum = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    colsum[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      colsum[j] = colsum[j] + a[i+j*m];
    }
  }
  return colsum;
}
/******************************************************************************/

void r8col_swap ( int m, int n, double a[], int j1, int j2 )

/******************************************************************************/
/*
  Purpose:

    R8COL_SWAP swaps columns J1 and J2 of an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Example:

    Input:

      M = 3, N = 4, J1 = 2, J2 = 4

      A = (
        1.  2.  3.  4.
        5.  6.  7.  8.
        9. 10. 11. 12. )

    Output:

      A = (
        1.  4.  3.  2.
        5.  8.  7.  6.
        9. 12. 11. 10. )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N], the M by N array.

    Input, int J1, J2, the columns to be swapped.
    These columns are 1-based.
*/
{
  int i;
  double temp;

  if ( j1 < 1 || n < j1 || j2 < 1 || n < j2 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8COL_SWAP - Fatal error!\n" );
    fprintf ( stderr, "  J1 or J2 is out of bounds.\n" );
    fprintf ( stderr, "  J1 =   %d\n", j1 );
    fprintf ( stderr, "  J2 =   %d\n", j2 );
    fprintf ( stderr, "  NCOL = %d\n", n );
    exit ( 1 );
  }

  if ( j1 == j2 )
  {
    return;
  }

  for ( i = 0; i < m; i++ )
  {
    temp          = a[i+(j1-1)*m];
    a[i+(j1-1)*m] = a[i+(j2-1)*m];
    a[i+(j2-1)*m] = temp;
  }

  return;
}
/******************************************************************************/

double *r8col_to_r8vec ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_TO_R8VEC converts an R8COL to an R8VEC.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    This routine is not really useful in our C++ implementation, since
    we actually store an M by N matrix exactly as a vector already.

  Example:

    M = 3, N = 4

    A =
      11 12 13 14
      21 22 23 24
      31 32 33 34

    R8COL_TO_R8VEC = ( 11, 21, 31, 12, 22, 32, 13, 23, 33, 14, 24, 34 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the M by N array.

    Output, double X[M*N], a vector containing the N columns of A.
*/
{
  int i;
  int j;
  int k;
  double *x;

  x = ( double * ) malloc ( m * n * sizeof ( double ) );

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[k] = a[i+j*m];
      k = k + 1;
    }
  }

  return x;
}
/******************************************************************************/

void r8col_tol_undex ( int m, int n, double a[], int unique_num, double tol,
  int undx[], int xdnu[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_TOL_UNDEX indexes tolerably unique entries in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The goal of this routine is to determine a vector UNDX,
    which points to the unique elements of A, in sorted order,
    and a vector XDNU, which identifies, for each entry of A, the index of
    the unique sorted element of A.

    This is all done with index vectors, so that the elements of
    A are never moved.

    The first step of the algorithm requires the indexed sorting
    of A, which creates arrays INDX and XDNI.  (If all the entries
    of A are unique, then these arrays are the same as UNDX and XDNU.)

    We then use INDX to examine the entries of A in sorted order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.

    Once this process has been completed, the vector A could be
    replaced by a compressed vector XU, containing the unique entries
    of A in sorted order, using the formula

      XU(*) = A(UNDX(*)).

    We could then, if we wished, reconstruct the entire vector A, or
    any element of it, by index, as follows:

      A(I) = XU(XDNU(I)).

    We could then replace A by the combination of XU and XDNU.

    Later, when we need the I-th entry of A, we can locate it as
    the XDNU(I)-th entry of XU.

    Here is an example of a vector A, the sort and inverse sort
    index vectors, and the unique sort and inverse unique sort vectors
    and the compressed unique sorted vector.

      I     A  Indx  Xdni       XU  Undx  Xdnu
    ----+-----+-----+-----+--------+-----+-----+
      0 | 11.     0     0 |    11.     0     0
      1 | 22.     2     4 |    22.     1     1
      2 | 11.     5     1 |    33.     3     0
      3 | 33.     8     7 |    55.     4     2
      4 | 55.     1     8 |                  3
      5 | 11.     6     2 |                  0
      6 | 22.     7     5 |                  1
      7 | 22.     3     6 |                  1
      8 | 11.     4     3 |                  0

    INDX(2) = 3 means that sorted item(2) is A(3).
    XDNI(2) = 5 means that A(2) is sorted item(5).

    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
    XDNU(8) = 2 means that A(8) is at unique sorted item(2).

    XU(XDNU(I))) = A(I).
    XU(I)        = A(UNDX(I)).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the dimension of the data values.

    Input, int N, the number of data values,

    Input, double A[M*N], the data values.

    Input, int UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.

    Input, double TOL, a tolerance for equality.

    Output, int UNDX[UNIQUE_NUM], the UNDX vector.

    Output, int XDNU[N], the XDNU vector.
*/
{
  int base = 0;
  double diff;
  int i;
  int i2;
  int *indx;
  int j;
  int k;
  int unique;
/*
  Implicitly sort the array.
*/
  indx = r8col_sort_heap_index_a ( m, n, base, a );
/*
  Consider entry I = 0.
  It is unique, so set the number of unique items to K.
  Set the K-th unique item to I.
  Set the representative of item I to the K-th unique item.
*/
  i = 0;
  k = 0;
  undx[k] = indx[i];
  xdnu[indx[i]] = k;
/*
  Consider entry I.

  If it is unique, increase the unique count K, set the
  K-th unique item to I, and set the representative of I to K.

  If it is not unique, set the representative of item I to a
  previously determined unique item that is close to it.
*/
  for ( i = 1; i < n; i++ )
  {
    unique = 1;
    for ( j = 0; j <= k; j++ )
    {
      diff = 0.0;
      for ( i2 = 0; i2 < m; i2++ )
      {
        diff = r8_max ( diff, r8_abs ( a[i2+indx[i]*m] - a[i2+undx[j]*m] ) );
      }
      if ( diff <= tol )
      {
        unique = 0;
        xdnu[indx[i]] = j;
        break;
      }
    }
    if ( unique )
    {
      k = k + 1;
      undx[k] = indx[i];
      xdnu[indx[i]] = k;
    }
  }
  free ( indx );

  return;
}
/******************************************************************************/

double *r8col_uniform_abvec_new ( int m, int n, double a[], double b[], int *seed )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNIFORM_ABVEC_NEW fills an R8COL with scaled pseudorandom numbers.

  Discussion:

    An R8COL is an array of R8 values, regarded as a set of column vectors.

    The user specifies a minimum and maximum value for each row.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 December 2011

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M], B[M], the upper and lower limits.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has been updated.

    Output, double R8COL_UNIFORM_ABVEC_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int j;
  int k;
  double *r;

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = a[i] 
        + ( b[i] - a[i] ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

int r8col_tol_unique_count ( int m, int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8COL_TOL_UNIQUE_COUNT counts tolerably unique entries in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The columns of the array may be ascending or descending sorted.

    If the tolerance is large enough, then the concept of uniqueness
    can become ambiguous.  If we have a tolerance of 1.5, then in the
    list ( 1, 2, 3, 4, 5, 6, 7, 8, 9 ) is it fair to say we have only
    one unique entry?  That would be because 1 may be regarded as unique,
    and then 2 is too close to 1 to be unique, and 3 is too close to 2 to
    be unique and so on.

    This seems wrongheaded.  So I prefer the idea that an item is not
    unique under a tolerance only if it is close to something that IS unique.
    Thus, the unique items are guaranteed to cover the space if we include
    a disk of radius TOL around each one.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the array of N columns of data.

    Input, double TOL, a tolerance for equality.

    Output, int R8COL_TOL_UNIQUE_COUNT, the number of unique columns.
*/
{
  int base = 0;
  double diff;
  int i;
  int i2;
  int *indx;
  int j;
  int k;
  int *undx;
  int unique;

  undx = ( int * ) malloc ( n * sizeof ( int ) );
/*
  Implicitly sort the array.
*/
  indx = r8col_sort_heap_index_a ( m, n, base, a );
/*
  Consider entry I = 0.
  It is unique, so set the number of unique items to K.
  Set the K-th unique item to I.
  Set the representative of item I to the K-th unique item.
*/
  i = 0;
  k = 0;
  undx[k] = indx[i];
/*
  Consider entry I.

  If it is unique, increase the unique count K, set the
  K-th unique item to I, and set the representative of I to K.

  If it is not unique, set the representative of item I to a
  previously determined unique item that is close to it.
*/
  for ( i = 1; i < n; i++ )
  {
    unique = 1;
    for ( j = 0; j <= k; j++ )
    {
      diff = 0.0;
      for ( i2 = 0; i2 < m; i2++ )
      {
        diff = r8_max ( diff, r8_abs ( a[i2+indx[i]*m] - a[i2+undx[j]*m] ) );
      }
      if ( diff <= tol )
      {
        unique = 0;
        break;
      }
    }
    if ( unique )
    {
      k = k + 1;
      undx[k] = indx[i];
    }
  }
  free ( indx );
  free ( undx );

  k = k + 1;

  return k;
}
/******************************************************************************/

int *r8col_tol_unique_index ( int m, int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8COL_TOL_UNIQUE_INDEX indexes tolerably unique entries in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A,
    gathered in order, then

      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of A.
    The length of an "element" of A, and the number of "elements".

    Input, double A[M*N], the array.

    Input, double TOL, a tolerance for equality.

    Output, int R8COL_TOL_UNIQUE_INDEX[N], the unique index.
*/
{
  double diff;
  int i;
  int j1;
  int j2;
  int *unique_index;
  int unique_num;

  unique_index = ( int * ) malloc ( n * sizeof ( int ) );

  for ( j1 = 0; j1 < n; j1++ )
  {
    unique_index[j1] = -1;
  }
  unique_num = 0;

  for ( j1 = 0; j1 < n; j1++ )
  {
    if ( unique_index[j1] == -1 )
    {
      unique_index[j1] = unique_num;

      for ( j2 = j1 + 1; j2 < n; j2++ )
      {
        diff = 0.0;
        for ( i = 0; i < m; i++ )
        {
          diff = r8_max ( diff, r8_abs ( a[i+j1*m] - a[i+j2*m] ) );
        }
        if ( diff <= tol )
        {
          unique_index[j2] = unique_num;
        }
      }
      unique_num = unique_num + 1;
    }
  }
  return unique_index;
}
/******************************************************************************/

void r8col_undex ( int m, int n, double a[], int unique_num, int undx[],
  int xdnu[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNDEX indexes unique entries in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The goal of this routine is to determine a vector UNDX,
    which points to the unique elements of A, in sorted order,
    and a vector XDNU, which identifies, for each entry of A, the index of
    the unique sorted element of A.

    This is all done with index vectors, so that the elements of
    A are never moved.

    The first step of the algorithm requires the indexed sorting
    of A, which creates arrays INDX and XDNI.  (If all the entries
    of A are unique, then these arrays are the same as UNDX and XDNU.)

    We then use INDX to examine the entries of A in sorted order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.

    Once this process has been completed, the vector A could be
    replaced by a compressed vector XU, containing the unique entries
    of A in sorted order, using the formula

      XU(*) = A(UNDX(*)).

    We could then, if we wished, reconstruct the entire vector A, or
    any element of it, by index, as follows:

      X(I) = XU(XDNU(I)).

    We could then replace A by the combination of XU and XDNU.

    Later, when we need the I-th entry of A, we can locate it as
    the XDNU(I)-th entry of XU.

    Here is an example of a vector A, the sort and inverse sort
    index vectors, and the unique sort and inverse unique sort vectors
    and the compressed unique sorted vector.

      I     A  Indx  Xdni       XU  Undx  Xdnu
    ----+-----+-----+-----+--------+-----+-----+
      0 | 11.     0     0 |    11.     0     0
      1 | 22.     2     4 |    22.     1     1
      2 | 11.     5     1 |    33.     3     0
      3 | 33.     8     7 |    55.     4     2
      4 | 55.     1     8 |                  3
      5 | 11.     6     2 |                  0
      6 | 22.     7     5 |                  1
      7 | 22.     3     6 |                  1
      8 | 11.     4     3 |                  0

    INDX(2) = 3 means that sorted item(2) is A(3).
    XDNI(2) = 5 means that A(2) is sorted item(5).

    UNDX(3) = 4 means that unique sorted item(3) is at A(4).
    XDNU(8) = 2 means that A(8) is at unique sorted item(2).

    XU(XDNU(I))) = A(I).
    XU(I)        = A(UNDX(I)).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the dimension of the data values.

    Input, int N, the number of data values,

    Input, double A[M*N], the data values.

    Input, int UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.

    Output, int UNDX[UNIQUE_NUM], the UNDX vector.

    Output, int XDNU[N], the XDNU vector.
*/
{
  int base = 0;
  double diff;
  int i;
  int *indx;
  int j;
  int k;
/*
  Implicitly sort the array.
*/
  indx = r8col_sort_heap_index_a ( m, n, base, a );
/*
  Walk through the implicitly sorted array X.
*/
  i = 0;

  j = 0;
  undx[j] = indx[i];

  xdnu[indx[i]] = j;

  for ( i = 1; i < n; i++ )
  {
    diff = 0.0;
    for ( k = 0; k < m; k++ )
    {
      diff = r8_max ( diff, r8_abs ( a[k+indx[i]*m] - a[k+undx[j]*m] ) );
    }
    if ( 0.0 < diff )
    {
      j = j + 1;
      undx[j] = indx[i];
    }
    xdnu[indx[i]] = j;
  }
  free ( indx );

  return;
}
/******************************************************************************/

int r8col_unique_count ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNIQUE_COUNT counts unique entries in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    The columns of the array may be ascending or descending sorted.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the array of N columns of data.

    Output, int R8COL_UNIQUE_COUNT, the number of unique columns.
*/
{
  double diff;
  int i;
  int j1;
  int j2;
  int *unique;
  int unique_num;

  unique_num = 0;

  unique = ( int * ) malloc ( n * sizeof ( int ) );

  for ( j1 = 0; j1 < n; j1++ )
  {
    unique_num = unique_num + 1;
    unique[j1] = 1;

    for ( j2 = 0; j2 < j1; j2++ )
    {
      diff = 0.0;
      for ( i = 0; i < m; i++ )
      {
        diff = r8_max ( diff, r8_abs ( a[i+j1*m] - a[i+j2*m] ) );
      }
      if ( diff == 0.0 )
      {
        unique_num = unique_num - 1;
        unique[j1] = 0;
        break;
      }
    }
  }
  free ( unique );

  return unique_num;
}
/******************************************************************************/

int *r8col_unique_index ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_UNIQUE_INDEX indexes unique entries in an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

    For element A(1:M,J) of the matrix, UNIQUE_INDEX(J) is the uniqueness index
    of A(1:M,J).  That is, if A_UNIQUE contains the unique elements of A,
    gathered in order, then

      A_UNIQUE ( 1:M, UNIQUE_INDEX(J) ) = A(1:M,J)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of A.
    The length of an "element" of A, and the number of "elements".

    Input, double A[M*N], the array.

    Output, int R8COL_UNIQUE_INDEX[N], the unique index.
*/
{
  double diff;
  int i;
  int j1;
  int j2;
  int *unique_index;
  int unique_num;

  unique_index = ( int * ) malloc ( n * sizeof ( int ) );

  for ( j1 = 0; j1 < n; j1++ )
  {
    unique_index[j1] = -1;
  }
  unique_num = 0;

  for ( j1 = 0; j1 < n; j1++ )
  {
    if ( unique_index[j1] == -1 )
    {
      unique_index[j1] = unique_num;

      for ( j2 = j1 + 1; j2 < n; j2++ )
      {
        diff = 0.0;
        for ( i = 0; i < m; i++ )
        {
          diff = r8_max ( diff, r8_abs ( a[i+j1*m] - a[i+j2*m] ) );
        }
        if ( diff == 0.0 )
        {
          unique_index[j2] = unique_num;
        }
      }
      unique_num = unique_num + 1;
    }
  }
  return unique_index;
}
/******************************************************************************/

double *r8col_variance ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8COL_VARIANCE returns the variances of an R8COL.

  Discussion:

    An R8COL is an M by N array of R8's, regarded as an array of N columns,
    each of length M.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the array.

    Input, double A[M*N], the array whose variances are desired.

    Output, double R8COL_VARIANCE[N], the variances of the rows.
*/
{
  int i;
  int j;
  double mean;
  double *variance;

  variance = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    mean = 0.0;
    for ( i = 0; i < m; i++ )
    {
      mean = mean + a[i+j*m];
    }
    mean = mean / ( double ) ( m );

    variance[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      variance[j] = variance[j] + pow ( a[i+j*m] - mean, 2 );
    }

    if ( 1 < m )
    {
      variance[j] = variance[j] / ( double ) ( m - 1 );
    }
    else
    {
      variance[j] = 0.0;
    }
  }

  return variance;
}
/******************************************************************************/

double r8int_to_r8int ( double rmin, double rmax, double r, double r2min,
  double r2max )

/******************************************************************************/
/*
  Purpose:

    R8INT_TO_R8INT maps one R8 interval to another.

  Discussion:

    The formula used is

      R2 := R2MIN + ( R2MAX - R2MIN ) * ( R - RMIN ) / ( RMAX - RMIN )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, double RMIN, RMAX, the first range.

    Input, double R, the number to be converted.

    Input, double R2MAX, R2MIN, the second range.

    Output, double R8INT_TO_R8INT, the corresponding value in
    the range [R2MIN,R2MAX].
*/
{
  double  r2;

  if ( rmax == rmin )
  {
    r2 = ( r2max + r2min ) / 2.0;
  }
  else
  {
    r2 = ( ( ( rmax - r        ) * r2min
           + (        r - rmin ) * r2max )
           / ( rmax     - rmin ) );
  }

  return r2;
}
/******************************************************************************/

int r8int_to_i4int ( double rmin, double rmax, double r, int imin, int imax )

/******************************************************************************/
/*
  Purpose:

    R8INT_TO_I4INT maps an R8 interval to an integer interval.

  Discussion:

    The formula used is

      I := IMIN + ( IMAX - IMIN ) * ( R - RMIN ) / ( RMAX - RMIN )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, double RMIN, RMAX, the range.

    Input, double R, the number to be converted.

    Input, int IMAX, IMIN, the integer range.

    Output, int R8INT_TO_I4INT, the corresponding value in the range [IMIN,IMAX].
*/
{
  int i;

  if ( rmax == rmin )
  {
    i = ( imax + imin ) / 2;
  }
  else
  {
    i = r8_nint (
      ( ( rmax - r        ) * ( double ) ( imin )
      + (        r - rmin ) * ( double ) ( imax ) )
      / ( rmax     - rmin ) );
  }

  return i;
}
/******************************************************************************/

void r8mat_add ( int m, int n, double alpha, double a[], double beta, 
  double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ADD computes C = alpha * A + beta * B for R8MAT's.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Thanks to Kjartan Halvorsen for pointing out an error, 09 November 2012.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double ALPHA, the multiplier for A.

    Input, double A[M*N], the first matrix.

    Input, double BETA, the multiplier for A.

    Input, double B[M*N], the second matrix.

    Output, double C[M*N], the sum of alpha*A+beta*B.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
    }
  }
  return;
}
/******************************************************************************/

double *r8mat_add_new ( int m, int n, double alpha, double a[], double beta, 
  double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ADD_NEW computes C = alpha * A + beta * B for R8MAT's.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double ALPHA, the multiplier for A.

    Input, double A[M*N], the first matrix.

    Input, double BETA, the multiplier for A.

    Input, double B[M*N], the second matrix.

    Output, double R8MAT_ADD_NEW[M*N], the sum of alpha*A+beta*B.
*/
{
  double *c;
  int i;
  int j;

  c = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = alpha * a[i+j*m] + beta * b[i+j*m];
    }
  }
  return c;
}
/******************************************************************************/

double r8mat_amax ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_AMAX returns the maximum absolute value entry of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Output, double R8MAT_AMAX, the maximum absolute value entry of A.
*/
{
  int i;
  int j;
  double value;

  value = r8_abs ( a[0+0*m] );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( value < r8_abs ( a[i+j*m] ) )
      {
        value = r8_abs ( a[i+j*m] );
      }
    }
  }
  return value;
}
/******************************************************************************/

double *r8mat_border_add ( int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_BORDER_ADD adds a "border" to an R8MAT.

  Discussion:

    We suppose the input data gives values of a quantity on nodes
    in the interior of a 2D grid, and we wish to create a new table
    with additional positions for the nodes that would be on the
    border of the 2D grid.

                  0 0 0 0 0 0
      * * * *     0 * * * * 0
      * * * * --> 0 * * * * 0
      * * * *     0 * * * * 0
                  0 0 0 0 0 0

    The illustration suggests the situation in which a 3 by 4 array
    is input, and a 5 by 6 array is to be output.

    The old data is shifted to its correct positions in the new array.

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the table data.

    Output, double TABLE2[(M+2)*(N+2)], the augmented table data.
*/
{
  int i;
  int j;
  double *table2;

  table2 = ( double * ) malloc ( ( m + 2 ) * ( n + 2 ) * sizeof ( double ) );

  for ( j = 0; j < n+2; j++ )
  {
    for ( i = 0; i < m+2; i++ )
    {
      if ( i == 0 || i == m+1 || j == 0 || j == n+1 )
      {
        table2[i+j*(m+2)] = 0.0;
      }
      else
      {
        table2[i+j*(m+2)] = table[(i-1)+(j-1)*m];
      }
    }
  }

  return table2;
}
/******************************************************************************/

double *r8mat_border_cut ( int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_BORDER_CUT cuts the "border" of an R8MAT.

  Discussion:

    We suppose the input data gives values of a quantity on nodes
    on a 2D grid, and we wish to create a new table corresponding only
    to those nodes in the interior of the 2D grid.

      0 0 0 0 0 0
      0 * * * * 0    * * * *
      0 * * * * 0 -> * * * *
      0 * * * * 0    * * * *
      0 0 0 0 0 0

    The illustration suggests the situation in which a 5 by 6 array
    is input, and a 3 by 4 array is to be output.

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the table data.

    Output, double TABLE2[(M-2)*(N-2)], the "interior" table data.
*/
{
  int i;
  int j;
  double *table2;

  if ( m <= 2 || n <= 2 )
  {
    return NULL;
  }

  table2 = ( double * ) malloc ( ( m - 2 ) * ( n - 2 ) * sizeof ( double ) );

  for ( j = 0; j < n-2; j++ )
  {
    for ( i = 0; i < m-2; i++ )
    {
      table2[i+j*(m-2)] = table[(i+1)+(j+1)*m];
    }
  }

  return table2;
}
/******************************************************************************/

double *r8mat_cholesky_factor ( int n, double a[], int *flag )

/******************************************************************************/
/*
  Purpose:

    R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The matrix must be symmetric and positive semidefinite.

    For a positive semidefinite symmetric matrix A, the Cholesky factorization
    is a lower triangular matrix L such that:

      A = L * L'

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix A.

    Input, double A[N*N], the N by N matrix.

    Output, int *FLAG, an error flag.
    0, no error was detected.
    1, the matrix is not positive definite.
    2, the matrix is not nonnegative definite.

    Output, double R8MAT_CHOLESKY_FACTOR[N*N], the N by N lower triangular
    Cholesky factor.
*/
{
  double *c;
  int i;
  int j;
  int k;
  double sum2;
  double tol;

  *flag = 0;
  tol = sqrt ( r8_epsilon ( ) );

  c = r8mat_copy_new ( n, n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      c[i+j*n] = 0.0;
    }
    for ( i = j; i < n; i++ )
    {
      sum2 = c[j+i*n];
      for ( k = 0; k < j; k++ )
      {
        sum2 = sum2 - c[j+k*n] * c[i+k*n];
      }
      if ( i == j )
      {
        if ( 0.0 < sum2 )
        {
          c[i+j*n] = sqrt ( sum2 );
        }
        else if ( sum2 < - tol )
        {
          *flag = 2;
          fprintf ( stderr, "\n" );
          fprintf ( stderr, "R8MAT_CHOLESKY_FACTOR - Fatal error!\n" );
          fprintf ( stderr, "  Matrix is not nonnegative definite.\n" );
          fprintf ( stderr, "  Diagonal I = %d\n", i );
          fprintf ( stderr, "  SUM2 = %g\n", sum2 );
          exit ( 1 );
        }
        else
        {
          *flag = 1;
          c[i+j*n] = 0.0;
        }
      }
      else
      {
        if ( c[j+j*n] != 0.0 )
        {
          c[i+j*n] = sum2 / c[j+j*n];
        }
        else
        {
          c[i+j*n] = 0.0;
        }
      }
    }
  }

  return c;
}
/******************************************************************************/

double *r8mat_cholesky_solve ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_CHOLESKY_SOLVE solves a Cholesky factored linear system A * x = b.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix A.

    Input, double A[N*N], the N by N Cholesky factor of the
    system matrix.

    Input, double B[N], the right hand side of the linear system.

    Output, double R8MAT_CHOLESKY_SOLVE[N], the solution of the linear system.
*/
{
  double *x;
  double *y;
/*
  Solve L * y = b.
*/
  y = r8mat_l_solve ( n, a, b );
/*
  Solve L' * x = y.
*/
  x = r8mat_lt_solve ( n, a, y );

  free ( y );

  return x;
}
/******************************************************************************/

double *r8mat_choresky_factor ( int n, double a[], int *flag )

/******************************************************************************/
/*
  Purpose:

    R8MAT_CHORESKY_FACTOR computes the "Choresky" factor of a symmetric R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The matrix must be symmetric and positive semidefinite.

    For a positive semidefinite symmetric matrix A, the Cholesky factorization
    is an upper triangular matrix R such that:

      A = R * R'

    Note that the usual Cholesky factor is a LOWER triangular matrix L
    such that

      A = L * L'

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 April 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix A.

    Input, double A[N*N], the N by N matrix.

    Output, int *FLAG, an error flag.
    0, no error was detected.
    1, the matrix was not positive definite.  A NULL factor was returned.

    Output, double R8MAT_CHORESKY_FACTOR[N*N], the N by N upper triangular
    "Choresky" factor.
*/
{
  double *c;
  int i;
  int j;
  int k;
  double sum2;

  *flag = 0;

  c = r8mat_copy_new ( n, n, a );

  r8mat_flip_rows ( n, n, c );
  r8mat_flip_cols ( n, n, c );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      c[i+j*n] = 0.0;
    }
    for ( i = j; i < n; i++ )
    {
      sum2 = c[j+i*n];
      for ( k = 0; k < j; k++ )
      {
        sum2 = sum2 - c[j+k*n] * c[i+k*n];
      }
      if ( i == j )
      {
        if ( sum2 <= 0.0 )
        {
          *flag = 1;
          return NULL;
        }
        c[i+j*n] = sqrt ( sum2 );
      }
      else
      {
        if ( c[j+j*n] != 0.0 )
        {
          c[i+j*n] = sum2 / c[j+j*n];
        }
        else
        {
          c[i+j*n] = 0.0;
        }
      }
    }
  }

  r8mat_flip_cols ( n, n, c );
  r8mat_flip_rows ( n, n, c );

  return c;
}
/******************************************************************************/

void r8mat_copy ( int m, int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_COPY copies one R8MAT to another.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A1[M*N], the matrix to be copied.

    Output, double A2[M*N], the copy of A1.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_copy_new ( int m, int n, double a1[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A1[M*N], the matrix to be copied.

    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
*/
{
  double *a2;
  int i;
  int j;

  a2 = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }

  return a2;
}
/******************************************************************************/

void r8mat_delete ( double **a, int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DELETE frees the memory set aside by R8MAT_NEW.

  Discussion:

    This function releases the memory associated with an array that was 
    created by a command like:

      double **a;
      a = r8mat_new ( m, n );

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2011

  Author:

    John Burkardt

  Parameters:

    Input, double **A, the array.

    Input, int M, N, the number of rows and columns.
*/
{
  int i;
/*
  Delete the pointers to rows.
*/
  for ( i = 0; i < m; i++ )
  {
    free ( a[i] );
  }
/*
  Delete the pointer to the pointers to rows.
*/
  free ( a );

  return;
}
/******************************************************************************/

double r8mat_det ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DET computes the determinant of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 May 2010

  Author:

    Original FORTRAN77 version by Helmut Spaeth.
    C version by John Burkardt.

  Reference:

    Helmut Spaeth,
    Cluster Analysis Algorithms
    for Data Reduction and Classification of Objects,
    Ellis Horwood, 1980, page 125-127.

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the matrix whose determinant is desired.

    Output, double R8MAT_DET, the determinant of the matrix.
*/
{
  double *b;
  double det;
  int i;
  int j;
  int k;
  int kk;
  int m;
  double temp;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i+j*n] = a[i+j*n];
    }
  }

  det = 1.0;

  for ( k = 1; k <= n; k++ )
  {
    m = k;
    for ( kk = k+1; kk <= n; kk++ )
    {
      if ( r8_abs ( b[m-1+(k-1)*n] ) < r8_abs ( b[kk-1+(k-1)*n] ) )
      {
        m = kk;
      }
    }

    if ( m != k )
    {
      det = -det;

      temp = b[m-1+(k-1)*n];
      b[m-1+(k-1)*n] = b[k-1+(k-1)*n];
      b[k-1+(k-1)*n] = temp;
    }

    det = det * b[k-1+(k-1)*n];

    if ( b[k-1+(k-1)*n] != 0.0 )
    {
      for ( i = k+1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] / b[k-1+(k-1)*n];
      }

      for ( j = k+1; j <= n; j++ )
      {
        if ( m != k )
        {
          temp = b[m-1+(j-1)*n];
          b[m-1+(j-1)*n] = b[k-1+(j-1)*n];
          b[k-1+(j-1)*n] = temp;
        }
        for ( i = k+1; i <= n; i++ )
        {
          b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + b[i-1+(k-1)*n] * b[k-1+(j-1)*n];
        }
      }
    }
  }

  free ( b );

  return det;
}
/******************************************************************************/

double r8mat_det_2d ( double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DET_2D computes the determinant of a 2 by 2 R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Discussion:

    The determinant of a 2 by 2 matrix is

      a11 * a22 - a12 * a21.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double A[2*2], the matrix whose determinant is desired.

    Output, double R8MAT_DET_2D, the determinant of the matrix.
*/
{
  double det;

  det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];

  return det;
}
/******************************************************************************/

double r8mat_det_3d ( double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DET_3D computes the determinant of a 3 by 3 R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The determinant of a 3 by 3 matrix is

        a11 * a22 * a33 - a11 * a23 * a32
      + a12 * a23 * a31 - a12 * a21 * a33
      + a13 * a21 * a32 - a13 * a22 * a31

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double A[3*3], the matrix whose determinant is desired.

    Output, double R8MAT_DET_3D, the determinant of the matrix.
*/
{
  double det;

  det =
      a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
    + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
    + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );

  return det;
}
/******************************************************************************/

double r8mat_det_4d ( double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double A[4*4], the matrix whose determinant is desired.

    Output, double R8MAT_DET_4D, the determinant of the matrix.
*/
{
  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
}
/******************************************************************************/

double r8mat_det_5d ( double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DET_5D computes the determinant of a 5 by 5 R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double A[5*5], the matrix whose determinant is desired.

    Output, double R8MAT_DET_5D, the determinant of the matrix.
*/
{
  double b[4*4];
  double det;
  int i;
  int inc;
  int j;
  int k;
  double sign;
/*
  Expand the determinant into the sum of the determinants of the
  five 4 by 4 matrices created by dropping row 1, and column k.
*/
  det = 0.0;
  sign = 1.0;

  for ( k = 0; k < 5; k++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      for ( j = 0; j < 4; j++ )
      {
        if ( j < k )
        {
          inc = 0;
        }
        else
        {
          inc = 1;
        }
        b[i+j*4] = a[i+1+(j+inc)*5];
      }
    }

    det = det + sign * a[0+k*5] * r8mat_det_4d ( b );

    sign = - sign;
  }

  return det;
}
/******************************************************************************/

void r8mat_diag_add_scalar ( int n, double a[], double s )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix.

    Input/output, double A[N*N], the N by N matrix to be modified.

    Input, double S, the value to be added to the diagonal
    of the matrix.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = a[i+i*n] + s;
  }

  return;
}
/******************************************************************************/

void r8mat_diag_add_vector ( int n, double a[], double v[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DIAG_ADD_VECTOR adds a vector to the diagonal of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix.

    Input/output, double A[N*N], the N by N matrix.

    Input, double V[N], the vector to be added to the diagonal of A.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = a[i+i*n] + v[i];
  }

  return;
}
/******************************************************************************/

double *r8mat_diag_get_vector ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix.

    Input, double A[N*N], the N by N matrix.

    Output, double R8MAT_DIAG_GET_VECTOR[N], the diagonal entries
    of the matrix.
*/
{
  int i;
  double *v;

  v = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    v[i] = a[i+i*n];
  }

  return v;
}
/******************************************************************************/

void r8mat_diag_set_scalar ( int n, double a[], double s )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DIAG_SET_SCALAR sets the diagonal of an R8MAT to a scalar value.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix.

    Input/output, double A[N*N], the N by N matrix to be modified.

    Input, double S, the value to be assigned to the diagonal
    of the matrix.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = s;
  }

  return;
}
/******************************************************************************/

void r8mat_diag_set_vector ( int n, double a[], double v[] )

/******************************************************************************/
/*
//  Purpose:
//
//    R8MAT_DIAG_SET_VECTOR sets the diagonal of an R8MAT to a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input/output, double A[N*N], the N by N matrix.
//
//    Input, double V[N], the vector to be assigned to the
//    diagonal of A.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i+i*n] = v[i];
  }

  return;
}
/******************************************************************************/

double r8mat_dif_fro ( int m, int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_DIF_FRO returns the Frobenius norm of the difference of R8MAT's.

  Discussion: 							    

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The Frobenius norm is defined as

      R8MAT_NORM_FRO = sqrt (
        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )

    The matrix Frobenius norm is not derived from a vector norm, but
    is compatible with the vector L2 norm, so that:

      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    19 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], double B[M*N], the matrices for which we
    want the Frobenius norm of the difference.

    Output, double R8MAT_DIF_FRO, the Frobenius norm of ( A - B ).
*/
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m] - b[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

double *r8mat_expand_linear ( int m, int n, double x[], int mfat, int nfat )

/******************************************************************************/
/*
  Purpose:

    R8MAT_EXPAND_LINEAR linearly interpolates new data into an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    In this routine, the expansion is specified by giving the number
    of intermediate values to generate between each pair of original
    data rows and columns.

    The interpolation is not actually linear.  It uses the functions

      1, x, y, and xy.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of input data.

    Input, double X[M*N], the original data.

    Input, int MFAT, NFAT, the number of data values to interpolate
    between each row, and each column, of original data values.

    Output, double XFAT[M2*N2], the fattened data, where
    M2 = (M-1)*(MFAT+1)+1,
    N2 = (N-1)*(NFAT+1)+1.
*/
{
  int i;
  int ihi;
  int ii;
  int iii;
  int ip1;
  int j;
  int jhi;
  int jj;
  int jjj;
  int jp1;
  int m2;
  int n2;
  double s;
  double t;
  double x00;
  double x01;
  double x10;
  double x11;
  double *xfat;

  m2 = ( m - 1 ) * ( mfat + 1 ) + 1;
  n2 = ( n - 1 ) * ( nfat + 1 ) + 1;

  xfat = ( double * ) malloc ( m2 * n2 * sizeof ( double ) );

  for ( i = 1; i <= m; i++ )
  {
    if ( i < m )
    {
      ihi = mfat;
    }
    else
    {
      ihi = 0;
    }

    for ( j = 1; j <= n; j++ )
    {
      if ( j < n )
      {
        jhi = nfat;
      }
      else
      {
        jhi = 0;
      }

      if ( i < m )
      {
        ip1 = i + 1;
      }
      else
      {
        ip1 = i;
      }

      if ( j < n )
      {
        jp1 = j + 1;
      }
      else
      {
        jp1 = j;
      }

      x00 = x[i-1+(j-1)*m];
      x10 = x[ip1-1+(j-1)*m];
      x01 = x[i-1+(jp1-1)*m];
      x11 = x[ip1-1+(jp1-1)*m];

      for ( ii = 0; ii <= ihi; ii++ )
      {
        s = ( double ) ( ii ) / ( double ) ( ihi + 1 );

        for ( jj = 0; jj <= jhi; jj++ )
        {
          t = ( double ) ( jj ) / ( double ) ( jhi + 1 );

          iii = 1 + ( i - 1 ) * ( mfat + 1 ) + ii;
          jjj = 1 + ( j - 1 ) * ( nfat + 1 ) + jj;

          xfat[iii-1+(jjj-1)*m2] =
                                            x00
              + s     * (       x10       - x00 )
              + t     * (             x01 - x00 )
              + s * t * ( x11 - x10 - x01 + x00 );
        }
      }
    }
  }

  return xfat;
}
/******************************************************************************/

double *r8mat_expand_linear2 ( int m, int n, double a[], int m2, int n2 )

/******************************************************************************/
/*
  Purpose:

    R8MAT_EXPAND_LINEAR2 expands an R8MAT by linear interpolation.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    In this version of the routine, the expansion is indicated
    by specifying the dimensions of the expanded array.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in A.

    Input, double A(M,N), a "small" M by N array.

    Input, int M2, N2, the number of rows and columns in A2.

    Output, double R8MAT_EXPAND_LINEAR2[M2*N2], the expanded array,
    which contains an interpolated version of the data in A.
*/
{
  double *a2;
  int i;
  int i1;
  int i2;
  int j;
  int j1;
  int j2;
  double r;
  double r1;
  double r2;
  double s;
  double s1;
  double s2;

  a2 = ( double * ) malloc ( m2 * n2 * sizeof ( double ) );

  for ( i = 1; i <= m2; i++ )
  {
    if ( m2 == 1 )
    {
      r = 0.5;
    }
    else
    {
      r = ( double ) ( i - 1 ) / ( double ) ( m2 - 1 );
    }

    i1 = 1 + ( int ) ( r * ( double ) ( m - 1 ) );
    i2 = i1 + 1;

    if ( m < i2 )
    {
      i1 = m - 1;
      i2 = m;
    }

    r1 = ( double ) ( i1 - 1 ) / ( double ) ( m - 1 );
    r2 = ( double ) ( i2 - 1 ) / ( double ) ( m - 1 );

    for ( j = 1; j <= n2; j++ )
    {
      if ( n2 == 1 )
      {
        s = 0.5;
      }
      else
      {
        s = ( double ) ( j - 1 ) / ( double ) ( n2 - 1 );
      }

      j1 = 1 + ( int ) ( s * ( double ) ( n - 1 ) );
      j2 = j1 + 1;

      if ( n < j2 )
      {
        j1 = n - 1;
        j2 = n;
      }

      s1 = ( double ) ( j1 - 1 ) / ( double ) ( n - 1 );
      s2 = ( double ) ( j2 - 1 ) / ( double ) ( n - 1 );

      a2[i-1+(j-1)*m2] =
        ( ( r2 - r ) * ( s2 - s ) * a[i1-1+(j1-1)*m]
        + ( r - r1 ) * ( s2 - s ) * a[i2-1+(j1-1)*m]
        + ( r2 - r ) * ( s - s1 ) * a[i1-1+(j2-1)*m]
        + ( r - r1 ) * ( s - s1 ) * a[i2-1+(j2-1)*m] )
        / ( ( r2 - r1 ) * ( s2 - s1 ) );
    }
  }

  return a2;
}
/******************************************************************************/

void r8mat_flip_cols ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_FLIP_COLS swaps the columns of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    To "flip" the columns of an R8MAT is to start with something like

      11 12 13 14 15
      21 22 23 24 25
      31 32 33 34 35
      41 42 43 44 45
      51 52 53 54 55

    and return

      15 14 13 12 11
      25 24 23 22 21
      35 34 33 32 31
      45 44 43 42 41
      55 54 53 52 51

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N], the matrix whose columns are to be flipped.
*/
{
  int i;
  int j;
  double t;

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < ( n / 2 ); j++ )
    {
      t              = a[i+     j *m];
      a[i+     j *m] = a[i+(n-1-j)*m];
      a[i+(n-1-j)*m] = t;
    }
  }
  return;
}
/******************************************************************************/

void r8mat_flip_rows ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_FLIP_ROWS swaps the rows of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    To "flip" the rows of an R8MAT is to start with something like

      11 12 13 14 15
      21 22 23 24 25
      31 32 33 34 35
      41 42 43 44 45
      51 52 53 54 55

    and return

      51 52 53 54 55
      41 42 43 44 45
      31 32 33 34 35
      21 22 23 24 25
      11 12 13 14 15

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N], the matrix whose rows are to be flipped.
*/
{
  int i;
  int j;
  double t;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < ( m / 2 ); i++ )
    {
      t            = a[    i+j*m];
      a[    i+j*m] = a[m-1-i+j*m];
      a[m-1-i+j*m] = t;
    }
  }
  return;
}
/******************************************************************************/

void r8mat_fss ( int n, double a[], int nb, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_FSS factors and solves a system with multiple right hand sides.

  Discussion:

    This routine uses partial pivoting, but no pivot vector is required.

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, double A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.
    Input, int NB, the number of right hand sides.

    Input/output, double X[N*NB], on input, the right hand sides of the
    linear systems.  On output, the solutions of the linear systems.
*/
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;

  for ( jcol = 1; jcol <= n; jcol++ )
  {
/*
  Find the maximum element in column I.
*/
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8MAT_FSS - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", jcol );
      exit ( 1 );
    }
/*
  Switch rows JCOL and IPIV, and X.
*/
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
/*
  Scale the pivot row.
*/
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
/*
  Use the pivot row to eliminate lower entries in that column.
*/
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
/*
  Back solve.
*/
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_fss_new ( int n, double a[], int nb, double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_FSS_NEW factors and solves a system with multiple right hand sides.

  Discussion:

    This routine uses partial pivoting, but no pivot vector is required.

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input/output, double A[N*N].
    On input, A is the coefficient matrix of the linear system.
    On output, A is in unit upper triangular form, and
    represents the U factor of an LU factorization of the
    original coefficient matrix.

    Input, int NB, the number of right hand sides.

    Input, double B[N*NB], the right hand sides of the linear systems.

    Output, double R8MAT_FSS_NEW[N*NB], the solutions of the linear systems.
*/
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;
  double *x;

  x = ( double * ) malloc ( n * nb * sizeof ( double ) );

  for ( j = 0; j < nb; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i+j*n] = b[i+j*n];
    }
  }
  for ( jcol = 1; jcol <= n; jcol++ )
  {
/*
  Find the maximum element in column I.
*/
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8MAT_FSS_NEW - Fatal error!\n" );
      fprintf ( stderr, "  Zero pivot on step %d\n", jcol );
      exit ( 1 );
    }
/*
  Switch rows JCOL and IPIV, and X.
*/
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t            = x[jcol-1+j*n];
        x[jcol-1+j*n] = x[ipiv-1+j*n];
        x[ipiv-1+j*n] = t;
      }
    }
/*
  Scale the pivot row.
*/
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      x[jcol-1+j*n] = x[jcol-1+j*n] / t;
    }
/*
  Use the pivot row to eliminate lower entries in that column.
*/
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          x[i-1+j*n] = x[i-1+j*n] + t * x[jcol-1+j*n];
        }
      }
    }
  }
/*
  Back solve.
*/
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        x[i-1+j*n] = x[i-1+j*n] - a[i-1+(jcol-1)*n] * x[jcol-1+j*n];
      }
    }
  }

  return x;
}
/******************************************************************************/

double *r8mat_givens_post ( int n, double a[], int row, int col )

/******************************************************************************/
/*
  Purpose:

    R8MAT_GIVENS_POST computes the Givens postmultiplier rotation matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The Givens post-multiplier matrix G(ROW,COL) has the property that
    the (ROW,COL)-th entry of A*G is zero.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrices A and G.

    Input, double A[N*N], the matrix to be operated upon.

    Input, int ROW, COL, the row and column of the
    entry of A*G which is to be zeroed out.

    Output, double R8MAT_GIVENS_POST[N*N], the Givens rotation matrix.
    G is an orthogonal matrix, that is, the inverse of
    G is the transpose of G.
*/
{
  double *g;
  double theta;

  g = r8mat_identity_new ( n );

  theta = atan2 ( a[row-1+(col-1)*n], a[row-1+(row-1)*n] );

  g[row-1+(row-1)*n] =  cos ( theta );
  g[row-1+(col-1)*n] = -sin ( theta );
  g[col-1+(row-1)*n] =  sin ( theta );
  g[col-1+(col-1)*n] =  cos ( theta );

  return g;
}
/******************************************************************************/

double *r8mat_givens_pre ( int n, double a[], int row, int col )

/******************************************************************************/
/*
  Purpose:

    R8MAT_GIVENS_PRE computes the Givens premultiplier rotation matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The Givens premultiplier rotation matrix G(ROW,COL) has the
    property that the (ROW,COL)-th entry of G*A is zero.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrices A and G.

    Input, double A[N*N], the matrix to be operated upon.

    Input, int ROW, COL, the row and column of the
    entry of the G*A which is to be zeroed out.

    Output, double R8MAT_GIVENS_PRE[N*N], the Givens rotation matrix.
    G is an orthogonal matrix, that is, the inverse of
    G is the transpose of G.
*/
{
  double *g;
  double theta;

  g = r8mat_identity_new ( n );

  theta = atan2 ( a[row-1+(col-1)*n], a[col-1+(col-1)*n] );

  g[row-1+(row-1)*n] =  cos ( theta );
  g[row-1+(col-1)*n] = -sin ( theta );
  g[col-1+(row-1)*n] =  sin ( theta );
  g[col-1+(col-1)*n] =  cos ( theta );

  return g;
}
/******************************************************************************/

double *r8mat_hess ( double (*fx) ( int n, double x[] ), int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HESS approximates a Hessian matrix via finite differences.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    H(I,J) = d2 F / d X(I) d X(J)

    The values returned by this routine will be only approximate.
    In some cases, they will be so poor that they are useless.
    However, one of the best applications of this routine is for
    checking your own Hessian calculations, since as Heraclitus
    said, you'll never get the same result twice when you differentiate
    a complicated expression by hand.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, double *FX ( int N, double X[] ), the name of the user
    function routine.

    Input, int N, the number of variables.

    Input, double X[N], the values of the variables.

    Output, double H[N*N], the approximated N by N Hessian matrix.
*/
{
  double eps;
  double f00;
  double fmm;
  double fmp;
  double fpm;
  double fpp;
  double *h;
  int i;
  int j;
  double *s;
  double xi;
  double xj;
/*
  Choose the stepsizes.
*/
  s = ( double * ) malloc ( n * sizeof ( double ) );

  eps = pow ( r8_epsilon ( ), 0.33 );

  for ( i = 0; i < n; i++ )
  {
    s[i] = eps * r8_max ( r8_abs ( x[i] ), 1.0 );
  }
/*
  Calculate the diagonal elements.
*/
  h = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    xi = x[i];

    f00 = fx ( n, x );

    x[i] = xi + s[i];
    fpp = fx ( n, x );

    x[i] = xi - s[i];
    fmm = fx ( n, x );

    h[i+i*n] = ( ( fpp - f00 ) + ( fmm - f00 ) ) / s[i] / s[i];

    x[i] = xi;
  }
/*
  Calculate the off diagonal elements.
*/
  for ( i = 0; i < n; i++ )
  {
    xi = x[i];

    for ( j = i+1; j < n; j++ )
    {
      xj = x[j];

      x[i] = xi + s[i];
      x[j] = xj + s[j];
      fpp = fx ( n, x );

      x[i] = xi + s[i];
      x[j] = xj - s[j];
      fpm = fx ( n, x );

      x[i] = xi - s[i];
      x[j] = xj + s[j];
      fmp = fx ( n, x );

      x[i] = xi - s[i];
      x[j] = xj - s[j];
      fmm = fx ( n, x );

      h[j+i*n] = ( ( fpp - fpm ) + ( fmm - fmp ) ) / ( 4.0 * s[i] * s[j] );

      h[i+j*n] = h[j+i*n];

      x[j] = xj;
    }
    x[i] = xi;
  }

  free ( s );

  return h;
}
/******************************************************************************/

void r8mat_house_axh ( int n, double a[], double v[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.

  Discussion: 							    

    An R8MAT is a doubly dimensioned array of double precision values, which
    may be stored as a vector in column-major order.

    The Householder matrix H(V) is defined by

      H(V) = I - 2 * v * v' / ( v' * v )

    This routine is not particularly efficient.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    07 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Input/output, double A[N*N], on input, the matrix to be postmultiplied.
    On output, A has been replaced by A*H.

    Input, double V[N], a vector defining a Householder matrix.
*/
{
  double *ah;
  int i;
  int j;
  int k;
  double v_normsq;

  v_normsq = 0.0;
  for ( i = 0; i < n; i++ )
  {
    v_normsq = v_normsq + v[i] * v[i];
  }
/*
  Compute A*H' = A*H
*/
  ah = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      ah[i+j*n] = a[i+j*n];
      for ( k = 0; k < n; k++ )
      {
        ah[i+j*n] = ah[i+j*n] - 2.0 * a[i+k*n] * v[k] * v[j] / v_normsq;
      }
    }
  }
/*
  Copy A = AH;
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i+j*n] = ah[i+j*n];
    }
  }
  free ( ah );

  return;
}
/******************************************************************************/

double *r8mat_house_axh_new ( int n, double a[], double v[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HOUSE_AXH_NEW computes A*H where H is a compact Householder matrix.

  Discussion: 							    

    An R8MAT is a doubly dimensioned array of double precision values, which
    may be stored as a vector in column-major order.

    The Householder matrix H(V) is defined by

      H(V) = I - 2 * v * v' / ( v' * v )

    This routine is not particularly efficient.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Input, double A[N*N], the matrix to be postmultiplied.

    Input, double V[N], a vector defining a Householder matrix.

    Output, double R8MAT_HOUSE_AXH[N*N], the product A*H.
*/
{
  double *ah;
  int i;
  int j;
  int k;
  double v_normsq;

  v_normsq = 0.0;
  for ( i = 0; i < n; i++ )
  {
    v_normsq = v_normsq + v[i] * v[i];
  }
/*
  Compute A*H' = A*H
*/
  ah = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      ah[i+j*n] = a[i+j*n];
      for ( k = 0; k < n; k++ )
      {
        ah[i+j*n] = ah[i+j*n] - 2.0 * a[i+k*n] * v[k] * v[j] / v_normsq;
      }
    }
  }
  return ah;
}
/******************************************************************************/

double *r8mat_house_form ( int n, double v[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HOUSE_FORM constructs a Householder matrix from its compact form.

  Discussion: 							    

    An R8MAT is a doubly dimensioned array of double precision values, which
    may be stored as a vector in column-major order.

    H(v) = I - 2 * v * v' / ( v' * v )

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, double V[N], the vector defining the Householder matrix.

    Output, double R8MAT_HOUSE_FORM[N*N], the Householder matrix.
*/
{
  double beta;
  double *h;
  int i;
  int j;
/*
  Compute the L2 norm of V.
*/
  beta = 0.0;
  for ( i = 0; i < n; i++ )
  {
    beta = beta + v[i] * v[i];
  }
/*
  Form the matrix H.
*/
  h = r8mat_identity_new ( n );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      h[i+j*n] = h[i+j*n] - 2.0 * v[i] * v[j] / beta;
    }
  }
  return h;
}
/******************************************************************************/

double *r8mat_house_hxa ( int n, double a[], double v[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HOUSE_HXA computes H*A where H is a compact Householder matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The Householder matrix H(V) is defined by

      H(V) = I - 2 * v * v' / ( v' * v )

    This routine is not particularly efficient.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Input, double A[N*N], the matrix to be premultiplied.

    Input, double V[N], a vector defining a Householder matrix.

    Output, double R8MAT_HOUSE_HXA[N*N], the product H*A.
*/
{
  double *ha;
  int i;
  int j;
  int k;
  double v_normsq;

  v_normsq = 0.0;
  for ( i = 0; i < n; i++ )
  {
    v_normsq = v_normsq + v[i] * v[i];
  }
/*
  Compute A*H' = A*H
*/
  ha = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      ha[i+j*n] = a[i+j*n];
      for ( k = 0; k < n; k++ )
      {
        ha[i+j*n] = ha[i+j*n] - 2.0 * v[i] * v[k] * a[k+j*n] / v_normsq;
      }
    }
  }

  return ha;
}
/******************************************************************************/

double *r8mat_house_post ( int n, double a[], int row, int col )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HOUSE_POST computes a Householder post-multiplier matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    H(ROW,COL) has the property that the ROW-th column of
    A*H(ROW,COL) is zero from entry COL+1 to the end.

    In the most common case, where a QR factorization is being computed,
    ROW = COL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrices.

    Input, double A[N*N], the matrix whose Householder matrix
    is to be computed.

    Input, int ROW, COL, specify the location of the
    entry of the matrix A which is to be preserved.  The entries in
    the same row, but higher column, will be zeroed out if
    A is postmultiplied by H.

    Output, double R8MAT_HOUSE_POST[N*N], the Householder matrix.
*/
{
  double *h;
  int j;
  double *v;
  double *w;
/*
  Set up the vector V.
*/
  w = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < col-1; j++ )
  {
    w[j] = 0.0;
  }
  for ( j = col-1; j < n; j++ )
  {
    w[j] = a[row+j*n];
  }

  v = r8vec_house_column ( n, w, col );
/*
  Form the matrix H(V).
*/
  h = r8mat_house_form ( n, v );

  free ( v );
  free ( w );

  return h;
}
/******************************************************************************/

double *r8mat_house_pre ( int n, double a[], int row, int col )

/******************************************************************************/
/*
  Purpose:

    R8MAT_HOUSE_PRE computes a Householder pre-multiplier matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    H(ROW,COL) has the property that the COL-th column of
    H(ROW,COL)*A is zero from entry ROW+1 to the end.

    In the most common case, where a QR factorization is being computed,
    ROW = COL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrices.

    Input, double A[N*N], the matrix whose Householder matrix
    is to be computed.

    Input, int ROW, COL, specify the location of the
    entry of the matrix A which is to be preserved.  The entries in
    the same column, but higher rows, will be zeroed out if A is
    premultiplied by H.

    Output, double R8MAT_HOUSE_PRE[N*N], the Householder matrix.
*/
{
  double *h;
  int i;
  double *v;
  double *w;

  w = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Set up the vector V.
*/
  for ( i = 0; i < row-1; i++ )
  {
    w[i] = 0.0;
  }
  for ( i = row-1; i < n; i++ )
  {
    w[i] = a[i+col*n];
  }

  v = r8vec_house_column ( n, w, row );
/*
  Form the matrix H(V).
*/
  h = r8mat_house_form ( n, v );

  free ( v );
  free ( w );

  return h;
}
/******************************************************************************/

void r8mat_identity  ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_IDENTITY sets an R8MAT to the identity matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Output, double A[N*N], the N by N identity matrix.
*/
{
  int i;
  int j;
  int k;

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_identity_new ( int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_IDENTITY_NEW returns an identity matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Output, double R8MAT_IDENTITY_NEW[N*N], the N by N identity matrix.
*/
{
  double *a;
  int i;
  int j;
  int k;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return a;
}
/******************************************************************************/

int r8mat_in_01 ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_IN_01 is TRUE if the entries of an R8MAT are in the range [0,1].

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix.

    Output, int R8MAT_IN_01, is TRUE if every entry of A is
    between 0 and 1.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < 0.0 || 1.0 < a[i+j*m] )
      {
        return 0;
      }
    }
  }

  return 1;
}
/******************************************************************************/

double *r8mat_indicator_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_INDICATOR_NEW sets up an "indicator" R8MAT.

  Discussion:

    The value of each entry suggests its location, as in:

      11  12  13  14
      21  22  23  24
      31  32  33  34

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 January 2005

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Output, double R8MAT_INDICATOR_NEW[M*N], the table.
*/
{
  double *a;
  int fac;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  fac = i4_power ( 10, i4_log_10 ( n ) + 1 );

  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      a[i-1+(j-1)*m] = ( double ) ( fac * i + j );
    }
  }
  return a;
}
/******************************************************************************/

int r8mat_insignificant ( int m, int n, double r[], double s[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_INSIGNIFICANT determines if an R8MAT is insignificant.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the dimension of the matrices.

    Input, double R[M*N], the vector to be compared against.

    Input, double S[M*N], the vector to be compared.

    Output, int R8MAT_INSIGNIFICANT, is TRUE if S is insignificant
    compared to R.
*/
{
  int i;
  int j;
  double t;
  double tol;
  int value;

  value = 1;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      t = r[i+j*m] + s[i+j*m];
      tol = r8_epsilon ( ) * r8_abs ( r[i+j*m] );

      if ( tol < r8_abs ( r[i+j*m] - t ) )
      {
        value = 0;
        break;
      }
    }
  }
  return value;
}
/******************************************************************************/

double *r8mat_inverse_2d ( double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_INVERSE_2D inverts a 2 by 2 R8MAT using Cramer's rule.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double A[2*2], the matrix to be inverted.

    Output, double R8MAT_INVERSE_2D[2*2], the inverse of the matrix A.
*/
{
  double *b;
  double det;
/*
  Compute the determinant of A.
*/
  det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];
/*
  If the determinant is zero, bail out.
*/
  if ( det == 0.0 )
  {
    return NULL;
  }
/*
  Compute the entries of the inverse matrix using an explicit formula.
*/
  b = ( double * ) malloc ( 2 * 2 * sizeof ( double ) );

  b[0+0*2] = + a[1+1*2] / det;
  b[0+1*2] = - a[0+1*2] / det;
  b[1+0*2] = - a[1+0*2] / det;
  b[1+1*2] = + a[0+0*2] / det;

  return b;
}
/******************************************************************************/

double *r8mat_inverse_3d ( double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_INVERSE_3D inverts a 3 by 3 R8MAT using Cramer's rule.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    If the determinant is zero, A is singular, and does not have an
    inverse.  In that case, the output is set to NULL.

    If the determinant is nonzero, its value is an estimate
    of how nonsingular the matrix A is.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double A[3*3], the matrix to be inverted.

    Output, double R8MAT3_INVERSE[3*3], the inverse of the matrix A.
*/
{
  double *b;
  double det;
/*
  Compute the determinant of A.
*/
  det =
     a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
   + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
   + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );

  if ( det == 0.0 )
  {
    return NULL;
  }

  b = ( double * ) malloc ( 3 * 3 * sizeof ( double ) );

  b[0+0*3] =   ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] ) / det;
  b[0+1*3] = - ( a[0+1*3] * a[2+2*3] - a[0+2*3] * a[2+1*3] ) / det;
  b[0+2*3] =   ( a[0+1*3] * a[1+2*3] - a[0+2*3] * a[1+1*3] ) / det;

  b[1+0*3] = - ( a[1+0*3] * a[2+2*3] - a[1+2*3] * a[2+0*3] ) / det;
  b[1+1*3] =   ( a[0+0*3] * a[2+2*3] - a[0+2*3] * a[2+0*3] ) / det;
  b[1+2*3] = - ( a[0+0*3] * a[1+2*3] - a[0+2*3] * a[1+0*3] ) / det;

  b[2+0*3] =   ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] ) / det;
  b[2+1*3] = - ( a[0+0*3] * a[2+1*3] - a[0+1*3] * a[2+0*3] ) / det;
  b[2+2*3] =   ( a[0+0*3] * a[1+1*3] - a[0+1*3] * a[1+0*3] ) / det;

  return b;
}
/******************************************************************************/

double *r8mat_inverse_4d ( double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_INVERSE_4D inverts a 4 by 4 matrix using Cramer's rule.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 September 2003

  Author:

    John Burkardt

  Parameters:

    Input, double A[4][4], the matrix to be inverted.

    Output, double R8MAT_INVERSE_4D[4][4], the inverse of the matrix A.
*/
{
  double *b;
  double det;
/*
  Compute the determinant of A.
*/
  det = r8mat_det_4d ( a );
/*
  If the determinant is zero, bail out.
*/
  if ( det == 0.0 )
  {
    return NULL;
  }
/*
  Compute the entries of the inverse matrix using an explicit formula.
*/
  b = ( double * ) malloc ( 4 * 4 * sizeof ( double ) );

  b[0+0*4] =
    +(
    + a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[1+2*4] * ( a[2+3*4] * a[3+1*4] - a[2+1*4] * a[3+3*4] )
    + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    ) / det;

  b[1+0*4] =
    -(
    + a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[1+2*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
    ) / det;

  b[2+0*4] =
    +(
    + a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
    + a[1+1*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[3+0*4] =
    -(
    + a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    + a[1+1*4] * ( a[2+2*4] * a[3+0*4] - a[2+0*4] * a[3+2*4] )
    + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[0+1*4] =
    -(
    + a[0+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[2+3*4] * a[3+1*4] - a[2+1*4] * a[3+3*4] )
    + a[0+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    ) / det;

  b[1+1*4] =
    +(
    + a[0+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
    ) / det;

  b[2+1*4] =
    -(
    + a[0+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
    + a[0+1*4] * ( a[2+3*4] * a[3+0*4] - a[2+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[3+1*4] =
    +(
    + a[0+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
    + a[0+1*4] * ( a[2+2*4] * a[3+0*4] - a[2+0*4] * a[3+2*4] )
    + a[0+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] )
    ) / det;

  b[0+2*4] =
    +(
    + a[0+1*4] * ( a[1+2*4] * a[3+3*4] - a[1+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[3+1*4] - a[1+1*4] * a[3+3*4] )
    + a[0+3*4] * ( a[1+1*4] * a[3+2*4] - a[1+2*4] * a[3+1*4] )
    ) / det;

  b[1+2*4] =
    -(
    + a[0+0*4] * ( a[1+2*4] * a[3+3*4] - a[1+3*4] * a[3+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[3+0*4] - a[1+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[3+2*4] - a[1+2*4] * a[3+0*4] )
    ) / det;

  b[2+2*4] =
    +(
    + a[0+0*4] * ( a[1+1*4] * a[3+3*4] - a[1+3*4] * a[3+1*4] )
    + a[0+1*4] * ( a[1+3*4] * a[3+0*4] - a[1+0*4] * a[3+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[3+1*4] - a[1+1*4] * a[3+0*4] )
    ) / det;

  b[3+2*4] =
    -(
    + a[0+0*4] * ( a[1+1*4] * a[3+2*4] - a[1+2*4] * a[3+1*4] )
    + a[0+1*4] * ( a[1+2*4] * a[3+0*4] - a[1+0*4] * a[3+2*4] )
    + a[0+2*4] * ( a[1+0*4] * a[3+1*4] - a[1+1*4] * a[3+0*4] )
    ) / det;

  b[0+3*4] =
    -(
    + a[0+1*4] * ( a[1+2*4] * a[2+3*4] - a[1+3*4] * a[2+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[2+1*4] - a[1+1*4] * a[2+3*4] )
    + a[0+3*4] * ( a[1+1*4] * a[2+2*4] - a[1+2*4] * a[2+1*4] )
    ) / det;

  b[1+3*4] =
    +(
    + a[0+0*4] * ( a[1+2*4] * a[2+3*4] - a[1+3*4] * a[2+2*4] )
    + a[0+2*4] * ( a[1+3*4] * a[2+0*4] - a[1+0*4] * a[2+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[2+2*4] - a[1+2*4] * a[2+0*4] )
    ) / det;

  b[2+3*4] =
    -(
    + a[0+0*4] * ( a[1+1*4] * a[2+3*4] - a[1+3*4] * a[2+1*4] )
    + a[0+1*4] * ( a[1+3*4] * a[2+0*4] - a[1+0*4] * a[2+3*4] )
    + a[0+3*4] * ( a[1+0*4] * a[2+1*4] - a[1+1*4] * a[2+0*4] )
    ) / det;

  b[3+3*4] =
    +(
    + a[0+0*4] * ( a[1+1*4] * a[2+2*4] - a[1+2*4] * a[2+1*4] )
    + a[0+1*4] * ( a[1+2*4] * a[2+0*4] - a[1+0*4] * a[2+2*4] )
    + a[0+2*4] * ( a[1+0*4] * a[2+1*4] - a[1+1*4] * a[2+0*4] )
    ) / det;

  return b;
}
/******************************************************************************/

double r8mat_is_identity ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_IS_IDENTITY determines if an R8MAT is the identity.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The routine returns the Frobenius norm of A - I.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 July 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, double A[N*N], the matrix.

    Output, double R8MAT_IS_IDENTITY, the Frobenius norm
    of the difference matrix A - I, which would be exactly zero
    if A were the identity matrix.
*/
{
  double error_frobenius;
  int i;
  int j;
  double t;

  error_frobenius = 0.0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i == j )
      {
        t = a[i+j*n] - 1.0;
      }
      else
      {
        t = a[i+j*n];
      }
      error_frobenius = error_frobenius + t * t;
    }
  }
  error_frobenius = sqrt ( error_frobenius );

  return error_frobenius;
}
/******************************************************************************/

double r8mat_is_symmetric ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_IS_SYMMETRIC checks an R8MAT for symmetry.

  Discussion:

    An R8MAT is a matrix of double precision real values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the order of the matrix.

    Input, double A[M*N], the matrix.

    Output, double RMAT_IS_SYMMETRIC, measures the 
    Frobenius norm of ( A - A' ), which would be zero if the matrix
    were exactly symmetric.
*/
{
  int i;
  int j;
  double value;

  if ( m != n )
  {
    value = r8_huge ( );
    return value;
  }

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m] - a[j+i*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

double *r8mat_jac ( int m, int n, double eps,
  double *(*fx) ( int m, int n, double x[] ), double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_JAC estimates a dense jacobian matrix of the function FX.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    FPRIME(I,J) = d F(I) / d X(J).

    The jacobian is assumed to be dense, and the LINPACK/LAPACK
    double precision general matrix storage mode ("DGE") is used.

    Forward differences are used, requiring N+1 function evaluations.

    Values of EPS have typically been chosen between
    sqrt ( EPSMCH ) and sqrt ( sqrt ( EPSMCH ) ) where EPSMCH is the
    machine tolerance.

    If EPS is too small, then F(X+EPS) will be the same as
    F(X), and the jacobian will be full of zero entries.

    If EPS is too large, the finite difference estimate will
    be inaccurate.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of functions.

    Input, int N, the number of variables.

    Input, double EPS, a tolerance to be used for shifting the
    X values during the finite differencing.  No single value
    of EPS will be reliable for all vectors X and functions FX.

    Input, double *(*FX) ( int m, int n, double x[] ), the name of
    the user written routine which evaluates the M-dimensional
    function at a given N-dimensional point X.

    Input, double X[N], the point where the jacobian
    is to be estimated.

    Output, double R8MAT_JAC[M*N], the estimated jacobian matrix.
*/
{
  double del;
  double *fprime;
  int i;
  int j;
  double xsave;
  double *work1;
  double *work2;

  fprime = ( double * ) malloc ( m * n * sizeof ( double ) );
/*
  Evaluate the function at the base point, X.
*/
  work2 = fx ( m, n, x );
/*
  Now, one by one, vary each component J of the base point X, and
  estimate DF(I)/DX(J) = ( F(X+) - F(X) )/ DEL.
*/
  for ( j = 0; j < n; j++ )
  {
    xsave = x[j];
    del = eps * ( 1.0 + r8_abs ( x[j] ) );
    x[j] = x[j] + del;
    work1 = fx ( m, n, x );
    x[j] = xsave;
    for ( i = 0; i < m; i++ )
    {
      fprime[i+j*m] = ( work1[i] - work2[i] ) / del;
    }
    free ( work1 );
  }
  free ( work2 );

  return fprime;
}
/******************************************************************************/

double *r8mat_l_inverse ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_L_INVERSE inverts a lower triangular R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    A lower triangular matrix is a matrix whose only nonzero entries
    occur on or below the diagonal.

    The inverse of a lower triangular matrix is a lower triangular matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, number of rows and columns in the matrix.

    Input, double A[N*N], the lower triangular matrix.

    Output, double R8MAT_L_INVERSE[N*N], the inverse matrix.
*/
{
  double *b;
  int i;
  int j;
  int k;
  double temp;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i < j )
      {
        b[i+j*n] = 0.0;
      }
      else if ( j == i )
      {
        b[i+j*n] = 1.0 / a[i+j*n];
      }
      else
      {
        temp = 0.0;
        for ( k = 0; k < i; k++ )
        {
          temp = temp + a[i+k*n] * b[k+j*n];
        }
        b[i+j*n] = -temp / a[i+i*n];
      }
    }
  }

  return b;
}
/******************************************************************************/

void r8mat_l_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_L_PRINT prints a lower triangular R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Example:

    M = 5, N = 5
    A = (/ 11, 21, 31, 41, 51, 22, 32, 42, 52, 33, 43, 53, 44, 54, 55 /)

    11
    21 22
    31 32 33
    41 42 43 44
    51 52 53 54 55

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[*], the M by N matrix.  Only the lower
    triangular elements are stored, in column major order.

    Input, char *TITLE, a title.
*/
{
  int i;
  int indx[10];
  int j;
  int jhi;
  int jlo;
  int jmax;
  int nn;
  int size;

  printf ( "\n" );
  printf ( "%s\n", title );

  jmax = i4_min ( n, m );

  if ( m <= n )
  {
    size = ( m * ( m + 1 ) ) / 2;
  }
  else if ( n < m )
  {
    size = ( n * ( n + 1 ) ) / 2 + ( m - n ) * n;
  }

  if ( r8vec_is_int ( size, a ) )
  {
    nn = 10;
    for ( jlo = 1; jlo <= jmax; jlo = jlo + nn )
    {
      jhi = i4_min ( jlo + nn - 1, i4_min ( m, jmax ) );
      printf ( "\n" );
      printf ( "  Col   " );
      for ( j = jlo; j <= jhi; j++ )
      {
        printf ( "%6d", j );
      }
      printf ( "\n" );
      printf ( "  Row  \n" );
      for ( i = jlo; i <= m; i++ )
      {
        jhi = i4_min ( jlo + nn - 1, i4_min ( i, jmax ) );
        for ( j = jlo; j <= jhi; j++ )
        {
          indx[j-jlo] = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2;
        }
        printf ( "  %6d", i );
        for ( j = 0; j <= jhi-jlo; j++ )
        {
          printf ( "%6g", a[indx[j]-1] );
        }
        printf ( "\n" );
      }
    }
  }
  else if ( r8vec_amax ( size, a ) < 1000000.0 )
  {
    nn = 5;
    for ( jlo = 1; jlo <= jmax; jlo = jlo + nn );
    {
      jhi = i4_min ( jlo + nn - 1, i4_min ( m - 1, jmax ) );
      printf ( "\n" );
      printf ( "  Col " );
      for ( j = jlo; j <= jhi; j++ )
      {
        printf ( "%6d", j );
      }
      printf ( "\n" );
      printf ( "  Row  \n" );
      for ( i = jlo; i <= m; i++ )
      {
        jhi = i4_min ( jlo + nn - 1, i4_min ( i, jmax ) );
        for ( j = jlo; j <= jhi; j++ )
        {
          indx[j-jlo] = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2;
        }
        printf ( "  %6d", i );
        for ( j = 0; j <= jhi-jlo; j++ )
        {
          printf ( "%14g", a[indx[j]-1] );
        }
        printf ( "\n" );
      }
    }
  }
  else
  {
    nn = 5;

    for ( jlo = 1; jlo <= jmax; jlo = jlo + nn )
    {
      jhi = i4_min ( jlo + nn - 1, i4_min ( m - 1, jmax ) );
      printf ( "\n" );
      printf ( "  Col " );
      for ( j = jlo; j <= jhi; j++ )
      {
        printf ( "%7d       ", j );
      }
      printf ( "\n" );
      printf ( "  Row \n");
      for ( i = jlo; i <= m; i++ )
      {
        jhi = i4_min ( jlo + nn - 1, i4_min ( i, jmax ) );
        for ( j = jlo; j <= jhi; j++ )
        {
          indx[j-jlo] = ( j - 1 ) * m + i - ( j * ( j - 1 ) ) / 2;
        }
        printf ( "  %6d", i );
        for ( j = 0; j <= jhi-jlo; j++ )
        {
          printf ( "%14g", a[indx[j]-1] );
        }
      }
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_l_solve ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_L_SOLVE solves a lower triangular linear system.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of
    the matrix A.

    Input, double A[N*N], the N by N lower triangular matrix.

    Input, double B[N], the right hand side of the linear system.

    Output, double R8MAT_L_SOLVE[N], the solution of the linear system.
*/
{
  double dot;
  int i;
  int j;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );
/*
  Solve L * x = b.
*/
  for ( i = 0; i < n; i++ )
  {
    dot = 0.0;
    for ( j = 0; j < i; j++ )
    {
      dot = dot + a[i+j*n] * x[j];
    }
    x[i] = ( b[i] - dot ) / a[i+i*n];
  }

  return x;
}
/******************************************************************************/

double *r8mat_l1_inverse ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_L1_INVERSE inverts a unit lower triangular R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    A unit lower triangular matrix is a matrix with only 1's on the main
    diagonal, and only 0's above the main diagonal.

    The inverse of a unit lower triangular matrix is also
    a unit lower triangular matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, number of rows and columns in the matrix.

    Input, double A[N*N], the unit lower triangular matrix.

    Output, double R8MAT_L1_INVERSE[N*N], the inverse matrix.
*/
{
  double *b;
  int i;
  int j;
  int k;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( i < j )
      {
        b[i+j*n] = 0.0;
      }
      else if ( j == i )
      {
        b[i+j*n] = 1.0;
      }
      else
      {
        b[i+j*n] = 0.0;
        for ( k = 0; k < i; k++ )
        {
          b[i+j*n] = b[i+j*n] - a[i+k*n] * b[k+j*n];
        }
      }
    }
  }

  return b;
}
/******************************************************************************/

double *r8mat_lt_solve ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_LT_SOLVE solves a transposed lower triangular linear system.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Given the lower triangular matrix A, the linear system to be solved is:

      A' * x = b

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix A.

    Input, double A[N*N], the N by N lower triangular matrix.

    Input, double B[N], the right hand side of the linear system.

    Output, double R8MAT_LT_SOLVE[N], the solution of the linear system.
*/
{
  int i;
  int j;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = n-1; 0 <= j; j-- )
  {
    x[j] = b[j];
    for ( i = j+1; i < n; i++ )
    {
      x[j] = x[j] - x[i] * a[i+j*n];
    }
    x[j] = x[j] / a[j+j*n];
  }

  return x;
}
/******************************************************************************/

void r8mat_lu ( int m, int n, double a[], double l[], double p[], double u[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_LU computes the LU factorization of a rectangular R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The routine is given an M by N matrix A, and produces

      L, an M by M unit lower triangular matrix,
      U, an M by N upper triangular matrix, and
      P, an M by M permutation matrix P,

    so that

      A = P' * L * U.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix to be factored.

    Output, double L[M*M], the M by M unit lower triangular factor.

    Output, double P[M*M], the M by M permutation matrix.

    Output, double U[M*N], the M by N upper triangular factor.
*/
{
  int i;
  int ipiv;
  int j;
  int k;
  double pivot;
/*
  Initialize:

    U:=A
    L:=Identity
    P:=Identity
*/
  r8mat_copy ( m, n, a, u );

  r8mat_zero ( m, m, l );
  r8mat_zero ( m, m, p );
  for ( i = 0; i < m; i++ )
  {
    l[i+i*m] = 1.0;
    p[i+i*m] = 1.0;
  }
/*
  On step J, find the pivot row, IPIV, and the pivot value PIVOT.
*/
  for ( j = 0; j < i4_min ( m - 1, n ); j++ )
  {
    pivot = 0.0;
    ipiv = -1;

    for ( i = j; i < m; i++ )
    {
      if ( pivot < r8_abs ( u[i+j*m] ) )
      {
        pivot = r8_abs ( u[i+j*m] );
        ipiv = i;
      }
    }
/*
  Unless IPIV is zero, swap rows J and IPIV.
*/
    if ( ipiv != -1 )
    {
      r8row_swap ( m, n, u, j+1, ipiv+1 );

      r8row_swap ( m, m, l, j+1, ipiv+1 );

      r8row_swap ( m, m, p, j+1, ipiv+1 );
/*
  Zero out the entries in column J, from row J+1 to M.
*/
      for ( i = j+1; i < m; i++ )
      {
        if ( u[i+j*m] != 0.0 )
        {
          l[i+j*m] = u[i+j*m] / u[j+j*m];

          u[i+j*m] = 0.0;

          for ( k = j+1; k < n; k++ )
          {
            u[i+k*m] = u[i+k*m] - l[i+j*m] * u[j+k*m];
          }
        }
      }
    }
  }

  return;
}
/******************************************************************************/

double r8mat_max ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MAX returns the maximum entry of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Output, double R8MAT_MAX, the maximum entry of A.
*/
{
  int i;
  int j;
  double value;

  value = a[0+0*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( value < a[i+j*m] )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
/******************************************************************************/

void r8mat_max_index ( int m, int n, double a[], int *i_max, int *j_max )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MAX_INDEX returns the location of the maximum entry of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Output, int *I_MAX, *J_MAX, the indices of the maximum entry of A.
*/
{
  int i;
  int i2;
  int j;
  int j2;

  i2 = -1;
  j2 = -1;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i2 == -1 && j2 == -1 )
      {
        i2 = i;
        j2 = j;
      }
      else if ( a[i2+j2*m] < a[i+j*m] )
      {
        i2 = i;
        j2 = j;
      }
    }
  }

  *i_max = i2 + 1;
  *j_max = j2 + 1;

  return;
}
/******************************************************************************/

double r8mat_maxcol_minrow ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MAXCOL_MINROW gets the maximum column minimum row of an M by N matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    R8MAT_MAXCOL_MINROW = max ( 1 <= I <= N ) ( min ( 1 <= J <= M ) A(I,J) )

    For a given matrix, R8MAT_MAXCOL_MINROW <= R8MAT_MINROW_MAXCOL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix.

    Output, double R8MAT_MAXCOL_MINROW, the maximum column
    minimum row entry of A.
*/
{
  int i;
  int j;
  double minrow;
  double value;

  value = - r8_huge ( );

  for ( i = 0; i < m; i++ )
  {
    minrow = r8_huge ( );

    for ( j = 0; j < n; j++ )
    {
      minrow = r8_min ( minrow, a[i+j*m] );
    }
    value = r8_max ( value, minrow );
  }

  return value;
}
/******************************************************************************/

double r8mat_maxrow_mincol ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MAXROW_MINCOL gets the maximum row minimum column of an M by N matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    R8MAT_MAXROW_MINCOL = max ( 1 <= J <= N ) ( min ( 1 <= I <= M ) A(I,J) )

    For a given matrix, R8MAT_MAXROW_MINCOL <= R8MAT_MINCOL_MAXROW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix.

    Output, double R8MAT_MAXROW_MINCOL, the maximum row
    minimum column entry of A.
*/
{
  int i;
  int j;
  double mincol;
  double value;

  value = - r8_huge ( );

  for ( j = 0; j < n; j++ )
  {
    mincol = r8_huge ( );
    for ( i = 0; i < m; i++ )
    {
      mincol = r8_min ( mincol, a[i+j*m] );
    }
    value = r8_max ( value, mincol );
  }
  return value;
}
/******************************************************************************/

double r8mat_min ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MIN returns the minimum entry of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 May 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Output, double R8MAT_MIN, the minimum entry of A.
*/
{
  int i;
  int j;
  double value;

  value = a[0+0*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < value )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
/******************************************************************************/

void r8mat_min_index ( int m, int n, double a[], int *i_min, int *j_min )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MIN_INDEX returns the location of the minimum entry of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Output, int *I_MIN, *J_MIN, the indices of the minimum entry of A.
*/
{
  int i;
  int i2;
  int j;
  int j2;

  i2 = -1;
  j2 = -1;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( i2 == -1 && j2 == -1 )
      {
        i2 = i;
        j2 = j;
      }
      else if ( a[i+j*m] < a[i2+j2*m] )
      {
        i2 = i;
        j2 = j;
      }
    }
  }

  *i_min = i2 + 1;
  *j_min = j2 + 1;

  return;
}
/******************************************************************************/

double r8mat_mincol_maxrow ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MINCOL_MAXROW gets the minimum column maximum row of an M by N matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    R8MAT_MINCOL_MAXROW = min ( 1 <= I <= N ) ( max ( 1 <= J <= M ) A(I,J) )

    For a given matrix, R8MAT_MAXROW_MINCOL <= R8MAT_MINCOL_MAXROW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A(M,N), the matrix.

    Output, double R8MAT_MINCOL_MAXROW, the minimum column
    maximum row entry of A.
*/
{
  int i;
  int j;
  double maxrow;
  double value;

  value = r8_huge ( );

  for ( i = 0; i < m; i++ )
  {
    maxrow = - r8_huge ( );
    for ( j = 0; j < n; j++ )
    {
      maxrow = r8_max ( maxrow, a[i+j*m] );
    }
    value = r8_min ( value, maxrow );
  }

  return value;
}
/******************************************************************************/

double r8mat_minrow_maxcol ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MINROW_MAXCOL gets the minimum row maximum column of an M by N matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    R8MAT_MINROW_MAXCOL = min ( 1 <= J <= N ) ( max ( 1 <= I <= M ) A(I,J) )

    For a given matrix, R8MAT_MAXCOL_MINROW <= R8MAT_MINROW_MAXCOL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix.

    Output, double R8MAT_MINROW_MAXCOL, the minimum row
    maximum column entry of A.
*/
{
  int i;
  int j;
  double maxcol;
  double value;;

  value = r8_huge ( );

  for ( j = 0; j < n; j++ )
  {
    maxcol = - r8_huge ( );
    for ( i = 0; i < m; i++ )
    {
      maxcol = r8_max ( maxcol, a[i+j*m] );
    }
    value = r8_min ( value, maxcol );
  }

  return value;
}
/******************************************************************************/

double r8mat_minvm ( int n1, int n2, double a[], double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MINVM returns C = inverse(A) * B for R8MAT's.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, the order of the matrices.

    Input, double A[N1*N1], B[N1*N2], the matrices.

    Output, double C[N1*N2], the result, C = inverse(A) * B.
*/
{
  double *alu;
  double *d;

  alu = r8mat_copy_new ( n1, n1, a );

  d = r8mat_fss_new ( n1, alu, n2, b );
 
  r8mat_copy ( n1, n2, d, c );

  free ( alu );
  free ( d );

  return 0.0;
}
/******************************************************************************/

double *r8mat_minvm_new ( int n1, int n2, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MINVM_NEW returns inverse(A) * B for R8MAT's.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, the order of the matrices.

    Input, double A[N1*N1], B[N1*N2], the matrices.

    Output, double R8MAT_MINVM_NEW[N1*N2], the result, C = inverse(A) * B.
*/
{
  double *alu;
  double *c;

  alu = r8mat_copy_new ( n1, n1, a );

  c = r8mat_fss_new ( n1, alu, n2, b );
 
  free ( alu );

  return c;
}
/******************************************************************************/

void r8mat_mm ( int n1, int n2, int n3, double a[], double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MM multiplies two matrices.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.

    Output, double C[N1*N3], the product matrix C = A * B.
*/
{
  double *d;
  int i;
  int j;
  int k;

  d = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      d[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        d[i+j*n1] = d[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

 for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = d[i+j*n1];
    }
  }

  free ( d );

  return;
}
/******************************************************************************/

double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MM_NEW multiplies two matrices.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.

    Output, double R8MAT_MM[N1*N3], the product matrix C = A * B.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
/******************************************************************************/

double *r8mat_mmt_new ( int n1, int n2, int n3, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MMT_NEW computes C = A * B'.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    13 November 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N3*N2], the matrices to multiply.

    Output, double R8MAT_MMT[N1*N3], the product matrix C = A * B'.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[j+k*n3];
      }
    }
  }

  return c;
}
/******************************************************************************/

double *r8mat_mtm_new ( int n1, int n2, int n3, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MTM_NEW computes C = A' * B.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N2*N1], double B[N2*N3], the matrices to multiply.

    Output, double R8MAT_MM[N1*N3], the product matrix C = A' * B.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[k+i*n2] * b[k+j*n2];
      }
    }
  }

  return c;
}
/******************************************************************************/

void r8mat_mtv ( int m, int n, double a[], double x[], double atx[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MTV multiplies a transposed matrix times a vector.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as an argument.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix.

    Input, double A[M,N], the M by N matrix.

    Input, double X[M], the vector to be multiplied by A.

    Output, double ATX[N], the product A'*X.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    atx[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      atx[j] = atx[j] + a[i+j*m] * x[i];
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_mtv_new ( int m, int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MTV_NEW multiplies a transposed matrix times a vector.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix.

    Input, double A[M,N], the M by N matrix.

    Input, double X[M], the vector to be multiplied by A.

    Output, double R8MAT_MTV_NEW[N], the product A'*X.
*/
{
  int i;
  int j;
  double *y;

  y = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    y[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      y[j] = y[j] + a[i+j*m] * x[i];
    }
  }

  return y;
}
/******************************************************************************/

void r8mat_mv ( int m, int n, double a[], double x[], double ax[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MV multiplies a matrix times a vector.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as an argument.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix.

    Input, double A[M,N], the M by N matrix.

    Input, double X[N], the vector to be multiplied by A.

    Output, double AX[M], the product A*X.
*/
{
  int i;
  int j;

  for ( i = 0; i < m; i++ )
  {
    ax[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      ax[i] = ax[i] + a[i+j*m] * x[j];
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_mv_new ( int m, int n, double a[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MV_NEW multiplies a matrix times a vector.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as the function value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 April 2007

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix.

    Input, double A[M,N], the M by N matrix.

    Input, double X[N], the vector to be multiplied by A.

    Output, double R8MAT_MV[M], the product A*X.
*/
{
  int i;
  int j;
  double *y;

  y = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
/******************************************************************************/

void r8mat_mxm ( int n1, int n2, int n3, double a[], double b[], double c[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MXM multiplies two matrices.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    For this routine, the result is returned as an argument.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.

    Output, double C[N1*N3], the product matrix C = A * B.
*/
{
  int i;
  int j;
  int k;

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_mxm_new ( int n1, int n2, int n3, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_MXM_NEW multiplies two matrices.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N1, N2, N3, the order of the matrices.

    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.

    Output, double R8MAT_MXM_NEW[N1*N3], the product matrix C = A * B.
*/
{
  double *c;
  int i;
  int j;
  int k;

  c = ( double * ) malloc ( n1 * n3 * sizeof ( double ) );

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
/******************************************************************************/

double **r8mat_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NEW sets up an R8MAT of the desired dimensions.

  Discussion:

    A declaration of the form
      double **a;
    is necesary.  Then an assignment of the form:
      a = r8mat_new ( m, n );
    allows the user to assign entries to the matrix using typical
    2D array notation:
      a[2][3] = 17.0;
      y = a[1][0];
    and so on.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 September 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, double **R8MAT_NEW, the array.
*/
{
  double **a;
  int i;

  a = ( double ** ) malloc ( m * sizeof ( double * ) );

  for ( i = 0; i < m; i++ )
  {
    a[i] = ( double * ) malloc ( n * sizeof ( double ) );
  }
  return a;
}
/******************************************************************************/

void r8mat_nint ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NINT rounds the entries of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of A.

    Input/output, double A[M*N], the matrix to be NINT'ed.
*/
{
  int i;
  int j;
  int s;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( a[i+j*m] < 0.0 )
      {
        s = -1;
      }
      else
      {
        s = 1;
      }
      a[i+j*m] = s * ( int ) ( r8_abs ( a[i+j*m] ) + 0.5 );
    }
  }

  return;
}
/******************************************************************************/

double r8mat_norm_eis ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_EIS returns the EISPACK norm of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The EISPACK norm is defined as:

      R8MAT_NORM_EIS =
        sum ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix whose EISPACK norm is desired.

    Output, double R8MAT_NORM_EIS, the EISPACK norm of A.
*/
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + r8_abs ( a[i+j*m] );
    }
  }

  return value;
}
/******************************************************************************/

double r8mat_norm_fro ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The Frobenius norm is defined as

      R8MAT_NORM_FRO = sqrt (
        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
    The matrix Frobenius norm is not derived from a vector norm, but
    is compatible with the vector L2 norm, so that:

      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix whose Frobenius
    norm is desired.

    Output, double R8MAT_NORM_FRO, the Frobenius norm of A.
*/
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

double r8mat_norm_fro_affine ( int m, int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_FRO_AFFINE returns the Frobenius norm of an R8MAT difference.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The Frobenius norm is defined as

      R8MAT_NORM_FRO = sqrt (
        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
    The matrix Frobenius norm is not derived from a vector norm, but
    is compatible with the vector L2 norm, so that:

      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows.

    Input, int N, the number of columns.

    Input, double A1[M*N], the matrices for whose difference the Frobenius
    norm is desired.

    Output, double R8MAT_NORM_FRO_AFFINE, the Frobenius norm of A1 - A2.
*/
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a1[i+j*m] - a2[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

double r8mat_norm_l1 ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_L1 returns the matrix L1 norm of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The matrix L1 norm is defined as:

      R8MAT_NORM_L1 = max ( 1 <= J <= N )
        sum ( 1 <= I <= M ) abs ( A(I,J) ).

    The matrix L1 norm is derived from the vector L1 norm, and
    satisifies:

      r8vec_norm_l1 ( A * x ) <= r8mat_norm_l1 ( A ) * r8vec_norm_l1 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A(M,N), the matrix whose L1 norm is desired.

    Output, double R8MAT_NORM_L1, the L1 norm of A.
*/
{
  double col_sum;
  int i;
  int j;
  double value;

  value = 0.0;

  for ( j = 0; j < n; j++ )
  {
    col_sum = 0.0;
    for ( i = 0; i < m; i++ )
    {
      col_sum = col_sum + r8_abs ( a[i+j*m] );
    }
    value = r8_max ( value, col_sum );
  }
  return value;
}
/******************************************************************************/

double r8mat_norm_l2 ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_L2 returns the matrix L2 norm of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The matrix L2 norm is defined as:

      R8MAT_NORM_L2 = sqrt ( max ( 1 <= I <= M ) LAMBDA(I) )

    where LAMBDA contains the eigenvalues of A * A'.

    The matrix L2 norm is derived from the vector L2 norm, and
    satisifies:

      r8vec_norm_l2 ( A * x ) <= r8mat_norm_l2 ( A ) * r8vec_norm_l2 ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A(M,N), the matrix whose L2 norm is desired.

    Output, double R8MAT_NORM_L2, the L2 norm of A.
*/
{
  double *at;
  double *b;
  double *diag;
  double value;

  at = r8mat_transpose_new ( m, n, a );
/*
  Compute B = A * A'.
*/
  b = r8mat_mm_new ( m, n, m, a, at );
/*
  Diagonalize B.
*/
  r8mat_symm_jacobi ( m, b );
/*
  Find the maximum eigenvalue, and take its square root.
*/
  diag = r8mat_diag_get_vector ( m, b );

  value = sqrt ( r8vec_max ( m, diag ) );

  free ( at );
  free ( b );
  free ( diag );

  return value;
}
/******************************************************************************/

double r8mat_norm_li ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORM_LI returns the matrix L-oo norm of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The matrix L-oo norm is defined as:

      R8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).

    The matrix L-oo norm is derived from the vector L-oo norm,
    and satisifies:

      r8vec_norm_li ( A * x ) <= r8mat_norm_li ( A ) * r8vec_norm_li ( x ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix whose L-oo
    norm is desired.

    Output, double R8MAT_NORM_LI, the L-oo norm of A.
*/
{
  int i;
  int j;
  double row_sum;
  double value;

  value = 0.0;

  for ( i = 0; i < m; i++ )
  {
    row_sum = 0.0;
    for ( j = 0; j < n; j++ )
    {
      row_sum = row_sum + r8_abs ( a[i+j*m] );
    }
    value = r8_max ( value, row_sum );
  }
  return value;
}
/******************************************************************************/

double *r8mat_normal_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NORMAL_01_NEW returns a unit pseudonormal R8MAT.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns in the array.

    Input/output, int *SEED, the "seed" value, which should NOT be 0.
    On output, SEED has been updated.

    Output, double R8MAT_NORMAL_01_NEW[M*N], the array of pseudonormal values.
*/
{
  double *r;

  r = r8vec_normal_01_new ( m * n, seed );

  return r;
}
/******************************************************************************/

double *r8mat_nullspace ( int m, int n, double a[], int nullspace_size )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NULLSPACE computes the nullspace of a matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Let A be an MxN matrix.

    If X is an N-vector, and A*X = 0, then X is a null vector of A.

    The set of all null vectors of A is called the nullspace of A.

    The 0 vector is always in the null space.

    If the 0 vector is the only vector in the nullspace of A, then A
    is said to have maximum column rank.  (Because A*X=0 can be regarded
    as a linear combination of the columns of A).  In particular, if A
    is square, and has maximum column rank, it is nonsingular.

    The dimension of the nullspace is the number of linearly independent
    vectors that span the nullspace.  If A has maximum column rank,
    its nullspace has dimension 0.

    This routine uses the reduced row echelon form of A to determine
    a set of NULLSPACE_SIZE independent null vectors.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 October 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of
    the matrix A.

    Input, double A[M*N], the matrix to be analyzed.

    Input, int NULLSPACE_SIZE, the size of the nullspace.

    Output, double R8MAT_NULLSPACE[N*NULLSPACE_SIZE], vectors that
    span the nullspace.
*/
{
  int *col;
  int i;
  int i2;
  int j;
  int j2;
  double *nullspace;
  int *row;
  double *rref;
/*
  Make a copy of A.
*/
  rref = r8mat_copy_new ( m, n, a );
/*
  Get the reduced row echelon form of A.
*/
  r8mat_rref ( m, n, rref );
/*
  Note in ROW the columns of the leading nonzeros.
  COL(J) = +J if there is a leading 1 in that column, and -J otherwise.
*/
  row = ( int * ) malloc ( m * sizeof ( int ) );

  for ( i = 0; i < m; i++ )
  {
    row[i] = 0;
  }

  col = ( int * ) malloc ( n * sizeof ( int ) );

  for ( j = 0; j < n; j++ )
  {
    col[j] = - ( j + 1 );
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( rref[i+j*m] == 1.0 )
      {
        row[i] = ( j + 1 );
        col[j] = ( j + 1 );
        break;
      }
    }
  }

  nullspace = r8mat_zero_new ( n, nullspace_size );

  j2 = 0;
/*
  If column J does not contain a leading 1, then it contains
  information about a null vector.
*/
  for ( j = 0; j < n; j++ )
  {
    if ( col[j] < 0 )
    {
      for ( i = 0; i < m; i++ )
      {
        if ( rref[i+j*m] != 0.0 )
        {
          i2 = row[i] - 1;
          nullspace[i2+j2*n] = - rref[i+j*m];
        }
      }
      nullspace[j+j2*n] = 1.0;
      j2 = j2 + 1;
    }
  }
  free ( col );
  free ( row );
  free ( rref );

  return nullspace;
}
/******************************************************************************/

int r8mat_nullspace_size ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_NULLSPACE_SIZE computes the size of the nullspace of a matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Let A be an MxN matrix.

    If X is an N-vector, and A*X = 0, then X is a null vector of A.

    The set of all null vectors of A is called the nullspace of A.

    The 0 vector is always in the null space.

    If the 0 vector is the only vector in the nullspace of A, then A
    is said to have maximum column rank.  (Because A*X=0 can be regarded
    as a linear combination of the columns of A).  In particular, if A
    is square, and has maximum column rank, it is nonsingular.

    The dimension of the nullspace is the number of linearly independent
    vectors that span the nullspace.  If A has maximum column rank,
    its nullspace has dimension 0.

    This routine ESTIMATES the dimension of the nullspace.  Cases of
    singularity that depend on exact arithmetic will probably be missed.

    The nullspace will be estimated by counting the leading 1's in the
    reduced row echelon form of A, and subtracting this from N.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of
    the matrix A.

    Input, double A[M*N], the matrix to be analyzed.

    Output, int R8MAT_NULLSPACE_SIZE, the estimated size
    of the nullspace.
*/
{
  int i;
  int j;
  int leading;
  int nullspace_size;
  double *rref;
/*
  Make a copy of A.
*/
  rref = r8mat_copy_new ( m, n, a );
/*
  Get the reduced row echelon form of A.
*/
  r8mat_rref ( m, n, rref );
/*
  Count the leading 1's in A.
*/
  leading = 0;
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( rref[i+j*m] == 1.0 )
      {
        leading = leading + 1;
        break;
      }
    }
  }
  nullspace_size = n - leading;

  free ( rref );

  return nullspace_size;
}
/******************************************************************************/

void r8mat_orth_uniform ( int n, int *seed, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ORTH_UNIFORM returns a random orthogonal matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The inverse of A is equal to A'.

    A * A'  = A' * A = I.

    Columns and rows of A have unit Euclidean norm.

    Distinct pairs of columns of A are orthogonal.

    Distinct pairs of rows of A are orthogonal.

    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.

    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.

    The determinant of A is +1 or -1.

    All the eigenvalues of A have modulus 1.

    All singular values of A are 1.

    All entries of A are between -1 and 1.

  Discussion:

    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
    National Academy of Sciences of Belarus, for convincingly
    pointing out the severe deficiencies of an earlier version of
    this routine.

    Essentially, the computation involves saving the Q factor of the
    QR factorization of a matrix whose entries are normally distributed.
    However, it is only necessary to generate this matrix a column at
    a time, since it can be shown that when it comes time to annihilate
    the subdiagonal elements of column K, these (transformed) elements of
    column K are still normally distributed random values.  Hence, there
    is no need to generate them at the beginning of the process and
    transform them K-1 times.

    For computational efficiency, the individual Householder transformations
    could be saved, as recommended in the reference, instead of being
    accumulated into an explicit matrix format.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Reference:

    Pete Stewart,
    Efficient Generation of Random Orthogonal Matrices With an Application
    to Condition Estimators,
    SIAM Journal on Numerical Analysis,
    Volume 17, Number 3, June 1980, pages 403-409.

  Parameters:

    Input, int N, the order of A.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double A[N*N], the orthogonal matrix.
*/
{
  double *a2;
  int i;
  int j;
  double *v;
  double *x;
/*
  Start with A = the identity matrix.
*/
  r8mat_identity ( n, a );
/*
  Now behave as though we were computing the QR factorization of
  some other random matrix.  Generate the N elements of the first column,
  compute the Householder matrix H1 that annihilates the subdiagonal elements,
  and set A := A * H1' = A * H.

  On the second step, generate the lower N-1 elements of the second column,
  compute the Householder matrix H2 that annihilates them,
  and set A := A * H2' = A * H2 = H1 * H2.

  On the N-1 step, generate the lower 2 elements of column N-1,
  compute the Householder matrix HN-1 that annihilates them, and
  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
  This is our random orthogonal matrix.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 1; j < n; j++ )
  {
/*
  Set the vector that represents the J-th column to be annihilated.
*/
    for ( i = 1; i < j; i++ )
    {
      x[i-1] = 0.0;
    }
    for ( i = j; i <= n; i++ )
    {
      x[i-1] = r8_normal_01 ( seed );
    }
/*
  Compute the vector V that defines a Householder transformation matrix
  H(V) that annihilates the subdiagonal elements of X.
*/
    v = r8vec_house_column ( n, x, j );
/*
  Postmultiply the matrix A by H'(V) = H(V).
*/
    a2 = r8mat_house_axh_new ( n, a, v );

    free ( v );

    r8mat_copy ( n, n, a2, a );

    free ( a2 );
  }

  free ( x );

  return;
}
/******************************************************************************/

double *r8mat_orth_uniform_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ORTH_UNIFORM_NEW returns a random orthogonal matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The inverse of A is equal to A'.

    A * A'  = A' * A = I.

    Columns and rows of A have unit Euclidean norm.

    Distinct pairs of columns of A are orthogonal.

    Distinct pairs of rows of A are orthogonal.

    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.

    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.

    The determinant of A is +1 or -1.

    All the eigenvalues of A have modulus 1.

    All singular values of A are 1.

    All entries of A are between -1 and 1.

  Discussion:

    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
    National Academy of Sciences of Belarus, for convincingly
    pointing out the severe deficiencies of an earlier version of
    this routine.

    Essentially, the computation involves saving the Q factor of the
    QR factorization of a matrix whose entries are normally distributed.
    However, it is only necessary to generate this matrix a column at
    a time, since it can be shown that when it comes time to annihilate
    the subdiagonal elements of column K, these (transformed) elements of
    column K are still normally distributed random values.  Hence, there
    is no need to generate them at the beginning of the process and
    transform them K-1 times.

    For computational efficiency, the individual Householder transformations
    could be saved, as recommended in the reference, instead of being
    accumulated into an explicit matrix format.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Reference:

    Pete Stewart,
    Efficient Generation of Random Orthogonal Matrices With an Application
    to Condition Estimators,
    SIAM Journal on Numerical Analysis,
    Volume 17, Number 3, June 1980, pages 403-409.

  Parameters:

    Input, int N, the order of A.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8MAT_ORTH_UNIFORM_NEW[N*N], the orthogonal matrix.
*/
{
  double *a;
  double *a2;
  int i;
  int j;
  double *v;
  double *x;
/*
  Start with A = the identity matrix.
*/
  a = r8mat_identity_new ( n );
/*
  Now behave as though we were computing the QR factorization of
  some other random matrix.  Generate the N elements of the first column,
  compute the Householder matrix H1 that annihilates the subdiagonal elements,
  and set A := A * H1' = A * H.

  On the second step, generate the lower N-1 elements of the second column,
  compute the Householder matrix H2 that annihilates them,
  and set A := A * H2' = A * H2 = H1 * H2.

  On the N-1 step, generate the lower 2 elements of column N-1,
  compute the Householder matrix HN-1 that annihilates them, and
  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
  This is our random orthogonal matrix.
*/
  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 1; j < n; j++ )
  {
/*
  Set the vector that represents the J-th column to be annihilated.
*/
    for ( i = 1; i < j; i++ )
    {
      x[i-1] = 0.0;
    }
    for ( i = j; i <= n; i++ )
    {
      x[i-1] = r8_normal_01 ( seed );
    }
/*
  Compute the vector V that defines a Householder transformation matrix
  H(V) that annihilates the subdiagonal elements of X.
*/
    v = r8vec_house_column ( n, x, j );
/*
  Postmultiply the matrix A by H'(V) = H(V).
*/
    a2 = r8mat_house_axh_new ( n, a, v );

    free ( v );

    r8mat_copy ( n, n, a2, a );

    free ( a2 );
  }

  free ( x );

  return a;
}
/******************************************************************************/

void r8mat_plot ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PLOT "plots" an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the matrix.

    Input, char *TITLE, a title.
*/
{
  int i;
  int j;
  int jhi;
  int jlo;

  printf ( "\n" );
  printf ( "%s\n", title );

  for ( jlo = 1; jlo <= n; jlo = jlo + 70 )
  {
    jhi = i4_min ( jlo + 70-1, n );
    printf ( "\n" );
    printf ( "          " );
    for ( j = jlo; j <= jhi; j++ )
    {
      printf ( "%d", ( j % 10 ) );
    }
    printf ( "\n" );
    printf ( "\n" );

    for ( i = 1; i <= m; i++ )
    {
      printf ( "%6d    ", i );
      for ( j = jlo; j <= jhi; j++ )
      {
        printf ( "%c", r8mat_plot_symbol ( a[i-1+(j-1)*m] ) );
      }
      printf ( "\n" );
    }
  }

  return;
}
/******************************************************************************/

char r8mat_plot_symbol ( double r )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PLOT_SYMBOL returns a symbol for entries of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, double R, a value whose symbol is desired.

    Output, char R8MAT_PLOT_SYMBOL, is
    '-' if R is negative,
    '0' if R is zero,
    '+' if R is positive.
*/
{
  char c;

  if ( r < 0.0 )
  {
    c = '-';
  }
  else if ( r == 0.0 )
  {
    c = '0';
  }
  else if ( 0.0 < r )
  {
    c = '+';
  }

  return c;
}
/******************************************************************************/

double *r8mat_poly_char ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_POLY_CHAR computes the characteristic polynomial of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix A.

    Input, double A[N*N], the N by N matrix.

    Output, double R8MAT_POLY_CHAR[N+1], the coefficients of the characteristic
    polynomial of A.  P(N) contains the coefficient of X**N
    (which will be 1), P(I) contains the coefficient of X**I,
    and P(0) contains the constant term.
*/
{
  int i;
  int order;
  double *p;
  double trace;
  double *work1;
  double *work2;

  p = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );
/*
  Initialize WORK1 to the identity matrix.
*/
  work1 = r8mat_identity_new ( n );

  p[n] = 1.0;

  for ( order = n-1; 0 <= order; order-- )
  {
/*
  Work2 = A * WORK1.
*/
    work2 = r8mat_mm_new ( n, n, n, a, work1 );
/*
  Take the trace.
*/
    trace = r8mat_trace ( n, work2 );
/*
  P(ORDER) = -Trace ( WORK2 ) / ( N - ORDER )
*/
    p[order] = -trace / ( double ) ( n - order );
/*
  WORK1 := WORK2 + P(IORDER) * Identity.
*/
    free ( work1 );

    r8mat_copy ( n, n, work2, work1 );

    free ( work2 );

    for ( i = 0; i < n; i++ )
    {
      work1[i+i*n] = work1[i+i*n] + p[order];
    }
  }

  free ( work1 );

  return p;
}
/******************************************************************************/

double *r8mat_power ( int n, double a[], int npow )

/******************************************************************************/
/*
  Purpose:

    R8MAT_POWER computes a nonnegative power of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The algorithm is:

      B = I
      do NPOW times:
        B = A * B
      end

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Input, double A[N*N], the matrix to be raised to a power.

    Input, int NPOW, the power to which A is to be raised.
    NPOW must be nonnegative.

    Output, double B[N*N], the value of A**NPOW.
*/
{
  double *b;
  double *c;
  int ipow;

  if ( npow < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_POWER - Fatal error!\n" );
    fprintf ( stderr, "  Input value of NPOW < 0.\n" );
    fprintf ( stderr, "  NPOW = %d\n", npow );
    exit ( 1 );
  }

  b = r8mat_identity_new ( n );

  for ( ipow = 1; ipow <= npow; ipow++ )
  {
    c = r8mat_mm_new ( n, n, n, a, b );
    r8mat_copy ( n, n, c, b );
    free ( c );
  }

  return b;
}
/******************************************************************************/

void r8mat_power_method ( int n, double a[], double *r, double v[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_POWER_METHOD applies the power method to a matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    If the power method has not converged, then calling the routine
    again immediately with the output from the previous call will
    continue the iteration.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Input, double A[N*N], the matrix.

    Output, double *R, the estimated eigenvalue.

    Input/output, double V[N], on input, an estimate
    for the eigenvector.  On output, an improved estimate for the
    eigenvector.
*/
{
  double *av;
  double eps;
  int i;
  int it;
  double it_eps = 0.0001;
  int it_max = 100;
  int it_min = 10;
  int j;
  double r2;
  double r_old;

  eps = sqrt ( r8_epsilon ( ) );

  *r = r8vec_norm ( n, v );

  if ( *r == 0.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      v[i] = 1.0;
    }
    *r = sqrt ( ( double ) n );
  }

  for ( i = 0; i < n; i++ )
  {
    v[i] = v[i] / *r;
  }

  for ( it = 1; it <= it_max; it++ )
  {
    av = r8mat_mv_new ( n, n, a, v );

    r_old = *r;
    *r = r8vec_norm ( n, av );

    if ( it_min < it )
    {
      if ( r8_abs ( *r - r_old ) <= it_eps * ( 1.0 + r8_abs ( *r ) ) )
      {
        break;
      }
    }

    r8vec_copy ( n, av, v );

    free ( av );

    if ( *r != 0.0 )
    {
      for ( i = 0; i < n; i++ )
      {
        v[i] = v[i] / *r;
      }
    }
//
//  Perturb V a bit, to avoid cases where the initial guess is exactly
//  the eigenvector of a smaller eigenvalue.
//
    if ( it < it_max / 2 )
    {
      j = ( ( it - 1 ) % n );
      v[j] = v[j] + eps * ( 1.0 + r8_abs ( v[j] ) );
      r2 = r8vec_norm ( n, v );
      for ( i = 0; i < n; i++ )
      {
        v[i] = v[i] / r2;
      }
    }
  }
  return;
}
/******************************************************************************/

void r8mat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT prints an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*M]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[M*N], the M by N matrix.

    Input, char *TITLE, a title.
*/
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_PRINT_SOME prints some of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows of the matrix.
    M must be positive.

    Input, int N, the number of columns of the matrix.
    N must be positive.

    Input, double A[M*N], the matrix.

    Input, int ILO, JLO, IHI, JHI, designate the first row and
    column, and the last row and column to be printed.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }
/*
  Print the columns of the matrix, in strips of 5.
*/
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    fprintf ( stdout, "\n" );
/*
  For each column J in the current range...

  Write the header.
*/
    fprintf ( stdout, "  Col:  ");
    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "  %7d     ", j - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row\n" );
    fprintf ( stdout, "\n" );
/*
  Determine the range of the rows in this strip.
*/
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
/*
  Print out (up to) 5 entries in row I, that lie in the current strip.
*/
      fprintf ( stdout, "%5d:", i - 1 );
      for ( j = j2lo; j <= j2hi; j++ )
      {
        fprintf ( stdout, "  %14f", a[i-1+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

void r8mat_ref ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_REF computes the row echelon form of a matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    A matrix is in row echelon form if:

    * The first nonzero entry in each row is 1.

    * The leading 1 in a given row occurs in a column to
      the right of the leading 1 in the previous row.

    * Rows which are entirely zero must occur last.

  Example:

    Input matrix:

     1.0  3.0  0.0  2.0  6.0  3.0  1.0
    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
     3.0  9.0  0.0  0.0  6.0  6.0  2.0
    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0

    Output matrix:

     1.0  3.0  0.0  2.0  6.0  3.0  1.0
     0.0  0.0  0.0  1.0  2.0  4.5  1.5
     0.0  0.0  0.0  0.0  0.0  1.0  0.3
     0.0  0.0  0.0  0.0  0.0  0.0  0.0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of
    the matrix A.

    Input/output, double A[M*N].  On input, the matrix to be
    analyzed.  On output, the REF form of the matrix.
*/
{
  int i;
  int j;
  int lead;
  int r;
  double temp;

  lead = 0;

  for ( r = 0; r < m; r++ )
  {
    if ( n - 1 < lead )
    {
      break;
    }

    i = r;

    while ( a[i+lead*m] == 0.0 )
    {
      i = i + 1;

      if ( m - 1 < i )
      {
        i = r;
        lead = lead + 1;
        if ( n - 1 < lead )
        {
          lead = -1;
          break;
         }
      }
    }

    if ( lead < 0 )
    {
      break;
    }

    for ( j = 0; j < n; j++ )
    {
      temp     = a[i+j*m];
      a[i+j*m] = a[r+j*m];
      a[r+j*m] = temp;
    }

    temp = a[r+lead*m];

    for ( j = 0; j < n; j++ )
    {
      a[r+j*m] = a[r+j*m] / temp;
    }

    for ( i = r + 1; i < m; i++ )
    {
      temp = a[i+lead*m];
      for ( j = 0; j < n; j++ )
      {
        a[i+j*m] = a[i+j*m] - temp * a[r+j*m];
      }
    }
    lead = lead + 1;
  }
  return;
}
/******************************************************************************/

double r8mat_rms ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_RMS returns the RMS norm of an R8MAT.

  Discussion:

    The matrix RMS norm is defined as:

      R8MAT_RMS = 
        sqrt ( sum ( 1 <= J <= N ) sum ( 1 <= I <= M ) A(I,J)^2 / M / N )

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the dimensions of the array.

    Input, double A[M*N], the array.

    Output, double R8MAT_RMS, the RMS norm of A.
*/
{
  int i;
  int j;
  double value;

  value = 0.0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + a[i+j*m] * a[i+j*m];
    }
    value = sqrt ( value / ( double ) ( m ) / ( double ) ( n ) );
  }

  return value;
}
/******************************************************************************/

void r8mat_rref ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_RREF computes the reduced row echelon form of a matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    A matrix is in row echelon form if:

    * The first nonzero entry in each row is 1.

    * The leading 1 in a given row occurs in a column to
      the right of the leading 1 in the previous row.

    * Rows which are entirely zero must occur last.

    The matrix is in reduced row echelon form if, in addition to
    the first three conditions, it also satisfies:

    * Each column containing a leading 1 has no other nonzero entries.

  Example:

    Input matrix:

     1.0  3.0  0.0  2.0  6.0  3.0  1.0
    -2.0 -6.0  0.0 -2.0 -8.0  3.0  1.0
     3.0  9.0  0.0  0.0  6.0  6.0  2.0
    -1.0 -3.0  0.0  1.0  0.0  9.0  3.0

    Output matrix:

     1.0  3.0  0.0  0.0  2.0  0.0  0.0
     0.0  0.0  0.0  1.0  2.0  0.0  0.0
     0.0  0.0  0.0  0.0  0.0  1.0  0.3
     0.0  0.0  0.0  0.0  0.0  0.0  0.0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of
    the matrix A.

    Input/output, double A[M*N].  On input, the matrix to be
    analyzed.  On output, the RREF form of the matrix.
*/
{
  int i;
  int j;
  int lead;
  int r;
  double temp;

  lead = 0;

  for ( r = 0; r < m; r++ )
  {
    if ( n - 1 < lead )
    {
      break;
    }

    i = r;

    while ( a[i+lead*m] == 0.0 )
    {
      i = i + 1;

      if ( m - 1 < i )
      {
        i = r;
        lead = lead + 1;
        if ( n - 1 < lead )
        {
          lead = -1;
          break;
         }
      }
    }

    if ( lead < 0 )
    {
      break;
    }

    for ( j = 0; j < n; j++ )
    {
      temp     = a[i+j*m];
      a[i+j*m] = a[r+j*m];
      a[r+j*m] = temp;
    }

    temp = a[r+lead*m];

    for ( j = 0; j < n; j++ )
    {
      a[r+j*m] = a[r+j*m] / temp;
    }

    for ( i = 0; i < m; i++ )
    {
      if ( i != r )
      {
        temp = a[i+lead*m];
        for ( j = 0; j < n; j++ )
        {
          a[i+j*m] = a[i+j*m] - temp * a[r+j*m];
        }
      }
    }
    lead = lead + 1;

  }
  return;
}
/******************************************************************************/

void r8mat_scale ( int m, int n, double s, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SCALE multiplies an R8MAT by a scalar.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double S, the scale factor.

    Input/output, double A[M*N], the matrix to be scaled.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] * s;
    }
  }
  return;
}
/******************************************************************************/

int r8mat_significant ( int m, int n, double r[], double s[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SIGNIFICANT determines if an R8MAT is significant compared to another.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the dimension of the matrices.

    Input, double R[M*N], the vector to be compared against.

    Input, double S[M*N], the vector to be compared.

    Output, int R8MAT_SIGNIFICANT, is TRUE if S is significant
    compared to R.
*/
{
  int i;
  int j;
  double t;
  double tol;
  int value;

  value = 0;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      t = r[i+j*m] + s[i+j*m];
      tol = r8_epsilon ( ) * r8_abs ( r[i+j*m] );

      if ( tol < r8_abs ( r[i+j*m] - t ) )
      {
        value = 1;
        break;
      }
    }
  }
  return value;
}
/******************************************************************************/

int r8mat_solve ( int n, int rhs_num, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    Entry A(I,J) is stored as A[I+J*N]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 October 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, int RHS_NUM, the number of right hand sides.  RHS_NUM
    must be at least 0.

    Input/output, double A[N*(N+RHS_NUM)], contains in rows and columns 1
    to N the coefficient matrix, and in columns N+1 through
    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
    area has been destroyed, while the right hand sides have
    been overwritten with the corresponding solutions.

    Output, int R8MAT_SOLVE, singularity flag.
    0, the matrix was not singular, the solutions were computed;
    J, factorization failed on step J, and the solutions could not
    be computed.
*/
{
  double apivot;
  double factor;
  int i;
  int ipivot;
  int j;
  int k;
  double temp;

  for ( j = 0; j < n; j++ )
  {
/*
  Choose a pivot row.
*/
    ipivot = j;
    apivot = a[j+j*n];

    for ( i = j; i < n; i++ )
    {
      if ( r8_abs ( apivot ) < r8_abs ( a[i+j*n] ) )
      {
        apivot = a[i+j*n];
        ipivot = i;
      }
    }

    if ( apivot == 0.0 )
    {
      return j;
    }
/*
  Interchange.
*/
    for ( i = 0; i < n + rhs_num; i++ )
    {
      temp          = a[ipivot+i*n];
      a[ipivot+i*n] = a[j+i*n];
      a[j+i*n]      = temp;
    }
/*
  A(J,J) becomes 1.
*/
    a[j+j*n] = 1.0;
    for ( k = j; k < n + rhs_num; k++ )
    {
      a[j+k*n] = a[j+k*n] / apivot;
    }
/*
  A(I,J) becomes 0.
*/
    for ( i = 0; i < n; i++ )
    {
      if ( i != j )
      {
        factor = a[i+j*n];
        a[i+j*n] = 0.0;
        for ( k = j; k < n + rhs_num; k++ )
        {
          a[i+k*n] = a[i+k*n] - factor * a[j+k*n];
        }
      }
    }
  }

  return 0;
}
/******************************************************************************/

double *r8mat_solve_2d ( double a[], double b[], double *det )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SOLVE_2D solves a 2 by 2 linear system using Cramer's rule.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    If the determinant DET is returned as zero, then the matrix A is
    singular, and does not have an inverse.  In that case, X is
    returned as the NULL vector.

    If DET is nonzero, then its value is roughly an estimate
    of how nonsingular the matrix A is.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double A[2*2], the matrix.

    Input, double B[2], the right hand side.

    Output, double *DET, the determinant of the system.

    Output, double R8MAT_SOLVE_2D[2], the solution of the system,
    if DET is nonzero.  Otherwise, the NULL vector.
*/
{
  double *x;
/*
  Compute the determinant.
*/
  *det = a[0+0*2] * a[1+1*2] - a[0+1*2] * a[1+0*2];
/*
  If the determinant is zero, bail out.
*/
  if ( *det == 0.0 )
  {
    return NULL;
  }
/*
  Compute the solution.
*/
  x = ( double * ) malloc ( 2 * sizeof ( double ) );

  x[0] = (  a[1+1*2] * b[0] - a[0+1*2] * b[1] ) / ( *det );
  x[1] = ( -a[1+0*2] * b[0] + a[0+0*2] * b[1] ) / ( *det );

  return x;
}
/******************************************************************************/

double *r8mat_solve_3d ( double a[], double b[], double *det )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SOLVE_3D solves a 3 by 3 linear system using Cramer's rule.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    If the determinant DET is returned as zero, then the matrix A is
    singular, and does not have an inverse.  In that case, X is
    returned as the NULL vector.

    If DET is nonzero, then its value is roughly an estimate
    of how nonsingular the matrix A is.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double A[3*3], the matrix.

    Input, double B[3], the right hand side.

    Output, double *DET, the determinant of the system.

    Output, double R8MAT_SOLVE_3D[3], the solution of the system,
    if DET is nonzero.  Otherwise, the NULL vector.
*/
{
  double *x;
/*
  Compute the determinant.
*/
  *det =  a[0+0*3] * ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] )
        + a[0+1*3] * ( a[1+2*3] * a[2+0*3] - a[1+0*3] * a[2+2*3] )
        + a[0+2*3] * ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] );
/*
  If the determinant is zero, bail out.
*/
  if ( *det == 0.0 )
  {
    return NULL;
  }
/*
  Compute the solution.
*/
  x = ( double * ) malloc ( 3 * sizeof ( double ) );

  x[0] = (   ( a[1+1*3] * a[2+2*3] - a[1+2*3] * a[2+1*3] ) * b[0]
           - ( a[0+1*3] * a[2+2*3] - a[0+2*3] * a[2+1*3] ) * b[1]
           + ( a[0+1*3] * a[1+2*3] - a[0+2*3] * a[1+1*3] ) * b[2] ) / ( *det );

  x[1] = ( - ( a[1+0*3] * a[2+2*3] - a[1+2*3] * a[2+0*3] ) * b[0]
           + ( a[0+0*3] * a[2+2*3] - a[0+2*3] * a[2+0*3] ) * b[1]
           - ( a[0+0*3] * a[1+2*3] - a[0+2*3] * a[1+0*3] ) * b[2] ) / ( *det );

  x[2] = (   ( a[1+0*3] * a[2+1*3] - a[1+1*3] * a[2+0*3] ) * b[0]
           - ( a[0+0*3] * a[2+1*3] - a[0+1*3] * a[2+0*3] ) * b[1]
           + ( a[0+0*3] * a[1+1*3] - a[0+1*3] * a[1+0*3] ) * b[2] ) / ( *det );

  return x;
}
/******************************************************************************/

double *r8mat_solve2 ( int n, double a[], double b[], int *ierror )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SOLVE2 computes the solution of an N by N linear system.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The linear system may be represented as

      A*X = B

    If the linear system is singular, but consistent, then the routine will
    still produce a solution.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of equations.

    Input/output, double A[N*N].
    On input, A is the coefficient matrix to be inverted.
    On output, A has been overwritten.

    Input/output, double B[N].
    On input, B is the right hand side of the system.
    On output, B has been overwritten.

    Output, double R8MAT_SOLVE2[N], the solution of the linear system.

    Output, int *IERROR.
    0, no error detected.
    1, consistent singularity.
    2, inconsistent singularity.
*/
{
  double amax;
  int i;
  int imax;
  int j;
  int k;
  int *piv;
  double *x;

  *ierror = 0;

  piv = i4vec_zero_new ( n );
  x = r8vec_zero_new ( n );
/*
  Process the matrix.
*/
  for ( k = 1; k <= n; k++ )
  {
/*
  In column K:
    Seek the row IMAX with the properties that:
      IMAX has not already been used as a pivot;
      A(IMAX,K) is larger in magnitude than any other candidate.
*/
    amax = 0.0;
    imax = 0;
    for ( i = 1; i <= n; i++ )
    {
      if ( piv[i-1] == 0 )
      {
        if ( amax < r8_abs ( a[i-1+(k-1)*n] ) )
        {
          imax = i;
          amax = r8_abs ( a[i-1+(k-1)*n] );
        }
      }
    }
/*
  If you found a pivot row IMAX, then,
    eliminate the K-th entry in all rows that have not been used for pivoting.
*/
    if ( imax != 0 )
    {
      piv[imax-1] = k;
      for ( j = k+1; j <= n; j++ )
      {
        a[imax-1+(j-1)*n] = a[imax-1+(j-1)*n] / a[imax-1+(k-1)*n];
      }
      b[imax-1] = b[imax-1] / a[imax-1+(k-1)*n];
      a[imax-1+(k-1)*n] = 1.0;

      for ( i = 1; i <= n; i++ )
      {
        if ( piv[i-1] == 0 )
        {
          for ( j = k+1; j <= n; j++ )
          {
            a[i-1+(j-1)*n] = a[i-1+(j-1)*n] - a[i-1+(k-1)*n] * a[imax-1+(j-1)*n];
          }
          b[i-1] = b[i-1] - a[i-1+(k-1)*n] * b[imax-1];
          a[i-1+(k-1)*n] = 0.0;
        }
      }
    }
  }
/*
  Now, every row with nonzero PIV begins with a 1, and
  all other rows are all zero.  Begin solution.
*/
  for ( j = n; 1 <= j; j-- )
  {
    imax = 0;
    for ( k = 1; k <= n; k++ )
    {
      if ( piv[k-1] == j )
      {
        imax = k;
      }
    }

    if ( imax == 0 )
    {
      x[j-1] = 0.0;

      if ( b[j-1] == 0.0 )
      {
        *ierror = 1;
        printf ( "\n" );
        printf ( "R8MAT_SOLVE2 - Warning:\n" );
        printf ( "  Consistent singularity, equation = %d\n", j );
      }
      else
      {
        *ierror = 2;
        printf ( "\n" );
        printf ( "R8MAT_SOLVE2 - Warning:\n" );
        printf ( "  Inconsistent singularity, equation = %d\n", j );
      }
    }
    else
    {
      x[j-1] = b[imax-1];

      for ( i = 1; i <= n; i++ )
      {
        if ( i != imax )
        {
          b[i-1] = b[i-1] - a[i-1+(j-1)*n] * x[j-1];
        }
      }
    }
  }

  free ( piv );

  return x;
}
/******************************************************************************/

double r8mat_sum ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SUM returns the sum of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    04 January 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the array.

    Output, double R8MAT_SUM, the sum of the entries.
*/
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + a[i+j*m];
    }
  }
  return value;
}
/******************************************************************************/

double *r8mat_symm_eigen ( int n, double x[], double q[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SYMM_EIGEN returns a symmetric matrix with given eigensystem.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The user must supply the desired eigenvalue vector, and the desired
    eigenvector matrix.  The eigenvector matrix must be orthogonal.  A
    suitable random orthogonal matrix can be generated by
    R8MAT_ORTH_UNIFORM_NEW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Input, double X[N], the desired eigenvalues for the matrix.

    Input, double Q[N*N], the eigenvector matrix of A.

    Output, double R8MAT_SYMM_EIGEN[N*N], a symmetric N by N matrix with
    eigenvalues X and eigenvectors the columns of Q.
*/
{
  double *a;
  int i;
  int j;
  int k;
/*
  Set A = Q * Lambda * Q'.
*/
  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      a[i+j*n] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        a[i+j*n] = a[i+j*n] + q[i+k*n] * x[k] * q[j+k*n];
      }
    }
  }

  return a;
}
/******************************************************************************/

void r8mat_symm_jacobi ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_SYMM_JACOBI applies Jacobi eigenvalue iteration to a symmetric matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    This code was modified so that it treats as zero the off-diagonal
    elements that are sufficiently close to, but not exactly, zero.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 October 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of A.

    Input/output, double A[N*N], a symmetric N by N matrix.
    On output, the matrix has been overwritten by an approximately
    diagonal matrix, with the eigenvalues on the diagonal.
*/
{
  double c;
  double eps = 0.00001;
  int i;
  int it;
  int it_max = 100;
  int j;
  int k;
  double norm_fro;
  double s;
  double sum2;
  double t;
  double t1;
  double t2;
  double u;

  norm_fro = r8mat_norm_fro ( n, n, a );

  it = 0;

  for ( ; ; )
  {
    it = it + 1;

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < i; j++ )
      {
        if ( eps * norm_fro < r8_abs ( a[i+j*n] ) + r8_abs ( a[j+i*n] ) )
        {
          u = ( a[j+j*n] - a[i+i*n] ) / ( a[i+j*n] + a[j+i*n] );

          t = r8_sign ( u ) / ( r8_abs ( u ) + sqrt ( u * u + 1.0 ) );
          c = 1.0 / sqrt ( t * t + 1.0 );
          s = t * c;
/*
  A -> A * Q.
*/
          for ( k = 0; k < n; k++ )
          {
            t1 = a[i+k*n];
            t2 = a[j+k*n];
            a[i+k*n] = t1 * c - t2 * s;
            a[j+k*n] = t1 * s + t2 * c;
          }
/*
  A -> QT * A
*/
          for ( k = 0; k < n; k++ )
          {
            t1 = a[k+i*n];
            t2 = a[k+j*n];
            a[k+i*n] = c * t1 - s * t2;
            a[k+j*n] = s * t1 + c * t2;
          }
        }
      }
    }
/*
  Test the size of the off-diagonal elements.
*/
    sum2 = 0.0;
    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < i; j++ )
      {
        sum2 = sum2 + r8_abs ( a[i+j*n] );
      }
    }

    if ( sum2 <= eps * ( norm_fro + 1.0 ) )
    {
      break;
    }

    if ( it_max <= it )
    {
      break;
    }

  }

  return;
}
/******************************************************************************/

int r8mat_to_r8plu ( int n, double a[], int pivot[], double lu[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TO_R8PLU factors a general matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    This routine is a simplified version of the LINPACK routine DGEFA.
    Fortran conventions are used to index doubly-dimensioned arrays.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 May 2010

  Author:

    John Burkardt

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, double A[N*N], the matrix to be factored.

    Output, int PIVOT[N], a vector of pivot indices.

    Output, double LU[N*N], an upper triangular matrix U and the multipliers
    L which were used to obtain it.  The factorization can be written
    A = L * U, where L is a product of permutation and unit lower
    triangular matrices and U is upper triangular.

    Output, int R8MAT_TO_R8PLU, singularity flag.
    0, no singularity detected.
    nonzero, the factorization failed on the R8MAT_TO_R8PLU-th step.
*/
{
  int i;
  int info;
  int j;
  int k;
  int l;
  double temp;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      lu[i+j*n] = a[i+j*n];
    }
  }
  info = 0;

  for ( k = 1; k <= n-1; k++ )
  {
/*
  Find L, the index of the pivot row.
*/
    l = k;
    for ( i = k+1; i <= n; i++ )
    {
      if ( r8_abs ( lu[l-1+(k-1)*n] ) < r8_abs ( lu[i-1+(k-1)*n] ) )
      {
        l = i;
      }
    }

    pivot[k-1] = l;
/*
  If the pivot index is zero, the algorithm has failed.
*/
    if ( lu[l-1+(k-1)*n] == 0.0 )
    {
      info = k;
      return info;
    }
/*
  Interchange rows L and K if necessary.
*/
    if ( l != k )
    {
      temp            = lu[l-1+(k-1)*n];
      lu[l-1+(k-1)*n] = lu[k-1+(k-1)*n];
      lu[k-1+(k-1)*n] = temp;
    }
/*
  Normalize the values that lie below the pivot entry A(K,K).
*/
    for ( i = k+1; i <= n; i++ )
    {
      lu[i-1+(k-1)*n] = -lu[i-1+(k-1)*n] / lu[k-1+(k-1)*n];
    }
/*
  Row elimination with column indexing.
*/
    for ( j = k+1; j <= n; j++ )
    {
      if ( l != k )
      {
        temp            = lu[l-1+(j-1)*n];
        lu[l-1+(j-1)*n] = lu[k-1+(j-1)*n];
        lu[k-1+(j-1)*n] = temp;
      }

      for ( i = k+1; i <= n; i++ )
      {
        lu[i-1+(j-1)*n] = lu[i-1+(j-1)*n] + lu[i-1+(k-1)*n] * lu[k-1+(j-1)*n];
      }
    }
  }

  pivot[n-1] = n;

  if ( lu[n-1+(n-1)*n] == 0.0 )
  {
    info = n;
  }

  return info;
}
/******************************************************************************/

double r8mat_trace ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRACE computes the trace of an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The trace of a square matrix is the sum of the diagonal elements.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix A.

    Input, double A[N*N], the matrix whose trace is desired.

    Output, double R8MAT_TRACE, the trace of the matrix.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i+i*n];
  }

  return value;
}
/******************************************************************************/

void r8mat_transpose_in_place ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_IN_PLACE transposes a square matrix in place.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of rows and columns of the matrix A.

    Input/output, double A[N*N], the matrix to be transposed.
*/
{
  int i;
  int j;
  double t;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      t        = a[i+j*n];
      a[i+j*n] = a[j+i*n];
      a[j+i*n] = t;
    }
  }
  return;
}
/******************************************************************************/

double *r8mat_transpose_new ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_NEW returns the transpose of a matrix.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 June 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns of the matrix A.

    Input, double A[M*N], the matrix whose transpose is desired.

    Output, double R8MAT_TRANSPOSE_NEW[N*M], the transposed matrix.
*/
{
  double *b;
  int i;
  int j;

  b = ( double * ) malloc ( n * m * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      b[j+i*n] = a[i+j*m];
    }
  }
  return b;
}
/******************************************************************************/

void r8mat_transpose_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 May 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, char *TITLE, a title.
*/
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
/******************************************************************************/

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], an M by N matrix to be printed.

    Input, int ILO, JLO, the first row and column to print.

    Input, int IHI, JHI, the last row and column to print.

    Input, char *TITLE, a title.
*/
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );

  if ( m <= 0 || n <= 0 )
  {
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  (None)\n" );
    return;
  }

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Row:" );
    for ( i = i2lo; i <= i2hi; i++ )
    {
      fprintf ( stdout, "  %7d     ", i - 1 );
    }
    fprintf ( stdout, "\n" );
    fprintf ( stdout, "  Col\n" );
    fprintf ( stdout, "\n" );

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      fprintf ( stdout, "%5d:", j - 1 );
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        fprintf ( stdout, "  %14f", a[(i-1)+(j-1)*m] );
      }
      fprintf ( stdout, "\n" );
    }
  }

  return;
# undef INCX
}
/******************************************************************************/

double *r8mat_u_inverse ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_U_INVERSE inverts an upper triangular R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    An upper triangular matrix is a matrix whose only nonzero entries
    occur on or above the diagonal.

    The inverse of an upper triangular matrix is an upper triangular matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, number of rows and columns in the matrix.

    Input, double A[N*N], the upper triangular matrix.

    Output, double R8MAT_U_INVERSE[N*N], the inverse matrix.
*/
{
  double *b;
  int i;
  int j;
  int k;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = n-1; 0 <= j; j-- )
  {
    for ( i = n-1; 0 <= i; i-- )
    {
      if ( j < i )
      {
        b[i+j*n] = 0.0;
      }
      else if ( i == j )
      {
        b[i+j*n] = 1.0 / a[i+j*n];
      }
      else
      {
        b[i+j*n] = 0.0;
        for ( k = i+1; k <= j; k++ )
        {
          b[i+j*n] = b[i+j*n] - a[i+k*n] * b[k+j*n];
        }
       b[i+j*n] = b[i+j*n] / a[i+i*n];
      }
    }
  }

  return b;
}
/******************************************************************************/

double *r8mat_u1_inverse ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_U1_INVERSE inverts a unit upper triangular R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    A unit upper triangular matrix is a matrix with only 1's on the main
    diagonal, and only 0's below the main diagonal.

    The inverse of a unit upper triangular matrix is also
    a unit upper triangular matrix.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    C translation by John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, number of rows and columns in the matrix.

    Input, double A[N*N], the unit upper triangular matrix.

    Output, double R8MAT_U1_INVERSE[N*N), the inverse matrix.
*/
{
  double *b;
  int i;
  int j;
  int k;

  b = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( j = n-1; 0 <= j; j-- )
  {
    for ( i = n-1; 0 <= i; i-- )
    {
      if ( j < i )
      {
        b[i+j*n] = 0.0;
      }
      else if ( i == j )
      {
        b[i+j*n] = 1.0;
      }
      else
      {
        b[i+j*n] = 0.0;
        for ( k = i+1; k <= j; k++ )
        {
          b[i+j*n] = b[i+j*n] - a[i+k*n] * b[k+j*n];
        }
       b[i+j*n] = b[i+j*n] / a[i+i*n];
      }
    }
  }

  return b;
}
/******************************************************************************/

void r8mat_uniform_01 ( int m, int n, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01 fills an R8MAT with pseudorandom values scaled to [0,1].

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0, otherwise the output value of SEED
    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
    been updated.

    Output, double R[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int j;
  int k;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }

      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_uniform_01_new ( int m, int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_01_NEW fills an R8MAT with pseudorandom values scaled to [0,1].

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Springer Verlag, pages 201-202, 1983.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, pages 362-376, 1986.

    Philip Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, pages 136-143, 1969.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0, otherwise the output value of SEED
    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
    been updated.

    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int j;
  int k;
  double *r;

  r = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

void r8mat_uniform_ab ( int m, int n, double a, double b, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_AB returns a pseudorandom R8MAT scaled to [A,B].

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A, B, the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, double R[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_uniform_ab_new ( int m, int n, double a, double b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_AB_NEW returns a scaled pseudorandom R8MAT.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A, B, the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, double R8MAT_UNIFORM_AB_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = (double*) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

void r8mat_uniform_abvec ( int m, int n, double a[], double b[], int *seed, 
  double r[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_ABVEC returns a scaled pseudorandom R8MAT.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M], B[M], the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, double R[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_ABVEC - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = a[i] + ( b[i] - a[i] ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
/******************************************************************************/

double *r8mat_uniform_abvec_new ( int m, int n, double a[], double b[], int *seed )

/******************************************************************************/
/*
  Purpose:

    R8MAT_UNIFORM_ABVEC_NEW returns a scaled pseudorandom R8MAT.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 October 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M], B[M], the limits of the pseudorandom values.

    Input/output, int *SEED, the "seed" value.  Normally, this
    value should not be 0.  On output, SEED has 
    been updated.

    Output, double R8MAT_UNIFORM_ABVEC_NEW[M*N], a matrix of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_UNIFORM_ABVEC_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = (double*) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }
      r[i+j*m] = a[i] + ( b[i] - a[i] ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
/******************************************************************************/

double *r8mat_vand2 ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_VAND2 returns the N by N row Vandermonde matrix A.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

    The row Vandermonde matrix returned by this routine reads "across"
    rather than down.  In particular, each row begins with a 1, followed by
    some value X, followed by successive powers of X.

    The formula for the matrix entries is:

      A(I,J) = X(I)**(J-1)

  Properties:

    A is nonsingular if, and only if, the X values are distinct.

    The determinant of A is

      det(A) = product ( 2 <= I <= N ) (
        product ( 1 <= J <= I-1 ) ( ( X(I) - X(J) ) ) ).

    The matrix A is generally ill-conditioned.

  Example:

    N = 5, X = (2, 3, 4, 5, 6)

    1 2  4   8   16
    1 3  9  27   81
    1 4 16  64  256
    1 5 25 125  625
    1 6 36 216 1296

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix desired.

    Input, double X[N], the values that define A.

    Output, double R8MAT_VAND2[N*N], the N by N row Vandermonde matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( n * n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( j == 0 && x[i] == 0.0 )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = pow ( x[i], j );
      }
    }
  }

  return a;
}
/******************************************************************************/

void r8mat_zero ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ZERO zeroes an R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, double A[M*N], a matrix of zeroes.
*/
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return;
}
/******************************************************************************/

double *r8mat_zero_new ( int m, int n )

/******************************************************************************/
/*
  Purpose:

    R8MAT_ZERO_NEW returns a new zeroed R8MAT.

  Discussion:

    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
    in column-major order.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 September 2008

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Output, double R8MAT_ZERO[M*N], the new zeroed matrix.
*/
{
  double *a;
  int i;
  int j;

  a = ( double * ) malloc ( m * n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 0.0;
    }
  }
  return a;
}
/******************************************************************************/

double r8plu_det ( int n, int pivot[], double lu[] )

/******************************************************************************/
/*
  Purpose:

    R8PLU_DET computes the determinant of a real PLU matrix.

  Discussion:

    The matrix should have been factored by R8MAT_TO_R8PLU.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 May 2010

  Author:

    John Burkardt

  Reference:

    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    LINPACK User's Guide,
    SIAM, 1979

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, int PIVOT[N], the pivot vector computed by R8MAT_TO_R8PLU.

    Input, double LU[N*N], the LU factors computed by R8MAT_TO_R8PLU.

    Output, double R8PLU_DET, the determinant of the matrix.
*/
{
  double det;
  int i;

  det = 1.0;

  for ( i = 0; i < n; i++ )
  {
    det = det * lu[i+i*n];
    if ( pivot[i] != i+1 )
    {
      det = -det;
    }
  }

  return det;
}
/******************************************************************************/

void r8plu_inverse ( int n, int pivot[], double lu[], double a_inverse[] )

/******************************************************************************/
/*
  Purpose:

    R8PLU_INVERSE computes the inverse of a real PLU matrix.

  Discussion:

    The matrix should have been factored by R8MAT_TO_R8PLU.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix A.

    Input, int PIVOT[N], the pivot vector from R8MAT_TO_R8PLU.

    Input, double LU[N*N], the LU factors computed by R8MAT_TO_R8PLU.

    Output, double A_INVERSE[N*N], the inverse of the original matrix
    A that was factored by R8MAT_TO_R8PLU.
*/
{
  int i;
  int j;
  int k;
  double temp;
  double *work;
//
  work = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      a_inverse[i+j*n] = lu[i+j*n];
    }
  }
/*
  Compute Inverse(U).
*/
  for ( k = 1; k <= n; k++ )
  {
    a_inverse[k-1+(k-1)*n]     = 1.0 / a_inverse[k-1+(k-1)*n];
    for ( i = 1; i <= k-1; i++ )
    {
      a_inverse[i-1+(k-1)*n] = -a_inverse[i-1+(k-1)*n] * a_inverse[k-1+(k-1)*n];
    }

    for ( j = k + 1; j <= n; j++ )
    {
      temp                     = a_inverse[k-1+(j-1)*n];
      a_inverse[k-1+(j-1)*n]   = 0.0;
      for ( i = 1; i <= k; i++ )
      {
        a_inverse[i-1+(j-1)*n] = a_inverse[i-1+(j-1)*n]
          + temp * a_inverse[i-1+(k-1)*n];
      }
    }
  }
/*
  Form Inverse(U) * Inverse(L).
*/
  for ( k = n - 1; 1 <= k; k-- )
  {
    for ( i = k + 1; i <= n; i++ )
    {
      work[i-1] = a_inverse[i-1+(k-1)*n];
      a_inverse[i-1+(k-1)*n] = 0.0;
    }

    for ( j = k + 1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        a_inverse[i-1+(k-1)*n] = a_inverse[i-1+(k-1)*n]
          + a_inverse[i-1+(j-1)*n] * work[j-1];
      }
    }

    if ( pivot[k-1] != k )
    {
      for ( i = 1; i <= n; i++ )
      {
        temp                            = a_inverse[i-1+(k-1)*n];
        a_inverse[i-1+(k-1)*n]          = a_inverse[i-1+(pivot[k-1]-1)*n];
        a_inverse[i-1+(pivot[k-1]-1)*n] = temp;
      }
    }
  }

  free ( work );

  return;
}
/******************************************************************************/

void r8plu_mul ( int n, int pivot[], double lu[], double x[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8PLU_MUL computes A * x using the PLU factors of A.

  Discussion:

    It is assumed that R8MAT_TO_R8PLU has computed the PLU factors of
    the matrix A.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, int PIVOT[N], the pivot vector computed by R8MAT_TO_R8PLU.

    Input, double LU[N*N], the matrix factors computed by R8MAT_TO_R8PLU.

    Input, double X[N], the vector to be multiplied.

    Output, double B[N], the result of the multiplication.
*/
{
  int i;
  int j;
  int k;
  double temp;

  for ( i = 0; i < n; i++ )
  {
    b[i] = x[i];
  }
/*
  Y = U * X.
*/
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 0; i < j-1; i++ )
    {
      b[i] = b[i] + lu[i+(j-1)*n] * b[j-1];
    }
    b[j-1] = lu[j-1+(j-1)*n] * b[j-1];
  }
/*
  B = PL * Y = PL * U * X = A * x.
*/
  for ( j = n-1; 1 <= j; j-- )
  {
    for ( i = j; i < n; i++ )
    {
      b[i] = b[i] - lu[i+(j-1)*n] * b[j-1];
    }

    k = pivot[j-1];

    if ( k != j )
    {
      temp = b[k-1];
      b[k-1] = b[j-1];
      b[j-1] = temp;
    }
  }

  return;
}
/******************************************************************************/

void r8plu_sol ( int n, int pivot[], double lu[], double b[], double x[] )

/******************************************************************************/
/*
  Purpose:

    R8PLU_SOL solves a linear system A*x=b from the PLU factors.

  Discussion:

    The PLU factors should have been computed by R8MAT_TO_R8PLU.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, int PIVOT[N], the pivot vector from R8MAT_TO_R8PLU.

    Input, double LU[N*N], the LU factors from R8MAT_TO_R8PLU.

    Input, double B[N], the right hand side vector.

    Output, double X[N], the solution vector.
*/
{
  int i;
  int j;
  int k;
  double temp;
/*
  Solve PL * Y = B.
*/
  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  for ( k = 1; k <= n-1; k++ )
  {
    j = pivot[k-1];

    if ( j != k )
    {
      temp   = x[j-1];
      x[j-1] = x[k-1];
      x[k-1] = temp;
    }

    for ( i = k+1; i <= n; i++ )
    {
      x[i-1] = x[i-1] + lu[i-1+(k-1)*n] * x[k-1];
    }
  }
/*
  Solve U * X = Y.
*/
  for ( k = n; 1 <= k; k-- )
  {
    x[k-1] = x[k-1] / lu[k-1+(k-1)*n];
    for ( i = 1; i <= k-1; i++ )
    {
      x[i-1] = x[i-1] - lu[i-1+(k-1)*n] * x[k-1];
    }
  }

  return;
}
/******************************************************************************/

void r8plu_to_r8mat ( int n, int pivot[], double lu[], double a[] )

/******************************************************************************/
/*
  Purpose:

    R8PLU_TO_R8MAT recovers the matrix A that was factored by R8MAT_TO_R8PLU.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.
    N must be positive.

    Input, int PIVOT[N], the pivot vector computed by R8MAT_TO_R8PLU.

    Input, double LU[N*N], the matrix factors computed by R8MAT_TO_R8PLU.

    Output, double A[N*N], the matrix whose factors are represented by
    LU and PIVOT.
*/
{
  int i;
  int j;
  int k;
  double temp;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= n; i++ )
    {
      for ( k = 1; k <= i-1; k++ )
      {
        a[k-1+(j-1)*n] = a[k-1+(j-1)*n] + lu[k-1+(i-1)*n] * a[i-1+(j-1)*n];
      }
      a[i-1+(j-1)*n] = lu[i-1+(i-1)*n] * a[i-1+(j-1)*n];
    }
/*
  B = PL * Y = PL * U * X = A * x.
*/
    for ( i = n-1; 1 <= i; i-- )
    {
      for ( k = i+1; k <= n; k++ )
      {
        a[k-1+(j-1)*n] = a[k-1+(j-1)*n] - lu[k-1+(i-1)*n] * a[i-1+(j-1)*n];
      }

      k = pivot[i-1];

      if ( k != i )
      {
        temp           = a[k-1+(j-1)*n];
        a[k-1+(j-1)*n] = a[i-1+(j-1)*n];
        a[i-1+(j-1)*n] = temp;
      }
    }
  }

  return;
}
/******************************************************************************/

int r8poly_degree ( int na, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8POLY_DEGREE returns the degree of a polynomial.

  Discussion:

    The degree of a polynomial is the index of the highest power
    of X with a nonzero coefficient.

    The degree of a constant polynomial is 0.  The degree of the
    zero polynomial is debatable, but this routine returns the
    degree as 0.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int NA, the dimension of A.

    Input, double A[NA+1], the coefficients of the polynomials.

    Output, int R8POLY_DEGREE, the degree of A.
*/
{
  int degree;

  degree = na;

  while ( 0 < degree )
  {
    if ( a[degree] != 0.0 )
    {
      return degree;
    }
    degree = degree - 1;
  }

  return degree;
}
/******************************************************************************/

double *r8poly_deriv ( int n, double c[], int p )

/******************************************************************************/
/*
  Purpose:

    R8POLY_DERIV returns the derivative of a polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the degree of the polynomial.

    Input, double C[N+1], the polynomial coefficients.
    C[I] is the coefficient of X^I.

    Input, int P, the order of the derivative.
    0 means no derivative is taken.
    1 means first derivative,
    2 means second derivative and so on.
    Values of P less than 0 are meaningless.  Values of P greater
    than N are meaningful, but the code will behave as though the
    value of P was N+1.

    Output, double R8POLY_DERIV CP[N-P+1], the polynomial coefficients of
    the derivative.
*/
{
  double *cp;
  double *cp_temp;
  int d;
  int i;

  if ( n < p )
  {
    return NULL;
  }
  cp_temp = r8vec_copy_new ( n + 1, c );

  for ( d = 1; d <= p; d++ )
  {
    for ( i = 0; i <= n-d; i++ )
    {
      cp_temp[i] = ( double ) ( i + 1 ) * cp_temp[i+1];
    }
    cp_temp[n-d+1] = 0.0;
  }

  cp = r8vec_copy_new ( n - p + 1, cp_temp );

  free ( cp_temp );

  return cp;
}
/******************************************************************************/

double r8poly_lagrange_0 ( int npol, double xpol[], double xval )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_0 evaluates the Lagrange factor at a point.

  Discussion:

    W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )

  Discussion:

    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
    polynomial L(IPOL)(X), has the property:

      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )

    and may be expressed as:

      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int NPOL, the number of abscissas.
    NPOL must be at least 1.

    Input, double XPOL[NPOL], the abscissas, which should be distinct.

    Input, double XVAL, the point at which the Lagrange factor is to be
    evaluated.

    Output, double R8POLY_LAGRANGE_0, the value of the Lagrange factor at XVAL.
*/
{
  int i;
  double wval;

  wval = 1.0;
  for ( i = 0; i < npol; i++ )
  {
    wval = wval * ( xval - xpol[i] );
  }

  return wval;
}
/******************************************************************************/

double r8poly_lagrange_1 ( int npol, double xpol[], double xval )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_1 evaluates the first derivative of the Lagrange factor.

  Discussion:

    W(XPOL(1:NPOL))(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )

    W'(XPOL(1:NPOL))(X)
      = Sum ( 1 <= J <= NPOL ) Product ( I /= J ) ( X - XPOL(I) )

    We also have the recursion:

      W'(XPOL(1:NPOL))(X) = d/dX ( ( X - XPOL(NPOL) ) * W(XPOL(1:NPOL-1))(X) )
                    = W(XPOL(1:NPOL-1))(X)
                    + ( X - XPOL(NPOL) ) * W'(XPOL(1:NPOL-1))(X)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int NPOL, the number of abscissas.

    Input, double XPOL[NPOL], the abscissas, which should be distinct.

    Input, double XVAL, the point at which the Lagrange factor is to be
    evaluated.

    Output, double R8POLY_LAGRANGE_1, the derivative of W with respect to XVAL.
*/
{
  double dwdx;
  int i;
  double w;

  dwdx = 0.0;
  w = 1.0;

  for ( i = 0; i < npol; i++ )
  {
    dwdx = w + ( xval - xpol[i] ) * dwdx;
    w = w * ( xval - xpol[i] );
  }

  return dwdx;
}
/******************************************************************************/

double r8poly_lagrange_2 ( int npol, double xpol[], double xval )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_2 evaluates the second derivative of the Lagrange factor.

  Discussion:

    W(X)  = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )

    W'(X) = Sum ( 1 <= J <= NPOL )
            Product ( I /= J ) ( X - XPOL(I) )

    W"(X) = Sum ( 1 <= K <= NPOL )
            Sum ( J =/ K )
            Product ( I /= K, J ) ( X - XPOL(I) )

    For a set of points XPOL(I), 1 <= I <= NPOL, the IPOL-th Lagrange basis
    polynomial L(IPOL)(X), has the property:

      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )

    and may be expressed as:

      L(IPOL)(X) = W(X) / ( ( X - XPOL(IPOL) ) * W'(XPOL(IPOL)) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 January 2004

  Author:

    John Burkardt

  Parameters:

    Input, int NPOL, the number of abscissas.
    NPOL must be at least 1.

    Input, double XPOL[NPOL], the abscissas, which should be distinct.

    Input, double XVAL, the point at which the Lagrange factor is to be
    evaluated.

    Output, double R8POLY_LAGRANGE_2, the second derivative of W with respect to XVAL.
*/
{
  double dw2dx2;
  int i;
  int j;
  int k;
  double term;

  dw2dx2 = 0.0;

  for ( k = 0; k < npol; k++ )
  {
    for ( j = 0; j < npol; j++ )
    {
      if ( j != k )
      {
        term = 1.0;
        for ( i = 0; i < npol; i++ )
        {
          if ( i != j && i != k )
          {
            term = term * ( xval - xpol[i] );
          }
        }
        dw2dx2 = dw2dx2 + term;
      }
    }
  }

  return dw2dx2;
}
/******************************************************************************/

double *r8poly_lagrange_coef ( int npol, int ipol, double xpol[] )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_COEF returns the coefficients of a Lagrange polynomial.

  Discussion:

    Given NPOL distinct abscissas, XPOL(*), the IPOL-th Lagrange
    polynomial P(IPOL)(X) is defined as the polynomial of degree
    NPOL - 1 which is 1 at XPOL(IPOL) and 0 at the NPOL - 1 other
    abscissas.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int NPOL, the number of abscissas.
    NPOL must be at least 1.

    Input, int IPOL, the index of the polynomial to evaluate.
    IPOL must be between 1 and NPOL.

    Input, double XPOL[NPOL], the abscissas of the Lagrange polynomials.
    The entries in XPOL must be distinct.

    Output, double R8POLY_LAGRANGE_COEF[NPOL], the polynomial coefficients
    of the IPOL-th Lagrange polynomial.
*/
{
  int i;
  int index;
  int j;
  double *pcof;
/*
  Make sure IPOL is legal.
*/
  if ( ipol < 1 || npol < ipol )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8POLY_LAGRANGE_COEF - Fatal error!\n" );
    fprintf ( stderr, "  1 <= IPOL <= NPOL is required.\n" );
    fprintf ( stderr, "  but IPOL = %d\n", ipol );
    fprintf ( stderr, "  and NPOL = %d\n", npol );
    exit ( 1 );
  }
/*
  Check that the abscissas are distinct.
*/
  if ( !r8vec_distinct ( npol, xpol ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8POLY_LAGRANGE_COEF - Fatal error!\n" );
    fprintf ( stderr, "  Two entries of XPOL are equal:\n" );
    exit ( 1 );
  }

  pcof = ( double * ) malloc ( npol * sizeof ( double ) );

  pcof[0] = 1.0;
  for ( i = 1; i < npol; i++ )
  {
    pcof[i] = 0.0;
  }

  index = 0;

  for ( i = 1; i <= npol; i++ )
  {
    if ( i != ipol )
    {
      index = index + 1;

      for ( j = index; 0 <= j; j-- )
      {
        pcof[j] = - xpol[i-1] * pcof[j] / ( xpol[ipol-1] - xpol[i-1] );

        if ( 0 < j )
        {
          pcof[j] = pcof[j] + pcof[j-1] / ( xpol[ipol-1] - xpol[i-1] );
        }
      }
    }
  }

  return pcof;
}
/******************************************************************************/

void r8poly_lagrange_factor ( int npol, double xpol[], double xval,
  double *wval, double *dwdx )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_FACTOR evaluates the polynomial Lagrange factor at a point.

  Discussion:

    Suppose F(X) is at least N times continuously differentiable in the
    interval [A,B].  Pick NPOL distinct points XPOL(I) in [A,B] and compute
    the interpolating polynomial P(X) of order NPOL ( and degree NPOL-1)
    which passes through all the points ( XPOL(I), F(XPOL(I)) ).
    Then in the interval [A,B], the maximum error

      abs ( F(X) - P(X) )

    is bounded by:

      C * FNMAX * W(X)

    where

      C is a constant,
      FNMAX is the maximum value of the NPOL-th derivative of F in [A,B],
      W(X) is the Lagrange factor.

    Thus, the value of W(X) is useful as part of an estimated bound
    for the interpolation error.

    The formula is:

      W(X) = Product ( 1 <= I <= NPOL ) ( X - XPOL(I) )

    Note that the Chebyshev abscissas have the property that they minimize
    the value of W(X) over the interval [A,B].  Hence, if the abscissas may
    be chosen arbitrarily, the Chebyshev abscissas have this advantage over
    other choices.

    For a set of points XPOL[I], 0 <= I <= NPOL-1, the IPOL-th Lagrange basis
    polynomial L(IPOL)(X), has the property:

      L(IPOL)( XPOL(J) ) = delta ( IPOL, J )

    and may be expressed as:

      L(IPOL)(X) = W(X) / ( ( X - XPOL[IPOL] ) * W'(XPOL[IPOL]) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int NPOL, the number of abscissas.
    NPOL must be at least 1.

    Input, double XPOL[NPOL], the abscissas, which should be distinct.

    Input, double XVAL, the point at which the Lagrange factor is to be evaluated.

    Output, double *WVAL, the value of the Lagrange factor at XVAL.

    Output, double *DWDX, the derivative of W with respect to XVAL.
*/
{
  int i;
  int j;
  double term;

  *wval = 1.0;
  for ( i = 0; i < npol; i++ )
  {
    *wval = *wval * ( xval - xpol[i] );
  }

  *dwdx = 0.0;

  for ( i = 0; i < npol; i++ )
  {
    term = 1.0;

    for ( j = 0; j < npol; j++ )
    {
      if ( i != j )
      {
        term = term * ( xval - xpol[j] );
      }
    }
    *dwdx = *dwdx + term;
  }

  return;
}
/******************************************************************************/

int r8poly_lagrange_val ( int npol, int ipol, double xpol[], double xval,
  double *pval, double *dpdx )

/******************************************************************************/
/*
  Purpose:

    R8POLY_LAGRANGE_VAL evaluates the IPOL-th Lagrange polynomial.

  Discussion:

    Given NPOL distinct abscissas, XPOL[*], the IPOL-th Lagrange
    polynomial P(IPOL)(X) is defined as the polynomial of degree
    NPOL - 1 which is 1 at XPOL[IPOL] and 0 at the NPOL - 1 other
    abscissas.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int NPOL, the number of abscissas.
    NPOL must be at least 1.

    Input, int IPOL, the index of the polynomial to evaluate.
    IPOL must be between 0 and NPOL-1.

    Input, double XPOL[NPOL], the abscissas of the Lagrange polynomials.
    The entries in XPOL must be distinct.

    Input, double XVAL, the point at which the IPOL-th Lagrange polynomial
    is to be evaluated.

    Output, double *PVAL, the value of the IPOL-th Lagrange polynomial at XVAL.

    Output, double *DPDX, the derivative of the IPOL-th Lagrange polynomial at XVAL.

    Output, int R8POLY_LAGRANGE_VAL, 0 if no error.
*/
{
  int i;
  int j;
  double p2;
/*
  Make sure IPOL is legal.
*/
  if ( ipol < 0 || npol-1 < ipol )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8POLY_LAGRANGE_VAL - Fatal error!\n" );
    fprintf ( stderr, "  0 <= IPOL <= NPOL-1 is required.\n" );
    exit ( 1 );
  }
/*
  Check that the abscissas are distinct.
*/
  for ( i = 1; i < npol; i++ )
  {
    for ( j = 0; j < i; j++ )
    {
      if ( xpol[i] == xpol[j] )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "R8POLY_LAGRANGE_VAL - Fatal error!\n" );
        fprintf ( stderr, "  Two entries of XPOL are equal:\n" );
        fprintf ( stderr, "  XPOL(%d) = %g.\n", i, xpol[i] );
        fprintf ( stderr, "  XPOL(%d) = %g.\n", j, xpol[j] );
        exit ( 1 );
      }
    }
  }
/*
  Evaluate the polynomial.
*/
  *pval = 1.0;

  for ( i = 0; i < npol; i++ )
  {
    if ( i != ipol )
    {
      *pval = *pval * ( xval - xpol[i] ) / ( xpol[ipol] - xpol[i] );
    }
  }
/*
  Evaluate the derivative, which can be found by summing up the result
  of differentiating one factor at a time, successively.
*/
  *dpdx = 0.0;

  for ( i = 0; i < npol; i++ )
  {
    if ( i != ipol )
    {
      p2 = 1.0;

      for ( j = 0; j < npol; j++ )
      {
        if ( j == i )
        {
          p2 = p2                      / ( xpol[ipol] - xpol[j] );
        }
        else if ( j != ipol )
        {
          p2 = p2 * ( xval - xpol[j] ) / ( xpol[ipol] - xpol[j] );
        }
      }
      *dpdx = *dpdx + p2;
    }
  }

  return 0;
}
/******************************************************************************/

int r8poly_order ( int na, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8POLY_ORDER returns the order of a polynomial.

  Discussion:

    The order of a polynomial is one more than the degree.

    The order of a constant polynomial is 1.  The order of the
    zero polynomial is debatable, but this routine returns the
    order as 1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int NA, the dimension of A.

    Input, double A[NA+1], the coefficients of the polynomials.

    Output, int R8POLY_ORDER, the order of A.
*/
{
  int order;

  order = na + 1;

  while ( 1 < order )
  {
    if ( a[order-1] != 0.0 )
    {
      return order;
    }
    order = order - 1;
  }

  return order;
}
/******************************************************************************/

void r8poly_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8POLY_PRINT prints out a polynomial.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of A.

    Input, double A[N+1], the polynomial coefficients.
    A(0) is the constant term and
    A(N) is the coefficient of X**N.

    Input, char *TITLE, a title.
*/
{
  int i;
  double mag;
  int n2;
  char plus_minus;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );

  n2 = r8poly_degree ( n, a );

  if ( n2 <= 0 )
  {
    printf ( "  p(x) = 0\n" );
    return;
  }

  if ( a[n2] < 0.0 )
  {
    plus_minus = '-';
  }
  else
  {
    plus_minus = ' ';
  }

  mag = r8_abs ( a[n2] );

  if ( 2 <= n2 )
  {
    printf ( "  p(x) = %c%g * x^%d\n", plus_minus, mag, n2 );
  }
  else if ( n2 == 1 )
  {
    printf ( "  p(x) = %c%g * x\n", plus_minus, mag );
  }
  else if ( n2 == 0 )
  {
    printf ( "  p(x) = %c%g\n", plus_minus, mag );
  }

  for ( i = n2 - 1; 0 <= i; i-- )
  {
    if ( a[i] < 0.0 )
    {
      plus_minus = '-';
    }
    else
    {
      plus_minus = '+';
    }

    mag = r8_abs ( a[i] );

    if ( mag != 0.0 )
    {
      if ( 2 <= i )
      {
        printf ( "         %c%g * x^%d\n", plus_minus, mag, i );
      }
      else if ( i == 1 )
      {
        printf ( "         %c%g * x\n", plus_minus, mag );
      }
      else if ( i == 0 )
      {
        printf ( "         %c%g\n", plus_minus, mag );
      }
    }
  }
  return;
}
/******************************************************************************/

void r8poly_shift ( double scale, double shift, int n, double poly_cof[] )

/******************************************************************************/
/*
  Purpose:

    R8POLY_SHIFT adjusts the coefficients of a polynomial for a new argument.

  Discussion:

    Assuming P(X) is a polynomial in the argument X, of the form:

      P(X) =
          C(N) * X^N
        + ...
        + C(1) * X
        + C(0),

    and that Z is related to X by the formula:

      Z = SCALE * X + SHIFT

    then this routine computes coefficients C for the polynomial Q(Z):

      Q(Z) =
          C(N) * Z^N
        + ...
        + C(1) * Z
        + C(0)

    so that:

      Q(Z(X)) = P(X)

  Example:

    P(X) = 2 * X^2 - X + 6

    Z = 2.0 * X + 3.0

    Q(Z) = 0.5 *         Z^2 -  3.5 * Z + 12

    Q(Z(X)) = 0.5 * ( 4.0 * X^2 + 12.0 * X +  9 )
            - 3.5 * (              2.0 * X +  3 )
                                           + 12

            = 2.0         * X^2 -  1.0 * X +  6

            = P(X)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 August 2011

  Reference:

    Press, Flannery, Teukolsky, Vetterling,
    Numerical Recipes: The Art of Scientific Computing,
    Cambridge University Press.

  Parameters:

    Input, double SHIFT, SCALE, the shift and scale applied to X,
    so that Z = SCALE * X + SHIFT.

    Input, int N, the number of coefficients.

    Input/output, double POLY_COF[N+1].
    On input, the coefficient array in terms of the X variable.
    On output, the coefficient array in terms of the Z variable.
*/
{
  int i;
  int j;

  for ( i = 1; i <= n; i++ )
  {
    for ( j = i; j <= n; j++ )
    {
      poly_cof[j] = poly_cof[j] / scale;
    }
  }

  for ( i = 0; i <= n - 1; i++ )
  {
    for ( j = n - 1; i <= j; j-- )
    {
      poly_cof[j] = poly_cof[j] - shift * poly_cof[j+1];
    }
  }

  return;
}
/******************************************************************************/

double *r8poly_value ( int m, double c[], int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8POLY_VALUE evaluates a polynomial.

  Discussion:

    The polynomial 

      p(x) = c1 + c2 * x + c3 * x^2 + ... + cm * x^(m-1)

    is to be evaluated at the vector of values X.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the degree.

    Input, double C[M+1], the polynomial coefficients.  
    C[0] is the constant term.

    Input, int N, the number of evaluation points.

    Input, double X[N], the evaluation points.

    Output, double R8POLY_VALUE[N], the value of the polynomial at the 
    evaluation points.
*/
{
  int i;
  int j;
  double *p;

  p = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    p[j] = c[m];
  }

  for ( i = m - 1; 0 <= i; i-- )
  {
    for ( j = 0; j < n; j++ )
    {
      p[j] = p[j] * x[j] + c[i];
    }
  }
  return p;
}
/******************************************************************************/

double r8poly_value_horner ( int n, double c[], double x )

/******************************************************************************/
/*
  Purpose:

    R8POLY_VALUE_HORNER evaluates a polynomial using Horner's method.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 September 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of C.

    Input, double C[N+1], the polynomial coefficients.
    C(I) is the coefficient of X^I.

    Input, double X, the point at which the polynomial is
    to be evaluated.

    Output, double R8POLY_VALUE_HORNER, the value of the polynomial at X./
*/
{
  int i;
  double value;

  value = c[n];
  for ( i = n-1; 0 <= i; i-- )
  {
    value = value * x + c[i];
  }

  return value;
}
/******************************************************************************/

double r8poly_value_old ( int n, double a[], double x )

/******************************************************************************/
/*
  Purpose:

    R8POLY_VALUE_OLD evaluates an R8POLY.

  Discussion:

    For sanity's sake, the value of N indicates the NUMBER of
    coefficients, or more precisely, the ORDER of the polynomial,
    rather than the DEGREE of the polynomial.  The two quantities
    differ by 1, but cause a great deal of confusion.

    Given N and A, the form of the polynomial is:

      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 March 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the polynomial.

    Input, double A[N], the coefficients of the polynomial.
    A[0] is the constant term.

    Input, double X, the point at which the polynomial is to be evaluated.

    Output, double R8POLY_VALUE, the value of the polynomial at X.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    value = value * x + a[i];
  }

  return value;
}
/******************************************************************************/

double *r8poly_value_2d ( int m, double c[], int n, double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    R8POLY_VALUE_2D evaluates a polynomial in 2 variables, X and Y.

  Discussion:

    We assume the polynomial is of total degree M, and has the form:

      p(x,y) = c00 
             + c10 * x                + c01 * y
             + c20 * x^2   + c11 * xy + c02 * y^2
             + ...
             + cm0 * x^(m) + ...      + c0m * y^m.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the degree of the polynomial.

    Input, double C[T(M+1)], the polynomial coefficients.  
    C[0] is the constant term.  T(M+1) is the M+1-th triangular number.
    The coefficients are stored consistent with the following ordering
    of monomials: 1, X, Y, X^2, XY, Y^2, X^3, X^2Y, XY^2, Y^3, X^4, ...

    Input, int N, the number of evaluation points.

    Input, double X[N], Y[N], the evaluation points.

    Output, double R8POLY_VALUE_2D[N], the value of the polynomial at the 
    evaluation points.
*/
{
  int ex;
  int ey;
  int i;
  int j;
  double *p;
  int s;

  p = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    p[i] = 0.0;
  }

  j = 0;
  for ( s = 0; s <= m; s++ )
  {
    for ( ex = s; 0 <= ex; ex-- )
    {
      ey = s - ex;
      for ( i = 0; i < n; i++ )
      {
        p[i] = p[i] + c[j] * pow ( x[i], ex ) * pow ( y[i], ey );
      }
      j = j + 1;
    }
  }
  return p;
}
/******************************************************************************/

int r8poly2_ex ( double x1, double y1, double x2, double y2, double x3,
  double y3, double *x, double *y )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_EX finds the extremal point of a parabola determined by three points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
    on the parabola.  X1, X2 and X3 must be distinct.

    Output, double *X, *Y, the X coordinate of the extremal point of the
    parabola, and the value of the parabola at that point.

    Output, int R8POLY2_EX, error flag.
    0, no error.
    1, two of the X values are equal.
    2, the data lies on a straight line; there is no finite extremal
    point.
    3, the data lies on a horizontal line; every point is "extremal".
*/
{
  double bot;

  *x = 0.0;
  *y = 0.0;

  if ( x1 == x2 || x2 == x3 || x3 == x1 )
  {
    return 1;
  }

  if ( y1 == y2 && y2 == y3 && y3 == y1 )
  {
    *x = x1;
    *y = y1;
    return 3;
  }

  bot = ( x2 - x3 ) * y1 + ( x3 - x1 ) * y2 + ( x1 - x2 ) * y3;

  if ( bot == 0.0 )
  {
    return 2;
  }

  *x = 0.5 * (
      x1 * x1 * ( y3 - y2 )
    + x2 * x2 * ( y1 - y3 )
    + x3 * x3 * ( y2 - y1 ) ) /
    ( ( x2 - x3 ) * y1 + ( x3 - x1 ) * y2 + ( x1 - x2 ) * y3 );

  *y = - (
      ( *x - x2 ) * ( *x - x3 ) * ( x2 - x3 ) * y1
    + ( *x - x1 ) * ( *x - x3 ) * ( x3 - x1 ) * y2
    + ( *x - x1 ) * ( *x - x2 ) * ( x1 - x2 ) * y3 ) /
    ( ( x1 - x2 ) * ( x2 - x3 ) * ( x3 - x1 ) );

  return 0;
}
/******************************************************************************/

int r8poly2_ex2 ( double x1, double y1, double x2, double y2, double x3,
  double y3, double *x, double *y, double *a, double *b, double *c )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_EX2 finds the extremal point of a parabola determined by three points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
    on the parabola.  X1, X2 and X3 must be distinct.

    Output, double *X, *Y, the X coordinate of the extremal point of the
    parabola, and the value of the parabola at that point.

    Output, double *A, *B, *C, the coefficients that define the parabola:
    P(X) = A * X**2 + B * X + C.

    Output, int R8POLY2_EX2, error flag.
    0, no error.
    1, two of the X values are equal.
    2, the data lies on a straight line; there is no finite extremal
    point.
    3, the data lies on a horizontal line; any point is an "extremal point".
*/
{
  double v[3*3];
  double *w;

  *a = 0.0;
  *b = 0.0;
  *c = 0.0;
  *x = 0.0;
  *y = 0.0;

  if ( x1 == x2 || x2 == x3 || x3 == x1 )
  {
    return 1;
  }

  if ( y1 == y2 && y2 == y3 && y3 == y1 )
  {
    *x = x1;
    *y = y1;
    return 3;
  }
/*
  Set up the Vandermonde matrix.
*/
  v[0+0*3] = 1.0;
  v[0+1*3] = x1;
  v[0+2*3] = x1 * x1;

  v[1+0*3] = 1.0;
  v[1+1*3] = x2;
  v[1+2*3] = x2 * x2;

  v[2+0*3] = 1.0;
  v[2+1*3] = x3;
  v[2+2*3] = x3 * x3;
/*
  Get the inverse.
*/
  w = r8mat_inverse_3d ( v );
/*
  Compute the parabolic coefficients.
*/
  *c = w[0+0*3] * y1 + w[0+1*3] * y2 + w[0+2*3] * y3;
  *b = w[1+0*3] * y1 + w[1+1*3] * y2 + w[1+2*3] * y3;
  *a = w[2+0*3] * y1 + w[2+1*3] * y2 + w[2+2*3] * y3;
/*
  Determine the extremal point.
*/
  if ( *a == 0.0 )
  {
    return 2;
  }

  *x = - *b / ( 2.0 * *a );
  *y = *a * *x * *x + *b * *x + *c;

  return 0;
}
/******************************************************************************/

void r8poly2_rroot ( double a, double b, double c, double *r1, double *r2 )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_RROOT returns the real parts of the roots of a quadratic polynomial.

  Example:

    A    B    C       roots              R1   R2
   --   --   --     ------------------   --   --
    1   -4    3     1          3          1    3
    1    0    4     2*i      - 2*i        0    0
    2   -6    5     3 +   i    3 -   i    3    3

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, double A, B, C, the coefficients of the quadratic
    polynomial A * X**2 + B * X + C = 0 whose roots are desired.
    A must not be zero.

    Output, double *R1, *R2, the real parts of the roots
    of the polynomial.
*/
{
  double disc;
  double q;

  if ( a == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8POLY2_RROOT - Fatal error!\n" );
    fprintf ( stderr, "  The coefficient A is zero.\n" );
    exit ( 1 );
  }

  disc = b * b - 4.0 * a * c;
  disc = r8_max ( disc, 0.0 );

  q = ( b + r8_sign ( b ) * sqrt ( disc ) );
  *r1 = -0.5 * q / a;
  *r2 = -2.0 * c / q;

  return;
}
/******************************************************************************/

void r8poly2_val ( double x1, double y1, double x2, double y2,
  double x3, double y3, double x, double *y, double *yp, double *ypp )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_VAL evaluates a parabola defined by three data values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, double X1, Y1, X2, Y2, X3, Y3, three pairs of data values.
    If the X values are distinct, then all the Y values represent
    actual values of the parabola.

    Three special cases are allowed:

      X1 = X2 =/= X3: Y2 is the derivative at X1;
      X1 =/= X2 = X3: Y3 is the derivative at X3;
      X1 = X2 = X3:   Y2 is the derivative at X1, and
                      Y3 is the second derivative at X1.

    Input, double X, an abscissa at which the parabola is to be
    evaluated.

    Output, double *Y, *YP, *YPP, the values of the parabola and
    its first and second derivatives at X.
*/
{
  int distinct;
  double dif1;
  double dif2;
  double temp;
/*
  If any X's are equal, put them and the Y data first.
*/
  if ( x1 == x2 && x2 == x3 )
  {
    distinct = 1;
  }
  else if ( x1 == x2 )
  {
    distinct = 2;
  }
  else if ( x1 == x3 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8POLY2_VAL - Fatal error!\n" );
    fprintf ( stderr, "  X1 = X3 =/= X2.\n" );
    return;
  }
  else if ( x2 == x3 )
  {
    distinct = 2;
    temp = x1;
    x1 = x3;
    x3 = temp;
    temp = y1;
    y1 = y2;
    y2 = y3;
    y3 = y1;
  }
  else
  {
    distinct = 3;
  }
/*
  Set up the coefficients.
*/
  if ( distinct == 1 )
  {
    dif1 = y2;
    dif2 = 0.5 * y3;
  }
  else if ( distinct == 2 )
  {
    dif1 = y2;
    dif2 = ( ( y3 - y1 ) / ( x3 - x1 )
             - y2 ) / ( x3 - x2 );
  }
  else if ( distinct == 3 )
  {
    dif1 = ( y2 - y1 ) / ( x2 - x1 );
    dif2 = ( ( y3 - y1 ) / ( x3 - x1 )
           - ( y2 - y1 ) / ( x2 - x1 ) ) / ( x3 - x2 );
  }
/*
  Evaluate.
*/
  *y = y1 + ( x - x1 ) * dif1 + ( x - x1 ) * ( x - x2 ) * dif2;
  *yp = dif1 + ( 2.0 * x - x1 - x2 ) * dif2;
  *ypp = 2.0 * dif2;

  return;
}
/******************************************************************************/

void r8poly2_val2 ( int ndata, double tdata[],
  double ydata[], int left, double tval, double *yval )

/******************************************************************************/
/*
  Purpose:

    R8POLY2_VAL2 evaluates a parabolic function through 3 points in a table.

  Discussion:

    This routine is a utility routine used by OVERHAUSER_SPLINE_VAL.
    It constructs the parabolic interpolant through the data in
    3 consecutive entries of a table and evaluates this interpolant
    at a given abscissa value.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int NDATA, the number of data points.
    NDATA must be at least 3.

    Input, double TDATA[NDATA], the abscissas of the data points.  The
    values in TDATA must be in strictly ascending order.

    Input, double YDATA[NDATA], the data points corresponding to
    the abscissas.

    Input, int LEFT, the location of the first of the three
    consecutive data points through which the parabolic interpolant
    must pass.  0 <= LEFT <= NDATA - 3.

    Input, double TVAL, the value of T at which the parabolic interpolant
    is to be evaluated.  Normally, TDATA[0] <= TVAL <= T[NDATA-1], and
    the data will be interpolated.  For TVAL outside this range,
    extrapolation will be used.

    Output, double *YVAL, the value of the parabolic interpolant
    at TVAL.
*/
{
  double dif1;
  double dif2;
  double t1;
  double t2;
  double t3;
  double y1;
  double y2;
  double y3;
/*
  Check.
*/
  if ( left < 0 || ndata-3 < left )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "RPOLY2_VAL2 - Fatal error!\n" );
    fprintf ( stderr, "  LEFT < 0 or NDATA-3 < LEFT.\n" );
    exit ( 1 );
  }
/*
  Copy out the three abscissas.
*/
  t1 = tdata[left];
  t2 = tdata[left+1];
  t3 = tdata[left+2];

  if ( t2 <= t1 || t3 <= t2 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "RPOLY2_VAL2 - Fatal error!\n" );
    fprintf ( stderr, "  T2 <= T1 or T3 <= T2.\n" );
    fprintf ( stderr, "  T1 = %g\n", t1 );
    fprintf ( stderr, "  T2 = %g\n", t2 );
    fprintf ( stderr, "  T3 = %g\n", t3 );
    exit ( 1 );
  }
/*
  Construct and evaluate a parabolic interpolant for the data.
*/
  y1 = ydata[left];
  y2 = ydata[left+1];
  y3 = ydata[left+2];

  dif1 = ( y2 - y1 ) / ( t2 - t1 );
  dif2 =
    ( ( y3 - y1 ) / ( t3 - t1 )
    - ( y2 - y1 ) / ( t2 - t1 ) ) / ( t3 - t2 );

  *yval = y1 + ( tval - t1 ) * ( dif1 + ( tval - t2 ) * dif2 );

  return;
}
/******************************************************************************/

int r8r8_compare ( double x1, double y1, double x2, double y2 )

/******************************************************************************/
/*
  Purpose:

    R8R8_COMPARE compares two R8R8's.

  Discussion:

    An R8R8 is simply a pair of R8 values, stored separately.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X1, Y1, the first vector.

    Input, double X2, Y2, the second vector.

    Output, int R8R8_COMPARE:
    -1, (X1,Y1) < (X2,Y2);
     0, (X1,Y1) = (X2,Y2);
    +1, (X1,Y1) > (X2,Y2).
*/
{
  int value;

  if ( x1 < x2 )
  {
    value = -1;
  }
  else if ( x2 < x1 )
  {
    value = +1;
  }
  else if ( y1 < y2 )
  {
    value = -1;
  }
  else if ( y2 < y1 )
  {
    value = +1;
  }
  else
  {
    value = 0;
  }

  return value;
}
/******************************************************************************/

void r8r8_print ( double a1, double a2, char *title )

/******************************************************************************/
/*
  Purpose:

    R8R8_PRINT prints an R8R8.

  Discussion:

    An R8R8 is a pair of R8 values, regarded as a single item.

    A format is used which suggests a coordinate pair:

  Example:

    Center : ( 1.23, 7.45 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double A1, A2, the coordinates of the vector.

    Input, char *TITLE, a title.
*/
{
  fprintf ( stdout, "%s: ( %g, %g )\n", title, a1, a2 );

  return;
}
/******************************************************************************/

int r8r8r8_compare ( double x1, double y1, double z1, double x2, double y2,
  double z2 )

/******************************************************************************/
/*
  Purpose:

    R8R8R8_COMPARE compares two R8R8R8's.

  Discussion:

    An R8R8R8 is simply 3 R8 values, stored as scalars.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double X1, Y1, Z1, the first vector.

    Input, double X2, Y2, Z2, the second vector.

    Output, int R8R8R8_COMPARE:
    -1, (X1,Y1,Z1) < (X2,Y2,Z2);
     0, (X1,Y1,Z1) = (X2,Y2,Z2);
    +1, (X1,Y1,Z1) > (X2,Y2,Z2).
*/
{
  int value;

  if ( x1 < x2 )
  {
    value = -1;
  }
  else if ( x2 < x1 )
  {
    value = +1;
  }
  else if ( y1 < y2 )
  {
    value = -1;
  }
  else if ( y2 < y1 )
  {
    value = +1;
  }
  else if ( z1 < z2 )
  {
    value = -1;
  }
  else if ( z2 < z1 )
  {
    value = +1;
  }
  else
  {
    value = 0;
  }

  return value;
}
/******************************************************************************/

void r8r8r8vec_index_insert_unique ( int maxn, int *n, double x[], double y[],
  double z[], int indx[], double xval, double yval, double zval, int *ival,
  int *ierror )

/******************************************************************************/
/*
  Purpose:

    R8R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8R8 value in an indexed sorted list.

  Discussion:

    If the input value does not occur in the current list, it is added,
    and N, X, Y, Z and INDX are updated.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int MAXN, the maximum size of the list.

    Input/output, int *N, the size of the list.

    Input/output, double X[N], Y[N], Z[N], the R8R8R8 vector.

    Input/output, int INDX[N], the sort index of the list.

    Input, double XVAL, YVAL, ZVAL, the value to be inserted
    if it is not already in the list.

    Output, int *IVAL, the index in X, Y, Z corresponding to the
    value XVAL, YVAL, ZVAL.

    Output, int *IERROR, 0 for no error, 1 if an error occurred.
*/
{
  int equal;
  int i;
  int less;
  int more;

  *ierror = 0;

  if ( *n <= 0 )
  {
    if ( maxn <= 0 )
    {
      *ierror = 1;
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!\n" );
      fprintf ( stderr, "  Not enough space to store new data.\n" );
      return;
    }
    *n = 1;
    x[0] = xval;
    y[0] = yval;
    z[0] = zval;
    indx[0] = 1;
    *ival = 1;
    return;
  }
/*
  Does ( XVAL, YVAL, ZVAL ) already occur in ( X, Y, Z)?
*/
  r8r8r8vec_index_search ( *n, x, y, z, indx, xval, yval, zval,
    &less, &equal, &more );

  if ( equal == 0 )
  {
    if ( maxn <= *n )
    {
      *ierror = 1;
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!\n" );
      fprintf ( stderr, "  Not enough space to store new data.\n" );
      return;
    }

    x[*n] = xval;
    y[*n] = yval;
    z[*n] = zval;
    *ival = *n + 1;
    for ( i = *n-1; more-1 <= i; i-- )
    {
      indx[i+1] = indx[i];
    }
    indx[more-1] = *n + 1;
    *n = *n + 1;
  }
  else
  {
    *ival = indx[equal-1];
  }

  return;
}
/******************************************************************************/

void r8r8r8vec_index_search ( int n, double x[], double y[], double z[],
  int indx[], double xval, double yval, double zval, int *less, int *equal,
  int *more )

/******************************************************************************/
/*
  Purpose:

    R8R8R8VEC_INDEX_SEARCH searches for an R8R8R8 value in an indexed sorted list.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the list.

    Input, double X[N], Y[N], Z[N], the list.

    Input, int INDX[N], the sort index of the list.

    Input, double XVAL, YVAL, ZVAL, the value to be sought.

    Output, int *LESS, *EQUAL, *MORE, the indexes in INDX of the
    entries of X that are just less than, equal to, and just greater
    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
    is the greatest entry of X, then MORE is N+1.
*/
{
  int compare;
  int hi;
  int lo;
  int mid;
  double xhi;
  double xlo;
  double xmid;
  double yhi;
  double ylo;
  double ymid;
  double zhi;
  double zlo;
  double zmid;

  if ( n <= 0 )
  {
    *less = 0;
    *equal = 0;
    *more = 0;
    return;
  }

  lo = 1;
  hi = n;

  xlo = x[indx[lo-1]-1];
  ylo = y[indx[lo-1]-1];
  zlo = z[indx[lo-1]-1];

  xhi = x[indx[hi-1]-1];
  yhi = y[indx[hi-1]-1];
  zhi = z[indx[hi-1]-1];

  compare = r8r8r8_compare ( xval, yval, zval, xlo, ylo, zlo );

  if ( compare == -1 )
  {
    *less = 0;
    *equal = 0;
    *more = 1;
    return;
  }
  else if ( compare == 0 )
  {
    *less = 0;
    *equal = 1;
    *more = 2;
    return;
  }

  compare = r8r8r8_compare ( xval, yval, zval, xhi, yhi, zhi );

  if ( compare == 1 )
  {
    *less = n;
    *equal = 0;
    *more = n + 1;
    return;
  }
  else if ( compare == 0 )
  {
    *less = n - 1;
    *equal = n;
    *more = n + 1;
    return;
  }

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      *less = lo;
      *equal = 0;
      *more = hi;
      return;
    }

    mid = ( lo + hi ) / 2;
    xmid = x[indx[mid-1]-1];
    ymid = y[indx[mid-1]-1];
    zmid = z[indx[mid-1]-1];

    compare = r8r8r8_compare ( xval, yval, zval, xmid, ymid, zmid );

    if ( compare == 0 )
    {
      *equal = mid;
      *less = mid - 1;
      *more = mid + 1;
      return;
    }
    else if ( compare == -1 )
    {
      hi = mid;
    }
    else if ( compare == +1 )
    {
      lo = mid;
    }
  }

  return;
}
/******************************************************************************/

void r8r8vec_index_insert_unique ( int maxn, int *n, double x[], double y[],
  int indx[], double xval, double yval, int *ival, int *ierror )

/******************************************************************************/
/*
  Purpose:

    R8R8VEC_INDEX_INSERT_UNIQUE inserts a unique R8R8 value in an indexed sorted list.

  Discussion:

    If the input value does not occur in the current list, it is added,
    and N, X, Y and INDX are updated.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int MAXN, the maximum size of the list.

    Input/output, int *N, the size of the list.

    Input/output, double X[N], Y[N], the list of R8R8 vectors.

    Input/output, int INDX[N], the sort index of the list.

    Input, double XVAL, YVAL, the value to be inserted if it is
    not already in the list.

    Output, int *IVAL, the index in X, Y corresponding to the
    value XVAL, YVAL.

    Output, int *IERROR, 0 for no error, 1 if an error occurred.
*/
{
  int equal;
  int i;
  int less;
  int more;

  *ierror = 0;

  if ( *n <= 0 )
  {
    if ( maxn <= 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!\n" );
      fprintf ( stderr, "  Not enough space to store new data.\n" );
      exit ( 1 );
    }

    *n = 1;
    x[0] = xval;
    y[0] = yval;
    indx[0] = 1;
    *ival = 1;
    return;
  }
//
//  Does ( XVAL, YVAL ) already occur in ( X, Y )?
//
  r8r8vec_index_search ( *n, x, y, indx, xval, yval, &less, &equal, &more );

  if ( equal == 0 )
  {
    if ( maxn <= *n )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8R8VEC_INDEX_INSERT_UNIQUE - Fatal error!\n" );
      fprintf ( stderr, "  Not enough space to store new data.\n" );
      exit ( 1 );
    }

    x[*n] = xval;
    y[*n] = yval;
    *ival = *n + 1;
    for ( i = *n-1; more-1 <= i; i-- )
    {
      indx[i+1] = indx[i];
    }
    indx[more-1] = *n + 1;
    *n = *n + 1;
  }
  else
  {
    *ival = indx[equal-1];
  }

  return;
}
/******************************************************************************/

void r8r8vec_index_search ( int n, double x[], double y[], int indx[],
  double xval, double yval, int *less, int *equal, int *more )

/******************************************************************************/
/*
  Purpose:

    R8R8VEC_INDEX_SEARCH searches for an R8R8 value in an indexed sorted list.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the current list.

    Input, double X[N], Y[N], the list.

    Input, int INDX[N], the sort index of the list.

    Input, double XVAL, YVAL, the value to be sought.

    Output, int *LESS, *EQUAL, *MORE, the indexes in INDX of the
    entries of X that are just less than, equal to, and just greater
    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
    is the greatest entry of X, then MORE is N+1.
*/
{
  int compare;
  int hi;
  int lo;
  int mid;
  double xhi;
  double xlo;
  double xmid;
  double yhi;
  double ylo;
  double ymid;

  if ( n <= 0 )
  {
    *less = 0;
    *equal = 0;
    *more = 0;
    return;
  }

  lo = 1;
  hi = n;

  xlo = x[indx[lo-1]-1];
  ylo = y[indx[lo-1]-1];

  xhi = x[indx[hi-1]-1];
  yhi = y[indx[hi-1]-1];

  compare = r8r8_compare ( xval, yval, xlo, ylo );

  if ( compare == -1 )
  {
    *less = 0;
    *equal = 0;
    *more = 1;
    return;
  }
  else if ( compare == 0 )
  {
    *less = 0;
    *equal = 1;
    *more = 2;
    return;
  }

  compare = r8r8_compare ( xval, yval, xhi, yhi );

  if ( compare == 1 )
  {
    *less = n;
    *equal = 0;
    *more = n + 1;
    return;
  }
  else if ( compare == 0 )
  {
    *less = n - 1;
    *equal = n;
    *more = n + 1;
    return;
  }

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      *less = lo;
      *equal = 0;
      *more = hi;
      return;
    }

    mid = ( lo + hi ) / 2;
    xmid = x[indx[mid-1]-1];
    ymid = y[indx[mid-1]-1];

    compare = r8r8_compare ( xval, yval, xmid, ymid );

    if ( compare == 0 )
    {
      *equal = mid;
      *less = mid - 1;
      *more = mid + 1;
      return;
    }
    else if ( compare == -1 )
    {
      hi = mid;
    }
    else if ( compare == +1 )
    {
      lo = mid;
    }
  }

  return;
}
/******************************************************************************/

int r8row_compare ( int m, int n, double a[], int i, int j )

/******************************************************************************/
/*
  Purpose:

    R8ROW_COMPARE compares rows in an R8ROW.

  Discussion:

    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
    each of length N.

  Example:

    Input:

      M = 4, N = 3, I = 2, J = 4

      A = (
        1  5  9
        2  6 10
        3  7 11
        4  8 12 )

    Output:

      R8ROW_COMPARE = -1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the M by N array.

    Input, int I, J, the rows to be compared.
    I and J must be between 1 and M.

    Output, int R8ROW_COMPARE, the results of the comparison:
    -1, row I < row J,
     0, row I = row J,
    +1, row J < row I.
*/
{
  int k;
  int value;
/*
  Check.
*/
  if ( i < 1 || m < i )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8ROW_COMPARE - Fatal error!\n" );
    fprintf ( stderr, "  Row index I is out of bounds.\n" );
    fprintf ( stderr, "  I = %d\n", i );
    exit ( 1 );
  }

  if ( j < 1 || m < j )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8ROW_COMPARE - Fatal error!\n" );
    fprintf ( stderr, "  Row index J is out of bounds.\n" );
    fprintf ( stderr, "  J = %d\n", j );
    exit ( 1 );
  }

  value = 0;

  if ( i == j )
  {
    return value;
  }

  k = 0;

  while ( k < n )
  {
    if ( a[i-1+k*m] < a[j-1+k*m] )
    {
      value = -1;
      return value;
    }
    else if ( a[j-1+k*m] < a[i-1+k*m] )
    {
      value = +1;
      return value;
    }
    k = k + 1;
  }

  return value;
}
/******************************************************************************/

double *r8row_max ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8ROW_MAX returns the row maximums of an R8ROW.

  Example:

    A =
      1  2  3
      2  6  7

    MAX =
      3
      7

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the array.

    Input, double A[M*N], the array to be examined.

    Output, double R8ROW_MAX[M], the maximums of the rows.
*/
{
  int i;
  int j;
  double *amax;

  amax = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    amax[i] = a[i+0*m];

    for ( j = 1; j < n; j++ )
    {
      if ( amax[i] < a[i+j*m] )
      {
        amax[i] = a[i+j*m];
      }
    }
  }

  return amax;
}
/******************************************************************************/

double *r8row_mean ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8ROW_MEAN returns the row means of an R8ROW.

  Example:

    A =
      1  2  3
      2  6  7

    MEAN =
      2
      5

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the array.

    Input, double A[M*N], the array to be examined.

    Output, double R8ROW_MEAN[M], the means, or averages, of the rows.
*/
{
  int i;
  int j;
  double *mean;

  mean = ( double * ) malloc ( m * sizeof ( m ) );

  for ( i = 0; i < m; i++ )
  {
    mean[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      mean[i] = mean[i] + a[i+j*m];
    }
    mean[i] = mean[i] / ( double ) ( n );
  }

  return mean;
}
/******************************************************************************/

double *r8row_min ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8ROW_MIN returns the row minimums of an R8ROW.

  Example:

    A =
      1  2  3
      2  6  7

    MIN =
      1
      2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the array.

    Input, double A[M*N], the array to be examined.

    Output, double R8ROW_MIN[M], the minimums of the rows.
*/
{
  int i;
  int j;
  double *amin;

  amin = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    amin[i] = a[i+0*m];
    for ( j = 1; j < n; j++ )
    {
      if ( a[i+j*m] < amin[i] )
      {
        amin[i] = a[i+j*m];
      }
    }
  }

  return amin;
}
/******************************************************************************/

void r8row_sort_heap_a ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8ROW_SORT_HEAP_A ascending heapsorts an R8ROW.

  Discussion:

    An R8ROW is an M by N array of R8's, regarded as an array of M rows,
    each of length N.

    In lexicographic order, the statement "X < Y", applied to two real
    vectors X and Y of length M, means that there is some index I, with
    1 <= I <= M, with the property that

      X(J) = Y(J) for J < I,
    and
      X(I) < Y(I).

    In other words, the first time they differ, X is smaller.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N].
    On input, the array of N columns of M-vectors.
    On output, the rows of A have been sorted in lexicographic order.
*/
{
  int i;
  int indx;
  int isgn;
  int j;

  if ( m <= 0 )
  {
    return;
  }

  if ( n <= 1 )
  {
    return;
  }
/*
  Initialize.
*/
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
/*
  Call the external heap sorter.
*/
  for ( ; ; )
  {
    sort_heap_external ( m, &indx, &i, &j, isgn );
/*
  Interchange the I and J objects.
*/
    if ( 0 < indx )
    {
      r8row_swap ( m, n, a, i, j );
    }
/*
  Compare the I and J objects.
*/
    else if ( indx < 0 )
    {
      isgn = r8row_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

double *r8row_sum ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8ROW_SUM returns the sums of the rows of an R8ROW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the M by N array.

    Output, double ROWSUM[M], the sum of the entries of
    each row.
*/
{
  int i;
  int j;
  double *rowsum;

  rowsum = ( double * ) malloc ( m * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    rowsum[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      rowsum[i] = rowsum[i] + a[i+j*m];
    }
  }

  return rowsum;
}
/******************************************************************************/

void r8row_swap ( int m, int n, double a[], int irow1, int irow2 )

/******************************************************************************/
/*
  Purpose:

    R8ROW_SWAP swaps two rows of an R8ROW.

  Discussion:

    The two dimensional information is stored as a one dimensional
    array, by columns.

    The row indices are 1 based, NOT 0 based.  However, a preprocessor
    variable, called OFFSET, can be reset from 1 to 0 if you wish to
    use 0-based indices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 September 2003

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input/output, double A[M*N], an array of data.

    Input, int IROW1, IROW2, the two rows to swap.
    These indices should be between 1 and M.
*/
{
# define OFFSET 1
  int j;
  double t;
/*
  Check.
*/
  if ( irow1 < 1 || m < irow1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8ROW_SWAP - Fatal error!\n" );
    fprintf ( stderr, "  IROW1 is out of range.\n" );
    exit ( 1 );
  }

  if ( irow2 < 1 || m < irow2 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8ROW_SWAP - Fatal error!\n" );
    fprintf ( stderr, "  IROW2 is out of range.\n" );
    exit ( 1 );
  }

  if ( irow1 == irow2 )
  {
    return;
  }

  for ( j = 0; j < n; j++ )
  {
    t              = a[irow1-1+j*m];
    a[irow1-1+j*m] = a[irow2-1+j*m];
    a[irow2-1+j*m] = t;
  }

  return;
# undef OFFSET
}
/******************************************************************************/

double *r8row_to_r8vec ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8ROW_TO_R8VEC converts an R8ROW into an R8VEC.

  Example:

    M = 3, N = 4

    A =
      11 12 13 14
      21 22 23 24
      31 32 33 34

    R8ROW_TO_R8VEC = ( 11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns.

    Input, double A[M*N], the M by N array.

    Output, double R8ROW_TO_R8VEC[M*N], a vector containing the M rows of A.
*/
{
  int i;
  int j;
  int k;
  double *x;

  x = ( double * ) malloc ( m * n * sizeof ( double ) );

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[k] = a[i+j*m];
      k = k + 1;
    }
  }

  return x;
}
/******************************************************************************/

double *r8row_variance ( int m, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8ROW_VARIANCE returns the variances of an R8ROW.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns in the array.

    Input, double A[M*N], the array whose variances are desired.

    Output, double R8ROW_VARIANCE[M], the variances of the rows.
*/
{
  int i;
  int j;
  double mean;
  double *variance;

  variance = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < m; i++ )
  {
    mean = 0.0;
    for ( j = 0; j < n; j++ )
    {
      mean = mean + a[i+j*m];
    }
    mean = mean / ( double ) ( n );

    variance[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      variance[i] = variance[i] + pow ( ( a[i+j*m] - mean ), 2 );
    }

    if ( 1 < n )
    {
      variance[i] = variance[i] / ( double ) ( n - 1 );
    }
    else
    {
      variance[i] = 0.0;
    }

  }

  return variance;
}
/******************************************************************************/

void r8slmat_print ( int m, int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8SLMAT_PRINT prints a strict lower triangular R8MAT.

  Example:

    M = 5, N = 5
    A = (/ 21, 31, 41, 51, 32, 42, 52, 43, 53, 54 /)

    21
    31 32
    41 42 43
    51 52 53 54

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int M, the number of rows in A.

    Input, int N, the number of columns in A.

    Input, double A[*], the M by N matrix.  Only the strict
    lower triangular elements are stored, in column major order.

    Input, char *TITLE, a title.
*/
{
  int i;
  int indx;
  int j;
  int jhi;
  int jlo;
  int jmax;
  int nn;

  printf ( "\n" );
  printf ( "%s\n", title );

  jmax = i4_min ( n, m - 1 );

  nn = 5;

  for ( jlo = 1; jlo <= jmax; jlo = jlo + nn )
  {
    jhi = i4_min ( jlo + nn - 1, i4_min ( m - 1, jmax ) );
    printf ( "\n" );
    printf ( "  Col   " );
    for ( j = jlo; j <= jhi; j++ )
    {
      printf ( "%d       ", j );
    }
    printf ( "\n" );
    printf ( "  Row\n" );
    for ( i = jlo + 1; i <= m; i++ )
    {
      printf ( "%5d:", i );
      jhi = i4_min ( jlo + nn - 1, i4_min ( i - 1, jmax ) );
      for ( j = jlo; j <= jhi; j++ )
      {
        indx = ( j - 1 ) * m + i - ( j * ( j + 1 ) ) / 2;
        printf ( "  %12f", a[indx-1] );
      }
      printf ( "\n" );
    }
  }

  return;
}
/******************************************************************************/

void r8vec_01_to_ab ( int n, double a[], double amax, double amin )

/******************************************************************************/
/*
  Purpose:

    R8VEC_01_TO_AB shifts and rescales data to lie within given bounds.

  Discussion:

    An R8VEC is a vector of R8's.

    On input, A contains the original data, which is presumed to lie
    between 0 and 1.  However, it is not necessary that this be so.

    On output, A has been shifted and rescaled so that all entries which
    on input lay in [0,1] now lie between AMIN and AMAX.  Other entries will
    be mapped in a corresponding way.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of data values.

    Input/output, double A[N], the vector to be rescaled.

    Input, double AMAX, AMIN, the maximum and minimum values
    allowed for A.
*/
{
  double amax2;
  double amax3;
  double amin2;
  double amin3;
  int i;

  if ( amax == amin )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = amin;
    }
    return;
  }

  amax2 = r8_max ( amax, amin );
  amin2 = r8_min ( amax, amin );

  amin3 = r8vec_min ( n, a );
  amax3 = r8vec_max ( n, a );

  if ( amax3 != amin3 )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( amax3 - a[i]         ) * amin2
             + (         a[i] - amin3 ) * amax2 )
             / ( amax3          - amin3 );
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0.5 * ( amax2 + amin2 );
    }
  }

  return;
}
/******************************************************************************/

void r8vec_ab_to_01 ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AB_TO_01 shifts and rescales data to lie within [0,1].

  Discussion:

    An R8VEC is a vector of R8's.

    On input, A contains the original data.  On output, A has been shifted
    and scaled so that all entries lie between 0 and 1.

    The formula is:

      A(I) := ( A(I) - AMIN ) / ( AMAX - AMIN )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of data values.

    Input/output, double A[N], the data to be rescaled.
*/
{
  double amax;
  double amin;
  int i;

  amax = r8vec_max ( n, a );
  amin = r8vec_min ( n, a );

  if ( amin == amax )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0.5;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( a[i] - amin ) / ( amax - amin );
    }
  }

  return;
}
/******************************************************************************/

double *r8vec_ab_to_cd ( int n, double a[], double bmin, double bmax )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AB_TO_CD shifts and rescales data to lie within a given pair of bounds.

  Discussion:

    An R8VEC is a vector of R8's.

    The mininum entry of A is mapped to BMIN, the maximum entry
    to BMAX, and values in between are mapped linearly.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of data values.

    Input, double A[N], the data to be remapped.

    Input, double BMIN, BMAX, the values to which min(A) and max(A)
    are to be assigned.

    Output, double R8VEC_AB_TO_CD[N], the remapped data.
*/
{
  double amax;
  double amin;
  double *b;
  int i;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  if ( bmax == bmin )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i] = bmin;
    }
    return b;
  }

  amax = r8vec_max ( n, a );
  amin = r8vec_min ( n, a );

  if ( amin == amax )
  {
    for ( i = 0; i < n; i++ )
    {
      b[i] = 0.5 * ( bmax + bmin );
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      b[i] = ( ( amax - a[i]        ) * bmin
             + (        a[i] - amin ) * bmax )
             / ( amax        - amin );
    }
  }

  return b;
}
/******************************************************************************/

void r8vec_add ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ADD adds one R8VEC to another.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 September 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], the vector to be added.

    Input/output, double A2[N], the vector to be increased.
    On output, A2 = A2 + A1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a2[i] + a1[i];
  }
  return;
}
/******************************************************************************/

int r8vec_all_nonpositive ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ALL_NONPOSITIVE: ( all ( A <= 0 ) ) for R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries.

    Input, double A[N], the vector to check.

    Output, int R8VEC_ALL_NONPOSITIVE is 1 if all entries
    of A are less than or equal to zero.
*/
{
  int i;
  int value;

  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < a[i] )
    {
      value = 0;
      return value;
    }
  }
  value = 1;

  return value;
}
/******************************************************************************/

double r8vec_amax ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AMAX returns the maximum absolute value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, double AMAX, the value of the entry
    of largest magnitude.
*/
{
  double amax;
  int i;

  amax = 0.0;
  for ( i = 0; i < n; i++ )
  {
    if ( amax < r8_abs ( a[i] ) )
    {
      amax = r8_abs ( a[i] );
    }
  }
  return amax;
}
/******************************************************************************/

int r8vec_amax_index ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AMAX_INDEX returns the index of the maximum absolute value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, int R8VEC_AMAX_INDEX, the index of the entry of largest magnitude.
*/
{
  double amax;
  int amax_index;
  int i;

  if ( n <= 0 )
  {
    amax_index = -1;
  }
  else
  {
    amax_index = 1;
    amax = r8_abs ( a[0] );

    for ( i = 2; i <= n; i++ )
    {
      if ( amax < r8_abs ( a[i-1] ) )
      {
        amax_index = i;
        amax = r8_abs ( a[i-1] );
      }
    }
  }

  return amax_index;
}
/******************************************************************************/

double r8vec_amin ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AMIN returns the minimum absolute value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, double R8VEC_AMIN, the value of the entry
    of smallest magnitude.
*/
{
  double amin;
  int i;

  amin = r8_huge ( );
  for ( i = 0; i < n; i++ )
  {
    if ( r8_abs ( a[i] ) < amin )
    {
      amin = r8_abs ( a[i] );
    }
  }

  return amin;
}
/******************************************************************************/

int r8vec_amin_index ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_AMIN_INDEX returns the index of the minimum absolute value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, int R8VEC_AMIN_INDEX, the index of the entry of smallest magnitude.
*/
{
  double amin;
  int amin_index;
  int i;

  if ( n <= 0 )
  {
    amin_index = -1;
  }
  else
  {
    amin_index = 1;
    amin = r8_abs ( a[0] );

    for ( i = 2; i <= n; i++ )
    {
      if ( r8_abs ( a[i-1] ) < amin )
      {
        amin_index = i;
        amin = r8_abs ( a[i-1] );
      }
    }
  }

  return amin_index;
}
/******************************************************************************/

int r8vec_any_negative ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ANY_NEGATIVE: ( any A < 0 ) for R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    09 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries.

    Input, double A[N], the vector to check.

    Output, int R8VEC_ANY_NEGATIVE is 1 if any entry
    of A is less than zero.
*/
{
  int i;
  int value;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < 0.0 )
    {
      value = 1;
      return value;
    }
  }
  value = 0;

  return value;
}
/******************************************************************************/

int r8vec_any_nonzero ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ANY_NONZERO: ( any A nonzero ) for R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 December 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries.

    Input, double A[N], the vector to check.

    Output, int R8VEC_ANY_NONZERO is 1 if any entry is nonzero.
*/
{
  int i;
  int value;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] != 0.0 )
    {
      value = 1;
      return value;
    }
  }
  value = 0;

  return value;
}
/******************************************************************************/

double *r8vec_any_normal ( int dim_num, double v1[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ANY_NORMAL returns some normal vector to V1.

  Discussion:

    An R8VEC is a vector of R8's.

    If DIM_NUM < 2, then no normal vector can be returned.

    If V1 is the zero vector, then any unit vector will do.

    No doubt, there are better, more robust algorithms.  But I will take
    just about ANY reasonable unit vector that is normal to V1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double V1[DIM_NUM], the vector.

    Output, double R8VEC_ANY_NORMAL[DIM_NUM], a vector that is
    normal to V2, and has unit Euclidean length.
*/
{
  int i;
  int j;
  int k;
  double *v2;
  double vj;
  double vk;

  if ( dim_num < 2 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_ANY_NORMAL - Fatal error!\n" );
    fprintf ( stderr, "  Called with DIM_NUM < 2.\n" );
    exit ( 1 );
  }

  v2 = ( double * ) malloc ( dim_num * sizeof ( double ) );

  if ( r8vec_norm ( dim_num, v1 ) == 0.0 )
  {
    r8vec_zero ( dim_num, v2 );
    v2[0] = 1.0;
    return v2;
  }
/*
  Seek the largest entry in V1, VJ = V1(J), and the
  second largest, VK = V1(K).

  Since V1 does not have zero norm, we are guaranteed that
  VJ, at least, is not zero.
*/
  j = -1;
  vj = 0.0;

  k = -1;
  vk = 0.0;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( r8_abs ( vk ) < r8_abs ( v1[i] ) || k == -1 )
    {
      if ( r8_abs ( vj ) < r8_abs ( v1[i] ) || j == -1 )
      {
        k = j;
        vk = vj;
        j = i;
        vj = v1[i];
      }
      else
      {
        k = i;
        vk = v1[i];
      }
    }
  }
/*
  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
  will just about do the trick.
*/
  r8vec_zero ( dim_num, v2 );

  v2[j] = -vk / sqrt ( vk * vk + vj * vj );
  v2[k] =  vj / sqrt ( vk * vk + vj * vj );

  return v2;
}
/******************************************************************************/

int r8vec_ascends ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ASCENDS determines if an R8VEC is (weakly) ascending.

  Discussion:

    An R8VEC is a vector of R8's.

    For example, if:

      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )

    then

      R8VEC_ASCENDS = TRUE

    The sequence is not required to be strictly ascending.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 July 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the array.

    Input, double X[N], the array to be examined.

    Output, int R8VEC_ASCENDS, is TRUE if the
    entries of X ascend.
*/
{
  int i;
  int value;

  for ( i = 0; i < n - 1; i++ )
  {
    if ( x[i+1] < x[i] )
    {
      value = 0;
      return value;
    }
  }
  value = 1;

  return value;
}
/******************************************************************************/

int r8vec_ascends_strictly ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ASCENDS_STRICTLY determines if an R8VEC is strictly ascending.

  Discussion:

    An R8VEC is a vector of R8's.

    Notice the effect of entry number 6 in the following results:

      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.4, 9.8 )
      Y = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
      Z = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.6, 9.8 )

      R8VEC_ASCENDS_STRICTLY ( X ) = FALSE
      R8VEC_ASCENDS_STRICTLY ( Y ) = FALSE
      R8VEC_ASCENDS_STRICTLY ( Z ) = TRUE

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 July 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the array.

    Input, double X[N], the array to be examined.

    Output, int R8VEC_ASCENDS_STRICTLY, is TRUE if the
    entries of X strictly ascend.
*/
{
  int i;
  int value;

  for ( i = 0; i < n - 1; i++ )
  {
    if ( x[i+1] <= x[i] )
    {
      value = 0;
      return value;
    }
  }
  value = 1;

  return value;
}
/******************************************************************************/

void r8vec_bin ( int n, double x[], int bin_num, double bin_min, double bin_max,
  int bin[], double bin_limit[] )

/******************************************************************************/
/*
//  Purpose:
//
//    R8VEC_BIN computes bins based on a given R8VEC.
//
//  Discussion:
//
//    The user specifies minimum and maximum bin values, BIN_MIN and
//    BIN_MAX, and the number of bins, BIN_NUM.  This determines a
//    "bin width":
//
//      H = ( BIN_MAX - BIN_MIN ) / BIN_NUM
//
//    so that bin I will count all entries X(J) such that
//
//      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
//
//    The array X does NOT have to be sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of X.
//
//    Input, double X[N], an (unsorted) array to be binned.
//
//    Input, int BIN_NUM, the number of bins.  Two extra bins,
//    #0 and #BIN_NUM+1, count extreme values.
//
//    Input, double BIN_MIN, BIN_MAX, define the range and size
//    of the bins.  BIN_MIN and BIN_MAX must be distinct.
//    Normally, BIN_MIN < BIN_MAX, and the documentation will assume
//    this, but proper results will be computed if BIN_MIN > BIN_MAX.
//
//    Output, int BIN[BIN_NUM+2].
//    BIN(0) counts entries of X less than BIN_MIN.
//    BIN(BIN_NUM+1) counts entries greater than or equal to BIN_MAX.
//    For 1 <= I <= BIN_NUM, BIN(I) counts the entries X(J) such that
//      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
//    where H is the bin spacing.
//
//    Output, double BIN_LIMIT[BIN_NUM+1], the "limits" of the bins.
//    BIN(I) counts the number of entries X(J) such that
//      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
*/
{
  int i;
  int j;
  double t;

  if ( bin_max == bin_min )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_BIN - Fatal error!\n" );
    fprintf ( stderr, "  BIN_MIN = BIN_MAX = %g\n", bin_max );
    exit ( 1 );
  }

  for ( i = 0; i <= bin_num + 1; i++ )
  {
    bin[i] = 0;
  }

  for ( i = 0; i < n; i++ )
  {
    t = ( x[i] - bin_min ) / ( bin_max - bin_min );

    if ( t < 0.0 )
    {
      j = 0;
    }
    else if ( 1.0 <= t )
    {
      j = bin_num + 1;
    }
    else
    {
      j = 1 + ( int ) ( ( double ) ( bin_num ) * t );
    }
    bin[j] = bin[j] + 1;
  }
/*
//  Compute the bin limits.
*/
  for ( i = 0; i <= bin_num; i++ )
  {
    bin_limit[i] = (   ( double ) ( bin_num - i ) * bin_min   
                     + ( double ) (           i ) * bin_max ) 
                     / ( double ) ( bin_num     );
  }

  return;
}
/******************************************************************************/

void r8vec_bracket ( int n, double x[], double xval, int *left,
  int *right )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET searches a sorted array for successive brackets of a value.

  Discussion:

    An R8VEC is a vector of R8's.

    If the values in the vector are thought of as defining intervals
    on the real line, then this routine searches for the interval
    nearest to or containing the given value.

    It is always true that RIGHT = LEFT+1.

    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
      XVAL   < X[0] < X[1];
    If X(1) <= XVAL < X[N-1], then
      X[LEFT-1] <= XVAL < X[RIGHT-1];
    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
      X[LEFT-1] <= X[RIGHT-1] <= XVAL.

    For consistency, this routine computes indices RIGHT and LEFT
    that are 1-based, although it would be more natural in C and
    C++ to use 0-based values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, length of input array.

    Input, double X[N], an array that has been sorted into ascending order.

    Input, double XVAL, a value to be bracketed.

    Output, int *LEFT, *RIGHT, the results of the search.
*/
{
  int i;

  for ( i = 2; i <= n - 1; i++ )
  {
    if ( xval < x[i-1] )
    {
      *left = i - 1;
      *right = i;
      return;
    }

   }

  *left = n - 1;
  *right = n;

  return;
}
/******************************************************************************/

void r8vec_bracket2 ( int n, double x[], double xval, int start, int *left,
  int *right )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET2 searches a sorted array for successive brackets of a value.

  Discussion:

    An R8VEC is a vector of R8's.

    If the values in the vector are thought of as defining intervals
    on the real line, then this routine searches for the interval
    containing the given value.

    R8VEC_BRACKET2 is a variation on R8VEC_BRACKET.  It seeks to reduce
    the search time by allowing the user to suggest an interval that
    probably contains the value.  The routine will look in that interval
    and the intervals to the immediate left and right.  If this does
    not locate the point, a binary search will be carried out on
    appropriate subportion of the sorted array.

    In the most common case, 1 <= LEFT < LEFT + 1 = RIGHT <= N,
    and X(LEFT) <= XVAL <= X(RIGHT).

    Special cases:
      Value is less than all data values:
    LEFT = -1, RIGHT = 1, and XVAL < X(RIGHT).
      Value is greater than all data values:
    LEFT = N, RIGHT = -1, and X(LEFT) < XVAL.
      Value is equal to a data value:
    LEFT = RIGHT, and X(LEFT) = X(RIGHT) = XVAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, length of the input array.

    Input, double X[N], an array that has been sorted into
    ascending order.

    Input, double XVAL, a value to be bracketed by entries of X.

    Input, int START, between 1 and N, specifies that XVAL
    is likely to be in the interval:
      [ X(START), X(START+1) ]
    or, if not in that interval, then either
      [ X(START+1), X(START+2) ]
    or
      [ X(START-1), X(START) ].

    Output, int *LEFT, *RIGHT, the results of the search.
*/
{
  int high;
  int low;
/*
  Check.
*/
  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_BRACKET2 - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }

  if ( start < 1 || n < start )
  {
    start = ( n + 1 ) / 2;
  }
/*
  XVAL = X(START)?
*/
  if ( x[start-1] == xval )
  {
    *left = start;
    *right = start;
    return;
  }
/*
  X(START) < XVAL?
*/
  else if ( x[start-1] < xval )
  {
/*
  X(START) = X(N) < XVAL < Infinity?
*/
    if ( n < start + 1 )
    {
      *left = start;
      *right = -1;
      return;
    }
/*
  XVAL = X(START+1)?
*/
    else if ( xval == x[start] )
    {
      *left = start + 1;
      *right = start + 1;
      return;
    }
/*
  X(START) < XVAL < X(START+1)?
*/
    else if ( xval < x[start] )
    {
      *left = start;
      *right = start + 1;
      return;
    }
/*
  X(START+1) = X(N) < XVAL < Infinity?
*/
    else if ( n < start + 2 )
    {
      *left = start + 1;
      *right = -1;
      return;
    }
/*
  XVAL = X(START+2)?
*/
    else if ( xval == x[start+1] )
    {
      *left = start + 2;
      *right = start + 2;
      return;
    }
/*
  X(START+1) < XVAL < X(START+2)?
*/
    else if ( xval < x[start+1] )
    {
      *left = start + 1;
      *right = start + 2;
      return;
    }
/*
  Binary search for XVAL in [ X(START+2), X(N) ],
  where XVAL is guaranteed to be greater than X(START+2).
*/
    else
    {
      low = start + 2;
      high = n;

      r8vec_bracket ( high + 1 - low, x+low-1, xval, left, right );

      *left = *left + low - 1;
      *right = *right + low - 1;
    }
  }
/*
  -Infinity < XVAL < X(START) = X(1).
*/
  else if ( start == 1 )
  {
    *left = -1;
    *right = start;
    return;
  }
/*
  XVAL = X(START-1)?
*/
  else if ( xval == x[start-2] )
  {
    *left = start - 1;
    *right = start - 1;
    return;
  }
/*
  X(START-1) < XVAL < X(START)?
*/
  else if ( x[start-2] <= xval )
  {
    *left = start - 1;
    *right = start;
    return;
  }
/*
  Binary search for XVAL in [ X(1), X(START-1) ],
  where XVAL is guaranteed to be less than X(START-1).
*/
  else
  {
    low = 1;
    high = start - 1;
    r8vec_bracket ( high + 1 - low, x, xval, left, right );
  }

  return;
}
/******************************************************************************/

void r8vec_bracket3 ( int n, double t[], double tval, int *left )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET3 finds the interval containing or nearest a given value.

  Discussion:

    An R8VEC is a vector of R8's.

    The routine always returns the index LEFT of the sorted array
    T with the property that either
    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
    *  T < T[LEFT] = T[0], or
    *  T > T[LEFT+1] = T[N-1].

    The routine is useful for interpolation problems, where
    the abscissa must be located within an interval of data
    abscissas for interpolation, or the "nearest" interval
    to the (extreme) abscissa must be found so that extrapolation
    can be carried out.

    This version of the function has been revised so that the value of
    LEFT that is returned uses the 0-based indexing natural to C++.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, length of the input array.

    Input, double T[N], an array that has been sorted into ascending order.

    Input, double TVAL, a value to be bracketed by entries of T.

    Input/output, int *LEFT.
    On input, if 0 <= LEFT <= N-2, LEFT is taken as a suggestion for the
    interval [ T[LEFT-1] T[LEFT] ] in which TVAL lies.  This interval
    is searched first, followed by the appropriate interval to the left
    or right.  After that, a binary search is used.
    On output, LEFT is set so that the interval [ T[LEFT], T[LEFT+1] ]
    is the closest to TVAL; it either contains TVAL, or else TVAL
    lies outside the interval [ T[0], T[N-1] ].
*/
{
  int high;
  int low;
  int mid;
/*
  Check the input data.
*/
  if ( n < 2 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_BRACKET3 - Fatal error\n" );
    fprintf ( stderr, "  N must be at least 2.\n" );
    exit ( 1 );
  }
/*
  If *LEFT is not between 0 and N-2, set it to the middle value.
*/
  if ( *left < 0 || n - 2 < *left )
  {
    *left = ( n - 1 ) / 2;
  }
/*
  CASE 1: TVAL < T[*LEFT]:
  Search for TVAL in (T[I],T[I+1]), for I = 0 to *LEFT-1.
*/
  if ( tval < t[*left] )
  {
    if ( *left == 0 )
    {
      return;
    }
    else if ( *left == 1 )
    {
      *left = 0;
      return;
    }
    else if ( t[*left-1] <= tval )
    {
      *left = *left - 1;
      return;
    }
    else if ( tval <= t[1] )
    {
      *left = 0;
      return;
    }
/*
  ...Binary search for TVAL in (T[I],T[I+1]), for I = 1 to *LEFT-2.
*/
    low = 1;
    high = *left - 2;

    for ( ; ; )
    {
      if ( low == high )
      {
        *left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid] <= tval )
      {
        low = mid;
      }
      else
      {
        high = mid - 1;
      }
    }
  }
/*
  CASE 2: T[*LEFT+1] < TVAL:
  Search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+1 to N-2.
*/
  else if ( t[*left+1] < tval )
  {
    if ( *left == n - 2 )
    {
      return;
    }
    else if ( *left == n - 3 )
    {
      *left = *left + 1;
      return;
    }
    else if ( tval <= t[*left+2] )
    {
      *left = *left + 1;
      return;
    }
    else if ( t[n-2] <= tval )
    {
      *left = n - 2;
      return;
    }
/*
  ...Binary search for TVAL in (T[I],T[I+1]) for intervals I = *LEFT+2 to N-3.
*/
    low = *left + 2;
    high = n - 3;

    for ( ; ; )
    {

      if ( low == high )
      {
        *left = low;
        return;
      }

      mid = ( low + high + 1 ) / 2;

      if ( t[mid] <= tval )
      {
        low = mid;
      }
      else
      {
        high = mid - 1;
      }
    }
  }
/*
  CASE 3: T[*LEFT] <= TVAL <= T[*LEFT+1]:
  T is just where the user said it might be.
*/
  else
  {
  }

  return;
}
/******************************************************************************/

void r8vec_bracket4 ( int nt, double t[], int ns, double s[], int left[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET4 finds the interval containing or nearest a given value.

  Discussion:

    An R8VEC is a vector of R8's.

    The routine always returns the index LEFT of the sorted array
    T with the property that either
    *  T is contained in the interval [ T[LEFT], T[LEFT+1] ], or
    *  T < T[LEFT] = T[0], or
    *  T > T[LEFT+1] = T[NT-1].

    The routine is useful for interpolation problems, where
    the abscissa must be located within an interval of data
    abscissas for interpolation, or the "nearest" interval
    to the (extreme) abscissa must be found so that extrapolation
    can be carried out.

    This version of the function has been revised so that the value of
    LEFT that is returned uses the 0-based indexing natural to C++.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int NT, length of the input array.

    Input, double T[NT], an array that has been sorted
    into ascending order.

    Input, int NS, the number of points to be bracketed.

    Input, double S[NS], values to be bracketed by entries of T.

    Output, int LEFT[NS].
    LEFT[I] is set so that the interval [ T[LEFT[I]], T[LEFT[I]+1] ]
    is the closest to S[I]; it either contains S[I], or else S[I]
    lies outside the interval [ T[0], T[NT-1] ].
*/
{
  int high;
  int i;
  int low;
  int mid;
/*
  Check the input data.
*/
  if ( nt < 2 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_BRACKET4 - Fatal error!\n" );
    fprintf ( stderr, "  NT must be at least 2.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < ns; i++ )
  {
    left[i] = ( nt - 1 ) / 2;
/*
  CASE 1: S[I] < T[LEFT]:
  Search for S[I] in (T[I],T[I+1]), for I = 0 to LEFT-1.
*/
    if ( s[i] < t[left[i]] )
    {
      if ( left[i] == 0 )
      {
        continue;
      }
      else if ( left[i] == 1 )
      {
        left[i] = 0;
        continue;
      }
      else if ( t[left[i]-1] <= s[i] )
      {
        left[i] = left[i] - 1;
        continue;
      }
      else if ( s[i] <= t[1] )
      {
        left[i] = 0;
        continue;
      }
/*
  ...Binary search for S[I] in (T[I],T[I+1]), for I = 1 to *LEFT-2.
*/
      low = 1;
      high = left[i] - 2;

      for ( ; ; )
      {
        if ( low == high )
        {
          left[i] = low;
          break;
        }

        mid = ( low + high + 1 ) / 2;

        if ( t[mid] <= s[i] )
        {
          low = mid;
        }
        else
        {
          high = mid - 1;
        }
      }
    }
/*
  CASE 2: T[LEFT+1] < S[I]:
  Search for S[I] in (T[I],T[I+1]) for intervals I = LEFT+1 to NT-2.
*/
    else if ( t[left[i]+1] < s[i] )
    {
      if ( left[i] == nt - 2 )
      {
        continue;
      }
      else if ( left[i] == nt - 3 )
      {
        left[i] = left[i] + 1;
        continue;
      }
      else if ( s[i] <= t[left[i]+2] )
      {
        left[i] = left[i] + 1;
        continue;
      }
      else if ( t[nt-2] <= s[i] )
      {
        left[i] = nt - 2;
        continue;
      }
/*
  ...Binary search for S[I] in (T[I],T[I+1]) for intervals I = LEFT+2 to NT-3.
*/
      low = left[i] + 2;
      high = nt - 3;

      for ( ; ; )
      {

        if ( low == high )
        {
          left[i] = low;
          break;
        }

        mid = ( low + high + 1 ) / 2;

        if ( t[mid] <= s[i] )
        {
          low = mid;
        }
        else
        {
          high = mid - 1;
        }
      }
    }
/*
  CASE 3: T[LEFT] <= S[I] <= T[LEFT+1]:
*/
    else
    {
    }
  }
  return;
}
/******************************************************************************/

int r8vec_bracket5 ( int nd, double xd[], double xi )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET5 brackets data between successive entries of a sorted R8VEC.

  Discussion:

    We assume XD is sorted.

    If XI is contained in the interval [XD(1),XD(N)], then the returned 
    value B indicates that XI is contained in [ XD(B), XD(B+1) ].

    If XI is not contained in the interval [XD(1),XD(N)], then B = -1.

    This code implements a version of binary search which is perhaps more
    understandable than the usual ones.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ND, the number of data values.

    Input, double XD[N], the sorted data.

    Input, double XD, the query value.

    Output, int R8VEC_BRACKET5, the bracket information.
*/
{
  int b;
  int l;
  int m;
  int r;

  if ( xi < xd[0] || xd[nd-1] < xi )
  {
    b = -1;
  }
  else
  {
    l = 0;
    r = nd - 1;

    while ( l + 1 < r )
    {
      m = ( l + r ) / 2;
      if ( xi < xd[m] )
      {
        r = m;
      }
      else
      {
        l = m;
      }
    }
    b = l;
  }

  return b;
}
/******************************************************************************/

int *r8vec_bracket6 ( int nd, double xd[], int ni, double xi[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET6 brackets data between successive entries of a sorted R8VEC.

  Discussion:

    We assume XD is sorted.

    If XI(I) is contained in the interval [XD(1),XD(N)], then the value of
    B(I) indicates that XI(I) is contained in [ XD(B(I)), XD(B(I)+1) ].

    If XI(I) is not contained in the interval [XD(1),XD(N)], then B(I) = -1.

    This code implements a version of binary search which is perhaps more
    understandable than the usual ones.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ND, the number of data values.

    Input, double XD[N], the sorted data.

    Input, int NI, the number of inquiry values.

    Input, double XD[NI], the query values.

    Output, int R8VEC_BRACKET6[NI], the bracket information.
*/
{
  int *b;
  int i;
  int l;
  int m;
  int r;

  b = ( int * ) malloc ( ni * sizeof ( int ) );

  for ( i = 0; i < ni; i++ )
  {
    if ( xi[i] < xd[0] || xd[nd-1] < xi[i] )
    {
      b[i] = -1;
    }
    else
    {
      l = 0;
      r = nd - 1;

      while ( l + 1 < r )
      {
        m = ( l + r ) / 2;
        if ( xi[i] < xd[m] )
        {
          r = m;
        }
        else
        {
          l = m;
        }
      }

      b[i] = l;
    }
  }

  return b;
}
/******************************************************************************/

double *r8vec_chebyspace_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CHEBYSPACE_NEW creates a vector of Chebyshev spaced values in [A,B].

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 June 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the interval.

    Output, double R8VEC_CHEBYSPACE_NEW[N], a vector of Chebyshev spaced data.
*/
{
  double c;
  int i;
  double pi = 3.141592653589793;
  double theta;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      theta = ( double ) ( n - i - 1 ) * pi / ( double ) ( n - 1 );

      c = cos ( theta );

      if ( ( n % 2 ) == 1 )
      {
        if ( 2 * i + 1 == n )
        {
          c = 0.0;
        }
      }

      x[i] = ( ( 1.0 - c ) * a  
             + ( 1.0 + c ) * b ) 
             /   2.0;
    }
  }

  return x;
}
/******************************************************************************/

double *r8vec_cheby1space_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CHEBY1SPACE_NEW creates Type 1 Chebyshev spaced values in [A,B].

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the interval.

    Output, double R8VEC_CHEBY1SPACE_NEW[N], a vector of Chebyshev spaced data.
*/
{
  double c;
  int i;
  double pi = 3.141592653589793;
  double theta;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      theta = ( double ) ( 2 * ( n - i ) - 1 ) * pi / ( double ) ( 2 * n );

      c = cos ( theta );

      if ( ( n % 2 ) == 1 )
      {
        if ( 2 * i + 1 == n )
        {
          c = 0.0;
        }
      }

      x[i] = ( ( 1.0 - c ) * a  
             + ( 1.0 + c ) * b ) 
             /   2.0;
    }
  }

  return x;
}
/******************************************************************************/

double *r8vec_cheby2space_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CHEBY2SPACE_NEW creates Type 2 Chebyshev spaced values in [A,B].

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the interval.

    Output, double R8VEC_CHEBY2SPACE_NEW[N], a vector of Chebyshev spaced data.
*/
{
  double c;
  int i;
  double pi = 3.141592653589793;
  double theta;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      theta = ( double ) ( n - i - 1 ) * pi / ( double ) ( n - 1 );

      c = cos ( theta );

      if ( ( n % 2 ) == 1 )
      {
        if ( 2 * i + 1 == n )
        {
          c = 0.0;
        }
      }

      x[i] = ( ( 1.0 - c ) * a  
             + ( 1.0 + c ) * b ) 
             /   2.0;
    }
  }

  return x;
}
/******************************************************************************/

double r8vec_circular_variance ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CIRCULAR_VARIANCE returns the circular variance of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double X(N), the vector whose variance is desired.

    Output, double R8VEC_CIRCULAR VARIANCE, the circular variance
    of the vector entries.
*/
{
  int i;
  double mean;
  double sum_c;
  double sum_s;
  double value;

  mean = r8vec_mean ( n, x );

  sum_c = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum_c = sum_c + cos ( x[i] - mean );
  }

  sum_s = 0.0;
  for ( i = 0; i < n; i++ )
  {
    sum_s = sum_s + sin ( x[i] - mean );
  }

  value = sqrt ( sum_c * sum_c + sum_s * sum_s ) / ( double ) n;

  value = 1.0 - value;

  return value;
}
/******************************************************************************/

int r8vec_compare ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COMPARE compares two R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The lexicographic ordering is used.

  Example:

    Input:

      A1 = ( 2.0, 6.0, 2.0 )
      A2 = ( 2.0, 8.0, 12.0 )

    Output:

      ISGN = -1

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A[N], B[N], the vectors to be compared.

    Output, int R8VEC_COMPARE, the results of the comparison:
    -1, A is lexicographically less than B,
     0, A is equal to B,
    +1, A is lexicographically greater than B.
*/
{
  int isgn;
  int k;

  isgn = 0;

  for ( k = 0; k < n; k++ )
  {
    if ( a[k] < b[k] )
    {
      isgn = -1;
      return isgn;
    }
    else if ( b[k] < a[k] )
    {
      isgn = +1;
      return isgn;
    }
  }
  return isgn;
}
/******************************************************************************/

double *r8vec_convolution ( int m, double x[], int n, double y[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CONVOLUTION returns the convolution of two R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The I-th entry of the convolution can be formed by summing the products 
    that lie along the I-th diagonal of the following table:

    Y3 | 3   4   5   6   7
    Y2 | 2   3   4   5   6
    Y1 | 1   2   3   4   5
       +------------------
        X1  X2  X3  X4  X5

    which will result in:

    Z = ( X1 * Y1,
          X1 * Y2 + X2 * Y1,
          X1 * Y3 + X2 * Y2 + X3 * Y1,
                    X2 * Y3 + X3 * Y2 + X4 * Y1,
                              X3 * Y3 + X4 * Y2 + X5 * Y1,
                                        X4 * Y3 + X5 * Y2,
                                                  X5 * Y3 )
            
  Example:

    Input:

      X = (/ 1, 2, 3, 4 /)
      Y = (/ -1, 5, 3 /)

    Output:

      Z = (/ -1, 3, 10, 17, 29, 12 /)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int M, the dimension of X.

    Input, double X[M], the first vector to be convolved.

    Input, int N, the dimension of Y.

    Input, double Y[N], the second vector to be convolved.

    Output, double R8VEC_CONVOLUTION[M+N-1], the convolution of X and Y.
*/
{
  int i;
  int j;
  double *z;

  z = ( double * ) malloc ( ( m + n - 1 ) * sizeof ( double ) );

  for ( i = 0; i < m + n - 1; i++ )
  {
    z[i] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      z[j+i] = z[j+i] + x[i] * y[j];
    }
  }

  return z;
}
/******************************************************************************/

double *r8vec_convolution_circ ( int n, double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CONVOLUTION_CIRC: discrete circular convolution of two R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    z(1+m) = xCCy(m) = sum ( 0 <= k <= n-1 ) x(1+k) * y(1+m-k)

    Here, if the index of Y becomes nonpositive, it is "wrapped around"
    by having N added to it.

    The circular convolution is equivalent to multiplication of Y by a
    circulant matrix formed from the vector X.

  Example:

    Input:

      X = (/ 1, 2, 3, 4 /)
      Y = (/ 1, 2, 4, 8 /)

    Output:

      Circulant form:

      Z = ( 1 4 3 2 )   ( 1 )
          ( 2 1 4 3 )   ( 2 )
          ( 3 2 1 4 ) * ( 4 )
          ( 4 3 2 1 )   ( 8 )

      The formula:

      Z = (/ 1*1 + 2*8 + 3*4 + 4*2,
             1*2 + 2*1 + 3*8 + 4*4,
             1*4 + 2*2 + 3*1 + 4*8,
             1*8 + 2*4 + 3*2 + 4*1 /)

      Result:

      Z = (/ 37, 44, 43, 26 /)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vectors.

    Input, double X[N], Y[N], the vectors to be convolved.

    Output, double R8VEC_CONVOLVE_CIRC[N], the circular convolution of X and Y.
*/
{
  int i;
  int m;
  double *z;

  z = ( double * ) malloc ( n * sizeof ( double ) );

  for ( m = 1; m <= n; m++ )
  {
    z[m-1] = 0.0;
    for ( i = 1; i <= m; i++ )
    {
      z[m-1] = z[m-1] + x[i-1] * y[m-i];
    }
    for ( i = m+1; i <= n; i++ )
    {
      z[m-1] = z[m-1] + x[i-1] * y[n+m-i];
    }
  }

  return z;
}
/******************************************************************************/

void r8vec_copy ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COPY copies an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 July 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], the vector to be copied.

    Input, double A2[N], the copy of A1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
/******************************************************************************/

double *r8vec_copy_new ( int n, double a1[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_COPY_NEW copies an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], the vector to be copied.

    Output, double R8VEC_COPY_NEW[N], the copy of A1.
*/
{
  double *a2;
  int i;

  a2 = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
/******************************************************************************/

double r8vec_correlation ( int n, double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CORRELATION returns the correlation of two R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    If X and Y are two nonzero vectors of length N, then

      correlation = (x/||x||)' (y/||y||)

    It is the cosine of the angle between the two vectors.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vectors.

    Input, double X[N], Y[N], the vectors to be convolved.

    Output, double R8VEC_CORRELATION, the correlation of X and Y.
*/
{
  double correlation;
  double x_norm;
  double xy_dot;
  double y_norm;

  x_norm = r8vec_norm ( n, x );
  y_norm = r8vec_norm ( n, y );
  xy_dot = r8vec_dot_product ( n, x, y );

  if ( x_norm == 0.0 || y_norm == 0.0 )
  {
    correlation = 0.0;
  }
  else
  {
    correlation = xy_dot / x_norm / y_norm;
  }

  return correlation;
}
/******************************************************************************/

double r8vec_cross_product_2d ( double v1[2], double v2[2] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CROSS_PRODUCT_2D finds the cross product of a pair of R8VEC's in 2D.

  Discussion:

    Strictly speaking, the vectors lie in the (X,Y) plane, and
    the cross product here is a vector in the Z direction.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, double V1[2], V2[2], the vectors.

    Output, double R8VEC_CROSS_PRODUCT_2D, the Z component of the cross product
    of V1 and V2.
*/
{
  double value;

  value = v1[0] * v2[1] - v1[1] * v2[0];

  return value;
}
/******************************************************************************/

double r8vec_cross_product_affine_2d ( double v0[2], double v1[2],
  double v2[2] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CROSS_PRODUCT_AFFINE_2D finds the affine cross product in 2D.

  Discussion:

    Strictly speaking, the vectors lie in the (X,Y) plane, and
    the cross product here is a vector in the Z direction.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double V0[2], the base vector.

    Input, double V1[2], V2[2], the vectors.

    Output, double R8VEC_CROSS_PRODUCT_AFFINE_2D, the Z component of the
    cross product of V1 and V2.
*/
{
  double value;

  value =
      ( v1[0] - v0[0] ) * ( v2[1] - v0[1] )
    - ( v2[0] - v0[0] ) * ( v1[1] - v0[1] );

  return value;
}
/******************************************************************************/

double *r8vec_cross_product_3d ( double v1[3], double v2[3] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CROSS_PRODUCT_3D computes the cross product of two R8VEC's in 3D.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, double V1[3], V2[3], the coordinates of the vectors.

    Output, double R8VEC_CROSS_PRODUCT_3D[3], the cross product vector.
*/
{
  double *v3;

  v3 = ( double * ) malloc ( 3 * sizeof ( double ) );

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}
/******************************************************************************/

double *r8vec_cross_product_affine_3d ( double v0[3], double v1[3],
  double v2[3] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CROSS_PRODUCT_AFFINE_3D computes the affine cross product in 3D.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double V0[3], the base vector.

    Input, double V1[3], V2[3], the coordinates of the vectors.

    Output, double R8VEC_CROSS_PRODUCT_AFFINE_3D[3], the cross product vector.
*/
{
  double *v3;

  v3 = ( double * ) malloc ( 3 * sizeof ( double ) );

  v3[0] =
      ( v1[1] - v0[1] ) * ( v2[2] - v0[2] )
    - ( v2[1] - v0[1] ) * ( v1[2] - v0[2] );

  v3[1] =
      ( v1[2] - v0[2] ) * ( v2[0] - v0[0] )
    - ( v2[2] - v0[2] ) * ( v1[0] - v0[0] );

  v3[2] =
      ( v1[0] - v0[0] ) * ( v2[1] - v0[1] )
    - ( v2[0] - v0[0] ) * ( v1[1] - v0[1] );

  return v3;
}
/******************************************************************************/

double *r8vec_cum_new ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CUM_NEW computes the cumulutive sums of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Input:

      A = { 1.0, 2.0, 3.0, 4.0 }

    Output:

      A_CUM = { 1.0, 3.0, 6.0, 10.0 }

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], the vector to be summed.

    Output, double R8VEC_CUM_NEW[N], the cumulative sums.
*/
{
  double *a_cum;
  int i;

  a_cum = ( double * ) malloc ( n * sizeof ( double ) );

  a_cum[0] = a[0];

  for ( i = 1; i < n; i++ )
  {
    a_cum[i] = a_cum[i-1] + a[i];
  }

  return a_cum;
}
/******************************************************************************/

double *r8vec_cum0_new ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_CUM0_NEW computes the cumulutive sums of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Input:

      A = { 1.0, 2.0, 3.0, 4.0 }

    Output:

      A_CUM = { 0.0, 1.0, 3.0, 6.0, 10.0 }

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], the vector to be summed.

    Output, double R8VEC_CUM0_NEW[N+1], the cumulative sums.
*/
{
  double *a_cum;
  int i;

  a_cum = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

  a_cum[0] = 0.0;

  for ( i = 1; i <= n; i++ )
  {
    a_cum[i] = a_cum[i-1] + a[i-1];
  }

  return a_cum;
}
/******************************************************************************/

double *r8vec_dif ( int n, double h )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIF computes coefficients for estimating the N-th derivative.

  Discussion:

    An R8VEC is a vector of R8's.

    The routine computes the N+1 coefficients for a centered finite difference
    estimate of the N-th derivative of a function.

    The estimate has the form

      FDIF(N,X) = Sum ( 0 <= I <= N ) COF(I) * F ( X(I) )

    To understand the computation of the coefficients, it is enough
    to realize that the first difference approximation is

      FDIF(1,X) = F(X+DX) - F(X-DX) ) / (2*DX)

    and that the second difference approximation can be regarded as
    the first difference approximation repeated:

      FDIF(2,X) = FDIF(1,X+DX) - FDIF(1,X-DX) / (2*DX)
         = F(X+2*DX) - 2 F(X) + F(X-2*DX) / (4*DX)

    and so on for higher order differences.

    Thus, the next thing to consider is the integer coefficients of
    the sampled values of F, which are clearly the Pascal coefficients,
    but with an alternating negative sign.  In particular, if we
    consider row I of Pascal's triangle to have entries j = 0 through I,
    then P(I,J) = P(I-1,J-1) - P(I-1,J), where P(*,-1) is taken to be 0,
    and P(0,0) = 1.

       1
      -1  1
       1 -2   1
      -1  3  -3   1
       1 -4   6  -4   1
      -1  5 -10  10  -5  1
       1 -6  15 -20  15 -6 1

    Next, note that the denominator of the approximation for the
    N-th derivative will be (2*DX)^N.

    And finally, consider the location of the N+1 sampling
    points for F:

      X-N*DX, X-(N-2)*DX, X-(N-4)*DX, ..., X+(N-4)*DX, X+(N-2*DX), X+N*DX.

    Thus, a formula for evaluating FDIF(N,X) is

      fdif = 0.0
      do i = 0, n
        xi = x + (2*i-n) * h
        fdif = fdif + cof(i) * f(xi)
      end do

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the derivative to be approximated.
    N must be 0 or greater.

    Input, double H, the half spacing between points.
    H must be positive.

    Output, double R8VEC_DIF[N+1], the coefficients needed to approximate
    the N-th derivative of a function F.
*/
{
  double *cof;
  int i;
  int j;

  if ( n < 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_DIF - Fatal error!\n" );
    fprintf ( stderr, "  Derivative order N = %d\n", n );
    fprintf ( stderr, "  but N must be at least 0.\n" );
    exit ( 1 );
  }

  if ( h <= 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_DIF - Fatal error!\n" );
    fprintf ( stderr, "  The half sampling spacing is H = %g\n", h );
    fprintf ( stderr, "  but H must be positive.\n" );
    exit ( 1 );
  }

  cof = ( double * ) malloc ( ( n + 1 ) * sizeof ( double ) );

  for ( i = 0; i <= n; i++ )
  {
    cof[i] = 1.0;

    for ( j = i - 1; 1 <= j; j-- )
    {
      cof[j] = - cof[j] + cof[j-1];
    }

    if ( 0 < i )
    {
      cof[0] = -cof[0];
    }
  }

  for ( i = 0; i <= n; i++ )
  {
    cof[i] = cof[i] / pow ( 2.0 * h, n );
  }

  return cof;
}
/******************************************************************************/

double r8vec_diff_norm ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIFF_NORM returns the L2 norm of the difference of R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The vector L2 norm is defined as:

      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], B[N], the vectors.

    Output, double R8VEC_DIFF_NORM, the L2 norm of A - B.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

double r8vec_diff_norm_l1 ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIFF_NORM_L1 returns the L1 norm of the difference of R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The vector L1 norm is defined as:

      R8VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], B[N], the vectors.

    Output, double R8VEC_DIFF_NORM_L1, the L1 norm of A - B.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + r8_abs ( a[i] - b[i] );
  }
  return value;
}
/******************************************************************************/

double r8vec_diff_norm_l2 ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIFF_NORM_L2 returns the L2 norm of the difference of R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The vector L2 norm is defined as:

      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], B[N], the vectors.

    Output, double R8VEC_DIFF_NORM_L2, the L2 norm of A - B.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

double r8vec_diff_norm_li ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIFF_NORM_LI returns the L-oo norm of the difference of R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The vector L-oo norm is defined as:

      R8VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 April 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], B[N], the vectors.

    Output, double R8VEC_DIFF_NORM_LI, the L-oo norm of A - B.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = r8_max ( value, r8_abs ( a[i] - b[i] ) );
  }
  return value;
}
/******************************************************************************/

double r8vec_diff_norm_squared ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIFF_NORM_SQUARED: square of the L2 norm of the difference of R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The square of the L2 norm of the difference of A and B is:

      R8VEC_DIFF_NORM_SQUARED = sum ( 1 <= I <= N ) ( A[I] - B[I] )^2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], B[N], the vectors.

    Output, double R8VEC_DIFF_NORM_SQUARED, the square of the L2 norm of A - B.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( a[i] - b[i] ) * ( a[i] - b[i] );
  }

  return value;
}
/******************************************************************************/

void r8vec_direct_product ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    To explain what is going on here, suppose we had to construct
    a multidimensional quadrature rule as the product of K rules
    for 1D quadrature.

    The product rule will be represented as a list of points and weights.

    The J-th item in the product rule will be associated with
      item J1 of 1D rule 1,
      item J2 of 1D rule 2,
      ...,
      item JK of 1D rule K.

    In particular,
      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
    and
      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)

    So we can construct the quadrature rule if we can properly
    distribute the information in the 1D quadrature rules.

    This routine carries out that task.

    Another way to do this would be to compute, one by one, the
    set of all possible indices (J1,J2,...,JK), and then index
    the appropriate information.  An advantage of the method shown
    here is that you can process the K-th set of information and
    then discard it.

  Example:

    Rule 1:
      Order = 4
      X(1:4) = ( 1, 2, 3, 4 )

    Rule 2:
      Order = 3
      X(1:3) = ( 10, 20, 30 )

    Rule 3:
      Order = 2
      X(1:2) = ( 100, 200 )

    Product Rule:
      Order = 24
      X(1:24) =
        ( 1, 10, 100 )
        ( 2, 10, 100 )
        ( 3, 10, 100 )
        ( 4, 10, 100 )
        ( 1, 20, 100 )
        ( 2, 20, 100 )
        ( 3, 20, 100 )
        ( 4, 20, 100 )
        ( 1, 30, 100 )
        ( 2, 30, 100 )
        ( 3, 30, 100 )
        ( 4, 30, 100 )
        ( 1, 10, 200 )
        ( 2, 10, 200 )
        ( 3, 10, 200 )
        ( 4, 10, 200 )
        ( 1, 20, 200 )
        ( 2, 20, 200 )
        ( 3, 20, 200 )
        ( 4, 20, 200 )
        ( 1, 30, 200 )
        ( 2, 30, 200 )
        ( 3, 30, 200 )
        ( 4, 30, 200 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int FACTOR_INDEX, the index of the factor being processed.
    The first factor processed must be factor 0.

    Input, int FACTOR_ORDER, the order of the factor.

    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values
    for factor FACTOR_INDEX.

    Input, int FACTOR_NUM, the number of factors.

    Input, int POINT_NUM, the number of elements in the direct product.

    Input/output, double X[FACTOR_NUM*POINT_NUM], the elements of the
    direct product, which are built up gradually.

  Local Parameters:

    Local, int START, the first location of a block of values to set.

    Local, int CONTIG, the number of consecutive values to set.

    Local, int SKIP, the distance from the current value of START
    to the next location of a block of values to set.

    Local, int REP, the number of blocks of values to set.
*/
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( j = 0; j < point_num; j++ )
    {
      for ( i = 0; i < factor_num; i++ )
      {
        x[i+j*factor_num] = 0.0;
      }
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( i = 0; i < factor_order; i++ )
  {
    start = 0 + i * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( j = start; j < start + contig; j++ )
      {
        x[factor_index+j*factor_num] = factor_value[i];
      }
      start = start + skip;
    }
  }
  contig = contig * factor_order;

  return;
}
/******************************************************************************/

void r8vec_direct_product2 ( int factor_index, int factor_order,
  double factor_value[], int factor_num, int point_num, double w[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    To explain what is going on here, suppose we had to construct
    a multidimensional quadrature rule as the product of K rules
    for 1D quadrature.

    The product rule will be represented as a list of points and weights.

    The J-th item in the product rule will be associated with
      item J1 of 1D rule 1,
      item J2 of 1D rule 2,
      ...,
      item JK of 1D rule K.

    In particular,
      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
    and
      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)

    So we can construct the quadrature rule if we can properly
    distribute the information in the 1D quadrature rules.

    This routine carries out that task for the weights W.

    Another way to do this would be to compute, one by one, the
    set of all possible indices (J1,J2,...,JK), and then index
    the appropriate information.  An advantage of the method shown
    here is that you can process the K-th set of information and
    then discard it.

  Example:

    Rule 1:
      Order = 4
      W(1:4) = ( 2, 3, 5, 7 )

    Rule 2:
      Order = 3
      W(1:3) = ( 11, 13, 17 )

    Rule 3:
      Order = 2
      W(1:2) = ( 19, 23 )

    Product Rule:
      Order = 24
      W(1:24) =
        ( 2 * 11 * 19 )
        ( 3 * 11 * 19 )
        ( 4 * 11 * 19 )
        ( 7 * 11 * 19 )
        ( 2 * 13 * 19 )
        ( 3 * 13 * 19 )
        ( 5 * 13 * 19 )
        ( 7 * 13 * 19 )
        ( 2 * 17 * 19 )
        ( 3 * 17 * 19 )
        ( 5 * 17 * 19 )
        ( 7 * 17 * 19 )
        ( 2 * 11 * 23 )
        ( 3 * 11 * 23 )
        ( 5 * 11 * 23 )
        ( 7 * 11 * 23 )
        ( 2 * 13 * 23 )
        ( 3 * 13 * 23 )
        ( 5 * 13 * 23 )
        ( 7 * 13 * 23 )
        ( 2 * 17 * 23 )
        ( 3 * 17 * 23 )
        ( 5 * 17 * 23 )
        ( 7 * 17 * 23 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int FACTOR_INDEX, the index of the factor being processed.
    The first factor processed must be factor 0.

    Input, int FACTOR_ORDER, the order of the factor.

    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values for
    factor FACTOR_INDEX.

    Input, int FACTOR_NUM, the number of factors.

    Input, int POINT_NUM, the number of elements in the direct product.

    Input/output, double W[POINT_NUM], the elements of the
    direct product, which are built up gradually.

  Local Parameters:

    Local, integer START, the first location of a block of values to set.

    Local, integer CONTIG, the number of consecutive values to set.

    Local, integer SKIP, the distance from the current value of START
    to the next location of a block of values to set.

    Local, integer REP, the number of blocks of values to set.
*/
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 1.0;
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( j = 0; j < factor_order; j++ )
  {
    start = 0 + j * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( i = start; i < start + contig; i++ )
      {
        w[i] = w[i] * factor_value[j];
      }
      start = start + skip;
    }
  }

  contig = contig * factor_order;

  return;
}
/******************************************************************************/

double r8vec_distance ( int dim_num, double v1[], double v2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DISTANCE returns the Euclidean distance between two R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int DIM_NUM, the spatial dimension.

    Input, double V1[DIM_NUM], V2[DIM_NUM], the vectors.

    Output, double R8VEC_DISTANCE, the Euclidean distance
    between the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < dim_num; i++ )
  {
    value = pow ( v1[i] - v2[i], 2 );
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

int r8vec_distinct ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double X[N], the vector to be checked.

    Output, int R8VEC_DISTINCT is true if all N elements of X
    are distinct.
*/
{
  int i;
  int j;

  for ( i = 1; i <= n-1; i++ )
  {
    for ( j = 1; j <= i - 1; j++ )
    {
      if ( x[i] == x[j] )
      {
        return 0;
      }
    }
  }
  return 1;
}
/******************************************************************************/

void r8vec_divide ( int n, double a[], double s )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DIVIDE divides an R8VEC by a nonzero scalar.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, double A[N].  On input, the vector to be scaled.
    On output, each entry has been divided by S.

    Input, double S, the divisor.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / s;
  }
  return;
}
/******************************************************************************/

double r8vec_dot_product ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 July 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], A2[N], the two vectors to be considered.

    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
/******************************************************************************/

double r8vec_dot_product_affine ( int n, double v0[], double v1[], double v2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DOT_PRODUCT_AFFINE computes the affine dot product.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double V0[N], the base vector.

    Input, double V1[N], V2[N], the two vectors to be considered.

    Output, double R8VEC_DOT_PRODUCT_AFFINE, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v2[i] - v0[i] );
  }
  return value;
}
/******************************************************************************/

int r8vec_eq ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EQ is true if two R8VEC's are equal.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 August 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], A2[N], two vectors to compare.

    Output, int R8VEC_EQ, is TRUE if every pair of elements A1(I) and A2(I) are equal,
    and FALSE otherwise.
*/
{
  int i;
  int value;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] != a2[i] )
    {
      value = 0;
      return value;
    }
  }
  value = 1;
  return value;
}
/******************************************************************************/

void r8vec_even ( int n, double alo, double ahi, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EVEN returns an R8VEC of values evenly spaced between ALO and AHI.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 February 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values.

    Input, double ALO, AHI, the low and high values.

    Output, double A[N], N evenly spaced values.
*/
{
  int i;

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else
  {
    for ( i = 1; i <= n; i++ )
    {
      a[i-1] = ( ( double ) ( n - i     ) * alo
               + ( double ) (     i - 1 ) * ahi )
               / ( double ) ( n     - 1 );
    }
  }

  return;
}
/******************************************************************************/

double *r8vec_even_new ( int n, double alo, double ahi )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EVEN_NEW returns an R8VEC of values evenly spaced between ALO and AHI.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 February 2004

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values.

    Input, double ALO, AHI, the low and high values.

    Output, double R8VEC_EVEN_NEW[N], N evenly spaced values.
*/
{
  double *a;
  int i;

  a = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    a[0] = 0.5 * ( alo + ahi );
  }
  else
  {
    for ( i = 1; i <= n; i++ )
    {
      a[i-1] = ( ( double ) ( n - i     ) * alo
               + ( double ) (     i - 1 ) * ahi )
               / ( double ) ( n     - 1 );
    }
  }

  return a;
}
/******************************************************************************/

double r8vec_even_select ( int n, double xlo, double xhi, int ival )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].

  Discussion:

    An R8VEC is a vector of R8's.

    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / ( N - 1 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values.

    Input, double XLO, XHI, the low and high values.

    Input, int IVAL, the index of the desired point.
    IVAL is normally between 1 and N, but may be any integer value.

    Output, double R8VEC_EVEN_SELECT, the IVAL-th of N evenly spaced values
    between XLO and XHI.
    Unless N = 1, X(1) = XLO and X(N) = XHI.
    If N = 1, then X(1) = 0.5*(XLO+XHI).
*/
{
  double xval;

  if ( n == 1 )
  {
    xval = 0.5 * ( xlo + xhi );
  }
  else
  {
    xval = ( ( double ) ( n - ival     ) * xlo
           + ( double ) (     ival - 1 ) * xhi )
           / ( double ) ( n        - 1 );
  }

  return xval;
}
/******************************************************************************/

void r8vec_even2 ( int maxval, int nfill[], int nold, double xold[],
  int *nval, double xval[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EVEN2 linearly interpolates new numbers into an R8VECa.

  Discussion:

    An R8VEC is a vector of R8's.

    The number of values created between two old values can vary from
    one pair of values to the next.

    The interpolated values are evenly spaced.

    This routine is a generalization of R8VEC_EVEN.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int MAXVAL, the size of the XVAL array, as declared by the
    user.  MAXVAL must be large enough to hold the NVAL values computed by
    this routine.  In other words, MAXVAL must be at least equal to
    NOLD + SUM (1 <= I <= NOLD-1) NFILL(I).

    Input, int NFILL[NOLD-1], the number of values
    to be interpolated between XOLD(I) and XOLD(I+1).
    NFILL(I) does not count the endpoints.  Thus, if
    NFILL(I) is 1, there will be one new point generated
    between XOLD(I) and XOLD(I+1).
    NFILL(I) must be nonnegative.

    Input, int NOLD, the number of values XOLD,
    between which extra values are to be interpolated.

    Input, double XOLD[NOLD], the original vector of numbers
    between which new values are to be interpolated.

    Output, int *NVAL, the number of values computed
    in the XVAL array.
    NVAL = NOLD + SUM ( 1 <= I <= NOLD-1 ) NFILL(I)

    Output, double XVAL[MAXVAL].  On output, XVAL contains the
    NOLD values of XOLD, as well as the interpolated
    values, making a total of NVAL values.
*/
{
  int i;
  int j;
  int nadd;

  *nval = 1;

  for ( i = 1; i <= nold - 1; i++ )
  {

    if ( nfill[i-1] < 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8VEC_EVEN2 - Fatal error!\n" );
      fprintf ( stderr, "  NFILL[I-1] is negative for I = %d\n", i );
      fprintf ( stderr, "  NFILL[I-1] = %d\n", nfill[i-1] );
      exit ( 1 );
    }

    if ( maxval < *nval + nfill[i-1] + 1 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "R8VEC_EVEN2 - Fatal error!\n" );
      fprintf ( stderr, "  MAXVAL = %d is not large enough.\n", maxval );
      fprintf ( stderr, "  for the storage for interval I = %d\n", i );
      exit ( 1 );
    }

    nadd = nfill[i-1] + 2;

    for ( j = 1; j <= nadd; j++ )
    {
      xval[*nval+j-2] = ( ( double ) ( nadd - j     ) * xold[i-1]
                        + ( double ) (        j - 1 ) * xold[i] )
                        / ( double ) ( nadd     - 1 );
    }

    *nval = *nval + nfill[i-1] + 1;
  }

  return;
}
/******************************************************************************/

double r8vec_even2_select ( int n, double xlo, double xhi, int ival )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EVEN2_SELECT returns the I-th of N evenly spaced midpoint values.

  Discussion:

    An R8VEC is a vector of R8's.

    This function returns the I-th of N evenly spaced midpoints of N
    equal subintervals of [XLO,XHI].

    XVAL = ( ( 2 * N - 2 * IVAL + 1 ) * XLO 
           + (         2 * IVAL - 1 ) * XHI ) 
           / ( 2 * N                ) 

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 July 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values.

    Input, double XLO, XHI, the low and high values.

    Input, int IVAL, the index of the desired point.
    IVAL is normally between 1 and N, but may be any integer value.

    Output, double R8VEC_EVEN2_SELECT, the IVAL-th of N evenly spaced midpoints
    between XLO and XHI.
*/
{
  double xval;

  xval = ( ( double ) ( 2 * n - 2 * ival + 1 ) * xlo
         + ( double ) (         2 * ival - 1 ) * xhi )
         / ( double ) ( 2 * n                );

  return xval;
}
/******************************************************************************/

void r8vec_even3 ( int nold, int nval, double xold[], double xval[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EVEN3 evenly interpolates new data into an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    This routine accepts a short vector of numbers, and returns a longer
    vector of numbers, created by interpolating new values between
    the given values.

    Between any two original values, new values are evenly interpolated.

    Over the whole vector, the new numbers are interpolated in
    such a way as to try to minimize the largest distance interval size.

    The algorithm employed is not "perfect".

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int NOLD, the number of values XOLD, between which extra
    values are to be interpolated.

    Input, int NVAL, the number of values to be computed
    in the XVAL array.  NVAL should be at least NOLD.

    Input, double XOLD[NOLD], the original vector of numbers
    between which new values are to be interpolated.

    Output, double XVAL[NVAL].  On output, XVAL contains the
    NOLD values of XOLD, as well as interpolated
    values, making a total of NVAL values.
*/
{
  double density;
  int i;
  int ival;
  int j;
  int nmaybe;
  int npts;
  int ntemp;
  int ntot;
  double xlen;
  double xleni;
  double xlentot;

  xlen = 0.0;
  for ( i = 1; i <= nold - 1; i++ )
  {
    xlen = xlen + r8_abs ( xold[i] - xold[i-1] );
  }

  ntemp = nval - nold;

  density = ( double ) ( ntemp ) / xlen;

  ival = 1;
  ntot = 0;
  xlentot = 0.0;

  for ( i = 1; i <= nold - 1; i++ )
  {
    xleni = r8_abs ( xold[i] - xold[i-1] );
    npts = ( int ) ( density * xleni );
    ntot = ntot + npts;
/*
  Determine if we have enough left-over density that it should
  be changed into a point.  A better algorithm would agonize
  more over where that point should go.
*/
    xlentot = xlentot + xleni;
    nmaybe = r8_nint ( xlentot * density );

    if ( ntot < nmaybe )
    {
      npts = npts + nmaybe - ntot;
      ntot = nmaybe;
    }
    for ( j = 1; j <= npts + 2; j++ )
    {
      xval[ival+j-2] = ( ( double ) ( npts+2 - j     ) * xold[i-1]
                       + ( double ) (          j - 1 ) * xold[i] )
                       / ( double ) ( npts+2     - 1 );
    }
    ival = ival + npts + 1;
  }

  return;
}
/******************************************************************************/

double *r8vec_expand_linear ( int n, double x[], int fat )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EXPAND_LINEAR linearly interpolates new data into an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of input data values.

    Input, double X[N], the original data.

    Input, int FAT, the number of data values to interpolate
    between each pair of original data values.

    Output, double R8VEC_EXPAND_LINEAR[(N-1)*(FAT+1)+1], the "fattened" data.
*/
{
  int i;
  int j;
  int k;
  double *xfat;

  xfat = ( double * ) malloc ( ( (n-1) * (fat+1) + 1 ) * sizeof ( double ) );

  k = 0;

  for ( i = 0; i < n-1; i++ )
  {
    xfat[k] = x[i];
    k = k + 1;

    for ( j = 1; j <= fat; j++ )
    {
      xfat[k] = ( ( double ) ( fat - j + 1 ) * x[i]
                + ( double ) (       j     ) * x[i+1] )
                / ( double ) ( fat     + 1 );
      k = k + 1;
    }
  }

  xfat[k] = x[n-1];
  k = k + 1;

  return xfat;
}
/******************************************************************************/

double *r8vec_expand_linear2 ( int n, double x[], int before, int fat, 
  int after )

/******************************************************************************/
/*
  Purpose:

    R8VEC_EXPAND_LINEAR2 linearly interpolates new data into an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    This routine starts with a vector of data.

    The intent is to "fatten" the data, that is, to insert more points
    between successive values of the original data.

    There will also be extra points placed BEFORE the first original
    value and AFTER that last original value.

    The "fattened" data is equally spaced between the original points.

    The BEFORE data uses the spacing of the first original interval,
    and the AFTER data uses the spacing of the last original interval.

  Example:

    N = 3
    BEFORE = 3
    FAT = 2
    AFTER = 1

    X    = (/                   0.0,           6.0,             7.0       /)
    XFAT = (/ -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 6.33, 6.66, 7.0, 7.66 /)
            3 "BEFORE's"        Old  2 "FATS"  Old    2 "FATS"  Old  1 "AFTER"

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 July 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of input data values.
    N must be at least 2.

    Input, double X[N], the original data.

    Input, int BEFORE, the number of "before" values.

    Input, int FAT, the number of data values to interpolate
    between each pair of original data values.

    Input, int AFTER, the number of "after" values.

    Output, double R8VEC_EXPAND_LINEAR2[BEFORE+(N-1)*(FAT+1)+1+AFTER], the
    "fattened" data.
*/
{
  int i;
  int j;
  int k;
  double *xfat;

  xfat = ( double * ) malloc ( (before+(n-1)*(fat+1)+1+after) * sizeof ( double ) );

  k = 0;
/*
  Points BEFORE.
*/
  for ( j = 1 - before + fat; j <= fat; j++ )
  {
    xfat[k] = ( ( double ) ( fat - j + 1 ) * ( x[0] - ( x[1] - x[0] ) ) 
              + ( double ) (       j     ) *   x[0]          ) 
              / ( double ) ( fat     + 1 );
    k = k + 1;
  }
/*
  Original points and FAT points.
*/
  for ( i = 0; i < n - 1; i++ )
  {
    xfat[k] = x[0];
    k = k + 1;
    for ( j = 1; j <= fat; j++ )
    {
      xfat[k] = ( ( double ) ( fat - j + 1 ) * x[i]
                + ( double ) (       j     ) * x[i+1] ) 
                / ( double ) ( fat     + 1 );
      k = k + 1;
    }
  }

  xfat[k] = x[n-1];
  k = k + 1;
/*
  Points AFTER.
*/
  for ( j = 1; j <= after; j++ )
  {
    xfat[k] = ( ( double ) ( fat - j + 1 ) * x[n-1]
              + ( double ) (       j     ) * ( x[n-1] + ( x[n-1] - x[n-2] ) ) ) 
              / ( double ) ( fat     + 1 );
    k = k + 1;
  }

  return xfat;
}
/******************************************************************************/

int *r8vec_first_index ( int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8VEC_FIRST_INDEX indexes the first occurrence of values in an R8VEC.

  Discussion:

    For element A(I) of the vector, FIRST_INDEX(I) is the index in A of
    the first occurrence of the value A(I).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input, double A[N], the unsorted array to examine.

    Input, double TOL, a tolerance for equality.

    Output, int R8VEC_FIRST_INDEX[N], the first occurrence index.
*/
{
  int *first_index;
  int i;
  int j;

  first_index = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    first_index[i] = -1;
  }
  for ( i = 0; i < n; i++ )
  {
    if ( first_index[i] == -1 )
    {
      first_index[i] = i;
      for ( j = i + 1; j < n; j++ )
      {
        if ( r8_abs ( a[i] - a[j] ) <= tol )
        {
          first_index[j] = i;
        }
      }
    }
  }
  return first_index;
}
/******************************************************************************/

double r8vec_frac ( int n, double a[], int k )

/******************************************************************************/
/*
  Purpose:

    R8VEC_FRAC searches for the K-th smallest entry in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Hoare's algorithm is used.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Parameters:

    Input, int N, the number of elements of A.

    Input/output, double A[N].
    On input, A is the array to search.
    On output, the elements of A have been somewhat rearranged.

    Input, int K, the fractile to be sought.  If K = 1, the minimum
    entry is sought.  If K = N, the maximum is sought.  Other values
    of K search for the entry which is K-th in size.  K must be at
    least 1, and no greater than N.

    Output, double R8VEC_FRAC, the value of the K-th fractile of A.
*/
{
  double frac;
  int i;
  int iryt;
  int j;
  int left;
  double temp;
  double x;

  if ( n <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_FRAC - Fatal error!\n" );
    fprintf ( stderr, "  Illegal nonpositive value of N = %d\n", n );
    exit ( 1 );
  }

  if ( k <= 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_FRAC - Fatal error!\n" );
    fprintf ( stderr, "  Illegal nonpositive value of K = %d\n", k );
    exit ( 1 );
  }

  if ( n < k )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_FRAC - Fatal error!\n" );
    fprintf ( stderr, "  Illegal N < K, K = %d\n", k );
    exit ( 1 );
  }

  left = 1;
  iryt = n;

  for ( ; ; )
  {
    if ( iryt <= left )
    {
      frac = a[k-1];
      break;
    }

    x = a[k-1];
    i = left;
    j = iryt;

    for ( ; ; )
    {
      if ( j < i )
      {
        if ( j < k )
        {
          left = i;
        }
        if ( k < i )
        {
          iryt = j;
        }
        break;
      }
/*
  Find I so that X <= A(I).
*/
      while ( a[i-1] < x )
      {
        i = i + 1;
      }
/*
  Find J so that A(J) <= X.
*/
      while ( x < a[j-1] )
      {
        j = j - 1;
      }

      if ( i <= j )
      {
        temp   = a[i-1];
        a[i-1] = a[j-1];
        a[j-1] = temp;
        i = i + 1;
        j = j - 1;
      }
    }
  }

  return frac;
}
/******************************************************************************/

double *r8vec_fraction ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_FRACTION returns the fraction parts of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    If we regard a real number as

      R8 = SIGN * ( WHOLE + FRACTION )

    where

      SIGN is +1 or -1,
      WHOLE is a nonnegative integer
      FRACTION is a nonnegative real number strictly less than 1,

    then this routine returns the value of FRACTION.

  Example:

     R8    R8_FRACTION

    0.00      0.00
    1.01      0.01
    2.02      0.02
   19.73      0.73
   -4.34      0.34

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    20 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of arguments.

    Input, double X[N], the arguments.

    Output, double R8_FRACTION[N], the fraction parts.
*/
{
  double *fraction;
  int i;

  fraction = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    fraction[i] = r8_abs ( x[i] ) - ( double ) ( ( int ) ( r8_abs ( x[i] ) ) );
  }

  return fraction;
}
/******************************************************************************/

int r8vec_gt ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_GT == ( A1 > A2 ) for two R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The comparison is lexicographic.

    A1 > A2  <=>                              A1(1) > A2(1) or
                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
                 ...
                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 January 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vectors.

    Input, double A1[N], A2[N], the vectors to be compared.

    Output, int R8VEC_GT, is TRUE if and only if A1 > A2.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {

    if ( a2[i] < a1[i] )
    {
       return 1;
    }
    else if ( a1[i] < a2[i] )
    {
      return 0;
    }

  }

  return 0;
}
/******************************************************************************/

void r8vec_heap_a ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_HEAP_A reorders an R8VEC into a ascending heap.

  Discussion:

    An R8VEC is a vector of R8's.

    An ascending heap is an array A with the property that, for every index J,
    A[J] <= A[2*J+1] and A[J] <= A[2*J+2], (as long as the indices
    2*J+1 and 2*J+2 are legal).

  Diagram:

                  A(0)
                /      \
            A(1)         A(2)
          /     \        /  \
      A(3)       A(4)  A(5) A(6)
      /  \       /   \
    A(7) A(8)  A(9) A(10)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the size of the input array.

    Input/output, double A[N].
    On input, an unsorted array.
    On output, the array has been reordered into a heap.
*/
{
  int i;
  int ifree;
  double key;
  int m;
/*
  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
*/
  for ( i = (n/2)-1; 0 <= i; i-- )
  {
/*
  Copy the value out of the parent node.
  Position IFREE is now "open".
*/
    key = a[i];
    ifree = i;

    for ( ; ; )
    {
/*
  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
  IFREE.  (One or both may not exist because they equal or exceed N.)
*/
      m = 2 * ifree + 1;
/*
  Does the first position exist?
*/
      if ( n <= m )
      {
        break;
      }
      else
      {
/*
  Does the second position exist?
*/
        if ( m + 1 < n )
        {
/*
  If both positions exist, take the larger of the two values,
  and update M if necessary.
*/
          if ( a[m+1] < a[m] )
          {
            m = m + 1;
          }
        }
/*
  If the large descendant is larger than KEY, move it up,
  and update IFREE, the location of the free position, and
  consider the descendants of THIS position.
*/
        if ( a[m] <= key )
        {
          break;
        }
        a[ifree] = a[m];
        ifree = m;
      }
    }
/*
  When you have stopped shifting items up, return the item you
  pulled out back to the heap.
*/
    a[ifree] = key;
  }

  return;
}
/******************************************************************************/

void r8vec_heap_d ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_HEAP_D reorders an R8VEC into a descending heap.

  Discussion:

    An R8VEC is a vector of R8's.

    A heap is an array A with the property that, for every index J,
    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
    2*J+1 and 2*J+2 are legal).

  Diagram:

                  A(0)
                /      \
            A(1)         A(2)
          /     \        /  \
      A(3)       A(4)  A(5) A(6)
      /  \       /   \
    A(7) A(8)  A(9) A(10)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the size of the input array.

    Input/output, double A[N].
    On input, an unsorted array.
    On output, the array has been reordered into a heap.
*/
{
  int i;
  int ifree;
  double key;
  int m;
/*
  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
*/
  for ( i = (n/2)-1; 0 <= i; i-- )
  {
/*
  Copy the value out of the parent node.
  Position IFREE is now "open".
*/
    key = a[i];
    ifree = i;

    for ( ; ; )
    {
/*
  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
  IFREE.  (One or both may not exist because they equal or exceed N.)
*/
      m = 2 * ifree + 1;
/*
  Does the first position exist?
*/
      if ( n <= m )
      {
        break;
      }
      else
      {
/*
  Does the second position exist?
*/
        if ( m + 1 < n )
        {
/*
  If both positions exist, take the larger of the two values,
  and update M if necessary.
*/
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
/*
  If the large descendant is larger than KEY, move it up,
  and update IFREE, the location of the free position, and
  consider the descendants of THIS position.
*/
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }
      }
    }
/*
  When you have stopped shifting items up, return the item you
  pulled out back to the heap.
*/
    a[ifree] = key;
  }

  return;
}
/******************************************************************************/

int *r8vec_histogram ( int n, double a[], double a_lo, double a_hi,
  int histo_num )

/******************************************************************************/
/*
  Purpose:

    R8VEC_HISTOGRAM histograms an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Values between A_LO and A_HI will be histogrammed into the bins
    1 through HISTO_NUM.  Values below A_LO are counted in bin 0,
    and values greater than A_HI are counted in bin HISTO_NUM+1.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input, double A[N], the array to examine.

    Input, double A_LO, A_HI, the lowest and highest
    values to be histogrammed.  These values will also define the bins.

    Input, int HISTO_NUM, the number of bins to use.

    Output, int HISTO_GRAM[HISTO_NUM+2], contains the number of
    entries of A in each bin.
*/
{
  double delta;
  int *histo_gram;
  int i;
  int j;

  histo_gram = ( int * ) malloc ( ( histo_num + 2 ) * sizeof ( int ) );

  i4vec_zero ( histo_num+2, histo_gram );

  delta = ( a_hi - a_lo ) / ( double ) ( 2 * histo_num );

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < a_lo )
    {
      histo_gram[0] = histo_gram[0] + 1;
    }
    else if ( a[i] <= a_hi )
    {
      j = r8_nint (
        ( ( a_hi -       delta - a[i]        ) * ( double ) ( 1         )
        + (      -       delta + a[i] - a_lo ) * ( double ) ( histo_num ) )
        / ( a_hi - 2.0 * delta        - a_lo ) );

      histo_gram[j] = histo_gram[j] + 1;
    }
    else if ( a_hi < a[i] )
    {
      histo_gram[histo_num+1] = histo_gram[histo_num+1] + 1;
    }
  }

  return histo_gram;
}
/******************************************************************************/

double *r8vec_house_column ( int n, double a[], int k )

/******************************************************************************/
/*
  Purpose:

    R8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.

  Discussion:

    An R8VEC is a vector of R8's.

    The routine returns a vector V that defines a Householder
    premultiplier matrix H(V) that zeros out the subdiagonal entries of
    column K of the matrix A.

       H(V) = I - 2 * v * v'

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix A.

    Input, double A[N], column K of the matrix A.

    Input, int K, the column of the matrix to be modified.

    Output, double R8VEC_HOUSE_COLUMN[N], a vector of unit L2 norm which
    defines an orthogonal Householder premultiplier matrix H with the property
    that the K-th column of H*A is zero below the diagonal.
*/
{
  int i;
  double s;
  double *v;

  v = r8vec_zero_new ( n );

  if ( k < 1 || n <= k )
  {
    return v;
  }

  s = r8vec_norm_l2 ( n+1-k, a+k-1 );

  if ( s == 0.0 )
  {
    return v;
  }

  v[k-1] = a[k-1] + r8_abs ( s ) * r8_sign ( a[k-1] );

  r8vec_copy ( n-k, a+k, v+k );

  s = r8vec_norm_l2 ( n-k+1, v+k-1 );

  for ( i = k-1; i < n; i++ )
  {
    v[i] = v[i] / s;
  }

  return v;
}
/******************************************************************************/

double r8vec_i4vec_dot_product ( int n, double r8vec[], int i4vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_I4VEC_DOT_PRODUCT computes the dot product of an R8VEC and an I4VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    An I4VEC is a vector of I4's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    30 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double R8VEC[N], the first vector.

    Input, int I4VEC[N], the second vector.

    Output, double R8VEC_I4VEC_DOT_PRODUCT, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + r8vec[i] * ( double ) ( i4vec[i] );
  }
  return value;
}
/******************************************************************************/

int r8vec_in_01 ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_IN_01 is TRUE if the entries of an R8VEC are in the range [0,1].

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    06 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries.

    Input, double X[N], the vector

    Output, int R8VEC_IN_01, is TRUE if every entry is
    between 0 and 1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < 0.0 || 1.0 < x[i] )
    {
      return 0;
    }
  }

  return 1;
}
/******************************************************************************/

int r8vec_in_ab ( int n, double x[], double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_IN_AB is TRUE if the entries of an R8VEC are in the range [A,B].

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 April 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries.

    Input, double X[N], the vector

    Input, double A, B, the limits of the range.

    Output, int R8VEC_IN_AB, is TRUE if every entry is
    between A and B.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < a || b < x[i] )
    {
      return 0;
    }
  }

  return 1;
}
/******************************************************************************/

void r8vec_index_delete_all ( int n, double x[], int indx[], double xval,
  int *n2, double x2[], int indx2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_DELETE_ALL deletes all occurrences of a value from an indexed sorted list.

  Discussion:

    An R8VEC is a vector of R8's.

    Note that the value of N is adjusted because of the deletions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 October 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the current list.

    Input, double X[N], the list.

    Input, int INDX[N], the sort index of the list.

    Input, double XVAL, the value to be sought.

    Output, int *N2, the size of the current list.

    Output, double X2[N2], the list.

    Output, int INDX2[N2], the sort index of the list.
*/
{
  int equal;
  int equal1;
  int equal2;
  int get;
  int i;
  int less;
  int more;
  int put;

  if ( n < 1 )
  {
    *n2 = 0;
    return;
  }

  i4vec_copy ( n, indx, indx2 );
  r8vec_copy ( n, x, x2 );
  *n2 = n;

  r8vec_index_search ( *n2, x2, indx2, xval, &less, &equal, &more );

  if ( equal == 0 )
  {
    return;
  }

  equal1 = equal;

  for ( ; ; )
  {
    if ( equal1 <= 1 )
    {
      break;
    }

    if ( x2[indx2[equal1-2]-1] != xval )
    {
      break;
    }
    equal1 = equal1 - 1;
  }

  equal2 = equal;

  for ( ; ; )
  {
    if ( *n2 <= equal2 )
    {
      break;
    }

    if ( x2[indx2[equal2]-1] != xval )
    {
      break;
    }
    equal2 = equal2 + 1;
  }
/*
  Discard certain X values.
*/
  put = 0;

  for ( get = 1; get <= *n2; get++ )
  {
    if ( x2[get-1] != xval )
    {
      put = put + 1;
      x2[put-1] = x2[get-1];
    }
  }
/*
  Adjust the INDX values.
*/
  for ( equal = equal1; equal <= equal2; equal++ )
  {
    for ( i = 1; i <= *n2; i++ )
    {
      if ( indx2[equal-1] < indx2[i-1] )
      {
        indx2[i-1] = indx2[i-1] - 1;
      }
    }
  }
/*
  Discard certain INDX values.
*/
  for ( i = 0; i <= *n2 - equal2 - 1; i++ )
  {
    indx2[equal1+i-1] = indx2[equal2+i];
  }
  for ( i = *n2 + equal1 - equal2; i <= *n2; i++ )
  {
    indx2[i-1] = 0;
  }
/*
  Adjust N.
*/
  *n2 = put;

  return;
}
/******************************************************************************/

void r8vec_index_delete_dupes ( int n, double x[], int indx[],
  int *n2, double x2[], int indx2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_DELETE_DUPES deletes duplicates from an indexed sorted list.

  Discussion:

    An R8VEC is a vector of R8's.

    The output quantities N2, X2, and INDX2 are computed from the
    input quantities by sorting, and eliminating duplicates.

    The output arrays should be dimensioned of size N, unless the user
    knows in advance what the value of N2 will be.

    The output arrays may be identified with the input arrays.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 October 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the input list.

    Input, double X[N], the list.

    Input, int INDX[N], the sort index of the list.

    Output, int *N2, the number of unique entries in X.

    Output, double X2[N2], a copy of the list which has
    been sorted, and made unique.

    Output, int INDX2[N2], the sort index of the new list.
*/
{
  int i;
  int n3;
  double *x3;

  i = 0;
  n3 = 0;
  x3 = ( double * ) malloc ( n * sizeof ( double ) );

  for ( ; ; )
  {
    i = i + 1;

    if ( n < i )
    {
      break;
    }

    if ( 1 < i )
    {
      if ( x[indx[i-1]-1] == x3[n3-1] )
      {
        continue;
      }
    }
    n3 = n3 + 1;
    x3[n3-1] = x[indx[i-1]-1];
  }
/*
  Set the output data.
*/
  *n2 = n3;
  r8vec_copy ( n3, x3, x2 );
  for ( i = 0; i < n3; i++ )
  {
    indx2[i] = i + 1;
  }

  free ( x3 );

  return;
}
/******************************************************************************/

void r8vec_index_delete_one ( int n, double x[], int indx[], double xval,
  int *n2, double x2[], int indx2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_DELETE_ONE deletes one copy of a value from an indexed sorted list.

  Discussion:

    An R8VEC is a vector of R8's.

    If the value occurs in the list more than once, only one copy is deleted.

    Note that the value of N is adjusted because of the deletions.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 October 2000

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the current list.

    Input, double X[N], the list.

    Input, int INDX[N], the sort index of the list.

    Input, double XVAL, the value to be sought.

    Output, int *N2, the size of the current list.

    Output, double X2[N2], the list.

    Output, int INDX2[N2], the sort index of the list.
*/
{
  int equal;
  int i;
  int j;
  int less;
  int more;

  if ( n < 1 )
  {
    *n2 = 0;
    return;
  }

  *n2 = n;
  i4vec_copy ( *n2, indx, indx2 );
  r8vec_copy ( *n2, x, x2 );

  r8vec_index_search ( *n2, x2, indx2, xval, &less, &equal, &more );

  if ( equal != 0 )
  {
    j = indx2[equal-1];
    for ( i = j; i <= *n2-1; i++ )
    {
      x2[i-1] = x[i];
    }
    for ( i = equal; i <= *n2-1; i++ )
    {
      indx2[i-1] = indx2[i];
    }
    for ( i = 1; i <= *n2 - 1; i++ )
    {
      if ( j < indx2[i-1] )
      {
        indx2[i-1] = indx2[i-1] - 1;
      }
    }
    *n2 = *n2 - 1;
  }

  return;
}
/******************************************************************************/

void r8vec_index_insert ( int *n, double x[], int indx[], double xval )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_INSERT inserts a value in an indexed sorted list.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2005

  Author:

    John Burkardt

  Parameters:

    Input/output, int *N, the size of the current list.

    Input, double X[N], the list.

    Input, int INDX[N], the sort index of the list.

    Input, double XVAL, the value to be sought.
*/
{
  int equal;
  int i;
  int less;
  int more;

  if ( *n <= 0 )
  {
    *n = 1;
    x[0] = xval;
    indx[0] = 1;
    return;
  }

  r8vec_index_search ( *n, x, indx, xval, &less, &equal, &more );

  x[*n] = xval;
  for ( i = *n; more <= i; i-- )
  {
    indx[i] = indx[i-1];
  }
  indx[more-1] = *n + 1;
  *n = *n + 1;

  return;
}
/******************************************************************************/

void r8vec_index_insert_unique ( int *n, double x[], int indx[], double xval )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_INSERT_UNIQUE inserts a unique value in an indexed sorted list.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2005

  Author:

    John Burkardt

  Parameters:

    Input/output, int *N, the size of the current list.
    If the input value XVAL does not already occur in X, then N is increased.

    Input/output, double X[N], the list.
    If the input value XVAL does not already occur in X, then it is added
    to X.

    Input/output, int INDX[N], the sort index of the list.
    If the input value XVAL does not already occur in X, then INDX is updated.

    Input, double XVAL, the value which will be inserted into the X
    vector if it is not there already.
*/
{
  int equal;
  int i;
  int less;
  int more;

  if ( *n <= 0 )
  {
    *n = 1;
    x[0] = xval;
    indx[0] = 1;
    return;
  }
/*
  Does XVAL already occur in X?
*/
  r8vec_index_search ( *n, x, indx, xval, &less, &equal, &more );

  if ( equal == 0 )
  {
    x[*n] = xval;
    for ( i = *n; more <= i; i-- )
    {
      indx[i] = indx[i-1];
    }
    indx[more-1] = *n + 1;
    *n = *n + 1;
  }

  return;
}
/******************************************************************************/

void r8vec_index_order ( int n, double x[], int indx[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_ORDER sorts an R8VEC using an index vector.

  Discussion:

    An R8VEC is a vector of R8's.

    The index vector itself is not modified.  Therefore, the pair
    (X,INDX) no longer represents an index sorted vector.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the current list.

    Input/output, double X[N], the list.  On output, the list
    has been sorted.

    Input, int INDX[N], the sort index of the list.
*/
{
  int i;
  double *y;

  y = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[indx[i]-1];
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = y[i];
  }
  free ( y );

  return;
}
/******************************************************************************/

void r8vec_index_search ( int n, double x[], int indx[], double xval, int *less,
  int *equal, int *more )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_SEARCH searches for a value in an indexed sorted R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the current list.

    Input, double X[N], the list.

    Input, int INDX[N], the sort index of the list.

    Input, double XVAL, the value to be sought.

    Output, int *LESS, *EQUAL, *MORE, the indexes in INDX of the
    entries of X that are just less than, equal to, and just greater
    than XVAL.  If XVAL does not occur in X, then EQUAL is zero.
    If XVAL is the minimum entry of X, then LESS is 0.  If XVAL
    is the greatest entry of X, then MORE is N+1.
*/
{
  int hi;
  int lo;
  int mid;
  double xhi;
  double xlo;
  double xmid;

  if ( n <= 0 )
  {
    *less = 0;
    *equal = 0;
    *more = 0;
    return;
  }

  lo = 1;
  hi = n;
  xlo = x[indx[lo-1]-1];
  xhi = x[indx[hi-1]-1];

  if ( xval < xlo )
  {
    *less = 0;
    *equal = 0;
    *more = 1;
    return;
  }
  else if ( xval == xlo )
  {
    *less = 0;
    *equal = 1;
    *more = 2;
    return;
  }

  if ( xhi < xval )
  {
    *less = n;
    *equal = 0;
    *more = n + 1;
    return;
  }
  else if ( xval == xhi )
  {
    *less = n - 1;
    *equal = n;
    *more = n + 1;
    return;
  }

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      *less = lo;
      *equal = 0;
      *more = hi;
      return;
    }

    mid = ( lo + hi ) / 2;
    xmid = x[indx[mid-1]-1];

    if ( xval == xmid )
    {
      *equal = mid;
      *less = mid - 1;
      *more = mid + 1;
      return;
    }
    else if ( xval < xmid )
    {
      hi = mid;
    }
    else if ( xmid < xval )
    {
      lo = mid;
    }
  }
  return;
}
/******************************************************************************/

void r8vec_index_sort_unique ( int n, double x[], int *n2, double x2[],
  int indx2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_SORT_UNIQUE creates a sort index for an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the current list.

    Input, double X[N], the list.

    Output, int *N2, the number of unique elements in X.

    Output, double X2[N2], a list of the unique elements of X.

    Output, int INDX2[N2], the sort index of the list.
*/
{
  int i;

  *n2 = 0;

  for ( i = 0; i < n; i++ )
  {
    r8vec_index_insert_unique ( n2, x2, indx2, x[i] );
  }

  for ( i = *n2; i < n; i++ )
  {
    x2[i] = -1;
  }
  for ( i = *n2; i < n; i++ )
  {
    indx2[i] = -1;
  }

  return;
}
/******************************************************************************/

void r8vec_index_sorted_range ( int n, double r[], int indx[], double r_lo,
  double r_hi, int *i_lo, int *i_hi )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEX_SORTED_RANGE: search index sorted vector for elements in a range.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of items in the vector.

    Input, double R[N], the index sorted vector.

    Input, int INDX[N], the vector used to sort R.
    The vector R[INDX[*]] is sorted.

    Input, double R_LO, R_HI, the limits of the range.

    Output, int *I_LO, *I_HI, the range of indices
    so that I_LO <= I <= I_HI => R_LO <= R[INDX[I]] <= R_HI.  If no
    values in R lie in the range, then I_HI < I_LO will be returned.
*/
{
  int i1;
  int i2;
  int j1;
  int j2;
/*
  Cases we can handle immediately.
*/
  if ( r[indx[n-1]] < r_lo )
  {
    *i_lo = n;
    *i_hi = n - 1;
    return;
  }

  if ( r_hi < r[indx[0]] )
  {
    *i_lo = 0;
    *i_hi = -1;
    return;
  }
/*
  Are there are least two intervals?
*/
  if ( n == 1 )
  {
    if ( r_lo <= r[indx[0]] && r[indx[0]] <= r_hi )
    {
      *i_lo = 1;
      *i_hi = 1;
    }
    else
    {
      *i_lo = 0;
      *i_hi = -1;
    }
    return;
  }
/*
  Bracket R_LO.
*/
  if ( r_lo <= r[indx[0]] )
  {
    i_lo = 0;
  }
  else
  {
/*
  R_LO is in one of the intervals spanned by R(INDX(J1)) to R(INDX(J2)).
  Examine the intermediate interval [R(INDX(I1)), R(INDX(I1+1))].
  Does R_LO lie here, or below or above?
*/
    j1 = 0;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_lo < r[indx[i1]] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[indx[i2]] < r_lo )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        *i_lo = i1;
        break;
      }
    }
  }
/*
  Bracket R_HI.
*/
  if ( r[indx[n-1]] <= r_hi )
  {
    *i_hi = n - 1;
  }
  else
  {
    j1 = *i_lo;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_hi < r[indx[i1]] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[indx[i2]] < r_hi )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        *i_hi = i2;
        break;
      }
    }
  }
/*
  We expect to have computed the largest I_LO and smallest I_HI such that
    R(INDX(I_LO)) <= R_LO <= R_HI <= R(INDX(I_HI))
  but what we want is actually
    R_LO <= R(INDX(I_LO)) <= R(INDX(I_HI)) <= R_HI
  which we can usually get simply by incrementing I_LO and decrementing I_HI.
*/
  if ( r[indx[*i_lo]] < r_lo )
  {
    *i_lo = *i_lo + 1;
    if ( n - 1 < *i_lo )
    {
      *i_hi = *i_lo - 1;
    }
  }

  if ( r_hi < r[indx[*i_hi]] )
  {
    *i_hi = *i_hi - 1;
    if ( i_hi < 0 )
    {
      *i_lo = *i_hi + 1;
    }
  }

  return;
}
/******************************************************************************/

void r8vec_indexed_heap_d ( int n, double a[], int indx[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEXED_HEAP_D creates a descending heap from an indexed R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
    each referencing an entry of the data vector.

    The function adjusts the index vector INDX so that, for 1 <= J <= N/2,
    we have:
      A[INDX[2*J+1]]   <= A[INDX[J]]
    and
      A[INDX[2*J+2]] <= A[INDX[J]]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms for Computers and Calculators,
    Academic Press, 1978,
    ISBN: 0-12-519260-6,
    LC: QA164.N54.

  Parameters:

    Input, int N, the size of the index array.

    Input, double A[*], the data vector.

    Input/output, int INDX[N], the index array.
    Each entry of INDX must be a valid index for the array A.
    On output, the indices have been reordered into a descending heap.
*/
{
  int i;
  int ifree;
  int key;
  int m;
/*
  Only nodes N/2 - 1 down to 0 can be "parent" nodes.
*/
  for ( i = ( n / 2 ) - 1; 0 <= i; i-- )
  {
/*
  Copy the value out of the parent node.
  Position IFREE is now "open".
*/
    key = indx[i];
    ifree = i;

    for ( ; ; )
    {
/*
  Positions 2*IFREE+1 and 2*IFREE+2 are the descendants of position
  IFREE.  (One or both may not exist because they exceed N-1.)
*/
      m = 2 * ifree + 1;
/*
  Does the first position exist?
*/
      if ( n - 1 < m )
      {
        break;
      }
/*
  Does the second position exist?
*/
      if ( m + 1 <= n - 1 )
      {
/*
  If both positions exist, take the larger of the two values,
  and update M if necessary.
*/
        if ( a[indx[m]] < a[indx[m+1]] )
        {
          m = m + 1;
        }
      }
/*
  If the large descendant is larger than KEY, move it up,
  and update IFREE, the location of the free position, and
  consider the descendants of THIS position.
*/
      if ( a[indx[m]] <= a[key] )
      {
        break;
      }

      indx[ifree] = indx[m];
      ifree = m;
    }
/*
  Once there is no more shifting to do, KEY moves into the free spot IFREE.
*/
    indx[ifree] = key;
  }

  return;
}
/******************************************************************************/

int r8vec_indexed_heap_d_extract ( int *n, double a[], int indx[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEXED_HEAP_D_EXTRACT: extract from heap descending indexed R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
    each referencing an entry of the data vector.

    The routine finds the maximum value in the heap, returns that value to the
    user, deletes that value from the heap, and restores the heap to its
    proper form.

    Note that the argument N must be a variable, which will be decremented
    before return, and that INDX will hold one less value on output than it
    held on input.

    This is one of three functions needed to model a priority queue.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt

  Reference:

    Thomas Cormen, Charles Leiserson, Ronald Rivest,
    Introduction to Algorithms,
    MIT Press, 2001,
    ISBN: 0262032937,
    LC: QA76.C662.

  Parameters:

    Input/output, int *N, the number of items in the index vector.

    Input, double A[*], the data vector.

    Input/output, int INDX[N], the index vector.

    Output, int R8VEC_INDEXED_HEAP_D_EXTRACT, the index in A of the item of
    maximum value, which has now been removed from the heap.
*/
{
  int indx_extract;

  if ( *n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_INDEXED_HEAP_D_EXTRACT - Fatal error!\n" );
    fprintf ( stderr, "  The heap is empty.\n" );
    exit ( 1 );
  }
/*
  Get the index of the maximum value.
*/
  indx_extract = indx[0];

  if ( *n == 1 )
  {
    *n = 0;
    return indx_extract;
  }
/*
  Shift the last index down.
*/
  indx[0] = indx[*n-1];
/*
  Restore the heap structure.
*/
  *n = *n - 1;
  r8vec_indexed_heap_d ( *n, a, indx );

  return indx_extract;
}
/******************************************************************************/

void r8vec_indexed_heap_d_insert ( int *n, double a[], int indx[],
  int indx_insert )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEXED_HEAP_D_INSERT: insert value into heap descending indexed R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
    each referencing an entry of the data vector.

    Note that the argument N must be a variable, and will be incremented before
    return, and that INDX must be able to hold one more entry on output than
    it held on input.

    This is one of three functions needed to model a priority queue.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt

  Reference:

    Thomas Cormen, Charles Leiserson, Ronald Rivest,
    Introduction to Algorithms,
    MIT Press, 2001,
    ISBN: 0262032937,
    LC: QA76.C662.

  Parameters:

    Input/output, int *N, the number of items in the index vector.

    Input, double A[*], the data vector.

    Input/output, int INDX[N], the index vector.

    Input, int INDX_INSERT, the index in A of the value
    to be inserted into the heap.
*/
{
  int i;
  int parent;

  *n = *n + 1;
  i = *n - 1;

  while ( 0 < i )
  {
    parent = ( i - 1 ) / 2;

    if ( a[indx_insert] <= a[indx[parent]] )
    {
      break;
    }

    indx[i] = indx[parent];
    i = parent;
  }

  indx[i] = indx_insert;

  return;
}
/******************************************************************************/

int r8vec_indexed_heap_d_max ( int n, double a[], int indx[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDEXED_HEAP_D_MAX: maximum value in heap descending indexed R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    An indexed R8VEC is an R8VEC of data values, and an R8VEC of N indices,
    each referencing an entry of the data vector.

    This is one of three functions needed to model a priority queue.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 August 2010

  Author:

    John Burkardt

  Reference:

    Thomas Cormen, Charles Leiserson, Ronald Rivest,
    Introduction to Algorithms,
    MIT Press, 2001,
    ISBN: 0262032937,
    LC: QA76.C662.

  Parameters:

    Input, int N, the number of items in the index vector.

    Input, double A[*], the data vector.

    Input, int INDX[N], the index vector.

    Output, int R8VEC_INDEXED_HEAP_D_MAX, the index in A of the maximum value
    in the heap.
*/
{
  int indx_max;

  indx_max = indx[0];

  return indx_max;
}
/******************************************************************************/

void r8vec_indicator ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDICATOR sets an R8VEC to the indicator vector {1,2,3...}.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    07 May 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Output, double A[N], the array.
*/
{
  int i;

  for ( i = 0; i <= n-1; i++ )
  {
    a[i] = ( double ) ( i + 1 );
  }

  return;
}
/******************************************************************************/

double *r8vec_indicator_new ( int n )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INDICATOR_NEW sets an R8VEC to the indicator vector {1,2,3...}.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Output, double R8VEC_INDICATOR_NEW[N], the array.
*/
{
  double *a;
  int i;

  a = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i <= n-1; i++ )
  {
    a[i] = ( double ) ( i + 1 );
  }

  return a;
}
/******************************************************************************/

void r8vec_insert ( int n, double a[], int pos, double value )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INSERT inserts a value into an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the array on input.

    Input/output, double A[N+1], the array.  On input, A is
    assumed to contain only N entries, while on output, A actually
    contains N+1 entries.

    Input, int POS, the position to be assigned the new entry.
    1 <= POS <= N+1.

    Input, double VALUE, the value to be inserted.
*/
{
  int i;

  if ( pos < 1 || n + 1 < pos )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_INSERT - Fatal error!\n" );
    fprintf ( stderr, "  Illegal insertion position = %d\n", n );
    exit ( 1 );
  }
  else
  {
    for ( i = n+1; pos+1 <= i; i-- )
    {
      a[i-1] = a[i-2];
    }

    a[pos-1] = value;
  }

  return;
}
/******************************************************************************/

int r8vec_insignificant ( int n, double r[], double s[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_INSIGNIFICANT determines if an R8VEC is insignificant.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vectors.

    Input, double R[N], the vector to be compared against.

    Input, double S[N], the vector to be compared.

    Output, int R8VEC_INSIGNIFICANT, is TRUE if S is insignificant
    compared to R.
*/
{
  int i;
  double t;
  double tol;
  int value;

  value = 1;

  for ( i = 0; i < n; i++ )
  {
    t = r[i] + s[i];
    tol = r8_epsilon ( ) * r8_abs ( r[i] );

    if ( tol < r8_abs ( r[i] - t ) )
    {
      value = 0;
      break;
    }
  }
  
  return value;
}
/******************************************************************************/

int r8vec_is_int ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_IS_INT is TRUE if an R8VEC is integral.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector

    Output, int R8VEC_IS_INT, is TRUE if every entry of A is an integer.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] != ( double ) ( int ) a[i] )
    {
      return 0;
    }
  }
  return 1;
}
/******************************************************************************/

int r8vec_is_nonnegative ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_IS_NONNEGATIVE is true if all entries in an R8VEC are nonnegative.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double X[N], the vector to be checked.

    Output, int R8VEC_IS_NONNEGATIVE is true if all elements of X
    are nonnegative.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < 0.0 )
    {
      return 0;
    }
  }
  return 1;
}
/******************************************************************************/

int r8vec_is_zero ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_IS_ZERO is true if the entries in an R8VEC are all zero.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double X[N], the vector to be checked.

    Output, int R8VEC_IS_ZERO is true if all N elements of X
    are zero.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] != 0.0 )
    {
      return 0;
    }
  }
  return 1;
}
/******************************************************************************/

double *r8vec_legendre_new ( int n, double a_first, double a_last )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LEGENDRE_NEW creates a vector of Chebyshev spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    17 June 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A_FIRST, A_LAST, the first and last entries.

    Output, double R8VEC_LEGENDRE_NEW[N], a vector of Legendre spaced data.
*/
{
  double *a;
  int i;

  a = legendre_zeros ( n );

  for ( i = 0; i < n; i++ )
  {
    a[i] = ( ( 1.0 - a[i] ) * a_first  
           + ( 1.0 + a[i] ) * a_last ) 
           /   2.0;
  }
  return a;
}
/******************************************************************************/

double *r8vec_linspace_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
 
    In other words, the interval is divided into N-1 even subintervals,
    and the endpoints of intervals are used as the points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the first and last entries.

    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - 1 - i ) * a 
             + ( double ) (         i ) * b ) 
             / ( double ) ( n - 1     );
    }
  }
  return x;
}
/******************************************************************************/

double *r8vec_linspace2_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LINSPACE2_NEW creates a vector of linearly spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

    4 points evenly spaced between 0 and 12 will yield 2, 4, 6, 8, 10.

    In other words, the interval is divided into N+1 even subintervals,
    and the endpoints of internal intervals are used as the points.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 September 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the first and last entries.

    Output, double R8VEC_LINSPACE2_NEW[N], a vector of linearly spaced data.
*/
{
  int i;
  double *x;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  if ( n == 1 )
  {
    x[0] = ( a + b ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( double ) ( n - i     ) * a
             + ( double ) (     i + 1 ) * b ) 
             / ( double ) ( n     + 1 );
    }
  }
  return x;
}
/******************************************************************************/

int r8vec_lt ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_LT == ( A1 < A2 ) for two R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The comparison is lexicographic.

    A1 < A2  <=>                              A1(1) < A2(1) or
                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
                 ...
                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 January 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vectors.

    Input, double A1[N], A2[N], the vectors to be compared.

    Output, int R8VEC_LT, is TRUE if and only if A1 < A2.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] < a2[i] )
    {
      return 1;
    }
    else if ( a2[i] < a1[i] )
    {
      return 0;
    }

  }

  return 0;
}
/******************************************************************************/

void r8vec_mask_print ( int n, double a[], int mask_num, int mask[],
  char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MASK_PRINT prints a masked R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, int MASK_NUM, the number of masked elements.

    Input, int MASK[MASK_NUM], the indices of the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  printf ( "\n" );
  printf ( "  Masked vector printout:\n" );

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );
  for ( i = 0; i < mask_num; i++ )
  {
    printf ( "  %6d  %6d  %12f\n", i, mask[i], a[mask[i]-1] );
  }

  return;
}
/******************************************************************************/

double r8vec_max ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MAX returns the value of the maximum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], a pointer to the first entry of the array.

    Output, double R8VEC_MAX, the value of the maximum element.  This
    is set to 0.0 if N <= 0.
*/
{
  int i;
  //double *r8vec_pointer;
  double value;

  if ( n <= 0 )
  {
    value = 0.0;
    return value;
  }

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
/******************************************************************************/

int r8vec_max_abs_index ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MAX_ABS_INDEX returns the index of the maximum absolute value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, int R8VEC_MAX_INDEX, the index of the entry of largest
    absolute value.
*/
{
  int i;
  int max_index;

  if ( n <= 0 )
  {
    max_index = -1;
  }
  else
  {
    max_index = 0;

    for ( i = 1; i < n; i++ )
    {
      if ( r8_abs ( a[max_index] ) < r8_abs ( a[i] ) )
      {
        max_index = i;
      }
    }
  }

  return max_index;
}
/******************************************************************************/

int r8vec_max_index ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MAX_INDEX returns the index of the maximum value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, int R8VEC_MAX_INDEX, the index of the largest entry.
*/
{
  int i;
  int max_index;

  if ( n <= 0 )
  {
    max_index = -1;
  }
  else
  {
    max_index = 0;

    for ( i = 1; i < n; i++ )
    {
      if ( a[max_index] < a[i] )
      {
        max_index = i;
      }
    }
  }

  return max_index;
}
/******************************************************************************/

double r8vec_mean ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MEAN returns the mean of a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double X[N], the vector whose mean is desired.

    Output, double R8VEC_MEAN, the mean, or average, of the vector entries.
*/
{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ )
  {
    mean = mean + x[i];
  }

  mean = mean / ( double ) n;

  return mean;
}
/******************************************************************************/

double r8vec_median ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MEDIAN returns the median of an unsorted R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Hoare's algorithm is used.  The values of the vector are
    rearranged by this routine.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input/output, double A[N], the array to search.  On output,
    the order of the elements of A has been somewhat changed.

    Output, double R8VEC_MEDIAN, the value of the median of A.
*/
{
  int k;
  double median;

  k = ( n + 1 ) / 2;

  median = r8vec_frac ( n, a, k );

  return median;
}
/******************************************************************************/

double *r8vec_midspace_new ( int n, double a, double b )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIDSPACE_NEW creates a vector of linearly spaced values.

  Discussion:

    An R8VEC is a vector of R8's.

    This function divides the interval [a,b] into n subintervals, and then
    returns the midpoints of those subintervals.

  Example:

    N = 5, A = 10, B = 20
    X = [ 11, 13, 15, 17, 19 ]

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the endpoints of the interval.

    Output, double R8VEC_MIDSPACE_NEW[N], a vector of linearly spaced data.
*/
{
  double *x;
  int i;

  x = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( 2 * n - 2 * i - 1 ) * a 
           + ( double ) (         2 * i + 1 ) * b ) 
           / ( double ) ( 2 * n );
  }

  return x;
}
/******************************************************************************/

double r8vec_min ( int n, double r8vec[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN returns the value of the minimum element in a R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double R8VEC[N], the array to be checked.

    Output, double R8VEC_MIN, the value of the minimum element.
*/
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
/******************************************************************************/

int r8vec_min_index ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN_INDEX returns the index of the minimum value in an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Output, int R8VEC_MIN_INDEX, the index of the smallest entry.
*/
{
  int i;
  int value;

  if ( n <= 0 )
  {
    value = - 1;
  }
  else
  {
    value = 0;
    for ( i = 1; i < n; i++ )
    {
      if ( a[i] < a[value] )
      {
        value = i;
      }
    }
  }
  return value;
}
/******************************************************************************/

double r8vec_min_pos ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIN_POS returns the minimum positive value of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 November 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries.

    Input, double A[N], the array.

    Output, double R8VEC_MIN_POS, the smallest positive entry,
    or R8_HUGE if no entry is positive.
*/
{
  int i;
  double r8_huge = 1.0E+30;
  double value;

  value = r8_huge;

  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < a[i] )
    {
      if ( a[i] < value )
      {
        value = a[i];
      }
    }
  }
  return value;
}
/******************************************************************************/

int r8vec_mirror_next ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_MIRROR_NEXT steps through all sign variations of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    In normal use, the user would set every element of A to be positive.
    The routine will take the input value of A, and output a copy in
    which the signs of one or more entries have been changed.  Repeatedly
    calling the routine with the output from the previous call will generate
    every distinct "variation" of A; that is, all possible sign variations.

    When the output variable DONE is TRUE (or equal to 1), then the
    output value of A_NEW is the last in the series.

    Note that A may have some zero values.  The routine will essentially
    ignore such entries; more exactly, it will not stupidly assume that -0
    is a proper "variation" of 0.

    Also, it is possible to call this routine with the signs of A set
    in any way you like.  The routine will operate properly, but it
    will nonethess terminate when it reaches the value of A in which
    every nonzero entry has negative sign.

    More efficient algorithms using the Gray code seem to require internal
    memory in the routine, which is not one of MATLAB's strong points,
    or the passing back and forth of a "memory array", or the use of
    global variables, or unnatural demands on the user.  This form of
    the routine is about as clean as I can make it.

  Example:

      Input         Output
    ---------    --------------
    A            A         DONE
    ---------    --------  ----
     1  2  3     -1  2  3  false
    -1  2  3      1 -2  3  false
     1 -2  3     -1 -2  3  false
    -1 -2  3      1  2 -3  false
     1  2 -3     -1  2 -3  false
    -1  2 -3      1 -2 -3  false
     1 -2 -3     -1 -2 -3  false
    -1 -2 -3      1  2  3  true

     1  0  3     -1  0  3  false
    -1  0  3      1  0 -3  false
     1  0 -3     -1  0 -3  false
    -1  0 -3      1  0  3  true

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, double A[N], a vector of real numbers.  On
    output, some signs have been changed.

    Output, int R8VEC_MIRROR_NEXT, is TRUE if the input vector A was
    the last element
    in the series (every entry was nonpositive); the output vector is reset
    so that all entries are nonnegative, but presumably the ride is over.
*/
{
  int done;
  int i;
  int positive;
/*
  Seek the first strictly positive entry of A.
*/
  positive = -1;
  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 < a[i] )
    {
      positive = i;
      break;
    }
  }
/*
  If there is no strictly positive entry of A, there is no successor.
*/
  if ( positive == -1 )
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = - a[i];
    }
    done = 1;
    return done;
  }
/*
  Otherwise, negate A up to the positive entry.
*/
  for ( i = 0; i <= positive; i++ )
  {
    a[i] = - a[i];
  }
  done = 0;

  return done;
}
/******************************************************************************/

int r8vec_negative_strict ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NEGATIVE_STRICT: all entries of R8VEC are strictly negative.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vector.

    Input, double A[N], the vector.

    Output, int R8VEC_NEGATIVE_STRICT, is TRUE if every entry of
    A is strictly negative.
*/
{
  int i;
  int value;

  for ( i = 0; i < n; i++ )
  {
    if ( 0 <= a[i] )
    {
      value = 0;
      return value;
    }
  }
  value = 1;
  return value;
}
/******************************************************************************/

double *r8vec_nint ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NINT rounds the entries of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector to be rounded.

    Output, double B[N], the rounded values.
*/
{
  double *b;
  int i;
  int s;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] < 0.0 )
    {
      s = -1;
    }
    else
    {
      s = 1;
    }
    b[i] = ( double ) ( s * ( int ) ( r8_abs ( a[i] ) + 0.5 ) );
  }

  return b;
}
/******************************************************************************/

double r8vec_norm ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM returns the L2 norm of an R8VEC.

  Discussion:

    The vector L2 norm is defined as:

      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector whose L2 norm is desired.

    Output, double R8VEC_NORM, the L2 norm of A.
*/
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}
/******************************************************************************/

double r8vec_norm_affine ( int n, double v0[], double v1[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_AFFINE returns the affine L2 norm of an R8VEC.

  Discussion:

    The affine vector L2 norm is defined as:

      R8VEC_NORM_AFFINE(V0,V1)
        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vectors.

    Input, double V0[N], the base vector.

    Input, double V1[N], the vector whose affine L2 norm is desired.

    Output, double R8VEC_NORM_AFFINE, the affine L2 norm of V1.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
  }
  value = sqrt ( value );

  return value;
}
/******************************************************************************/

int r8vec_norm_l0 ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_L0 returns the l0 "norm" of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    The l0 "norm" simply counts the number of nonzero entries in the vector.
    It is not a true norm, but has some similarities to one.  It is useful
    in the study of compressive sensing.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 June 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], the vector.

    Output, int R8VEC_NORM_L0, the value of the norm.
*/
{
  int i;
  int value;

  value = 0;
  for ( i = 0; i < n; i++ )
  {
    if ( a[i] != 0.0 )
    {
      value = value + 1;
    }
  }
  return value;
}
/******************************************************************************/

double r8vec_norm_l1 ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_L1 returns the L1 norm of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    The vector L1 norm is defined as:

      R8VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector whose L1 norm is desired.

    Output, double R8VEC_NORM_L1, the L1 norm of A.
*/
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + r8_abs ( a[i] );
  }

  return v;
}
/******************************************************************************/

double r8vec_norm_l2 ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_L2 returns the L2 norm of an R8VEC.

  Discussion:

    The vector L2 norm is defined as:

      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 March 2003

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector whose L2 norm is desired.

    Output, double R8VEC_NORM_L2, the L2 norm of A.
*/
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}
/******************************************************************************/

double r8vec_norm_li ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_LI returns the L-oo norm of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    The vector L-oo norm is defined as:

      R8VEC_NORM_LI = max ( 1 <= I <= N ) abs ( A(I) ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector whose L-oo norm is desired.

    Output, double R8VEC_NORM_LI, the L-oo norm of A.
*/
{
  int i;
  double v1;
  double v2;

  v1 = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v2 = r8_abs ( a[i] );

    if ( v1 < v2 )
    {
      v1 = v2;
    }
  }

  return v1;
}
/******************************************************************************/

double r8vec_norm_lp ( int n, double a[], double p )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORM_LP returns the LP norm of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    The vector LP norm is defined as:

      R8VEC_NORM_LP = ( sum ( 1 <= I <= N ) ( abs ( A(I) ) )^P )^(1/P).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    21 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector whose LP norm is desired.

    Input, double P, the index of the norm.

    Output, double R8VEC_NORML_LP, the LP norm of A.
*/
{
  int i;
  double v;

  v = 0.0;

  if ( p == 1.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + r8_abs ( a[i] );
    }
  }
  else if ( p == 2.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + a[i] * a[i];
    }
    v = sqrt ( v );
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + pow ( r8_abs ( a[i] ), p );
    }
    v = pow (  ( double ) v, 1.0 / p );
  }

  return v;
}
/******************************************************************************/

void r8vec_normal ( int n, double b, double c, int *seed, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMAL returns a scaled pseudonormal R8VEC.

  Discussion:

    The scaled normal probability distribution function (PDF) has
    mean A and standard deviation B.

    This routine can generate a vector of values on one call.  It
    has the feature that it should provide the same results
    in the same order no matter how we break up the task.

    Before calling this routine, the user may call RANDOM_SEED
    in order to set the seed of the random number generator.

    The Box-Muller method is used, which is efficient, but
    generates an even number of values each time.  On any call
    to this routine, an even number of new values are generated.
    Depending on the situation, one value may be left over.
    In that case, it is saved for the next call.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    02 February 2005

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values desired.  If N is negative,
    then the code will flush its internal memory; in particular,
    if there is a saved value to be used on the next call, it is
    instead discarded.  This is useful if the user has reset the
    random number seed, for instance.

    Input, double B, C, the mean and standard deviation.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double X[N], a sample of the standard normal PDF.

  Local parameters:

    Local, int MADE, records the number of values that have
    been computed.  On input with negative N, this value overwrites
    the return value of N, so the user can get an accounting of
    how much work has been done.

    Local, double R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.

    Local, int SAVED, is 0 or 1 depending on whether there is a
    single saved value left over from the previous call.

    Local, int X_LO, X_HI, records the range of entries of
    X that we need to compute.  This starts off as 1:N, but is adjusted
    if we have a saved value that can be immediately stored in X(1),
    and so on.

    Local, double Y, the value saved from the previous call, if
    SAVED is 1.
*/
{
# define R8_PI 3.141592653589793

  int i;
  int m;
  static int made = 0;
  double *r;
  static int saved = 0;
  int x_hi;
  int x_lo;
  static double y = 0.0;
/*
  I'd like to allow the user to reset the internal data.
  But this won't work properly if we have a saved value Y.
  I'm making a crock option that allows the user to signal
  explicitly that any internal memory should be flushed,
  by passing in a negative value for N.
*/
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return;
  }
  else if ( n == 0 )
  {
    return;
  }
/*
  Record the range of X we need to fill in.
*/
  x_lo = 1;
  x_hi = n;
/*
  Use up the old value, if we have it.
*/
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
/*
  Maybe we don't need any more values.
*/
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
/*
  If we need just one new value, do that here to avoid null arrays.
*/
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * R8_PI * r[1] );
    y =         sqrt ( - 2.0 * log ( r[0] ) ) * sin ( 2.0 * R8_PI * r[1] );

    saved = 1;

    made = made + 2;

    free ( r );
  }
/*
  If we require an even number of values, that's easy.
*/
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    free ( r );
  }
/*
  If we require an odd number of values, we generate an even number,
  and handle the last pair specially, storing one in X(N), and
  saving the other for later.
*/
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );
    }

    i = 2 * m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
    y           = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    free ( r );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = b + c * x[i];
  }

  return;
# undef R8_PI
}
/******************************************************************************/

void r8vec_normal_01 ( int n, int *seed, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

    This routine can generate a vector of values on one call.  It
    has the feature that it should provide the same results
    in the same order no matter how we break up the task.

    Before calling this routine, the user may call RANDOM_SEED
    in order to set the seed of the random number generator.

    The Box-Muller method is used, which is efficient, but
    generates an even number of values each time.  On any call
    to this routine, an even number of new values are generated.
    Depending on the situation, one value may be left over.
    In that case, it is saved for the next call.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values desired.  If N is negative,
    then the code will flush its internal memory; in particular,
    if there is a saved value to be used on the next call, it is
    instead discarded.  This is useful if the user has reset the
    random number seed, for instance.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double X[N], a sample of the standard normal PDF.

  Local parameters:

    Local, int MADE, records the number of values that have
    been computed.  On input with negative N, this value overwrites
    the return value of N, so the user can get an accounting of
    how much work has been done.

    Local, double R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.

    Local, int SAVED, is 0 or 1 depending on whether there is a
    single saved value left over from the previous call.

    Local, int X_LO, X_HI, records the range of entries of
    X that we need to compute.  This starts off as 1:N, but is adjusted
    if we have a saved value that can be immediately stored in X(1),
    and so on.

    Local, double Y, the value saved from the previous call, if
    SAVED is 1.
*/
{
# define R8_PI 3.141592653589793

  int i;
  int m;
  static int made = 0;
  double *r;
  static int saved = 0;
  int x_hi;
  int x_lo;
  static double y = 0.0;
/*
  I'd like to allow the user to reset the internal data.
  But this won't work properly if we have a saved value Y.
  I'm making a crock option that allows the user to signal
  explicitly that any internal memory should be flushed,
  by passing in a negative value for N.
*/
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return;
  }
  else if ( n == 0 )
  {
    return;
  }
/*
  Record the range of X we need to fill in.
*/
  x_lo = 1;
  x_hi = n;
/*
  Use up the old value, if we have it.
*/
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
/*
  Maybe we don't need any more values.
*/
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
/*
  If we need just one new value, do that here to avoid null arrays.
*/
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * R8_PI * r[1] );
    y =         sqrt ( - 2.0 * log ( r[0] ) ) * sin ( 2.0 * R8_PI * r[1] );

    saved = 1;

    made = made + 2;

    free ( r );
  }
/*
  If we require an even number of values, that's easy.
*/
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    free ( r );
  }
/*
  If we require an odd number of values, we generate an even number,
  and handle the last pair specially, storing one in X(N), and
  saving the other for later.
*/
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
    y           = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    free ( r );
  }

  return;
# undef R8_PI
}
/******************************************************************************/

double *r8vec_normal_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.

  Discussion:

    The standard normal probability distribution function (PDF) has
    mean 0 and standard deviation 1.

    This routine can generate a vector of values on one call.  It
    has the feature that it should provide the same results
    in the same order no matter how we break up the task.

    Before calling this routine, the user may call RANDOM_SEED
    in order to set the seed of the random number generator.

    The Box-Muller method is used, which is efficient, but
    generates an even number of values each time.  On any call
    to this routine, an even number of new values are generated.
    Depending on the situation, one value may be left over.
    In that case, it is saved for the next call.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 February 2012

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of values desired.  If N is negative,
    then the code will flush its internal memory; in particular,
    if there is a saved value to be used on the next call, it is
    instead discarded.  This is useful if the user has reset the
    random number seed, for instance.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.

  Local parameters:

    Local, int MADE, records the number of values that have
    been computed.  On input with negative N, this value overwrites
    the return value of N, so the user can get an accounting of
    how much work has been done.

    Local, double R[N+1], is used to store some uniform random values.
    Its dimension is N+1, but really it is only needed to be the
    smallest even number greater than or equal to N.

    Local, int SAVED, is 0 or 1 depending on whether there is a
    single saved value left over from the previous call.

    Local, int X_LO, X_HI, records the range of entries of
    X that we need to compute.  This starts off as 1:N, but is adjusted
    if we have a saved value that can be immediately stored in X(1),
    and so on.

    Local, double Y, the value saved from the previous call, if
    SAVED is 1.
*/
{
# define R8_PI 3.141592653589793

  int i;
  int m;
  static int made = 0;
  double *r;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;

  x = ( double * ) malloc ( n * sizeof ( double ) );
/*
  I'd like to allow the user to reset the internal data.
  But this won't work properly if we have a saved value Y.
  I'm making a crock option that allows the user to signal
  explicitly that any internal memory should be flushed,
  by passing in a negative value for N.
*/
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }
/*
  Record the range of X we need to fill in.
*/
  x_lo = 1;
  x_hi = n;
/*
  Use up the old value, if we have it.
*/
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
/*
  Maybe we don't need any more values.
*/
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
/*
  If we need just one new value, do that here to avoid null arrays.
*/
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * R8_PI * r[1] );
    y =         sqrt ( - 2.0 * log ( r[0] ) ) * sin ( 2.0 * R8_PI * r[1] );

    saved = 1;

    made = made + 2;

    free ( r );
  }
/*
  If we require an even number of values, that's easy.
*/
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    free ( r );
  }
/*
  If we require an odd number of values, we generate an even number,
  and handle the last pair specially, storing one in X(N), and
  saving the other for later.
*/
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
    y           = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    free ( r );
  }

  return x;
# undef R8_PI
}
/******************************************************************************/

void r8vec_normalize ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMALIZE normalizes an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    11 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, double A[N], the vector to be normalized.
    On output, A should have unit Euclidean norm.
*/
{
  int i;
  double norm;

  norm = 0.0;
  for ( i = 0; i < n; i++ )
  {
    norm = norm + a[i] * a[i];
  }
  norm = sqrt ( norm );

  if ( norm == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_NORMALIZE - Fatal error!\n" );
    fprintf ( stderr, "  The vector norm is 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / norm;
  }

  return;
}
/******************************************************************************/

void r8vec_normalize_l1 ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMALIZE_L1 normalizes an R8VEC to have unit sum.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, double A[N], the vector to be normalized.
    On output, the entries of A should have unit sum.  However, if
    the input vector has zero sum, the routine halts.
*/
{
  double a_sum;
  int i;

  a_sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    a_sum = a_sum + a[i];
  }

  if ( a_sum == 0.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_NORMALIZE_L1 - Fatal error!\n" );
    fprintf ( stderr, "  The vector entries sum to 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = a[i] / a_sum;
  }

  return;
}
/******************************************************************************/

double r8vec_normsq ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMSQ returns the squared L2 norm of an R8VEC.

  Discussion:

    The squared vector L2 norm is defined as:

      R8VEC_NORMSQ =  sum ( 1 <= I <= N ) A(I)^2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the vector dimension.

    Input, double A[N], the vector.

    Output, double R8VEC_NORMSQ, the squared L2 norm.
*/
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  return v;
}
/******************************************************************************/

double r8vec_normsq_affine ( int n, double v0[], double v1[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_NORMSQ_AFFINE returns the sqaured affine L2 norm of an R8VEC.

  Discussion:

    The sqaured affine vector L2 norm is defined as:

      R8VEC_NORMSQ_AFFINE(V0,V1)
        = sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vectors.

    Input, double V0[N], the base vector.

    Input, double V1[N], the vector.

    Output, double R8VEC_NORMSQ_AFFINE, the squared affine L2 norm.
*/
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
  }
  return value;
}
/******************************************************************************/

double *r8vec_ones_new ( int n )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ONES_NEW creates a vector of 1's.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    14 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Output, double R8VEC_ONES_NEW[N], a vector of 1's.
*/
{
  double *a;
  int i;

  a = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    a[i] = 1.0;
  }
  return a;
}
/******************************************************************************/

int r8vec_order_type ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ORDER_TYPE determines if an R8VEC is (non)strictly ascending/descending.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the array.

    Input, double X[N], the array to be checked.

    Output, int R8VEC_ORDER_TYPE, order indicator:
    -1, no discernable order;
    0, all entries are equal;
    1, ascending order;
    2, strictly ascending order;
    3, descending order;
    4, strictly descending order.
*/
{
  int i;
  int order;
/*
  Search for the first value not equal to X(0).
*/
  i = 0;

  for ( ; ; )
  {
    i = i + 1;
    if ( n-1 < i )
    {
      order = 0;
      return order;
    }

    if ( x[0] < x[i] )
    {
      if ( i == 1 )
      {
        order = 2;
        break;
      }
      else
      {
        order = 1;
        break;
      }
    }
    else if ( x[i] < x[0] )
    {
      if ( i == 1 )
      {
        order = 4;
        break;
      }
      else
      {
        order = 3;
        break;
      }
    }
  }
/*
  Now we have a "direction".  Examine subsequent entries.
*/
  for ( ; ; )
  {
    i = i + 1;
    if ( n - 1 < i )
    {
      break;
    }

    if ( order == 1 )
    {
      if ( x[i] < x[i-1] )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 2 )
    {
      if ( x[i] < x[i-1] )
      {
        order = -1;
        break;
      }
      else if ( x[i] == x[i-1] )
      {
        order = 1;
      }
    }
    else if ( order == 3 )
    {
      if ( x[i-1] < x[i] )
      {
        order = -1;
        break;
      }
    }
    else if ( order == 4 )
    {
      if ( x[i-1] < x[i] )
      {
        order = -1;
        break;
      }
      else if ( x[i] == x[i-1] )
      {
        order = 3;
      }
    }
  }
  return order;
}
/******************************************************************************/

void r8vec_part_quick_a ( int n, double a[], int *l, int *r )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PART_QUICK_A reorders an R8VEC as part of a quick sort.

  Discussion:

    An R8VEC is a vector of R8's.

    The routine reorders the entries of A.  Using A[0] as a
    key, all entries of A that are less than or equal to A[0] will
    precede A[0] which precedes all entries that are greater than A[0].

  Example:

    Input:

  N = 8

  A = ( 6, 7, 3, 1, 6, 8, 2, 9 )

    Output:

  L = 3, R = 6

  A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
        -------        -------

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of A.

    Input/output, double A[N].  On input, the array to be checked.
    On output, A has been reordered as described above.

    Output, int L, R, the indices of A that define the three segments.
    Let KEY = the input value of A[0].  Then
    I <= L             A(I) < KEY;
     L < I < R         A(I) = KEY;
             R <= I    A(I) > KEY.
*/
{
  int i;
  double key;
  int m;
  double temp;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_PART_QUICK_A - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }
  else if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key = a[0];
  m = 1;
/*
  The elements of unknown size have indices between L+1 and R-1.
*/
  *l = 1;
  *r = n + 1;

  for ( i = 2; i <= n; i++ )
  {

    if ( key < a[*l] )
    {
      *r = *r - 1;
      temp = a[*r-1];
      a[*r-1] = a[*l];
      a[*l] = temp;
    }
    else if ( a[*l] == key )
    {
      m = m + 1;
      temp = a[m-1];
      a[m-1] = a[*l];
      a[*l] = temp;
      *l = *l + 1;
    }
    else if ( a[*l] < key )
    {
      *l = *l + 1;
    }

  }
/*
  Now shift small elements to the left, and KEY elements to center.
*/
  for ( i = 1; i <= *l -m; i++ )
  {
    a[i-1] = a[i+m-1];
  }

  *l = *l - m;

  for ( i = *l+1; i <= *l+m; i++ )
  {
    a[i-1] = key;
  }

  return;
}
/******************************************************************************/

void r8vec_permute ( int n, int p[], int base, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PERMUTE permutes an R8VEC in place.

  Discussion:

    An R8VEC is a vector of R8's.

    This routine permutes an array of real "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.

  Example:

    Input:

      N = 5
      P = (   2,   4,   5,   1,   3 )
      A = ( 1.0, 2.0, 3.0, 4.0, 5.0 )
      BASE = 1

    Output:

      A    = ( 2.0, 4.0, 5.0, 1.0, 3.0 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of objects.

    Input, int P[N], the permutation.

    Input, int BASE, is 0 for a 0-based permutation and 1 for a 1-based permutation.

    Input/output, double A[N], the array to be permuted.
*/
{
  double a_temp;
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm_check ( n, p, base ) )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_PERMUTE - Fatal error!\n" );
    fprintf ( stderr, "  PERM_CHECK rejects this permutation.\n" );
    exit ( 1 );
  }
/*
  In order for the sign negation trick to work, we need to assume that the
  entries of P are strictly positive.  Presumably, the lowest number is BASE.
  So temporarily add 1-BASE to each entry to force positivity.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }
/*
  Search for the next element of the permutation that has not been used.
*/
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp = a[istart-1];
      iget = istart;
/*
  Copy the new value into the vacated entry.
*/
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          fprintf ( stderr, "\n" );
          fprintf ( stderr, "R8VEC_PERMUTE - Fatal error!\n" );
          fprintf ( stderr, "  A permutation index is out of range.\n" );
          fprintf ( stderr, "  P(%d) = %d\n", iput, iget );
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[iput-1] = a_temp;
          break;
        }
        a[iput-1] = a[iget-1];
      }
    }
  }
/*
  Restore the signs of the entries.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
/*
  Restore the base of the entries.
*/
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 +  base;
  }
  return;
}
/******************************************************************************/

void r8vec_permute_cyclic ( int n, int k, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PERMUTE_CYCLIC performs a cyclic permutation of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    For 0 <= K < N, this function cyclically permutes the input vector
    to have the form

     ( A[K], A[K+1], ..., A[N-1], A[0], ..., A[K-1] )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of objects.

    Input, int K, the increment used.

    Input/output, double A[N], the array to be permuted.
*/
{
  double *b;
  int i;
  int ipk;

  b = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    ipk = i4_wrap ( i + k, 0, n - 1 );
    b[i] = a[ipk];
  }

  for ( i = 0; i < n; i++ )
  {
    a[i] = b[i];
  }

  free ( b );

  return;
}
/******************************************************************************/

void r8vec_permute_uniform ( int n, double a[], int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PERMUTE_UNIFORM randomly permutes an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of objects.

    Input/output, double A[N], the array to be permuted.

    Input/output, int *SEED, a seed for the random number generator.
*/
{
  int base = 0;
  int *p;

  p = perm_uniform_new ( n, base, seed );

  r8vec_permute ( n, p, base, a );

  free ( p );

  return;
}
/******************************************************************************/

void r8vec_polarize ( int n, double a[], double p[], double a_normal[],
  double a_parallel[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_POLARIZE decomposes an R8VEC into normal and parallel components.

  Discussion:

    An R8VEC is a vector of R8's.

    The (nonzero) vector P defines a direction.

    The vector A can be written as the sum

      A = A_normal + A_parallel

    where A_parallel is a linear multiple of P, and A_normal
    is perpendicular to P.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the vector to be polarized.

    Input, double P[N], the polarizing direction.

    Output, double A_NORMAL[N], A_PARALLEL[N], the normal
    and parallel components of A.
*/
{
  double a_dot_p;
  int i;
  double p_norm;

  p_norm = 0.0;
  for ( i = 0; i < n; i++ )
  {
    p_norm = p_norm + pow ( p[i], 2 );
  }
  p_norm = sqrt ( p_norm );

  if ( p_norm == 0.0 )
  {
    for ( i = 0; i < n; i++ )
    {
      a_normal[i] = a[i];
    }
    for ( i = 0; i < n; i++ )
    {
      a_parallel[i] = 0.0;
    }
    return;
  }
  a_dot_p = 0.0;
  for ( i = 0; i < n; i++ )
  {
    a_dot_p = a_dot_p + a[i] * p[i];
  }
  a_dot_p = a_dot_p / p_norm;

  for ( i = 0; i < n; i++ )
  {
    a_parallel[i] = a_dot_p * p[i] / p_norm;
  }

  for ( i = 0; i < n; i++ )
  {
    a_normal[i] = a[i] - a_parallel[i];
  }

  return;
}
/******************************************************************************/

int r8vec_positive_strict ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_POSITIVE_STRICT: all entries of R8VEC are strictly positive.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 June 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the dimension of the vector.

    Input, double A[N], the vector.

    Output, int R8VEC_POSITIVE_STRICT, is TRUE if every entry of
    A is strictly positive.
*/
{
  int i;
  int value;

  for ( i = 0; i < n; i++ )
  {
    if ( a[i] <= 0.0 )
    {
      value = 0;
      return value;
    }
  }
  value = 1;
  return value;
}
/******************************************************************************/

void r8vec_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT prints an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 April 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

void r8vec_print_part ( int n, double a[], int max_print, char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT_PART prints "part" of an R8VEC.

  Discussion:

    The user specifies MAX_PRINT, the maximum number of lines to print.

    If N, the size of the vector, is no more than MAX_PRINT, then
    the entire vector is printed, one entry per line.

    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    followed by a line of periods suggesting an omission,
    and the last entry.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 February 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double A[N], the vector to be printed.

    Input, int MAX_PRINT, the maximum number of lines
    to print.

    Input, char *TITLE, a title.
*/
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print - 2; i++ )
    {
      fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
    }
    fprintf ( stdout, "  ......  ..............\n" );
    i = n - 1;
    fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      fprintf ( stdout, "  %8d: %14f\n", i, a[i] );
    }
    i = max_print - 1;
    fprintf ( stdout, "  %8d: %14f  ...more entries...\n", i, a[i] );
  }

  return;
}
/******************************************************************************/

void r8vec_print_some ( int n, double a[], int i_lo, int i_hi, char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRINT_SOME prints "some" of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    16 October 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vector.

    Input, double A[N], the vector to be printed.

    Input, integer I_LO, I_HI, the first and last indices to print.
    The routine expects 1 <= I_LO <= I_HI <= N.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = i4_max ( 1, i_lo ); i <= i4_min ( n, i_hi ); i++ )
  {
    fprintf ( stdout, "  %8d: %14f\n", i, a[i-1] );
  }

  return;
}
/******************************************************************************/

double r8vec_product ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_PRODUCT returns the product of the entries of an R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], the vector.

    Output, double R8VEC_PRODUCT, the product of the vector.
*/
{
  int i;
  double product;

  product = 1.0;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}
/******************************************************************************/

void r8vec_range ( int n, double x[], double xmin, double xmax, double y[],
  double *ymin, double *ymax )

/******************************************************************************/
/*
  Purpose:

    R8VEC_RANGE finds the range of Y's within a restricted X range.

  Discussion:

    An R8VEC is a vector of R8's.

    The routine is given a set of pairs of points (X,Y), and a range
    XMIN to XMAX of valid X values.  Over this range, it seeks
    YMIN and YMAX, the minimum and maximum values of Y for
    valid X's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double X[N], the X array.

    Input, double XMIN, XMAX, the range of X values to check.

    Input, double Y[N], the Y array.

    Output, double *YMIN, *YMAX, the range of Y values whose
    X value is within the X range.
*/
{
  int i;

  *ymin =   r8_huge ( );
  *ymax = - r8_huge ( );

  for ( i = 0; i < n; i++ )
  {
    if ( xmin <= x[i] && x[i] <= xmax )
    {
      *ymin = r8_min ( *ymin, y[i] );
      *ymax = r8_max ( *ymax, y[i] );
    }
  }

  return;
}
/******************************************************************************/

void r8vec_range_2 ( int n, double a[], double *amin, double *amax )

/******************************************************************************/
/*
  Purpose:

    R8VEC_RANGE_2 updates a range to include a new R8VEC

  Discussion:

    An R8VEC is a vector of R8's.

    Given a range AMIN to AMAX, and an array A, the routine will
    decrease AMIN if necessary, or increase AMAX if necessary, so that
    every entry of A is between AMIN and AMAX.

    However, AMIN will not be increased, nor AMAX decreased.

    This routine may be used to compute the maximum and minimum of a
    collection of arrays one at a time.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], the array.

    Input/output, double *AMIN, *AMAX.  On input, the
    current legal range of values for A.  On output, AMIN and AMAX
    are either unchanged, or else "widened" so that all entries
    of A are within the range.
*/
{
  *amax = r8_max ( *amax, r8vec_max ( n, a ) );
  *amin = r8_min ( *amin, r8vec_min ( n, a ) );

  return;
}
/******************************************************************************/

void r8vec_reverse ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_REVERSE reverses the elements of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Example:

    Input:

      N = 5, A = ( 11.0, 12.0, 13.0, 14.0, 15.0 ).

    Output:

      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input/output, double A[N], the array to be reversed.
*/
{
  int i;
  double temp;

  for ( i = 1; i <= n/2; i++ )
  {
    temp   = a[i-1];
    a[i-1] = a[n-i];
    a[n-i] = temp;
  }

  return;
}
/******************************************************************************/

double r8vec_rms ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_RMS returns the RMS norm of an R8VEC.

  Discussion:

    The vector RMS norm is defined as:

      R8VEC_RMS = sqrt ( sum ( 1 <= I <= N ) A(I)^2 / N )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 October 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in A.

    Input, double A[N], the vector whose L2 norm is desired.

    Output, double R8VEC_RMS, the RMS norm of A.
*/
{
  int i;
  double v;

  v = 0.0;

  if ( 0 < n )
  {
    for ( i = 0; i < n; i++ )
    {
      v = v + a[i] * a[i];
    }
    v = sqrt ( v / ( double ) ( n ) );
  }

  return v;
}
/******************************************************************************/

void r8vec_rotate ( int n, double a[], int m )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ROTATE "rotates" the entries of an R8VEC in place.

  Discussion:

    An R8VEC is a vector of R8's.

    This routine rotates an array of real "objects", but the same
    logic can be used to permute an array of objects of any arithmetic
    type, or an array of objects of any complexity.  The only temporary
    storage required is enough to store a single object.  The number
    of data movements made is N + the number of cycles of order 2 or more,
    which is never more than N + N/2.

  Example:

    Input:

      N = 5, M = 2
      A    = ( 1.0, 2.0, 3.0, 4.0, 5.0 )

    Output:

      A    = ( 4.0, 5.0, 1.0, 2.0, 3.0 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of objects.

    Input, int M, the number of positions to the right that
    each element should be moved.  Elements that shift pass position
    N "wrap around" to the beginning of the array.

    Input/output, double A[N], the array to be rotated.
*/
{
  int iget;
  int iput;
  int istart;
  int mcopy;
  int nset;
  double temp;
/*
  Force M to be positive, between 0 and N-1.
*/
  mcopy = i4_modp ( m, n );

  if ( mcopy == 0 )
  {
    return;
  }

  istart = 0;
  nset = 0;

  for ( ; ; )
  {
    istart = istart + 1;

    if ( n < istart )
    {
      break;
    }

    temp = a[istart-1];
    iget = istart;
/*
  Copy the new value into the vacated entry.
*/
    for ( ; ; )
    {
      iput = iget;

      iget = iget - mcopy;
      if ( iget < 1 )
      {
        iget = iget + n;
      }

      if ( iget == istart )
      {
        break;
      }

      a[iput-1] = a[iget-1];
      nset = nset + 1;
    }

    a[iput-1] = temp;
    nset = nset + 1;

    if ( n <= nset )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

double r8vec_scalar_triple_product ( double v1[3], double v2[3], double v3[3] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SCALAR_TRIPLE_PRODUCT computes the scalar triple product.

  Discussion:

    STRIPLE = V1 dot ( V2 x V3 ).

    STRIPLE is the volume of the parallelogram whose sides are
    formed by V1, V2 and V3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double V1[3], V2[3], V3[3], the three vectors.

    Output, double R8VEC_SCALAR_TRIPLE_PRODUCT, the scalar
    triple product.
*/
{
  double value;

  value =
      v1[0] * ( v2[1] * v3[2] - v2[2] * v3[1] )
    + v1[1] * ( v2[2] * v3[0] - v2[0] * v3[2] )
    + v1[2] * ( v2[0] * v3[1] - v2[1] * v3[0] );

  return value;
}
/******************************************************************************/

void r8vec_scale ( double s, int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SCALE multiples an R8VEC by a scale factor.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 September 2011

  Author:

    John Burkardt

  Parameters:

    Input, double S, the scale factor.

    Input, int N, the number of entries in the vectors.

    Input/output, double A[N], the vector to be scaled.
    On output, A[] = S * A[].
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = s * a[i];
  }
  return;
}
/******************************************************************************/

int r8vec_search_binary_a ( int n, double a[], double aval )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SEARCH_BINARY_A searches an ascending sorted R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Binary search is used.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Algorithm 1.9,
    Combinatorial Algorithms,
    CRC Press, 1998, page 26.

  Parameters:

    Input, int N, the number of elements in the array.

    Input, double A[N], the array to be searched.  The array must
    be sorted in ascending order.

    Input, double AVAL, the value to be searched for.

    Output, int R8VEC_SEARCH_BINARY_A, the result of the search.
    -1, AVAL does not occur in the array.
    I, A(I) = AVAL.
*/
{
  int high;
  int indx;
  int low;
  int mid;

  indx = -1;

  low = 1;
  high = n;

  while ( low <= high )
  {
    mid = ( low + high ) / 2;

    if ( a[mid-1] == aval )
    {
      indx = mid;
      break;
    }
    else if ( a[mid-1] < aval )
    {
      low = mid + 1;
    }
    else if ( aval < a[mid-1] )
    {
      high = mid - 1;
    }
  }

  return indx;
}
/******************************************************************************/

void r8vec_shift ( int shift, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SHIFT performs a shift on an R8VEC.

  Discussion:

    An R8VEC is a vector of R8 values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int SHIFT, the amount by which each entry is to
    be shifted.

    Input, int N, the length of the vector.

    Input/output, double X[N], the vector to be shifted.
*/
{
  int i;
  int ihi;
  int ilo;
  //int j;
  double *y;

  y = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[i];
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  ilo = i4_max ( 0, shift );
  ihi = i4_min ( n, n + shift );

  for ( i = ilo; i < ihi; i++ )
  {
    x[i] = y[i-shift];
  }

  free ( y );

  return;
}
/******************************************************************************/

void r8vec_shift_circular ( int shift, int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SHIFT_CIRCULAR performs a circular shift on an R8VEC.

  Discussion:

    An R8VEC is a vector of R8 values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int SHIFT, the amount by which each entry is to
    be shifted.

    Input, int N, the length of the vector.

    Input/output, double X[N], the vector to be shifted.
*/
{
  int i;
  int j;
  double *y;

  y = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    y[i] = x[i];
  }

  for ( i = 0; i < n; i++ )
  {
    j = i4_wrap ( i - shift, 0, n - 1 );
    x[i] = y[j];
  }

  free ( y );

  return;
}
/******************************************************************************/

void r8vec_sort_bubble_a ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_BUBBLE_A ascending sorts an R8VEC using bubble sort.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, length of input array.

    Input/output, double A[N].
    On input, an unsorted array of doubles.
    On output, A has been sorted.
*/
{
  int i;
  int j;
  double temp;

  for ( i = 0; i < n-1; i++ )
  {
    for ( j = i+1; j < n; j++ )
    {
      if ( a[j] < a[i] )
      {
        temp = a[i];
        a[i] = a[j];
        a[j] = temp;
      }
    }
  }
  return;
}
/******************************************************************************/

void r8vec_sort_bubble_d ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_BUBBLE_D descending sorts an R8VEC using bubble sort.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, length of input array.

    Input/output, double A[N].
    On input, an unsorted array of doubles.
    On output, A has been sorted.
*/
{
  int i;
  int j;
  double temp;

  for ( i = 0; i < n-1; i++ )
  {
    for ( j = i+1; j < n; j++ )
    {
      if ( a[i] < a[j] )
      {
        temp = a[i];
        a[i] = a[j];
        a[j] = temp;
      }
    }
  }
  return;
}
/******************************************************************************/

void r8vec_sort_heap_a ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_A ascending sorts an R8VEC using heap sort.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the number of entries in the array.

    Input/output, double A[N].
    On input, the array to be sorted;
    On output, the array has been sorted.
*/
{
  int n1;
  double temp;

  if ( n <= 1 )
  {
    return;
  }
/*
  1: Put A into descending heap form.
*/
  r8vec_heap_d ( n, a );
/*
  2: Sort A.

  The largest object in the heap is in A[0].
  Move it to position A[N-1].
*/
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
/*
  Consider the diminished heap of size N1.
*/
  for ( n1 = n - 1; 2 <= n1; n1-- )
  {
/*
  Restore the heap structure of the initial N1 entries of A.
*/
    r8vec_heap_d ( n1, a );
/*
  Take the largest object from A[0] and move it to A[N1-1].
*/
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }

  return;
}
/******************************************************************************/

void r8vec_sort_heap_d ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_D descending sorts an R8VEC using heap sort.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the number of entries in the array.

    Input/output, double A[N].
    On input, the array to be sorted;
    On output, the array has been sorted.
*/
{
  int n1;
  double temp;

  if ( n <= 1 )
  {
    return;
  }
/*
  1: Put A into ascending heap form.
*/
  r8vec_heap_a ( n, a );
/*
  2: Sort A.

  The smallest object in the heap is in A[0].
  Move it to position A[N-1].
*/
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
/*
  Consider the diminished heap of size N1.
*/
  for ( n1 = n - 1; 2 <= n1; n1-- )
  {
/*
  Restore the heap structure of the initial N1 entries of A.
*/
    r8vec_heap_a ( n1, a );
/*
  Take the largest object from A[0] and move it to A[N1-1].
*/
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }

  return;
}
/******************************************************************************/

void r8vec_sort_heap_index_a ( int n, double a[], int indx[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC

  Discussion:

    An R8VEC is a vector of R8's.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    Once the index array is computed, the sorting can be carried out
    "implicitly:

      a(indx(*))

    or explicitly, by the call

      r8vec_permute ( n, indx, 0, a )

    after which a(*) is sorted.

    Note that the index vector is 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], an array to be index-sorted.

    Output, int INDX[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
  double aval;
  int i;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]] < a[indx[j]] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return;
}
/******************************************************************************/

int *r8vec_sort_heap_index_a_new ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_INDEX_A_NEW does an indexed heap ascending sort of an R8VEC

  Discussion:

    An R8VEC is a vector of R8's.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    Once the index array is computed, the sorting can be carried out
    "implicitly:

      a(indx(*))

    or explicitly, by the call

      r8vec_permute ( n, indx, 0, a )

    after which a(*) is sorted.

    Note that the index vector is 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], an array to be index-sorted.

    Output, int R8VEC_SORT_HEAP_INDEX_A_NEW[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]] < a[indx[j]] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}
/******************************************************************************/

void r8vec_sort_heap_index_d ( int n, double a[], int indx[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_INDEX_D does an indexed heap descending sort of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    Once the index array is computed, the sorting can be carried out
    "implicitly:

      a(indx(*))

    or explicitly, by the call

      r8vec_permute ( n, indx, 0, a )

    after which a(*) is sorted.

    Note that the index vector is 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], an array to be index-sorted.

    Output, int INDX[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
  double aval;
  int i;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j]] < a[indx[j-1]] )
        {
          j = j + 1;
        }
      }

      if ( a[indx[j-1]] < aval )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }

    indx[i-1] = indxt;
  }
  return;
}
/******************************************************************************/

int *r8vec_sort_heap_index_d_new ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_INDEX_D_NEW does an indexed heap descending sort of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    Once the index array is computed, the sorting can be carried out
    "implicitly:

      a(indx(*))

    or explicitly, by the call

      r8vec_permute ( n, indx, 0, a )

    after which a(*) is sorted.

    Note that the index vector is 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], an array to be index-sorted.

    Output, int R8VEC_SORT_HEAP_INDEX_D_NEW[N], contains the sort index.  The
    I-th element of the sorted array is A(INDX(I)).
*/
{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j]] < a[indx[j-1]] )
        {
          j = j + 1;
        }
      }

      if ( a[indx[j-1]] < aval )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }

    indx[i-1] = indxt;
  }
  return indx;
}
/******************************************************************************/

int *r8vec_sort_heap_mask_a ( int n, double a[], int mask_num, int mask[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_HEAP_MASK_A: indexed heap ascending sort of a masked R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    An array A is given.  An array MASK of indices into A is given.
    The routine produces a vector INDX, which is a permutation of the
    entries of MASK, so that:

      A(MASK(INDX(I)) <= A(MASK(INDX(J))

    whenever

      I <= J

    In other words, only the elements of A that are indexed by MASK
    are to be considered, and the only thing that happens is that
    a rearrangment of the indices in MASK is returned that orders the
    masked elements.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], an array to be index-sorted.

    Input, int MASK_NUM, the number of mask elements.

    Input, int MASK[MASK_NUM], the mask array.  This is
    simply a list of indices of A.  The entries of MASK should
    be unique, and each one should be between 1 and N.

    Output, int INDX[MASK_NUM], the sort index.  There are MASK_NUM
    elements of A selected by MASK.  If we want to list those elements
    in order, then the I-th element is A(MASK(INDX(I))).
*/
{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  if ( mask_num < 1 )
  {
    return NULL;
  }

  if ( mask_num == 1 )
  {
    indx = ( int * ) malloc ( 1 * sizeof ( int ) );
    indx[0] = 1;
    return indx;
  }

  indx = i4vec_indicator_new ( mask_num );

  l = mask_num / 2 + 1;
  ir = mask_num;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[mask[indxt-1]-1];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[mask[indxt-1]-1];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[mask[indx[j-1]-1]-1] < a[mask[indx[j]-1]-1] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[mask[indx[j-1]-1]-1] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}
/******************************************************************************/

void r8vec_sort_insert_a ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_INSERT_A ascending sorts an R8VEC using an insertion sort.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Algorithm 1.1,
    Combinatorial Algorithms,
    CRC Press, 1998, page 11.

  Parameters:

    Input, int N, the number of items in the vector.
    N must be positive.

    Input/output, double A[N].

    On input, A contains data to be sorted.
    On output, the entries of A have been sorted in ascending order.
*/
{
  int i;
  int j;
  double x;

  for ( i = 1; i < n; i++ )
  {
    x = a[i];

    j = i;

    while ( 1 <= j && x < a[j-1] )
    {
      a[j] = a[j-1];
      j = j - 1;
    }
    a[j] = x;
  }

  return;
}
/******************************************************************************/

int *r8vec_sort_insert_index_a ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_INSERT_INDEX_A ascending index sorts an R8VEC using insertion.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Reference:

    Donald Kreher, Douglas Simpson,
    Combinatorial Algorithms,
    CRC Press, 1998, page 11.

  Parameters:

    Input, int N, the number of items in the vector.
    N must be positive.

    Input, double A[N], the array to be sorted.

    Output, int R8VEC_SORT_INSET_INDEX_A[N], the sorted indices.  The array
    is sorted when listed from A(INDX(1)) through A(INDX(N)).
*/
{
  int i;
  int *indx;
  int j;
  double x;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = i4vec_indicator_new ( n );

  for ( i = 2; i <= n; i++ )
  {
    x = a[i-1];

    j = i - 1;

    while ( 1 <= j )
    {
      if ( a[indx[j-1]-1] <= x )
      {
        break;
      }

      indx[j] = indx[j-1];
      j = j - 1;
    }
    indx[j] = i;
  }

  return indx;
}
/******************************************************************************/

void r8vec_sort_quick_a ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_QUICK_A ascending sorts an R8VEC using quick sort.

  Discussion:

    An R8VEC is a vector of R8's.

  Example:

    Input:

      N = 7

      A = ( 6, 7, 3, 2, 9, 1, 8 )

    Output:

      A = ( 1, 2, 3, 6, 7, 8, 9 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of A.

    Input/output, double A[N].  On input, the array to be sorted.
    On output, A has been reordered into ascending order.
*/
{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_SORT_QUICK_A - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }
  else if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[0] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {
/*
  Partition the segment.
*/
    r8vec_part_quick_a ( n_segment, a+base-1, &l_segment, &r_segment );
/*
  If the left segment has more than one element, we need to partition it.
*/
    if ( 1 < l_segment )
    {

      if ( LEVEL_MAX < level )
      {
        fprintf ( stderr, "\n" );
        fprintf ( stderr, "R8VEC_SORT_QUICK_A - Fatal error!\n" );
        fprintf ( stderr, "  Exceeding recursion maximum of %d\n", LEVEL_MAX );
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
/*
  The left segment and the middle segment are sorted.
  Must the right segment be partitioned?
*/
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
/*
  Otherwise, we back up a level if there is an earlier one.
*/
    else
    {
      for ( ; ; )
      {
        if ( 1 < level )
        {
          base = rsave[level-1];
          n_segment = rsave[level-2] - rsave[level-1];
          level = level - 1;
          if ( 0 < n_segment )
          {
            break;
          }
        }
        else
        {
          n_segment = 0;
          break;
        }
      }
    }
  }

  return;
# undef LEVEL_MAX
}
/******************************************************************************/

void r8vec_sort_shell_a ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORT_SHELL_A ascending sorts an R8VEC using Shell's sort.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input/output, double A[N].
    On input, an array to be sorted.
    On output, the sorted array.
*/
{
  double asave;
  int i;
  int ifree;
  int inc;
  int ipow;
  int j;
  int k;
  int maxpow;
  int test;

  if ( n <= 1 )
  {
    return;
  }
/*
  Determine the smallest MAXPOW so that
    N <= ( 3**MAXPOW - 1 ) / 2
*/
  maxpow = 1;
  test = 3;

  while ( test < 2 * n + 1 )
  {
    maxpow = maxpow + 1;
    test = test * 3;
  }

  if ( 1 < maxpow )
  {
    maxpow = maxpow - 1;
    test = test / 3;
  }
/*
  Now sort groups of size ( 3**IPOW - 1 ) / 2.
*/
  for ( ipow = maxpow; 1 <= ipow; ipow-- )
  {
    inc = ( test - 1 ) / 2;
    test = test / 3;
/*
  Sort the values with indices equal to K mod INC.
*/
    for ( k = 1; k <= inc; k++ )
    {
/*
  Insertion sort of the items with index
  INC+K, 2*INC+K, 3*INC+K, ...
*/
      for ( i = inc+k; i <= n; i = i + inc )
      {
        asave = a[i-1];
        ifree = i;
        j = i - inc;

        for ( ; ; )
        {
          if ( j < 1 )
          {
            break;
          }

          if ( a[j-1] <= asave )
          {
            break;
          }

          ifree = j;
          a[j+inc-1] = a[j-1];
          j = j - inc;
        }
        a[ifree-1] = asave;
      }
    }
  }

  return;
}
/******************************************************************************/

double *r8vec_sorted_merge_a ( int na, double a[], int nb, double b[], int *nc )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_MERGE_A merges two ascending sorted R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

    The elements of A and B should be sorted in ascending order.

    The elements in the output array C will also be in ascending order,
    and unique.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int NA, the dimension of A.

    Input, double A[NA], the first sorted array.

    Input, int NB, the dimension of B.

    Input, double B[NB], the second sorted array.

    Output, int *NC, the number of entries in the merged vector.

    Output, double R8VEC_SORTED_MERGE_A[NC], the merged unique sorted array.
*/
{
  double *c;
  double *d;
  int j;
  int ja;
  int jb;
  int na2;
  int nb2;
  int nd;
  int order;

  na2 = na;
  nb2 = nb;

  ja = 0;
  jb = 0;
  *nc = 0;
  nd = 0;
  d = ( double * ) malloc ( ( na + nb ) * sizeof ( double ) );

  order = r8vec_order_type ( na2, a );

  if ( order < 0 || 2 < order )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_SORTED_MERGE_A - Fatal error!\n" );
    fprintf ( stderr, "  The input array A is not ascending sorted.\n" );
    return NULL;
  }

  order = r8vec_order_type ( nb2, b );

  if ( order < 0 || 2 < order )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_SORTED_MERGE_A - Fatal error!\n" );
    fprintf ( stderr, "  The input array B is not ascending sorted.\n" );
    return NULL;
  }

  for ( ; ; )
  {
/*
  If we've used up all the entries of A, stick the rest of B on the end.
*/
    if ( na2 <= ja )
    {
      for ( j = 1; j <= nb2 - jb; j++ )
      {
        jb = jb + 1;
        if ( nd == 0 )
        {
          nd = nd + 1;
          d[nd-1] = b[jb-1];
        }
        else if ( d[nd-1] < b[jb-1] )
        {
          nd = nd + 1;
          d[nd-1] = b[jb-1];
        }
      }
      break;
    }
/*
  If we've used up all the entries of B, stick the rest of A on the end.
*/
    else if ( nb2 <= jb )
    {
      for ( j = 1; j <= na2 - ja; j++ )
      {
        ja = ja + 1;
        if ( nd == 0 )
        {
          nd = nd + 1;
          d[nd-1] = a[ja-1];
        }
        else if ( d[nd-1] < a[ja-1] )
        {
          nd = nd + 1;
          d[nd-1] = a[ja-1];
        }
      }
      break;
    }
/*
  Otherwise, if the next entry of A is smaller, that's our candidate.
*/
    else if ( a[ja] <= b[jb] )
    {
      ja = ja + 1;
      if ( nd == 0 )
      {
        nd = nd + 1;
        d[nd-1] = a[ja-1];
      }
      else if ( d[nd-1] < a[ja-1] )
      {
        nd = nd + 1;
        d[nd-1] = a[ja-1];
      }
    }
/*
  ...or if the next entry of B is the smaller, consider that.
*/
    else
    {
      jb = jb + 1;
      if ( nd == 0 )
      {
        nd = nd + 1;
        d[nd-1] = b[jb-1];
      }
      else if ( d[nd-1] < b[jb-1] )
      {
        nd = nd + 1;
        d[nd-1] = b[jb-1];
      }
    }
  }

  *nc = nd;

  c = r8vec_copy_new ( nd, d );

  free ( d );

  return c;
}
/******************************************************************************/

int r8vec_sorted_nearest ( int n, double a[], double value )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_NEAREST returns the nearest element in a sorted R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input, double A[N], a sorted vector.

    Input, double VALUE, the value whose nearest vector entry is sought.

    Output, int R8VEC_SORTED_NEAREST, the index of the nearest
    entry in the vector.
*/
{
  int hi;
  int lo;
  int mid;

  if ( n < 1 )
  {
    return (-1);
  }

  if ( n == 1 )
  {
    return 1;
  }

  if ( a[0] < a[n-1] )
  {
    if ( value < a[0] )
    {
      return 1;
    }
    else if ( a[n-1] < value )
    {
      return n;
    }
/*
  Seek an interval containing the value.
*/
    lo = 1;
    hi = n;

    while ( lo < hi - 1 )
    {
      mid = ( lo + hi ) / 2;

      if ( value == a[mid-1] )
      {
        return mid;
      }
      else if ( value < a[mid-1] )
      {
        hi = mid;
      }
      else
      {
        lo = mid;
      }
    }
/*
  Take the nearest.
*/
    if ( r8_abs ( value - a[lo-1] ) < r8_abs ( value - a[hi-1] ) )
    {
      return lo;
    }
    else
    {
      return hi;
    }
  }
/*
  A descending sorted vector A.
*/
  else
  {
    if ( value < a[n-1] )
    {
      return n;
    }
    else if ( a[0] < value )
    {
      return 1;
    }
/*
  Seek an interval containing the value.
*/
    lo = n;
    hi = 1;

    while ( lo < hi - 1 )
    {
      mid = ( lo + hi ) / 2;

      if ( value == a[mid-1] )
      {
        return mid;
      }
      else if ( value < a[mid-1] )
      {
        hi = mid;
      }
      else
      {
        lo = mid;
      }
    }
/*
  Take the nearest.
*/
    if ( r8_abs ( value - a[lo-1] ) < r8_abs ( value - a[hi-1] ) )
    {
      return lo;
    }
    else
    {
      return hi;
    }
  }
}
/******************************************************************************/

void r8vec_sorted_range ( int n, double r[], double r_lo, double r_hi,
  int *i_lo, int *i_hi )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_RANGE searches a sorted vector for elements in a range.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 September 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of items in the vector.

    Input, double R[N], the sorted vector.

    Input, double R_LO, R_HI, the limits of the range.

    Output, int *I_LO, *I_HI, the range of indices
    so that I_LO <= I <= I_HI => R_LO <= R(I) <= R_HI.  If no
    values in R lie in the range, then I_HI < I_LO will be returned.
*/
{
  int i1;
  int i2;
  int j1;
  int j2;
/*
  Cases we can handle immediately.
*/
  if ( r[n-1] < r_lo )
  {
    *i_lo = - 1;
    *i_hi = - 2;
    return;
  }

  if ( r_hi < r[0] )
  {
    *i_lo = - 1;
    *i_hi = - 2;
    return;
  }
/*
  Are there are least two intervals?
*/
  if ( n == 1 )
  {
    if ( r_lo <= r[0] && r[0] <= r_hi )
    {
      *i_lo = 1;
      *i_hi = 1;
    }
    else
    {
      *i_lo = - 1;
      *i_hi = - 2;
    }
    return;
  }
/*
  Bracket R_LO.
*/
  if ( r_lo <= r[0] )
  {
    *i_lo = 0;
  }
  else
  {
/*
  R_LO is in one of the intervals spanned by R(J1) to R(J2).
  Examine the intermediate interval [R(I1), R(I1+1)].
  Does R_LO lie here, or below or above?
*/
    j1 = 0;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_lo < r[i1] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[i2] < r_lo )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        *i_lo = i1;
        break;
      }
    }
  }
/*
  Bracket R_HI
*/
  if ( r[n-1] <= r_hi )
  {
    *i_hi = n - 1;
  }
  else
  {
    j1 = *i_lo;
    j2 = n - 1;
    i1 = ( j1 + j2 - 1 ) / 2;
    i2 = i1 + 1;

    for ( ; ; )
    {
      if ( r_hi < r[i1] )
      {
        j2 = i1;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else if ( r[i2] < r_hi )
      {
        j1 = i2;
        i1 = ( j1 + j2 - 1 ) / 2;
        i2 = i1 + 1;
      }
      else
      {
        *i_hi = i2;
        break;
      }
    }
  }
/*
  We expect to have computed the largest I_LO and smallest I_HI such that
    R(I_LO) <= R_LO <= R_HI <= R(I_HI)
  but what we want is actually
    R_LO <= R(I_LO) <= R(I_HI) <= R_HI
  which we can usually get simply by incrementing I_LO and decrementing I_HI.
*/
  if ( r[*i_lo] < r_lo )
  {
    *i_lo = *i_lo + 1;
    if ( n - 1 < *i_lo )
    {
      *i_hi = *i_lo - 1;
    }
  }

  if ( r_hi < r[*i_hi] )
  {
    *i_hi = *i_hi - 1;
    if ( *i_hi < 0 )
    {
      *i_lo = *i_hi + 1;
    }
  }

  return;
}
/******************************************************************************/

void r8vec_sorted_split ( int n, double a[], double split, int *i_lt,
  int *i_gt )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_SPLIT "splits" a sorted R8VEC, given a splitting value.

  Discussion:

    An R8VEC is a vector of R8's.

    Given a splitting value SPLIT, the routine seeks indices
    I_LT and I_GT so that

      A(I_LT) < SPLIT < A(I_GT),

    and if there are intermediate index values between I_LT and
    I_GT, then those entries of A are exactly equal to SPLIT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters

    Input, int N, the number of entries in A.

    Input, double A[N], a sorted array.

    Input, double SPLIT, a value to which the entries in A are
    to be compared.

    Output, int *I_LT:
    0 if no entries are less than SPLIT;
    N if all entries are less than SPLIT;
    otherwise, the index of the last entry in A less than SPLIT.

    Output, int *I_GT:
    1 if all entries are greater than SPLIT;
    N+1 if no entries are greater than SPLIT;
    otherwise the index of the first entry in A greater than SPLIT.
*/
{
  int hi;
  int i;
  int lo;
  int mid;

  if ( n < 1 )
  {
    *i_lt = -1;
    *i_gt = -1;
    return;
  }

  if ( split < a[0] )
  {
    *i_lt = 0;
    *i_gt = 1;
    return;
  }

  if ( a[n-1] < split )
  {
    *i_lt = n;
    *i_gt = n + 1;
    return;
  }

  lo = 1;
  hi = n;

  for ( ; ; )
  {
    if ( lo + 1 == hi )
    {
      *i_lt = lo;
      break;
    }

    mid = ( lo + hi ) / 2;

    if ( split <= a[mid-1] )
    {
      hi = mid;
    }
    else
    {
      lo = mid;
    }
  }

  for ( i = *i_lt + 1; i <= n; i++ )
  {
    if ( split < a[i-1] )
    {
      *i_gt = i;
      return;
    }
  }

  *i_gt = n + 1;

  return;
}
/******************************************************************************/

void r8vec_sorted_undex ( int x_num, double x_val[], int x_unique_num,
  double tol, int undx[], int xdnu[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_UNDEX returns unique sorted indexes for a sorted R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    The goal of this routine is to determine a vector UNDX,
    which points, to the unique elements of X, in sorted order,
    and a vector XDNU, which identifies, for each entry of X, the index of
    the unique sorted element of X.

    This is all done with index vectors, so that the elements of
    X are never moved.

    Assuming X is already sorted, we examine the entries of X in order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.

    Once this process has been completed, the vector X could be
    replaced by a compressed vector XU, containing the unique entries
    of X in sorted order, using the formula

      XU(I) = X(UNDX(I)).

    We could then, if we wished, reconstruct the entire vector X, or
    any element of it, by index, as follows:

      X(I) = XU(XDNU(I)).

    We could then replace X by the combination of XU and XDNU.

    Later, when we need the I-th entry of X, we can locate it as
    the XDNU(I)-th entry of XU.

    Here is an example of a vector X, the sort and inverse sort
    index vectors, and the unique sort and inverse unique sort vectors
    and the compressed unique sorted vector.

      I      X      XU  Undx  Xdnu
    ----+------+------+-----+-----+
      0 | 11.0 |  11.0    0     0
      1 | 11.0 |  22.0    4     0
      2 | 11.0 |  33.0    7     0
      3 | 11.0 |  55.0    8     0
      4 | 22.0 |                1
      5 | 22.0 |                1
      6 | 22.0 |                1
      7 | 33.0 |                2
      8 | 55.0 |                3

    INDX(2) = 3 means that sorted item(2) is X(3).
    XDNI(2) = 5 means that X(2) is sorted item(5).

    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
    XDNU(8) = 2 means that X(8) is at unique sorted item(2).

    XU(XDNU(I))) = X(I).
    XU(I)        = X(UNDX(I)).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, the number of data values.

    Input, double X_VAL[X_NUM], the data values.

    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.

    Input, double TOL, a tolerance for equality.

    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.

    Output, int XDNU[X_NUM], the XDNU vector.
*/
{
  int i;
  int j;
/*
  Walk through the sorted array.
*/
  i = 0;

  j = 0;
  undx[j] = i;

  xdnu[i] = j;

  for ( i = 1; i < x_num; i++ )
  {
    if ( tol < r8_abs ( x_val[i] - x_val[undx[j]] ) )
    {
      j = j + 1;
      undx[j] = i;
    }
    xdnu[i] = j;
  }

  return;
}
/******************************************************************************/

double *r8vec_sorted_unique ( int n, double a[], double tol, int *unique_num )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_UNIQUE finds the unique elements in a sorted R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    If the data is not sorted, the results of the routine will
    be garbage.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input, double A[N], the sorted array of N elements;

    Input, double TOL, a tolerance for checking equality.

    Output, int *UNIQUE_NUM, the number of unique elements of A.

    Output, double R8VEC_SORTED_UNIQUE[UNIQUE_NUM], the unique elements of A.
*/
{
  double *a_unique;
  int i;
  int iuniq;

  *unique_num = 0;

  if ( n <= 0 )
  {
    return NULL;
  }
/*
  Determine the number of unique elements.
*/
  iuniq = 0;
  *unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( tol < r8_abs ( a[i] - a[iuniq] ) )
    {
       iuniq = i;
      *unique_num = *unique_num + 1;
    }
  }
/*
  Set aside space for the unique elements.
*/
  a_unique = ( double * ) malloc ( *unique_num * sizeof ( double ) );
/*
  Repeat the search, but now store the unique elements.
*/
  *unique_num = 0;

  a_unique[*unique_num] = a[0];
  *unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( tol < r8_abs ( a[i] - a_unique[*unique_num-1] ) )
    {
      a_unique[*unique_num] = a[i];
      *unique_num = *unique_num + 1;
    }
  }

  return a_unique;
}
/******************************************************************************/

int r8vec_sorted_unique_count ( int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_UNIQUE_COUNT counts unique elements in a sorted R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Because the array is sorted, this algorithm is O(N).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input, double A[N], the sorted array to examine.

    Input, double TOL, a tolerance for checking equality.

    Output, int R8VEC_SORTED_UNIQUE_COUNT, the number of unique elements of A.
*/
{
  int i;
  int unique_num;

  unique_num = 0;

  if ( n < 1 )
  {
    return unique_num;
  }

  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( tol < r8_abs ( a[i-1] - a[i] ) )
    {
      unique_num = unique_num + 1;
    }
  }

  return unique_num;
}
/******************************************************************************/

void r8vec_sorted_unique_hist ( int n, double a[], double tol, int maxuniq,
  int *unique_num, double auniq[], int acount[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SORTED_UNIQUE_HIST histograms unique elements of a sorted R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input, double A[N], the array to examine, which must have been
    sorted.

    Input, double TOL, a tolerance for checking equality.

    Input, int MAXUNIQ, the maximum number of unique elements
    that can be handled.  If there are more than MAXUNIQ unique
    elements in A, the excess will be ignored.

    Output, int *UNIQUE_NUM, the number of unique elements of A.

    Output, double AUNIQ[UNIQUE_NUM], the unique elements of A.

    Output, int ACOUNT[UNIQUE_NUM], the number of times each element
    of AUNIQ occurs in A.
*/
{
  int i;
  int index;
/*
  Start taking statistics.
*/
  index = -1;

  for ( i = 0; i < n; i++ )
  {

    if ( i == 0 )
    {
      index = 0;
      auniq[index] = a[0];
      acount[index] = 1;
    }
    else if ( r8_abs ( a[i] - auniq[index] ) <= tol )
    {
      acount[index] = acount[index] + 1;
    }
    else if ( index + 1 < maxuniq )
    {
      index = index + 1;
      auniq[index] = a[i];
      acount[index] = 1;
    }
  }

  *unique_num = index + 1;

  return;
}
/******************************************************************************/

int r8vec_split ( int n, double a[], double split )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SPLIT "splits" an unsorted R8VEC based on a splitting value.

  Discussion:

    An R8VEC is a vector of R8's.

    If the vector is already sorted, it is simpler to do a binary search
    on the data than to call this routine.

    The vector is not assumed to be sorted before input, and is not
    sorted during processing.  If sorting is not needed, then it is
    more efficient to use this routine.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input/output, double A[N], the array to split.  On output,
    all the entries of A that are less than or equal to SPLIT
    are in A(1:ISPLIT).

    Input, double SPLIT, the value used to split the vector.
    It is not necessary that any value of A actually equal SPLIT.

    Output, int R8VEC_SPLIT, indicates the position of the last
    entry of the split vector that is less than or equal to SPLIT.
*/
{
  int i;
  int i1;
  int i2;
  int i3;
  int isplit;
  int j1;
  int j2;
  int j3;
  double temp;
/*
  Partition the vector into A1, A2, A3, where
    A1 = A(I1:J1) holds values <= SPLIT,
    A2 = A(I2:J2) holds untested values,
    A3 = A(I3:J3) holds values > SPLIT.
*/
  i1 = 1;
  j1 = 0;

  i2 = 1;
  j2 = n;

  i3 = n + 1;
  j3 = n;
/*
  Pick the next item from A2, and move it into A1 or A3.
  Adjust indices appropriately.
*/
  for ( i = 1; i <= n; i++ )
  {
    if ( a[i2-1] <= split )
    {
      i2 = i2 + 1;
      j1 = j1 + 1;
    }
    else
    {
      temp = a[i2-1];
      a[i2-1] = a[i3-2];
      a[i3-2] = temp;
      i3 = i3 - 1;
      j2 = j2 - 1;
    }
  }

  isplit = j1;

  return isplit;
}
/******************************************************************************/

double r8vec_std ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_STD returns the standard deviation of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    The standard deviation of a vector X of length N is defined as

      mean ( X(1:n) ) = sum ( X(1:n) ) / n

      std ( X(1:n) ) = sqrt ( sum ( ( X(1:n) - mean )^2 ) / ( n - 1 ) )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.
    N should be at least 2.

    Input, double A[N], the vector.

    Output, double R8VEC_STD, the standard deviation of the vector.
*/
{
  int i;
  double mean;
  double std;

  if ( n < 2 )
  {
    std = 0.0;
  }
  else
  {
    mean = 0.0;
    for ( i = 0; i < n; i++ )
    {
      mean = mean + a[i];
    }
    mean = mean / ( ( double ) n );

    std = 0.0;
    for ( i = 0; i < n; i++ )
    {
      std = std + ( a[i] - mean ) * ( a[i] - mean );
    }
    std = sqrt ( std / ( ( double ) ( n - 1 ) ) );
  }

  return std;
}
/******************************************************************************/

void r8vec_stutter ( int n, double a[], int m, double am[] )

/******************************************************************************/
/*

  Purpose:

    R8VEC_STUTTER makes a "stuttering" copy of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Applying a stuttering factor M of 3, the vector A = ( 1, 5, 8 ) becomes
    AM = ( 1, 1, 1, 5, 5, 5, 8, 8, 8 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the input vector.

    Input, double A[N], the vector.

    Input, int M, the "stuttering factor".

    Output, double AM[M*N], the stuttering vector.
*/
{
  int i;
  int j;
  int k;

  k = 0;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      am[k] = a[i];
      k = k + 1;
    }
  }
  return;
}
/******************************************************************************/

double *r8vec_stutter_new ( int n, double a[], int m )

/******************************************************************************/
/*

  Purpose:

    R8VEC_STUTTER_NEW makes a "stuttering" copy of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Applying a stuttering factor M of 3, the vector A = ( 1, 5, 8 ) becomes
    AM = ( 1, 1, 1, 5, 5, 5, 8, 8, 8 ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 March 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the size of the input vector.

    Input, double A[N], the vector.

    Input, int M, the "stuttering factor".

    Output, double R8VEC_STUTTER_NEW[M*N], the stuttering vector.
*/
{
  double *am;
  int i;
  int j;
  int k;

  am = ( double * ) malloc ( m * n * sizeof ( double ) );

  k = 0;
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      am[k] = a[i];
      k = k + 1;
    }
  }
  return am;
}
/******************************************************************************/

double r8vec_sum ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SUM returns the sum of an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], the vector.

    Output, double R8VEC_SUM, the sum of the vector.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }

  return value;
}
/******************************************************************************/

void r8vec_swap ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_SWAP swaps the entries of two R8VEC's.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 January 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the arrays.

    Input/output, double A1[N], A2[N], the vectors to swap.
*/
{
  int i;
  double temp;

  for ( i = 0; i < n; i++ )
  {
    temp  = a1[i];
    a1[i] = a2[i];
    a2[i] = temp;
  }

  return;
}
/******************************************************************************/

void r8vec_transpose_print ( int n, double a[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".

  Discussion:

    An R8VEC is a vector of R8's.

  Example:

    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
    TITLE = 'My vector:  '

    My vector:

        1.0    2.1    3.2    4.3    5.4
        6.5    7.6    8.7    9.8   10.9
       11.0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    12 November 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A[N], the vector to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;
  int ihi;
  int ilo;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );

  if ( n <= 0 )
  {
    printf ( "  (Empty)\n" );
    return;
  }

  for ( ilo = 0; ilo < n; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 5, n );
    for ( i = ilo; i < ihi; i++ )
    {
      printf ( "  %12f", a[i] );
    }
    printf ( "\n" );
  }

  return;
}
/******************************************************************************/

void r8vec_undex ( int x_num, double x_val[], int x_unique_num, double tol,
  int undx[], int xdnu[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNDEX returns unique sorted indexes for an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    The goal of this routine is to determine a vector UNDX,
    which points, to the unique elements of X, in sorted order,
    and a vector XDNU, which identifies, for each entry of X, the index of
    the unique sorted element of X.

    This is all done with index vectors, so that the elements of
    X are never moved.

    The first step of the algorithm requires the indexed sorting
    of X, which creates arrays INDX and XDNI.  (If all the entries
    of X are unique, then these arrays are the same as UNDX and XDNU.)

    We then use INDX to examine the entries of X in sorted order,
    noting the unique entries, creating the entries of XDNU and
    UNDX as we go.

    Once this process has been completed, the vector X could be
    replaced by a compressed vector XU, containing the unique entries
    of X in sorted order, using the formula

      XU(*) = X(UNDX(*)).

    We could then, if we wished, reconstruct the entire vector X, or
    any element of it, by index, as follows:

      X(I) = XU(XDNU(I)).

    We could then replace X by the combination of XU and XDNU.

    Later, when we need the I-th entry of X, we can locate it as
    the XDNU(I)-th entry of XU.

    Here is an example of a vector X, the sort and inverse sort
    index vectors, and the unique sort and inverse unique sort vectors
    and the compressed unique sorted vector.

      I     X  Indx  Xdni       XU  Undx  Xdnu
    ----+-----+-----+-----+--------+-----+-----+
      0 | 11.     0     0 |    11.     0     0
      1 | 22.     2     4 |    22.     1     1
      2 | 11.     5     1 |    33.     3     0
      3 | 33.     8     7 |    55.     4     2
      4 | 55.     1     8 |                  3
      5 | 11.     6     2 |                  0
      6 | 22.     7     5 |                  1
      7 | 22.     3     6 |                  1
      8 | 11.     4     3 |                  0

    INDX(2) = 3 means that sorted item(2) is X(3).
    XDNI(2) = 5 means that X(2) is sorted item(5).

    UNDX(3) = 4 means that unique sorted item(3) is at X(4).
    XDNU(8) = 2 means that X(8) is at unique sorted item(2).

    XU(XDNU(I))) = X(I).
    XU(I)        = X(UNDX(I)).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    02 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, int X_NUM, the number of data values.

    Input, double X_VAL[X_NUM], the data values.

    Input, int X_UNIQUE_NUM, the number of unique values in X_VAL.
    This value is only required for languages in which the size of
    UNDX must be known in advance.

    Input, double TOL, a tolerance for equality.

    Output, int UNDX[X_UNIQUE_NUM], the UNDX vector.

    Output, int XDNU[X_NUM], the XDNU vector.
*/
{
  int i;
  int *indx;
  int j;
/*
  Implicitly sort the array.
*/
  indx = r8vec_sort_heap_index_a_new ( x_num, x_val );
/*
  Walk through the implicitly sorted array X.
*/
  i = 0;

  j = 0;
  undx[j] = indx[i];

  xdnu[indx[i]] = j;

  for ( i = 1; i < x_num; i++ )
  {
    if ( tol < r8_abs ( x_val[indx[i]] - x_val[undx[j]] ) )
    {
      j = j + 1;
      undx[j] = indx[i];
    }
    xdnu[indx[i]] = j;
  }
  free ( indx );

  return;
}
/******************************************************************************/

void r8vec_uniform_01 ( int n, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_01 - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return;
}
/******************************************************************************/

double *r8vec_uniform_01_new ( int n, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_01_NEW returns a unit pseudorandom R8VEC.

  Discussion:

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 August 2004

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_01_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = (double*) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

void r8vec_uniform_ab ( int n, double a, double b, int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_AB returns a scaled pseudorandom R8VEC.

  Discussion:

    Each dimension ranges from A to B.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_AB - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return;
}
/******************************************************************************/

double *r8vec_uniform_ab_new ( int n, double a, double b, int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R8VEC.

  Discussion:

    Each dimension ranges from A to B.

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A, B, the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_AB_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = (double*) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
/******************************************************************************/

void r8vec_uniform_abvec ( int n, double a[], double b[], int *seed, double r[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_ABVEC returns a scaled pseudorandom R8VEC.

  Discussion:

    Dimension I ranges from A[I] to B[I].

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], B[N], the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_ABVEC - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = a[i] + ( b[i] - a[i] ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return;
}
/******************************************************************************/

double *r8vec_uniform_abvec_new ( int n, double a[], double b[], int *seed )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIFORM_ABVEC_NEW returns a scaled pseudorandom R8VEC.

  Discussion:

    Dimension I ranges from A[I] to B[I].

    This routine implements the recursion

      seed = 16807 * seed mod ( 2^31 - 1 )
      unif = seed / ( 2^31 - 1 )

    The integer arithmetic never requires more than 32 bits,
    including a sign bit.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    30 January 2005

  Author:

    John Burkardt

  Reference:

    Paul Bratley, Bennett Fox, Linus Schrage,
    A Guide to Simulation,
    Second Edition,
    Springer, 1987,
    ISBN: 0387964673,
    LC: QA76.9.C65.B73.

    Bennett Fox,
    Algorithm 647:
    Implementation and Relative Efficiency of Quasirandom
    Sequence Generators,
    ACM Transactions on Mathematical Software,
    Volume 12, Number 4, December 1986, pages 362-376.

    Pierre L'Ecuyer,
    Random Number Generation,
    in Handbook of Simulation,
    edited by Jerry Banks,
    Wiley, 1998,
    ISBN: 0471134031,
    LC: T57.62.H37.

    Peter Lewis, Allen Goodman, James Miller,
    A Pseudo-Random Number Generator for the System/360,
    IBM Systems Journal,
    Volume 8, Number 2, 1969, pages 136-143.

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double A[N], B[N], the lower and upper limits of the pseudorandom values.

    Input/output, int *SEED, a seed for the random number generator.

    Output, double R8VEC_UNIFORM_ABVEC_NEW[N], the vector of pseudorandom values.
*/
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_UNIFORM_ABVEC_NEW - Fatal error!\n" );
    fprintf ( stderr, "  Input value of SEED = 0.\n" );
    exit ( 1 );
  }

  r = (double*) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = a[i] + ( b[i] - a[i] ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double *r8vec_uniform_unit_new ( int m, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_UNIT_NEW generates a random unit vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_UNIT_NEW[M], a random direction 
//    vector, with unit norm.
//
{
  double *a;
  int i;
  double norm;
//
//  Take M random samples from the normal distribution.
//
  a = r8vec_normal_01_new ( m, seed );
//
//  Compute the norm.
//
  norm = 0.0;
  for ( i = 0; i < m; i++ )
  {
    norm = norm + a[i] * a[i];
  }
  norm = sqrt ( norm );
//
//  Normalize.
//
  for ( i = 0; i < m; i++ )
  {
    a[i] = a[i] / norm;
  }

  return a;
}
/******************************************************************************/

int r8vec_unique_count ( int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIQUE_COUNT counts the unique elements in an unsorted R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

    Because the array is unsorted, this algorithm is O(N^2).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    29 April 2004

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input, double A[N], the array to examine, which does NOT have to
    be sorted.

    Input, double TOL, a tolerance for checking equality.

    Output, int R8VEC_UNIQUE_COUNT, the number of unique elements of A.
*/
{
  int i;
  int j;
  int unique_num;

  unique_num = 0;

  for ( i = 0; i < n; i++ )
  {
    unique_num = unique_num + 1;

    for ( j = 0; j < i; j++ )
    {
      if ( r8_abs ( a[i] - a[j] ) <= tol )
      {
        unique_num = unique_num - 1;
        break;
      }
    }
  }
  return unique_num;
}
/******************************************************************************/

int *r8vec_unique_index ( int n, double a[], double tol )

/******************************************************************************/
/*
  Purpose:

    R8VEC_UNIQUE_INDEX indexes the unique occurrence of values in an R8VEC.

  Discussion:

    For element A(I) of the vector, UNIQUE_INDEX(I) is the uniqueness index
    of A(I).  That is, if A_UNIQUE contains the unique elements of A,
    gathered in order, then

      A_UNIQUE ( UNIQUE_INDEX(I) ) = A(I)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of elements of A.

    Input, double A[N], the unsorted array to examine.

    Input, double TOL, a tolerance for equality.

    Output, int R8VEC_UNIQUE_INDEX[N], the unique index.
*/
{
  int i;
  int j;
  int *unique_index;
  int unique_num;

  unique_index = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    unique_index[i] = -1;
  }
  unique_num = 0;

  for ( i = 0; i < n; i++ )
  {
    if ( unique_index[i] == -1 )
    {
      unique_index[i] = unique_num;
      for ( j = i + 1; j < n; j++ )
      {
        if ( r8_abs ( a[i] - a[j] ) <= tol )
        {
          unique_index[j] = unique_num;
        }
      }
      unique_num = unique_num + 1;
    }
  }
  return unique_index;
}
/******************************************************************************/

double r8vec_variance ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_VARIANCE returns the variance of an R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 May 2006

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Input, double X[N], the vector whose variance is desired.

    Output, double R8VEC_VARIANCE, the variance of the vector entries.
*/
{
  int i;
  double mean;
  double variance;

  mean = r8vec_mean ( n, x );

  variance = 0.0;
  for ( i = 0; i < n; i++ )
  {
    variance = variance + ( x[i] - mean ) * ( x[i] - mean );
  }

  if ( 1 < n )
  {
    variance = variance / ( double ) ( n - 1 );
  }
  else
  {
    variance = 0.0;
  }

  return variance;
}
/******************************************************************************/

double *r8vec_vector_triple_product ( double v1[3], double v2[3], double v3[3] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_VECTOR_TRIPLE_PRODUCT computes the vector triple product.

  Discussion:

    VTRIPLE = V1 x (V2 x V3)

    VTRIPLE is a vector perpendicular to V1, lying in the plane
    spanned by V2 and V3.  The norm of VTRIPLE is the product
    of the norms of V1, V2 and V3.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 October 2010

  Author:

    John Burkardt

  Parameters:

    Input, double V1[3], V2[3], V3[3], the three vectors.

    Output, double R8VEC_VECTOR_TRIPLE_PRODUCT[3], the vector triple product.
*/
{
  double *v123;
  double *v23;

  v23 = r8vec_cross_product_3d ( v2, v3 );

  v123 = r8vec_cross_product_3d ( v1, v23 );

  free ( v23 );

  return v123;
}
/******************************************************************************/

void r8vec_write ( int n, double r[], char *output_file )

/******************************************************************************/
/*
  Purpose:

    R8VEC_WRITE writes an R8VEC to a file.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 August 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the order of the matrix.

    Input, double R[N], the vector to be written.

    Input, char *OUTPUT_FILE, the name of the file to which
    the information is to be written.
*/
{
  int i;
  FILE *output;

  output = fopen ( output_file, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8VEC_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    return;
  }

  for ( i = 0; i < n; i++ )
  {
    fprintf ( output, "  %16f\n", r[i] );
  }

  fclose ( output );

  return;
}
/******************************************************************************/

void r8vec_zero ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ZERO zeroes an R8VEC.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 August 2008

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Output, double A[N], a vector of zeroes.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}
/******************************************************************************/

double *r8vec_zero_new ( int n )

/******************************************************************************/
/*
  Purpose:

    R8VEC_ZERO_NEW creates and zeroes an R8VEC.

  Discussion:

    An R8VEC is a vector of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    25 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vector.

    Output, double R8VEC_ZERO_NEW[N], a vector of zeroes.
*/
{
  double *a;
  int i;

  a = ( double * ) malloc ( n * sizeof ( double ) );

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return a;
}
/******************************************************************************/

int r8vec2_compare ( int n, double a1[], double a2[], int i, int j )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_COMPARE compares two elements of an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of data items.

    Input, double A1[N], A2[N], contain the two components of each item.

    Input, int I, J, the items to be compared.  These values will be
    1-based indices for the arrays A1 and A2.

    Output, int R8VEC2_COMPARE, the results of the comparison:
    -1, item I < item J,
     0, item I = item J,
    +1, item J < item I.
*/
{
  int isgn;

  isgn = 0;

  if ( a1[i-1] < a1[j-1] )
  {
    isgn = -1;
  }
  else if ( a1[i-1] == a1[j-1] )
  {
    if ( a2[i-1] < a2[j-1] )
    {
      isgn = -1;
    }
    else if ( a2[i-1] < a2[j-1] )
    {
      isgn = 0;
    }
    else if ( a2[j-1] < a2[i-1] )
    {
      isgn = +1;
    }
  }
  else if ( a1[j-1] < a1[i-1] )
  {
    isgn = +1;
  }

  return isgn;
}
/******************************************************************************/

void r8vec2_print ( int n, double a1[], double a2[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_PRINT prints an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A1[N], double A2[N], the vectors to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );
  for ( i = 0; i < n; i++ )
  {
    fprintf ( stdout, "  %4d: %14f  %14f\n", i, a1[i], a2[i] );
  }

  return;
}
/******************************************************************************/

void r8vec2_print_some ( int n, double x1[], double x2[], int max_print,
  char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_PRINT_SOME prints "some" of an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

    The user specifies MAX_PRINT, the maximum number of lines to print.

    If N, the size of the vectors, is no more than MAX_PRINT, then
    the entire vectors are printed, one entry of each per line.

    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
    followed by a line of periods suggesting an omission,
    and the last entry.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries of the vectors.

    Input, double X1[N], X2[N], the vector to be printed.

    Input, int MAX_PRINT, the maximum number of lines to print.

    Input, char *TITLE, a title.
*/
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  fprintf ( stdout, "\n" );
  fprintf ( stdout, "%s\n", title );
  fprintf ( stdout, "\n" );

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      fprintf ( stdout, "  %4d: %14f  %14f\n", i, x1[i], x2[i] );
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print-2; i++ )
    {
      fprintf ( stdout, "  %4d: %14f  %14f\n", i, x1[i], x2[i] );
    }
    fprintf ( stdout, "......  ..............  ..............\n" );
    i = n - 1;
    fprintf ( stdout, "  %4d: %14f  %14f\n", i, x1[i], x2[i] );
  }
  else
  {
    for ( i = 0; i < max_print - 1; i++ )
    {
      fprintf ( stdout, "  %4d: %14f  %14f\n", i, x1[i], x2[i] );
    }
    i = max_print - 1;
    fprintf ( stdout, "  %4d: %14f  %14f  ...more entries...\n", i, x1[i], x2[i] );
  }

  return;
}
/******************************************************************************/

void r8vec2_sort_a ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SORT_A ascending sorts an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

    Each item to be sorted is a pair of reals (X,Y), with the X
    and Y values stored in separate vectors A1 and A2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of items of data.

    Input/output, double A1[N], A2[N], the data to be sorted.
*/
{
  int i;
  int indx;
  int isgn;
  int j;
  double temp;
/*
  Initialize.
*/
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
/*
  Call the external heap sorter.
*/
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
/*
  Interchange the I and J objects.
*/
    if ( 0 < indx )
    {
      temp    = a1[i-1];
      a1[i-1] = a1[j-1];
      a1[j-1] = temp;

      temp    = a2[i-1];
      a2[i-1] = a2[j-1];
      a2[j-1] = temp;
    }
/*
  Compare the I and J objects.
*/
    else if ( indx < 0 )
    {
      isgn = r8vec2_compare ( n, a1, a2, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

void r8vec2_sort_d ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SORT_D descending sorts an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

    Each item to be sorted is a pair of reals (X,Y), with the X
    and Y values stored in separate vectors A1 and A2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of items of data.

    Input/output, double A1[N], A2[N], the data to be sorted.
*/
{
  int i;
  int indx;
  int isgn;
  int j;
  double temp;
/*
  Initialize.
*/
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
/*
  Call the external heap sorter.
*/
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
/*
  Interchange the I and J objects.
*/
    if ( 0 < indx )
    {
      temp    = a1[i-1];
      a1[i-1] = a1[j-1];
      a1[j-1] = temp;

      temp    = a2[i-1];
      a2[i-1] = a2[j-1];
      a2[j-1] = temp;
    }
/*
  Compare the I and J objects.
*/
    else if ( indx < 0 )
    {
      isgn = - r8vec2_compare ( n, a1, a2, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

int *r8vec2_sort_heap_index_a ( int n, double x[], double y[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

    The sorting is not actually carried out.  Rather an index array is
    created which defines the sorting.  This array may be used to sort
    or index the array, or to sort or index related arrays keyed on the
    original array.

    ( X(I), Y(I) ) < ( X(J), Y(J) ) if:

    * X(I) < X(J), or

    * X(I) = X(J), and Y(I) < Y(J).

    Once the index array is computed, the sorting can be carried out
    implicitly:

      ( x(indx(*)), y(indx(*) )

    or explicitly, by the calls

      r8vec_permute ( n, indx, 0, x )
      r8vec_permute ( n, indx, 0, y )

    after which ( x(*), y(*) ), is sorted.

    Note that the index vector is 0-based.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double X[N], Y[N], pairs of X, Y coordinates of points.

    Output, int INDX[N], the sort index.  The
    I-th element of the sorted array has coordinates
    ( X(INDX(I)), Y(INDX(I) ).
*/
{
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;
  double xval;
  double yval;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = ( int * ) malloc ( n * sizeof ( int ) );

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      xval = x[indxt];
      yval = y[indxt];
    }
    else
    {
      indxt = indx[ir-1];
      xval = x[indxt];
      yval = y[indxt];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( x[indx[j-1]] < x[indx[j]] ||
          ( x[indx[j-1]] == x[indx[j]] && y[indx[j-1]] < y[indx[j]] ) )
        {
          j = j + 1;
        }
      }

      if ( xval < x[indx[j-1]] ||
         ( xval == x[indx[j-1]] && yval < y[indx[j-1]] ) )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
  return indx;
}
/******************************************************************************/

void r8vec2_sorted_unique ( int n, double a1[], double a2[], int *unique_num )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SORTED_UNIQUE keeps the unique elements in an R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

    Item I is stored as the pair A1(I), A2(I).

    The items must have been sorted, or at least it must be the
    case that equal items are stored in adjacent vector locations.

    If the items were not sorted, then this routine will only
    replace a string of equal values by a single representative.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of items.

    Input/output, double A1[N], A2[N].
    On input, the array of N items.
    On output, an array of UNIQUE_NUM unique items.

    Output, int *UNIQUE_NUM, the number of unique items.
*/
{
  int itest;

  *unique_num = 0;

  if ( n <= 0 )
  {
    return;
  }

  *unique_num = 1;

  for ( itest = 1; itest < n; itest++ )
  {
    if ( a1[itest] != a1[*unique_num-1] ||
         a2[itest] != a2[*unique_num-1] )
    {
      a1[*unique_num] = a1[itest];
      a2[*unique_num] = a2[itest];
      *unique_num = *unique_num + 1;
    }
  }

  return;
}
/******************************************************************************/

void r8vec2_sorted_unique_index ( int n, double a1[], double a2[],
  int *unique_num, int indx[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SORTED_UNIQUE_INDEX indexes unique elements in a sorted R8VEC2.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

    Item I is stored as the pair A1(I), A2(I).

    The items must have been sorted, or at least it should be the
    case that equal items are stored in adjacent vector locations.

    If the items are not sorted, then this routine will only
    replace a string of equal values by a single representative.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of items.

    Input/output, double A1[N], A2[N].
    On input, the array of N items.
    On output, an array of unique items.

    Output, int *UNIQUE_NUM, the number of unique items.

    Output, int INDX[N], contains in entries 1 through UNIQUE_NUM an index
    array of the unique items.  To build new arrays with no repeated elements:
      B1(*) = A1(INDX(*))
*/
{
  int itest;

  if ( n <= 0 )
  {
    *unique_num = 0;
    return;
  }
  i4vec_zero ( n, indx );

  *unique_num = 1;
  indx[0] = 1;

  for ( itest = 2; itest <= n; itest++ )
  {
    if ( a1[itest-2] != a1[itest-1] || a2[itest-2] != a2[itest-1] )
    {
      *unique_num = *unique_num + 1;
      indx[*unique_num-1] = itest;
    }
  }

  return;
}
/******************************************************************************/

int r8vec2_sum_max_index ( int n, double a[], double b[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC2_SUM_MAX_INDEX returns the index of the maximum sum of two R8VEC's.

  Discussion:

    An R8VEC2 is a dataset consisting of N pairs of real values, stored
    as two separate vectors A1 and A2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    26 March 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the array.

    Input, double A[N], B[N], two arrays whose sum
    is to be examined.

    Output, int R8VEC2_SUM_MAX_INDEX, the index of the largest entry in A+B.
*/
{
  int i;
  double sum_max;
  int sum_max_index;

  if ( n <= 0 )
  {
    sum_max_index = -1;
  }
  else
  {
    sum_max_index = 1;
    sum_max = a[0] + b[0];

    for ( i = 2; i <= n; i++ )
    {
      if ( sum_max < a[i-1] + b[i-1] )
      {
        sum_max = a[i-1] + b[i-1];
        sum_max_index = i;
      }
    }
  }
  return sum_max_index;
}
/******************************************************************************/

void r8vec3_print ( int n, double a1[], double a2[], double a3[], char *title )

/******************************************************************************/
/*
  Purpose:

    R8VEC3_PRINT prints a triple of real vectors.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    27 August 2011

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of components of the vector.

    Input, double A1[N], double A2[N], double A3[N], the vectors
    to be printed.

    Input, char *TITLE, a title.
*/
{
  int i;

  printf ( "\n" );
  printf ( "%s\n", title );
  printf ( "\n" );
  for ( i = 0; i < n; i++ )
  {
    printf ( "  %4d  %10g  %10g  %10g\n", i, a1[i], a2[i], a3[i] );
  }

  return;
}
/******************************************************************************/

double *roots_to_r8poly ( int n, double x[] )

/******************************************************************************/
/*
  Purpose:

    ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    19 July 2010

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of roots specified.

    Input, double X[N], the roots.

    Output, double ROOTS_TO_R8POLY[N+1], the coefficients of the polynomial.
*/
{
  double *c;
  int i;
  int j;

  c = r8vec_zero_new ( n + 1 );
/*
  Initialize C to (0, 0, ..., 0, 1).
  Essentially, we are setting up a divided difference table.
*/
  c[n] = 1.0;
/*
  Convert to standard polynomial form by shifting the abscissas
  of the divided difference table to 0.
*/
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= n+1-j; i++ )
    {
      c[n-i] = c[n-i] - x[n+1-i-j] * c[n-i+1];
    }
  }
  return c;
}
/******************************************************************************/

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

/******************************************************************************/
/*
  Purpose:

    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.

  Discussion:

    The actual list is not passed to the routine.  Hence it may
    consist of integers, reals, numbers, names, etc.  The user,
    after each return from the routine, will be asked to compare or
    interchange two items.

    The current version of this code mimics the FORTRAN version,
    so the values of I and J, in particular, are FORTRAN indices.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    05 February 2004

  Author:

    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    C version by John Burkardt.

  Reference:

    Albert Nijenhuis, Herbert Wilf,
    Combinatorial Algorithms,
    Academic Press, 1978, second edition,
    ISBN 0-12-519260-6.

  Parameters:

    Input, int N, the length of the input list.

    Input/output, int *INDX.
    The user must set INDX to 0 before the first call.
    On return,
      if INDX is greater than 0, the user must interchange
      items I and J and recall the routine.
      If INDX is less than 0, the user is to compare items I
      and J and return in ISGN a negative value if I is to
      precede J, and a positive value otherwise.
      If INDX is 0, the sorting is done.

    Output, int *I, *J.  On return with INDX positive,
    elements I and J of the user's list should be
    interchanged.  On return with INDX negative, elements I
    and J are to be compared by the user.

    Input, int ISGN. On return with INDX negative, the
    user should compare elements I and J of the list.  If
    item I is to precede item J, set ISGN negative,
    otherwise set ISGN positive.
*/
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
/*
  INDX = 0: This is the first call.
*/
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
/*
  INDX < 0: The user is returning the results of a comparison.
*/
  else if ( *indx < 0 )
  {
    if ( *indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      if ( n1 == 1 )
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
/*
  0 < INDX: the user was asked to make an interchange.
*/
  else if ( *indx == 1 )
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 )
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 )
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 )
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}

} // extern "C"
