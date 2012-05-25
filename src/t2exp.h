/* --------------------------------------------------------------
    Header file for two tables-driven exponent function.

    Home page: www.imach.uran.ru/exp

    Copyright 2001-2002 by Dr. Raul N.Shakirov, IMach of RAS(UB),
    Phillip S. Pang, Ph.D. Biochemistry and Molecular Biophysics.
    Columbia University. NYC. All Rights Reserved.

    Permission has been granted to copy, distribute and modify
    software in any context without fee, including a commercial
    application, provided that the aforesaid copyright statement
    is present here as well as exhaustive description of changes.

    THE SOFTWARE IS DISTRIBUTED "AS IS". NO WARRANTY OF ANY KIND
    IS EXPRESSED OR IMPLIED. YOU USE AT YOUR OWN RISK. THE AUTHOR
    WILL NOT BE LIABLE FOR DATA LOSS, DAMAGES, LOSS OF PROFITS OR
    ANY OTHER KIND OF LOSS WHILE USING OR MISUSING THIS SOFTWARE.
-------------------------------------------------------------- */

#ifndef T2EXP_H
#define T2EXP_H

#ifdef  __cplusplus
extern "C" {
#endif/*__cplusplus*/

/* --------------------------------------------------------------
    Name:       t2exp

    Purpose:    Fast two table-driven exponent algorithm,
                effective for arg <= 0.

    Usage:      t2exp (arg)

    Domain:     Same as for standard exp() function
                (approximately -709 <= arg <= 709).

    Result:     Approximate exp of arg; if arg is outside the
                exp() domain, results are same as for standard
                exp() function - that is either 0 or INF.
-------------------------------------------------------------- */

extern double t2exp  (double arg);


/* --------------------------------------------------------------
    Name:       t2expini

    Purpose:    Build tables for t2exp().

    Usage:      t2expini()

    Note:       Used for development purposes only!
-------------------------------------------------------------- */

extern void t2expini (void);

/* --------------------------------------------------------------
    Name:       t2expinl

    Purpose:    Print tables for t2exp() in format of file t2exp.inl.

    Usage:      t2expinl()

    Note:       Used for development purposes only!
-------------------------------------------------------------- */

extern void t2expinl (void);

#ifdef  __cplusplus
}
#endif/*__cplusplus*/

#endif/*T2EXP_H*/
