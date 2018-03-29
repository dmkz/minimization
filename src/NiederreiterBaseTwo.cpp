/*
 * Авторы: Бураханова А., Золкин А., Казарян М., Хохлов С., Малоенко С.
 * Генератор сетоr Нидеррайтера, использующий поле F2
 * Файл: NiederreiterBaseTwo.cpp
*/

#include "NiederreiterBaseTwo.hpp"

int NiederreiterBaseTwo::Init() {
    return Init(DIM_MAX, 0);
}

int NiederreiterBaseTwo::Init(uint32_t dimension) {
    return Init(dimension, 0);
}

int NiederreiterBaseTwo::Init(uint32_t dimension_, int64_t seed_) {
    // Инициализация неприводимых многочленов:
    irred = {
        {0, 1, 0, 0, 0, 0, 0, 0, 0}, //  0   x
        {1, 1, 0, 0, 0, 0, 0, 0, 0}, //  1   1+x
        {1, 1, 1, 0, 0, 0, 0, 0, 0}, //  2   1+x+x^2
        {1, 1, 0, 1, 0, 0, 0, 0, 0}, //  3   1+x+x^3
        {1, 0, 1, 1, 0, 0, 0, 0, 0}, //  4   1+x^2+x^3
        {1, 1, 0, 0, 1, 0, 0, 0, 0}, //  5   1+x+x^4
        {1, 0, 0, 1, 1, 0, 0, 0, 0}, //  6   1+x^3+x^4
        {1, 0, 0, 1, 0, 1, 0, 0, 0}, //  7   1+x^3+x^5 
        {1, 1, 1, 1, 0, 1, 0, 0, 0}, //  8   1+x+x^2+x^3+x^5
        {1, 1, 1, 0, 1, 1, 0, 0, 0}, //  9   1+x+x^2+x^4+x^5
        {1, 1, 0, 1, 1, 1, 0, 0, 0}, // 10   1+x+x^3+x^4+x^5
        {1, 0, 1, 1, 1, 1, 0, 0, 0}, // 11   1+x^2+x^3+x^4+x^5
        {1, 1, 0, 0, 0, 0, 1, 0, 0}, // 12   1+x+x^6
        {1, 0, 0, 1, 0, 0, 1, 0, 0}, // 13   1+x^3+x^6
        {1, 1, 1, 0, 1, 0, 1, 0, 0}, // 14   1+x+x^2+x^4+x^6
        {1, 1, 0, 1, 1, 0, 1, 0, 0}, // 15   1+x+x^3+x^4+x^6
        {1, 0, 0, 0, 0, 1, 1, 0, 0}, // 16   1+x^5+x^6
        {1, 1, 1, 0, 0, 1, 1, 0, 0}, // 17   1+x+x^2+x^5+x^6
        {1, 0, 1, 1, 0, 1, 1, 0, 0}, // 18   1+x^2+x^3+x^5+x^6
        {1, 1, 0, 0, 1, 1, 1, 0, 0}, // 19   1+x+x^4+x^5+x^6 */
        {1, 0, 1, 0, 1, 1, 1, 0, 0}, // 20   1+x^2+x^4+x^5+x^6
        {1, 1, 0, 0, 0, 0, 0, 1, 0}, // 21   1+x+x^7
        {1, 0, 0, 1, 0, 0, 0, 1, 0}, // 22   1+x^3+x^7
        {1, 1, 1, 1, 0, 0, 0, 1, 0}, // 23   1+x+x^2+x^3+x^7
        {1, 0, 0, 0, 1, 0, 0, 1, 0}, // 24   1+x^4+x^7
        {1, 0, 1, 1, 1, 0, 0, 1, 0}, // 25   1+x^2+x^3+x^4+x^7
        {1, 1, 1, 0, 0, 1, 0, 1, 0}, // 26   1+x+x^2+x^5+x^7
        {1, 1, 0, 1, 0, 1, 0, 1, 0}, // 27   1+x+x^3+x^5+x^7
        {1, 0, 0, 1, 1, 1, 0, 1, 0}, // 28   1+x^3+x^4+x^5+x^7
        {1, 1, 1, 1, 1, 1, 0, 1, 0}, // 29   1+x+x^2+x^3+x^4+x^5+x^7
        {1, 0, 0, 0, 0, 0, 1, 1, 0}, // 30   1+x^6+x^7
        {1, 1, 0, 1, 0, 0, 1, 1, 0}, // 31   1+x+x^3+x^6+x^7
        {1, 1, 0, 0, 1, 0, 1, 1, 0}, // 32   1+x+x^4+x^6+x^7
        {1, 0, 1, 0, 1, 0, 1, 1, 0}, // 33   1+x^2+x^4+x^6+x^7
        {1, 0, 1, 0, 0, 1, 1, 1, 0}, // 34   1+x^2+x^5+x^6+x^7
        {1, 0, 0, 0, 1, 1, 1, 1, 0}, // 35   1+x^4+x^5+x^6+x^7
        {1, 1, 1, 0, 1, 1, 1, 1, 0}, // 36   1+x+x^2+x^4+x^5+x^6+x^7
        {1, 0, 1, 1, 1, 1, 1, 1, 0}, // 37   1+x^2+x^3+x^4+x^5+x^6+x^7
        {1, 1, 0, 1, 1, 0, 0, 0, 1}, // 38   1+x+x^3+x^4+x^8
        {1, 0, 1, 1, 1, 0, 0, 0, 1}, // 39   1+x^2+x^3+x^4+x^8
        {1, 1, 0, 1, 0, 1, 0, 0, 1}, // 40   1+x+x^3+x^5+x^8
        {1, 0, 1, 1, 0, 1, 0, 0, 1}, // 41   1+x^2+x^3+x^5+x^8
        {1, 0, 0, 1, 1, 1, 0, 0, 1}, // 42   1+x^3+x^4+x^5+x^8
        {1, 1, 1, 1, 1, 1, 0, 0, 1}, // 43   1+x+x^2+x^3+x^4+x^5+x^8
        {1, 0, 1, 1, 0, 0, 1, 0, 1}, // 44   1+x^2+x^3+x^6+x^8
        {1, 1, 1, 1, 1, 0, 1, 0, 1}, // 45   1+x+x^2+x^3+x^4+x^6+x^8
        {1, 1, 0, 0, 0, 1, 1, 0, 1}, // 46   1+x+x^5+x^6+x^8
        {1, 0, 1, 0, 0, 1, 1, 0, 1}, // 47   1+x^2+x^5+x^6+x^8
        {1, 0, 0, 1, 0, 1, 1, 0, 1}, // 48   1+x^3+x^5+x^6+x^8
        {1, 0, 0, 0, 1, 1, 1, 0, 1}, // 49   1+x^4+x^5+x^6+x^8
        {1, 1, 1, 0, 1, 1, 1, 0, 1}, // 50   1+x+x^2+x^4+x^5+x^6+x^8
        {1, 1, 0, 1, 1, 1, 1, 0, 1}, // 51   1+x+x^3+x^4+x^5+x^6+x^8
        {1, 1, 1, 0, 0, 0, 0, 1, 1}, // 52   1+x+x^2+x^7+x^8
        {1, 1, 0, 1, 0, 0, 0, 1, 1}, // 53   1+x+x^3+x^7+x^8
        {1, 0, 1, 1, 0, 0, 0, 1, 1}, // 54   1+x^2+x^3+x^7+x^8
        {1, 1, 1, 1, 1, 0, 0, 1, 1}, // 55   1+x+x^2+x^3+x^4+x^7+x^8
        {1, 0, 1, 1, 1, 1, 0, 1, 1}, // 56   1+x^2+x^3+x^4+x^5+x^7+x^8
        {1, 1, 0, 0, 0, 0, 1, 1, 1}, // 57   1+x+x^6+x^7+x^8
        {1, 1, 1, 0, 1, 0, 1, 1, 1}, // 58   1+x+x^2+x^4+x^6+x^7+x^8
        {1, 0, 1, 1, 1, 0, 1, 1, 1}, // 59   1+x^2+x^3+x^4+x^6+x^7+x^8
        {1, 1, 1, 0, 0, 1, 1, 1, 1}, // 60   1+x+x^2+x^5+x^6+x^7+x^8
        {1, 1, 0, 0, 1, 1, 1, 1, 1}, // 61   1+x+x^4+x^5+x^6+x^7+x^8
        {1, 0, 1, 0, 1, 1, 1, 1, 1}  // 62   1+x^2+x^4+x^5+x^6+x^7+x^8
        // {1, 0, 0, 1, 1, 1, 1, 1, 1, 0}, // 63   1+x^3+x^4+x^5+x^6+x^7+x^8 
        // {1, 0, 0, 0, 1, 0, 0, 0, 0, 1} // 65   1+x^4+x^9
        
    };
    
    assert((int)irred.size() == DIM_MAX);
    MAXE = 0;
    
    // Подсчет степеней многочленов:
    for (int i = 0; i < (int)irred.size(); ++i) {
        int deg = (int)irred[i].size()-1;;
        while (deg >= 0 && irred[i][deg] == 0) {
            --deg;
        }
        assert(deg > 0 && irred[i][deg] == 1);
        irred_deg.push_back(deg);
        MAXE = std::max(MAXE, deg);
    }

    dimension = dimension_;
    seed = seed_;
    setfld2();
    dim_save = 0;
    seed_save = 0;
    calcc2();
    
    return 0;
}

//****************************************************************************80

void NiederreiterBaseTwo::setfld2()

//****************************************************************************80
//
//  Purpose:
//
//    SETFLD2 sets up arithmetic tables for the finite field of order 2.
//
//  Discussion:
//
//    SETFLD2 sets up addition, multiplication, and subtraction tables 
//    for the finite field of order QIN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 March 2003
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int64_t ADD[2][2], MUL[2][2], SUB[2][2], the addition, multiplication, 
//    and subtraction tables, mod 2.
//
{
  int64_t i;
  int64_t j;
  int64_t p = 2;
  int64_t q = 2;
//
  for ( i = 0; i < q; i++ )
  {
    for ( j = 0; j < q; j++ )
    {
      add[i][j] = ( i + j ) % p;
      mul[i][j] = ( i * j ) % p;
    }
  }
//
//  Use the addition table to set the subtraction table.
//
  for ( i = 0; i < q; i++ )
  {
    for ( j = 0; j < q; j++ )
    {
      sub[ add[i][j] ] [i] = j;
    }
  }

  return;
}
//****************************************************************************80

//****************************************************************************80

void NiederreiterBaseTwo::plymul2 (
    int64_t pa_deg, int64_t pa[MAXDEG+1], 
    int64_t pb_deg, int64_t pb[MAXDEG+1], 
    int64_t *pc_deg, int64_t pc[MAXDEG+1] ) const

//****************************************************************************80
//
//  Purpose:
//
//    PLYMUL2 multiplies two polynomials in the field of order 2
//
//  Discussion:
//
//    Polynomials stored as arrays have the coefficient of degree N in 
//    POLY(N), and the degree of the polynomial in POLY(-1).  
//
//    A polynomial which is identically 0 is given degree -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 March 2003
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int64_t ADD[2][2], MUL[2][2], 
//    the addition and multiplication tables, mod 2.
//
//    Input, int64_t PA_DEG, the degree of PA.
//
//    Input, int64_t PA[MAXDEG+1], the first polynomial factor.
//
//    Input, int64_t PB_DEG, the degree of PB.
//
//    Input, int64_t PB[MAXDEG+1], the second polynomial factor.
//
//    Output, int64_t *PC_DEG, the degree of the product.
//
//    Output, int64_t PC[MAXDEG+1], the product polynomial.
//
{
  int64_t i;
  int64_t j;
  int64_t jhi;
  int64_t jlo;
  int64_t pt[MAXDEG+1];
  int64_t term;

  if ( pa_deg == -1 || pb_deg == -1 )
  {
    *pc_deg = -1;
  }
  else
  {
    *pc_deg = pa_deg + pb_deg;
  }

  if ( MAXDEG < *pc_deg )
  {
    std::cout << "\n";
    std::cout << "PLYMUL2 - Fatal error!\n";
    std::cout << "  Degree of the product exceeds MAXDEG.\n";
    std::exit ( 1 );
  }

  for ( i = 0; i <= *pc_deg; i++ )
  {

    jlo = i - pa_deg;
    if ( jlo < 0 ) 
    {
      jlo = 0;
    }

    jhi = pb_deg;
    if ( i < jhi ) 
    {
      jhi = i;
    }

    term = 0;

    for ( j = jlo; j <= jhi; j++ ) 
    {
      term = add [ term ] [ mul [ pa[i-j] ] [ pb[j] ] ];
    }
    pt[i] = term;
  }

  for ( i = 0; i <= *pc_deg; i++ )
  {
    pc[i] = pt[i];
  }

  for ( i = *pc_deg + 1; i <= MAXDEG; i++ )
  {
    pc[i] = 0;
  }

  return;
}
//****************************************************************************80

//****************************************************************************80
void NiederreiterBaseTwo::calcc2 ()

//****************************************************************************80
//
//  Purpose:
//
//    CALCC2 computes values of the constants C(I,J,R).
//
//  Discussion:
//
//    This program calculates the values of the constants C(I,J,R).
//
//    As far as possible, Niederreiter's notation is used.
//
//    For each value of I, we first calculate all the corresponding
//    values of C.  These are held in the array CI.  All these
//    values are either 0 or 1.  
//
//    Next we pack the values into the
//    array CJ, in such a way that CJ(I,R) holds the values of C
//    for the indicated values of I and R and for every value of
//    J from 1 to NBITS.  The most significant bit of CJ(I,R)
//    (not counting the sign bit) is C(I,1,R) and the least
//    significant bit is C(I,NBITS,R).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 March 2003
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    R Lidl, Harald Niederreiter, 
//    Finite Fields,
//    Cambridge University Press, 1984, page 553.
//
//    Harald Niederreiter,
//    Low-discrepancy and low-dispersion sequences,
//    Journal of Number Theory,
//    Volume 30, 1988, pages 51-70.
//
//  Parameters:
//
//    Input, int64_t DIM_NUM, the dimension of the sequence to be generated.
//
//    Output, int64_t CJ[DIM_MAX][NBITS], the packed values of 
//    Niederreiter's C(I,J,R)
//
//  Local Parameters:
//
//    Local, int64_t MAXE; we need DIM_MAX irreducible polynomials over Z2.
//    MAXE is the highest degree among these.
//
//    Local, int64_t MAXV, the maximum possible index used in V.
//
{
    
    int64_t b[MAXDEG+1];
    int64_t b_deg;
    int64_t ci[NBITS][NBITS];
    // int64_t count; не используется
    int64_t e;
    
  
    int64_t maxv = NBITS + MAXE;
    // int64_t nextq[DIM_MAX]; // не используется
    // int64_t p; // не используется
    int64_t px[MAXDEG+1];
    int64_t px_deg;
    // int64_t q; // не используется
    int64_t term;
    int64_t u;
    int64_t v[NBITS+MAXE+1];

    for (int64_t i = 0; i < dimension; i++) {
        //
        //  For each dimension, we need to calculate powers of an
        //  appropriate irreducible polynomial:  see Niederreiter
        //  page 65, just below equation (19).
        //
        //  Copy the appropriate irreducible polynomial into PX,
        //  and its degree into E.  Set polynomial B = PX ** 0 = 1.
        //  M is the degree of B.  Subsequently B will hold higher
        //  powers of PX.
        //
        e = irred_deg[i];

        px_deg = irred_deg[i];

        for (int64_t j = 0; j <= px_deg; j++) {
            px[j] = irred[i][j];
        }

        b_deg = 0;
        b[0] = 1;
        //
        //  Niederreiter (page 56, after equation (7), defines two
        //  variables Q and U.  We do not need Q explicitly, but we do need U.
        //
        u = 0;

        for (int64_t j = 0; j < NBITS; j++) {
            //
            //  If U = 0, we need to set B to the next power of PX
            //  and recalculate V.  This is done by subroutine CALCV.
            //
            if (u == 0) {
                calcv2(maxv, px_deg, px, &b_deg, b, v);
            }
            //
            //  Now C is obtained from V.  Niederreiter obtains A from V (page 65, 
            //  near the bottom), and then gets C from A (page 56, equation (7)).  
            //  However this can be done in one step.  Here CI(J,R) corresponds to
            //  Niederreiter's C(I,J,R).
            //
            for (int64_t r = 0; r < NBITS; r++) {
                ci[j][r] = v[r+u];
            }
            //
            //  Increment U.  
            //
            //  If U = E, then U = 0 and in Niederreiter's
            //  paper Q = Q + 1.  Here, however, Q is not used explicitly.
            //
            u = u + 1;
            if (u == e) {
                u = 0;
            }
        }
        //
        //  The array CI now holds the values of C(I,J,R) for this value
        //  of I.  We pack them into array CJ so that CJ(I,R) holds all
        //  the values of C(I,J,R) for J from 1 to NBITS.
        //
        for (int64_t r = 0; r < NBITS; r++) {
            term = 0;
            for (int64_t j = 0; j < NBITS; j++) {
                term = 2 * term + ci[j][r];
            }
            cj[i][r] = term;
        }
    }

    return;
}
//****************************************************************************80

void NiederreiterBaseTwo::calcv2 ( int64_t maxv, int64_t px_deg, int64_t px[MAXDEG+1], 
    int64_t *b_deg, int64_t b[MAXDEG+1], int64_t v[] ) const

//****************************************************************************80
//
//  Purpose:
//
//    CALCV2 calculates the value of the constants V(J,R).
//
//  Discussion:
//
//    This program calculates the values of the constants V(J,R) as
//    described in the reference (BFN) section 3.3.  It is called from CALCC2.  
//
//    Polynomials stored as arrays have the coefficient of degree N 
//    in POLY(N).  
//
//    A polynomial which is identically 0 is given degree -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 March 2003
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Harald Niederreiter,
//    Algorithm 738: 
//    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 20, Number 4, pages 494-495, 1994.
//
//  Parameters:
//
//    Input, int64_t MAXV, the dimension of the array V.
//
//    Input, int64_t PX_DEG, the degree of PX.
//
//    Input, int64_t PX[MAXDEG+1], the appropriate irreducible polynomial 
//    for the dimension currently being considered.  
//
//    Input, int64_t ADD[2][2], MUL[2][2], SUB[2][2], the addition, multiplication, 
//    and subtraction tables, mod 2.
//
//    Input/output, int64_t *B_DEG, the degree of the polynomial B.
//
//    Input/output, int64_t B[MAXDEG+1].  On input, B is the polynomial 
//    defined in section 2.3 of BFN.  The degree of B implicitly defines 
//    the parameter J of section 3.3, by degree(B) = E*(J-1).  On output,
//    B has been multiplied by PX, so its degree is now E * J.
//
//    Output, int64_t V[MAXV+1], the computed V array.
//
//  Local Parameters:
//
//    Local, int64_t ARBIT, indicates where the user can place
//    an arbitrary element of the field of order 2.  This means 
//    0 <= ARBIT < 2.  
//
//    Local, int64_t BIGM, is the M used in section 3.3.
//    It differs from the [little] m used in section 2.3,
//    denoted here by M.
//
//    Local, int64_t NONZER, shows where the user must put an arbitrary 
//    non-zero element of the field.  For the code, this means 
//    0 < NONZER < 2.
//
{
  static int64_t arbit = 1;
  int64_t bigm;
  // int64_t e = px_deg; // не используется
  int64_t h[MAXDEG+1];
  int64_t h_deg;
  int64_t i;
  // int64_t j; // не используется
  int64_t kj;
  int64_t m;
  static int64_t nonzer = 1;
  // static int64_t p = 2; // не используется
  // static int64_t q = 2; // не используется
  int64_t pb_deg;
  int64_t r;
  int64_t term;
//
//
//  The polynomial H is PX**(J-1), which is the value of B on arrival.
//
//  In section 3.3, the values of Hi are defined with a minus sign:
//  don't forget this if you use them later!
//
  h_deg = *b_deg;

  for ( i = 0; i <= h_deg; i++ )
  {
    h[i] = b[i];
  }

  bigm = h_deg;
//
//  Multiply B by PX so B becomes PX**J.
//  In section 2.3, the values of Bi are defined with a minus sign:
//  don't forget this if you use them later!
//
  pb_deg = *b_deg;

  plymul2 ( px_deg, px, pb_deg, b, &pb_deg, b );

  *b_deg = pb_deg;
  m = *b_deg;
//
//  We don't use J explicitly anywhere, but here it is just in case.
//
//  j = m / e;
//
//  Now choose a value of Kj as defined in section 3.3.
//  We must have 0 <= Kj < E*J = M.
//  The limit condition on Kj does not seem very relevant
//  in this program.
//
  kj = bigm;
//
//  Choose values of V in accordance with the conditions in section 3.3.
//
  for ( r = 0; r < kj; r++ )
  {
    v[r] = 0;
  }
  v[kj] = 1;

  if ( kj < bigm )
  {
    term = sub [ 0 ] [ h[kj] ];

    for ( r = kj+1; r <= bigm-1; r++ )
    {
      v[r] = arbit;
//
//  Check the condition of section 3.3,
//  remembering that the H's have the opposite sign.
//
      term = sub [ term ] [ mul [ h[r] ] [ v[r] ] ];

    }
//
//  Now V(BIGM) is anything but TERM.
//
    v[bigm] = add [ nonzer] [ term ];

    for ( r = bigm+1; r <= m-1; r++ )
    {
      v[r] = arbit;
    }
  }
  else
  {
    for ( r = kj+1; r <= m-1; r++ )
    {
      v[r] = arbit;
    }

  }
//
//  Calculate the remaining V's using the recursion of section 2.3,
//  remembering that the B's have the opposite sign.
//
  for ( r = 0; r <= maxv - m; r++ )
  {
    term = 0;
    for ( i = 0; i <= m-1; i++ )
    {
      term = sub [ term] [ mul [ b[i] ] [ v[r+i] ] ];
    }
    v[r+m] = term;
  }

  return;
}
//****************************************************************************80

std::vector<Real> NiederreiterBaseTwo::GeneratePoint ()

//****************************************************************************80
//
//  Purpose:
//
//    NIEDERREITER2 returns an element of the Niederreiter sequence base 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 March 2003
//
//  Author:
//
//    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Harald Niederreiter,
//    Low-discrepancy and low-dispersion sequences,
//    Journal of Number Theory,
//    Volume 30, 1988, pages 51-70.
//
//  Parameters:
//
//    Input, int64_t DIM_NUM, the dimension of the sequence to be generated.
//
//    Input/output, int64_t *SEED, the index of the element entry to
//    compute.  On output, SEED is typically reset by this routine
//    to SEED+1.
//
//    Output, long double QUASI[DIM_NUM], the next quasirandom vector.
//
//  Local Parameters:
//
//    Local, int64_t CJ(DIM_MAX,0:NBITS-1), the packed values of 
//    Niederreiter's C(I,J,R).
//
//    Local, int64_t DIM_SAVE, the spatial dimension of the sequence
//    as specified on an initialization call.
//
//    Local, int64_t COUNT, the index of the current item in the sequence,
//    expressed as an array of bits.  COUNT(R) is the same as Niederreiter's
//    AR(N) (page 54) except that N is implicit.
//
//    Local, int64_t NEXTQ[DIM_MAX], the numerators of the next item in the
//    series.  These are like Niederreiter's XI(N) (page 54) except that
//    N is implicit, and the NEXTQ are integers.  To obtain
//    the values of XI(N), multiply by RECIP.
//
{  
  int64_t gray;
  int64_t i;
  int64_t r;
  // int64_t skip; // не используется
  
//
//  Initialization.
//
  if ( dim_save < 1LL || dimension != dim_save || seed <= 0LL )
  {
    if ( dimension <= 0LL || DIM_MAX < dimension )
    {
      std::cout << "\n";
      std::cout << "NIEDERREITER2 - Fatal error!\n";
      std::cout << "  Bad spatial dimension.\n";
      std::exit ( 1 );
    }

    dim_save = dimension;

    if ( seed < 0LL )
    {
      seed = 0LL;
    }

    seed_save = seed;
//
//  Calculate the C array.
//
    calcc2 ();
  }
//
//  Set up NEXTQ appropriately, depending on the Gray code of SEED.
//
//  You can do this every time, starting NEXTQ back at 0,
//  or you can do it once, and then carry the value of NEXTQ
//  around from the previous computation.
//
  if ( seed != seed_save + 1LL )
  {
    gray = ( seed ) ^ ( seed >> 1 );

    for ( i = 0; i < dim_save; i++ )
    {
      nextq[i] = 0;
    }

    r = 0;

    while ( gray != 0 )
    {
      if ( ( gray % 2 ) != 0 )
      {
        for ( i = 0; i < dim_save; i++ )
        {
          nextq[i] = ( nextq[i] ) ^ ( cj[i][r] );
        }
      }
      gray /= 2;
      r = r + 1;
    }
  }
//
//  Multiply the numerators in NEXTQ by RECIP to get the next
//  quasi-random vector.
//
    std::vector<long double> quasi(dim_save);
  for ( i = 0; i < dim_save; i++ )
  {
    quasi[i] = ( ( long double ) nextq[i] ) * RECIP;
  }
//
//  Find the position of the right-hand zero in SEED.  This
//  is the bit that changes in the Gray-code representation as
//  we go from SEED to SEED+1.
//
  r = 0;
  i = seed;

  while ( ( i % 2 ) != 0 )
  {
    r = r + 1;
    i = i / 2;
  }
//
//  Check that we have not passed 2**NBITS calls.
//
  if ( NBITS <= r )
  {
    std::cout << "\n";
    std::cout << "NIEDERREITER2 - Fatal error!\n";
    std::cout << "  Too many calls!\n";
    std::exit ( 1 );
  }
//
//  Compute the new numerators in vector NEXTQ.
//
  for ( i = 0; i < dim_save; i++ )
  {
    nextq[i] = ( nextq[i] ) ^ ( cj[i][r] );
  }

  seed_save = seed;
  seed = seed + 1;

  return quasi;
}