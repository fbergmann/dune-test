#ifndef SBMLHELPER_HH
#define SBMLHELPER_HH

#include <math.h>
#include <stdarg.h>


double piecewise(
    double a, double b,
    double z)
{
  if (b != 0) return a;
  return z;
}

double piecewise(
    double a1, double b1,
    double a2, double b2,
    double z)
{
  if (b1 != 0) return a1;
  if (b2 != 0) return a2;
  return z;
}

double piecewise(
    double a1, double b1,
    double a2, double b2,
    double a3, double b3,
    double z)
{
  if (b1 != 0) return a1;
  if (b2 != 0) return a2;
  if (b3 != 0) return a3;
  return z;
}

double piecewise(
    double a1, double b1,
    double a2, double b2,
    double a3, double b3,
    double a4, double b4,
    double z)
{
  if (b1 != 0) return a1;
  if (b2 != 0) return a2;
  if (b3 != 0) return a3;
  if (b4 != 0) return a4;
  return z;
}


double piecewise(
    double a1, double b1,
    double a2, double b2,
    double a3, double b3,
    double a4, double b4,
    double a5, double b5,
    double z)
{
  if (b1 != 0) return a1;
  if (b2 != 0) return a2;
  if (b3 != 0) return a3;
  if (b4 != 0) return a4;
  if (b5 != 0) return a5;
  return z;
}


double geq(double a, double b)
{
  return a >= b;
}

double leq(double a, double b)
{
  return a <= b;
}

double eq(double a, double b)
{
  return a == b;
}

double neq(double a, double b)
{
  return a != b;
}

double gt(double a, double b)
{
  return a > b;
}

double lt(double a, double b)
{
  return a < b;
}

/*double and(double a, double b)
{
    return (double)((bool)a && (bool)b);
}*/

#endif // SBMLHELPER_HH
