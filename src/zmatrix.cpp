#include "zmatrix.hpp"
#include "atom.hpp"
#include "geograd.hpp"
#include "htab.hpp"
#include "units.h"
#include "dmat.hpp"
#include "random.h"
#include "linalg.hpp"

void ZMatrix::check(int n) const
{
  Vec<bool> seen(n, false);
  for (int m = 0; m < size(); m++) {
    if ((*this)[m].i < 0 || (*this)[m].i >= n)
      die("ZMatrix::check: index i=%d is out of range", (*this)[m].i);
    if (m > 0 && ((*this)[m].a < 0 || (*this)[m].a >= n))
      die("ZMatrix::check: index a=%d is out of range", (*this)[m].a);
    if (m > 1 && ((*this)[m].b < 0 || (*this)[m].b >= n))
      die("ZMatrix::check: index b=%d is out of range", (*this)[m].b);
    if (m > 2 && ((*this)[m].c < 0 || (*this)[m].c >= n))
      die("ZMatrix::check: index c=%d is out of range", (*this)[m].c);
    if (seen[(*this)[m].i])
      die("ZMatrix::check: repeat index i=%d", (*this)[m].i);
    if (m > 0 && !seen[(*this)[m].a])
      die("ZMatrix::check: index a=%d has not been referenced", (*this)[m].a);
    if (m > 1 && !seen[(*this)[m].b])
      die("ZMatrix::check: index b=%d has not been referenced", (*this)[m].b);
    if (m > 2 && !seen[(*this)[m].c])
      die("ZMatrix::check: index c=%d has not been referenced", (*this)[m].c);
    seen[(*this)[m].i] = true;
  }
}

void ZMatrix::to_cartesian(Atom *at) const
{
  Cartesian a(0,0,0), b(0,0,1), c(0,1,0);
  double r = 0, theta = 0, phi = 0;
  for (int m = 0; m < size(); m++) {
    if (m > 0) {
      a = at[(*this)[m].a].position;
      r = (*this)[m].r;
    }
    if (m > 1) {
      b = at[(*this)[m].b].position;
      theta = degrees_to_radians((*this)[m].theta);
    }
    if (m > 2) {
      c = at[(*this)[m].c].position;
      phi = degrees_to_radians((*this)[m].phi);
    }
    at[(*this)[m].i].position = ZLocation(a, b, c, r, theta, phi);
  }
}

void ZMatrix::from_cartesian(const Atom *at)
{
  Cartesian a(0,0,0), b(0,0,1), c(0,1,0);
  for (int m = 0; m < size(); m++) {
    if (m > 0)
      a = at[(*this)[m].a].position;
    if (m > 1)
      b = at[(*this)[m].b].position;
    if (m > 2)
      c = at[(*this)[m].c].position;
    const Cartesian i = at[(*this)[m].i].position;
    (*this)[m].r = i.distance(a);
    (*this)[m].theta = radians_to_degrees(Angle(i,a,b));
    (*this)[m].phi = radians_to_degrees(Dihedral(i,a,b,c));
  }
}

void ZMatrix::check_ignore_dihedral(const HSet<int> &ignore_dihedral) const
{
  HSetItr<int> i;
  for (i.init(ignore_dihedral); i.ok(); i.next())
    insist(3 <= i() && i() < size());
}

void ZMatrix::to_internal(DVec &q, const HSet<int> &ignore_dihedral) const
{
  insist(q.size() == degrees_of_freedom() - ignore_dihedral.size());
  check_ignore_dihedral(ignore_dihedral);
  int i = 0;
  for (int m = 0; m < size(); m++) {
    if (m > 0)
      q[i++] = (*this)[m].r;
    if (m > 1)
      q[i++] = degrees_to_radians((*this)[m].theta);
    if (m > 2 && !ignore_dihedral.exists(m))
      q[i++] = degrees_to_radians((*this)[m].phi);
  }
  insist(i == q.size());
}

void ZMatrix::from_internal(const DVec &q, const HSet<int> &ignore_dihedral)
{
  insist(q.size() == degrees_of_freedom() - ignore_dihedral.size());
  check_ignore_dihedral(ignore_dihedral);
  int i = 0;
  for (int m = 0; m < size(); m++) {
    if (m > 0)
      (*this)[m].r = q[i++];
    if (m > 1)
      (*this)[m].theta = radians_to_degrees(q[i++]);
    if (m > 2 && !ignore_dihedral.exists(m))
      (*this)[m].phi = radians_to_degrees(q[i++]);
  }
  insist(i == q.size());
}

ostream & operator<<(ostream &s, const ZMatrix &z)
{
  s << "{\n";
  IndentPush();
  for (int m = 0; m < z.size(); m++) {
    Out() << Indent() << "{ " << "i " << z[m].i;
    if (m > 0)
      Out() << "  a " << z[m].a << "  r " << z[m].r;
    if (m > 1)
      Out() << "  b " << z[m].b << "  theta " << z[m].theta;
    if (m > 2)
      Out() << "  c " << z[m].c << "  phi " << z[m].phi;
    Out() << " }\n";
  }
  IndentPop();
  return s << Indent() << "}\n";
}

int ZMatrix::degrees_of_freedom() const
{
  return number_of_bonds() + number_of_angles() + number_of_dihedrals();
}

int ZMatrix::number_of_bonds() const
{
  const int n = size() - 1;
  return n > 0 ? n : 0;
}

int ZMatrix::number_of_angles() const
{
  const int n = size() - 2;
  return n > 0 ? n : 0;
}

int ZMatrix::number_of_dihedrals() const
{
  const int n = size() - 3;
  return n > 0 ? n : 0;
}

double ZMatrix::jacobian() const
{
  switch (size()) {
  case 0:
  case 1:
  case 2:
    return 1;
  case 3:
    return (*this)[1].r * (*this)[2].r;
  default:
    double j = sq((*this)[1].r);
    for (int m = 2; m < size(); m++)
      j *= sq((*this)[m].r) * sin(degrees_to_radians((*this)[m].theta));
    return j;
  }
}

/***
 * Calculate \int exp(-(1/2) q cov^-1 q) J(q) dq 
 * by a Taylor expansion of J(q) up to fourth order
 ***/
void ZMatrix::harmonic_partition_function(const DVec &avgpos,
					  const DMat &cov,
					  const HSet<int> &ignore_dihedral,
					  double &logz0,
					  double &logz2,
					  double &logz4)
{
  const int m = degrees_of_freedom() - ignore_dihedral.size();
  insist(avgpos.size() == m);
  insist(cov.rows() == m);
  insist(cov.is_square());
  insist(ignore_dihedral.size() <= number_of_dihedrals());
  from_internal(avgpos,ignore_dihedral);
  logz0 = log (jacobian()) + LogDeterminant(cov)/2 + (m/2.0+ignore_dihedral.size())*log (2*M_PI);
  logz2 = logz4 = 0;
  if (size() < 3)
    return;
  const int n = number_of_bonds() + number_of_angles();
  insist(n > 1);
  DVec d(n, 0.0), d2(n, 0.0), d3(n, 0.0), d4(n, 0.0);
  DMat c(n,n);
  Vec<int> in(n);
  int i, j;
  if (size() == 3) {
    insist(m == 3);
    insist(n == 3);
    d[0] = 1/(*this)[1].r;
    d[1] = 1/(*this)[2].r;
    for (i = 0; i < 3; i++)
      in[i] = i;
  } else {
    i = 0;
    j = 0;
    for (int ia = 1; ia < size(); ia++) {
      const double r = (*this)[ia].r;
      d[i] = 2/r;
      d2[i] = 2/sq(r);
      in[i] = j;
      i++;
      j++;
      if (ia == 1)
	continue;
      const double t = degrees_to_radians((*this)[ia].theta);
      const double cot = cos(t)/sin(t);
      d[i] = cot;
      d2[i] = -1;
      d3[i] = -cot;
      d4[i] = 1;
      in[i] = j;
      i++;
      j++;
      if (ia == 2 || ignore_dihedral.exists(ia))
	continue;
      j++;
    }
    insist(i == n);
    insist(j == m);
  }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      c(i,j) = cov(in[i],in[j]);
#define C(i,j,k,l) (c(i,j)*c(k,l)+c(i,k)*c(j,l)+c(i,l)*c(j,k))
  double a2 = 0, a4 = 0;
  for (i = 0; i < n; i++) {
    a2 += d2[i]*c(i,i);
    a4 += d4[i]*C(i,i,i,i);
    for (j = i+1; j < n; j++) {
      a2 += 2*d[i]*d[j]*c(i,j);
      a4 += 4*(d3[i]*d[j]*C(i,i,i,j) + d[i]*d3[j]*C(i,j,j,j)) +
	6*d2[i]*d2[j]*C(i,i,j,j);
      for (int k = j+1; k < n; k++) {
	a4 += 12*(d2[i]*d[j]*d[k]*C(i,i,j,k) +
		  d[i]*d2[j]*d[k]*C(i,j,j,k) +
		  d[i]*d[j]*d2[k]*C(i,j,k,k));
	for (int l = k+1; l < n; l++)
	  a4 += 24*d[i]*d[j]*d[k]*d[l]*C(i,j,k,l);
      }
    }
  }
#undef C
  logz2 = log(1+a2/2);
  logz4 = log(1+a2/2+a4/24);
  Out() << "Ignored dihedrals: " << ignore_dihedral <<"\n";
}
