#include "euler.hpp"
#include "coord.hpp"
#include "fns.h"

/***
 * We use the "x-convention" for Euler angles,
 * as described in Goldstein 4-4 (p. 145)
 * or Marion & Thornton, 11.7 (p. 431)
 ***/

static double ArcCos(double x)
{
  if (x < -1)
    return acos(-1.0);
  else if (x > 1)
    return acos(1.0);
  else
    return acos(x);
}

static Tensor Skew(const Cartesian &u)
{
  return Tensor(0,-u.z,u.y,
		u.z,0,-u.x,
		-u.y,u.x,0);
}

static double Score(const Tensor &a, const Tensor &b)
{
  return sq(a.xx-a.xx) + sq(a.xy-b.xy) + sq(a.xz-b.xz) +
    sq(a.yx-a.yx) + sq(a.yy-b.yy) + sq(a.yz-b.yz) +
    sq(a.zx-a.zx) + sq(a.zy-b.zy) + sq(a.zz-b.zz);
}

void EulerAnglesToRotationMatrix(Tensor &m, double phi, double theta, double psi)
{
  const double ch = cos(phi);
  const double sh = sin(phi);
  const double ct = cos(theta);
  const double st = sin(theta);
  const double cs = cos(psi);
  const double ss = sin(psi);
  m.xx = cs*ch - ct*sh*ss;
  m.yx = -ss*ch - ct*sh*cs;
  m.zx = st*sh;
  m.xy = cs*sh + ct*ch*ss;
  m.yy = -ss*sh + ct*ch*cs;
  m.zy = -st*ch;
  m.xz = ss*st;
  m.yz = cs*st;
  m.zz = ct;
}

void RotationMatrixToEulerAngles(const Tensor &m, double &phi, double &theta, double &psi)
{
  Tensor u;
  double t[2], h[2], s[2];
  t[0] = ArcCos(m.zz);
  if (is_almost_zero(t[0])) {
    double best = 1e8;
    theta = psi = 0;
    h[0] = ArcCos(m.xx);
    h[1] = -h[0];
    for (int j = 0; j < 2; j++) {
      EulerAnglesToRotationMatrix(u, h[j], 0, 0);
      const double sc = Score(m,u);
      if (sc < best) {
	phi = h[j];
	best = sc;
      }
      if (is_almost_zero(best))
	return;
    }
  } else {
    t[1] = -t[0];
    double best = 1e8;
    for (int i = 0; i < 2; i++) {
      const double st = sin(t[i]);
      h[0] = ArcCos(-m.zy/st);
      h[1] = -h[0];
      s[0] = ArcCos(m.yz/st);
      s[1] = -s[0];
      for (int j = 0; j < 2; j++)
	for(int k = 0; k < 2; k++) {
	  EulerAnglesToRotationMatrix(u, h[j], t[i], s[k]);
	  const double sc = Score(m,u);
	  if (sc < best) {
	    phi = h[j];
	    theta = t[i];
	    psi = s[k];
	    best = sc;
	  }
	  if (is_almost_zero(best))
	    return;
	}
    }
  }
}

void AxisAngleToRotationMatrix(Tensor &t, const Cartesian &u, double theta)
{
  const Tensor n = Skew(u.as_unit_vector());
  t = Tensor(1) + sin(theta)*n + (1-cos(theta))*(n*n);
}

void EulerAngleJacobian(double phi, double theta, double psi,
			Tensor &dphi, Tensor &dtheta, Tensor &dpsi)
{
  const double ch = cos(phi);
  const double sh = sin(phi);
  const double ct = cos(theta);
  const double st = sin(theta);
  const double cs = cos(psi);
  const double ss = sin(psi);
  dphi = Tensor(-cs*sh - ch*ct*ss, ch*cs - ct*sh*ss, 0,
		-ch*cs*ct + sh*ss, -cs*ct*sh - ch*ss, 0,
		ch*st, sh*st, 0);
  dtheta = Tensor(sh*ss*st, -ch*ss*st, ct*ss,
		  cs*sh*st, -ch*cs*st, cs*ct,
		  ct*sh, -ch*ct, -st);
  dpsi = Tensor(-cs*ct*sh - ch*ss, ch*cs*ct - sh*ss, cs*st,
		-ch*cs + ct*sh*ss, -cs*sh - ch*ct*ss, -ss*st,
		0,0,0);
}
