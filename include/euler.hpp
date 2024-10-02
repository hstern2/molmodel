#ifndef EULER_H
#define EULER_H

class Cartesian;
class Tensor;

/***
 * We use the "x-convention" for Euler angles,
 * as described in Goldstein 4-4 (p. 145)
 * or Marion & Thornton, 11.7 (p. 431)
 ***/
void EulerAnglesToRotationMatrix(Tensor &r, double phi, double theta, double psi);
void RotationMatrixToEulerAngles(const Tensor &r, double &phi, double &theta, double &psi);
void AxisAngleToRotationMatrix(Tensor &r, const Cartesian &u, double theta);

/***
 * Return the derivatives of a rotation matrix 
 * with respect to Euler angles
 ***/
void EulerAngleJacobian(double phi, double theta, double psi, 
			Tensor &dphi, Tensor &dtheta, Tensor &dpsi);

#endif /* EULER_H */
