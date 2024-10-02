#ifndef RANDOM_H
#define RANDOM_H

#ifdef __cplusplus
extern "C" {
#endif

/***
 * Seed the random number generator
 ***/
void RandomSeed(int);

/***
 * Return a random number sampled from uniform distribution on [0,1)
 ***/
double UniformRandom();

/*** 
 * Return a random number sampled from Gaussian distribution with unit variance 
 * That is, rho(x) = 1/sqrt(2 pi) exp(-x^2/2)
 *
 * To sample from rho(x) = sqrt(a/pi) exp(-a x^2),
 * take GaussianRandom()/sqrt(2a)
 ***/
double GaussianRandom();

/***
 * Return a random number sampled from exponential distribtion rho(x) = exp(-x)
 * To sample from rho(x) = k exp(-kx), divide result by k
 ***/
double ExponentialRandom();

#ifdef __cplusplus
}
#endif

#endif /* RANDOM_H */
