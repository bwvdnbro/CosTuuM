/**
 * @file CTM.cpp
 *
 * @brief Main program entry point.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Configuration.hpp"
#include "Matrix.hpp"
#include "SpecialFunctions.hpp"
#include "TMatrix.hpp"

#include <cinttypes>
#include <cmath>
#include <complex>
#include <vector>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Main program entry point.
 *
 * For now, does the entire T-matrix calculation for the Mishchenko default
 * parameters and outputs the scattering and extinction factors.
 *
 * @param argc Number of command line arguments (currently ignored).
 * @param argv Command line arguments (currently ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// input parameters (should become real parameters at some point
  // size of the particle (in same units as the wavelength)
  const float_type axi = 10.;
  // ratio between the equal surface area sphere radius and equal volume sphere
  // radius (is recomputed if not equal to 1)
  float_type ratio_of_radii = 0.1;
  // wavelength of incoming radiation (in same units as the particle size)
  const float_type wavelength = 2. * M_PI;
  // refractive index
  const std::complex<float_type> mr(1.5, 0.02);
  // ratio of horizontal and vertical axis of the spheroidal particle
  const float_type axis_ratio = 0.5;
  // tolerance for the calculation
  const float_type tolerance = 1.e-4;
  // number of Gauss-Legendre points to use as a multiplier of the maximum
  // order of spherical harmonics used to decompose the electric field, during
  // the first loop of the algorithm
  const uint_fast32_t ndgs = 2;

  /// hardcoded program parameters
  // maximum number of iterations during the first loop
  const uint_fast32_t maximum_order = 200;
  // maximum number of iterations during the second loop
  const uint_fast32_t maximum_ngauss = 500;

  // make sure 'rat' contains the right ratio if it is not 1
  if (abs(ratio_of_radii - 1.) > 1.e-8) {
    ratio_of_radii =
        SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
            axis_ratio);
  }

  // R_V is the equivalent sphere radius
  const float_type R_V = ratio_of_radii * axi;
  // we need a reasonable initial guess for the order of the spherical harmonics
  // the below expression provides this guess
  const float_type xev = 2. * M_PI * R_V / wavelength;
  uint_fast32_t nmax = static_cast<uint_fast32_t>(
      std::max(float_type(4.), xev + 4.05 * cbrt(xev)));

  // we need to find a maximum expansion order and number of Gauss-Legendre
  // quadrature points that leads to a converged T-matrix
  // To this end, we will start with a reasonable guess for the order, and then
  // keep track of how the accuracy of the scattering and extinction factors
  // changes as we first increase the order and then the number of quadrature
  // points
  // To simplify the calculation, we only compute the m=0 degree terms in the
  // T-matrix during the convergence loops
  TMatrix *active_Tmatrix = nullptr;

  // loop control variables
  float_type old_qext = 0.;
  float_type old_qsca = 0.;
  float_type dext = 1.;
  float_type dsca = 1.;
  while (nmax < maximum_order && (dext > tolerance || dsca > tolerance)) {
    // initially we assume that the number of quadrature points is a fixed
    // multiple of the order
    const uint_fast32_t ngauss = ndgs * nmax;

    // delete the old matrix (if it exists) and create a new one
    delete active_Tmatrix;
    active_Tmatrix = new TMatrix(wavelength, mr, R_V, axis_ratio, nmax, ngauss);
    TMatrix &T = *active_Tmatrix;

    // calculate the scattering and extinction factors for this iteration
    float_type qsca = 0.;
    float_type qext = 0.;
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      const float_type dn1 = 2. * n + 1.;
      qsca += dn1 *
              (std::norm(T(0, n, 0, 0, n, 0)) + std::norm(T(1, n, 0, 1, n, 0)));
      qext += dn1 * (T(0, n, 0, 0, n, 0).real() + T(1, n, 0, 1, n, 0).real());
    }
    // compute the relative difference w.r.t. the previous values
    dsca = abs((old_qsca - qsca) / qsca);
    dext = abs((old_qext - qext) / qext);
    old_qext = qext;
    old_qsca = qsca;

    // some (temporary) diagnostic output
    ctm_warning("dsca: %g, dext: %g", double(dsca), double(dext));

    // increase the order
    ++nmax;
  }

  // check if we managed to converge
  if (nmax == maximum_order) {
    ctm_error("Unable to converge!");
  } else {
    // correct for overshoot in last iteration
    --nmax;
    ctm_warning("Converged for nmax = %" PRIuFAST32, nmax);
  }

  // now start increasing the number of quadrature points
  uint_fast32_t ngauss = nmax * ndgs + 1;
  dext = tolerance + 1.;
  dsca = tolerance + 1.;
  while (ngauss < maximum_ngauss && (dext > tolerance || dsca > tolerance)) {

    // delete the old matrix and create a new one
    delete active_Tmatrix;
    active_Tmatrix = new TMatrix(wavelength, mr, R_V, axis_ratio, nmax, ngauss);
    TMatrix &T = *active_Tmatrix;

    // calculate the scattering and extinction factors
    float_type qsca = 0.;
    float_type qext = 0.;
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      const float_type dn1 = 2. * n + 1.;
      qsca += dn1 *
              (std::norm(T(0, n, 0, 0, n, 0)) + std::norm(T(1, n, 0, 1, n, 0)));
      qext += dn1 * (T(0, n, 0, 0, n, 0).real() + T(1, n, 0, 1, n, 0).real());
    }
    // compute the relative difference w.r.t. the old values
    dsca = abs((old_qsca - qsca) / qsca);
    dext = abs((old_qext - qext) / qext);
    old_qext = qext;
    old_qsca = qsca;

    // some diagnostic output
    ctm_warning("dsca: %g, dext: %g", double(dsca), double(dext));

    // increase the number of quadrature points
    ++ngauss;
  }

  // check if we reached convergence
  if (ngauss == maximum_ngauss) {
    ctm_error("Unable to converge!");
  } else {
    // correct for overshoot in final iteration
    --ngauss;
    ctm_warning("Converged for ngauss = %" PRIuFAST32, ngauss);
  }

  // ok: we found the right order and number of quadrature points
  // we now need to compute the additional elements for m=/=0
  active_Tmatrix->compute_additional_elements();

  // compute the actual scattering and extinction factors using the full matrix
  TMatrix &T = *active_Tmatrix;
  old_qsca = 0.;
  old_qext = 0.;
  for (uint_fast32_t n1 = 1; n1 < nmax + 1; ++n1) {
    for (uint_fast32_t n2 = 1; n2 < nmax + 1; ++n2) {
      for (int_fast32_t m1 = -n1; m1 < static_cast<int_fast32_t>(n1 + 1);
           ++m1) {
        for (int_fast32_t m2 = -n2; m2 < static_cast<int_fast32_t>(n2 + 1);
             ++m2) {
          old_qsca += std::norm(T(0, n1, m1, 0, n2, m2));
          old_qsca += std::norm(T(0, n1, m1, 1, n2, m2));
          old_qsca += std::norm(T(1, n1, m1, 0, n2, m2));
          old_qsca += std::norm(T(1, n1, m1, 1, n2, m2));
        }
      }
    }
  }
  for (uint_fast32_t n1 = 1; n1 < nmax + 1; ++n1) {
    for (int_fast32_t m1 = -n1; m1 < static_cast<int_fast32_t>(n1 + 1); ++m1) {
      old_qext += T(0, n1, m1, 0, n1, m1).real();
      old_qext += T(0, n1, m1, 1, n1, m1).real();
      old_qext += T(1, n1, m1, 0, n1, m1).real();
      old_qext += T(1, n1, m1, 1, n1, m1).real();
    }
  }
  // output the factors
  ctm_warning("qsca: %g", double(old_qsca));
  ctm_warning("qext: %g", double(old_qext));
  ctm_warning("walb: %g", double(-old_qsca / old_qext));

  // clean up
  delete active_Tmatrix;

  // done!
  return 0;
}
