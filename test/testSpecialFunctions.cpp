/**
 * @file testSpecialFunctions.cpp
 *
 * @brief Unit test for the special functions in SpecialFunctions.hpp.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "Configuration.hpp"
#include "SpecialFunctions.hpp"
#include <fstream>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/*! @brief Maximum order to test for Bessel and Wigner D functions,
 *  @f$n_{max}@f$. */
const uint_fast32_t TESTSPECIALFUNCTIONS_NMAX = 80;

/*! @brief Order of Gauss-Legendre quadrature to test. */
const uint_fast32_t TESTSPECIALFUNCTIONS_NGAUSS = 200;

/**
 * @brief Basic function to test the Gauss-Legendre quadrature routine.
 *
 * @param x X value.
 * @param args Additional arguments (ignored).
 * @return Test function value.
 */
static float_type function(const float_type x, void *args) {

  const float_type cosx = cos(x);
  const float_type px = 0.5 + 0.1 * (cosx * cosx - 1.);
  float_type dnx[10], ddnx[10];
  SpecialFunctions::wigner_dn_0m(cosx, 10, 0, dnx, ddnx);
  return sin(x) * px * dnx[9];
}

/**
 * @brief Unit test for the special functions in SpecialFunctions.hpp.
 *
 * We call the functions SpecialFunctions::spherical_j_jdj_array() and
 * SpecialFunctions::spherical_y_ydy_array() for a range of input values (real
 * and complex for the first kind), and write them to two text files, called
 * test_bessel_complex.txt and test_bessel_real.txt. The accompanying script
 * plot_test_special_functions.py reads these files and computes the same
 * functions for the same input values using the Bessel functions that are part
 * of scipy.special. It then plots both versions and the relative difference
 * between the two. Note that we need to test both the first order and
 * @f$n_max@f$th order derivatives, as they use different bits of code.
 *
 * We then test the Wigner D function by calling
 * SpecialFunctions::wigner_dn_0m() on a range of input values, and write the
 * functions and their derivatives to a text file called test_wigner_d.txt.
 * The accompanying script reads this file and computes the same functions for
 * the same input values using the associated Legendre polynomials in
 * scipy.special and the relation between Wigner D functions and associated
 * Legendre polynomials. It then plots both versions and the relative
 * difference between the two. We test both @f$m = 0@f$ and @f$m = n_{max}@f$
 * to cover all possible paths through the code.
 *
 * We then test
 * SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio() by
 * calling it with 3 different axis ratios (1, <1, >1) and comparing with
 * results obtained with Python.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  /// Bessel functions

  // open the output text files
  std::ofstream cfile("test_bessel_complex.txt");
  std::ofstream rfile("test_bessel_real.txt");
  // loop over a spatial range
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    // this expression should match the np.arange(0.05, 100., 0.1) in the
    // Python script
    const float_type x = 0.05 + 0.1 * i;
    // we create a complex argument to test the complex version
    const std::complex<float_type> z(x, x);
    // create output arrays for the complex version test
    std::complex<float_type> j[TESTSPECIALFUNCTIONS_NMAX],
        dj[TESTSPECIALFUNCTIONS_NMAX];
    // call the spherical Bessel function of the first kind for complex input
    // and output
    SpecialFunctions::spherical_j_jdj_array(TESTSPECIALFUNCTIONS_NMAX, z, j,
                                            dj);
    // write the 1st and nth order function to the complex output file
    cfile << j[0].real() << "\t" << j[0].imag() << "\t" << dj[0].real() << "\t"
          << dj[0].imag() << "\n";
    cfile << j[TESTSPECIALFUNCTIONS_NMAX - 1].real() << "\t"
          << j[TESTSPECIALFUNCTIONS_NMAX - 1].imag() << "\t"
          << dj[TESTSPECIALFUNCTIONS_NMAX - 1].real() << "\t"
          << dj[TESTSPECIALFUNCTIONS_NMAX - 1].imag() << "\n";
    // create output arrays for the real version tests
    float_type jr[TESTSPECIALFUNCTIONS_NMAX], djr[TESTSPECIALFUNCTIONS_NMAX],
        yr[TESTSPECIALFUNCTIONS_NMAX], dyr[TESTSPECIALFUNCTIONS_NMAX];
    // call the spherical Bessel function of the first kind for real input
    // and output
    SpecialFunctions::spherical_j_jdj_array(TESTSPECIALFUNCTIONS_NMAX, x, jr,
                                            djr);
    // call the spherical Bessel function of the second kind for real input
    // and output
    SpecialFunctions::spherical_y_ydy_array(TESTSPECIALFUNCTIONS_NMAX, x, yr,
                                            dyr);
    // write the 1st and nth order function to the real output file
    rfile << jr[0] << "\t" << djr[0] << "\t" << yr[0] << "\t" << dyr[0] << "\n";
    rfile << jr[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
          << djr[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
          << yr[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
          << dyr[TESTSPECIALFUNCTIONS_NMAX - 1] << "\n";
  }

  /// Wigner D function

  // open the output text file
  std::ofstream dfile("test_wigner_d.txt");
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    // this expression should match the np.arange(-0.999, 1., 0.002) expression
    // in the Python script
    const float_type cosx = -0.999 + 0.002 * i;
    // create output arrays
    float_type d[TESTSPECIALFUNCTIONS_NMAX], dd[TESTSPECIALFUNCTIONS_NMAX];
    // call the function with m=0
    SpecialFunctions::wigner_dn_0m(cosx, TESTSPECIALFUNCTIONS_NMAX, 0, d, dd);
    // write an output line
    dfile << d[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
          << dd[TESTSPECIALFUNCTIONS_NMAX - 1] << "\n";
    // call the function with m=nmax
    SpecialFunctions::wigner_dn_0m(cosx, TESTSPECIALFUNCTIONS_NMAX,
                                   TESTSPECIALFUNCTIONS_NMAX, d, dd);
    // write an output line
    dfile << d[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
          << dd[TESTSPECIALFUNCTIONS_NMAX - 1] << "\n";
  }

  // open the output text file
  std::ofstream dsfile("test_wigner_d_sinx.txt");
  for (uint_fast32_t i = 0; i < 1000; ++i) {
    // this expression should match the np.arange(-0.999, 1., 0.002) expression
    // in the Python script
    const float_type cosx = -0.999 + 0.002 * i;
    const float_type sinx = sqrt(1. - cosx * cosx);
    const float_type sinx_inv = 1. / sinx;
    // create output arrays
    float_type d[TESTSPECIALFUNCTIONS_NMAX], dd[TESTSPECIALFUNCTIONS_NMAX];
    // call the function with m=0
    SpecialFunctions::wigner_dn_0m_sinx(cosx, sinx, sinx_inv,
                                        TESTSPECIALFUNCTIONS_NMAX, 0, d, dd);
    // write an output line
    dsfile << d[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
           << dd[TESTSPECIALFUNCTIONS_NMAX - 1] << "\n";
    // call the function with m=nmax
    SpecialFunctions::wigner_dn_0m_sinx(cosx, sinx, sinx_inv,
                                        TESTSPECIALFUNCTIONS_NMAX,
                                        TESTSPECIALFUNCTIONS_NMAX, d, dd);
    // write an output line
    dsfile << d[TESTSPECIALFUNCTIONS_NMAX - 1] << "\t"
           << dd[TESTSPECIALFUNCTIONS_NMAX - 1] << "\n";
  }

  /// Equal volume sphere to equal surface area sphere radius ratio

  assert_condition(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
          float_type(1.)) == float_type(1.));
  const double prolate = double(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
          float_type(0.5)));
  assert_values_equal_rel(prolate, 0.9637112829756893, 1.e-10);
  const double oblate = double(
      SpecialFunctions::get_equal_volume_to_equal_surface_area_sphere_ratio(
          float_type(2.)));
  assert_values_equal_rel(oblate, 0.955443263377723, 1.e-10);

  /// Gauss-Legendre quadrature

  const double tolerance = 1.e-7;
  std::vector<float_type> x(TESTSPECIALFUNCTIONS_NGAUSS),
      w(TESTSPECIALFUNCTIONS_NGAUSS);
  // first test against the values on Wikipedia
  {
    // n = 1
    SpecialFunctions::get_gauss_legendre_points_and_weights(1, x, w);
    ctm_warning("x: [%g], w: [%g]", double(x[0]), double(w[0]));
    assert_condition(x[0] == 0.);
    assert_condition(w[0] == 2.);
  }
  {
    // n = 2
    SpecialFunctions::get_gauss_legendre_points_and_weights(2, x, w);
    ctm_warning("x: [%g, %g], w: [%g, %g]", double(x[0]), double(x[1]),
                double(w[0]), double(w[1]));
    assert_values_equal_rel(double(x[0]), -1. / std::sqrt(3.), tolerance);
    assert_values_equal_rel(double(w[0]), 1., tolerance);
    assert_values_equal_rel(double(x[1]), 1. / std::sqrt(3.), tolerance);
    assert_values_equal_rel(double(w[1]), 1., tolerance);
  }
  {
    // n = 3
    SpecialFunctions::get_gauss_legendre_points_and_weights(3, x, w);
    ctm_warning("x: [%g, %g, %g], w: [%g, %g, %g]", double(x[0]), double(x[1]),
                double(x[2]), double(w[0]), double(w[1]), double(w[2]));
    assert_values_equal_rel(double(x[0]), -std::sqrt(3. / 5.), tolerance);
    assert_values_equal_rel(double(w[0]), 5. / 9., tolerance);
    assert_values_equal_tol(double(x[1]), 0., tolerance);
    assert_values_equal_rel(double(w[1]), 8. / 9., tolerance);
    assert_values_equal_rel(double(x[2]), std::sqrt(3. / 5.), tolerance);
    assert_values_equal_rel(double(w[2]), 5. / 9., tolerance);
  }
  {
    // n = 4
    SpecialFunctions::get_gauss_legendre_points_and_weights(4, x, w);
    ctm_warning("x: [%g, %g, %g, %g], w: [%g, %g, %g, %g]", double(x[0]),
                double(x[1]), double(x[2]), double(x[3]), double(w[0]),
                double(w[1]), double(w[2]), double(w[3]));
    assert_values_equal_rel(double(x[0]),
                            -std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)),
                            tolerance);
    assert_values_equal_rel(double(w[0]), (18 - std::sqrt(30.)) / 36.,
                            tolerance);
    assert_values_equal_rel(double(x[1]),
                            -std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)),
                            tolerance);
    assert_values_equal_rel(double(w[1]), (18. + std::sqrt(30.)) / 36.,
                            tolerance);
    assert_values_equal_rel(double(x[2]),
                            std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.)),
                            tolerance);
    assert_values_equal_rel(double(w[2]), (18. + std::sqrt(30.)) / 36.,
                            tolerance);
    assert_values_equal_rel(double(x[3]),
                            std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)),
                            tolerance);
    assert_values_equal_rel(double(w[3]), (18 - std::sqrt(30.)) / 36.,
                            tolerance);
  }
  {
    // n = 5
    SpecialFunctions::get_gauss_legendre_points_and_weights(5, x, w);
    ctm_warning("x: [%g, %g, %g, %g, %g], w: [%g, %g, %g, %g, %g]",
                double(x[0]), double(x[1]), double(x[2]), double(x[3]),
                double(x[4]), double(w[0]), double(w[1]), double(w[2]),
                double(w[3]), double(w[4]));
    assert_values_equal_rel(double(x[0]),
                            -1. / 3. * std::sqrt(5. + 2. * std::sqrt(10. / 7.)),
                            tolerance);
    assert_values_equal_rel(double(w[0]), (322. - 13. * std::sqrt(70.)) / 900.,
                            tolerance);
    assert_values_equal_rel(double(x[1]),
                            -1. / 3. * std::sqrt(5. - 2. * std::sqrt(10. / 7.)),
                            tolerance);
    assert_values_equal_rel(double(w[1]), (322. + 13. * std::sqrt(70.)) / 900.,
                            tolerance);
    assert_values_equal_tol(double(x[2]), 0., tolerance);
    assert_values_equal_rel(double(w[2]), 128. / 225., tolerance);
    assert_values_equal_rel(double(x[3]),
                            1. / 3. * std::sqrt(5. - 2. * std::sqrt(10. / 7.)),
                            tolerance);
    assert_values_equal_rel(double(w[3]), (322. + 13. * std::sqrt(70.)) / 900.,
                            tolerance);
    assert_values_equal_rel(double(x[4]),
                            1. / 3. * std::sqrt(5. + 2. * std::sqrt(10. / 7.)),
                            tolerance);
    assert_values_equal_rel(double(w[4]), (322. - 13. * std::sqrt(70.)) / 900.,
                            tolerance);
  }
  {
    // interval [0, 1], n = 2
    SpecialFunctions::get_gauss_legendre_points_and_weights_ab<float_type>(
        2, 0., 1., x, w);
    ctm_warning("x: [%g, %g], w: [%g, %g]", double(x[0]), double(x[1]),
                double(w[0]), double(w[1]));
    assert_values_equal_rel(double(x[0]), 0.5 - 0.5 / std::sqrt(3.), tolerance);
    assert_values_equal_rel(double(w[0]), 0.5, tolerance);
    assert_values_equal_rel(double(x[1]), 0.5 + 0.5 / std::sqrt(3.), tolerance);
    assert_values_equal_rel(double(w[1]), 0.5, tolerance);
  }

  {
    // test integration of a simple function
    const float_type old_quad =
        SpecialFunctions::gauss_legendre_quadrature<float_type>(
            function, 0., M_PI, nullptr, 10, 1000, 1.e-10, 1.e-5);
    ctm_warning("Quad: %g", double(old_quad));
  }

  // now print additional values for a much higher order to analyse with the
  // script
  SpecialFunctions::get_gauss_legendre_points_and_weights(
      TESTSPECIALFUNCTIONS_NGAUSS, x, w);
  std::ofstream gfile("test_gauss_legendre_quadrature.txt");
  for (uint_fast32_t i = 0; i < TESTSPECIALFUNCTIONS_NGAUSS; ++i) {
    gfile << x[i] << "\t" << w[i] << "\n";
  }

  /// Clebsch-Gordan coefficients
  /// We test our implementation against the values for all coefficients
  /// with n1 = 1 and n2 = 1 provided on Wikipedia:
  /// https://en.wikipedia.org/wiki/Table_of_Clebsch%E2%80%93Gordan_coefficients
  {
    const float_type C111122 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 1, 1, 1,
                                                                     2, 2);
    assert_values_equal_rel(double(C111122), 1., 1.e-10);

    const float_type C111021 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 1, 1, 0,
                                                                     2, 1);
    assert_values_equal_rel(double(C111021), std::sqrt(0.5), 1.e-10);
    const float_type C101121 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 0, 1, 1,
                                                                     2, 1);
    assert_values_equal_rel(double(C101121), std::sqrt(0.5), 1.e-10);

    const float_type C111011 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 1, 1, 0,
                                                                     1, 1);
    assert_values_equal_rel(double(C111011), std::sqrt(0.5), 1.e-10);
    const float_type C101111 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 0, 1, 1,
                                                                     1, 1);
    assert_values_equal_rel(double(C101111), -std::sqrt(0.5), 1.e-10);

    const float_type C111m120 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 1, 1,
                                                                     -1, 2, 0);
    assert_values_equal_rel(double(C111m120), std::sqrt(1. / 6.), 1.e-10);
    const float_type C101020 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 0, 1, 0,
                                                                     2, 0);
    assert_values_equal_rel(double(C101020), std::sqrt(2. / 3.), 1.e-10);
    const float_type C1m11120 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, -1, 1,
                                                                     1, 2, 0);
    assert_values_equal_rel(double(C1m11120), std::sqrt(1. / 6.), 1.e-10);

    const float_type C111m110 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 1, 1,
                                                                     -1, 1, 0);
    assert_values_equal_rel(double(C111m110), std::sqrt(0.5), 1.e-10);
    const float_type C101010 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 0, 1, 0,
                                                                     1, 0);
    assert_condition(double(C101010) == 0.);
    const float_type C1m11110 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, -1, 1,
                                                                     1, 1, 0);
    assert_values_equal_rel(double(C1m11110), -std::sqrt(0.5), 1.e-10);

    const float_type C111m100 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 1, 1,
                                                                     -1, 0, 0);
    assert_values_equal_rel(double(C111m100), std::sqrt(1. / 3.), 1.e-10);
    const float_type C101000 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 0, 1, 0,
                                                                     0, 0);
    assert_values_equal_rel(double(C101000), -std::sqrt(1. / 3.), 1.e-10);
    const float_type C1m11100 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, -1, 1,
                                                                     1, 0, 0);
    assert_values_equal_rel(double(C1m11100), std::sqrt(1. / 3.), 1.e-10);

    const float_type C1m1102m1 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, -1, 1,
                                                                     0, 2, -1);
    assert_values_equal_rel(double(C1m1102m1), std::sqrt(0.5), 1.e-10);
    const float_type C101m12m1 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 0, 1,
                                                                     -1, 2, -1);
    assert_values_equal_rel(double(C101m12m1), std::sqrt(0.5), 1.e-10);

    const float_type C1m1101m1 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, -1, 1,
                                                                     0, 1, -1);
    assert_values_equal_rel(double(C1m1101m1), -std::sqrt(0.5), 1.e-10);
    const float_type C101m11m1 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 0, 1,
                                                                     -1, 1, -1);
    assert_values_equal_rel(double(C101m11m1), std::sqrt(0.5), 1.e-10);

    const float_type C1m11m12m2 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, -1, 1,
                                                                     -1, 2, -2);
    assert_values_equal_rel(double(C1m11m12m2), 1., 1.e-10);

    // some additional checks to check the n1<n2 cases
    const float_type C211122 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(2, 1, 1, 1,
                                                                     2, 2);
    assert_values_equal_rel(double(C211122), -std::sqrt(1. / 3.), 1.e-10);

    const float_type C112122 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 1, 2, 1,
                                                                     2, 2);
    assert_values_equal_rel(double(C112122), std::sqrt(1. / 3.), 1.e-10);

    const float_type C221032 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(2, 2, 1, 0,
                                                                     3, 2);
    assert_values_equal_rel(double(C221032), std::sqrt(1. / 3.), 1.e-10);
    const float_type C102232 =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(1, 0, 2, 2,
                                                                     3, 2);
    assert_values_equal_rel(double(C102232), std::sqrt(1. / 3.), 1.e-10);

    // check a large value, without reference, just to check nothing breaks
    const float_type large_coefficient =
        SpecialFunctions::get_clebsch_gordan_coefficient<float_type>(
            100, 40, 120, -40, 220, 0);
    ctm_warning("Random large value: %g", double(large_coefficient));
    assert_condition(double(large_coefficient) == double(large_coefficient));

    // test the other Clebsch-Gordan function for normal input values
    std::vector<float_type> C112 =
        SpecialFunctions::get_clebsch_gordan_coefficients<float_type>(1, 1, 2);
    assert_condition(C112.size() == 3);
    assert_condition(C112[0] == C1m11120);
    assert_condition(C112[1] == C101020);
    assert_condition(C112[2] == C111m120);

    // test the other Clebsch-Gordan function for unequal n1 and n2
    std::vector<float_type> C213 =
        SpecialFunctions::get_clebsch_gordan_coefficients<float_type>(2, 1, 3);
    assert_condition(C213.size() == 3);
    assert_values_equal_rel(double(C213[0]), std::sqrt(1. / 5.), 1.e-10);
    assert_values_equal_rel(double(C213[1]), std::sqrt(3. / 5.), 1.e-10);
    assert_values_equal_rel(double(C213[2]), std::sqrt(1. / 5.), 1.e-10);

    std::vector<float_type> C212 =
        SpecialFunctions::get_clebsch_gordan_coefficients<float_type>(2, 1, 2);
    assert_condition(C213.size() == 3);
    assert_values_equal_rel(double(C212[0]), -std::sqrt(0.5), 1.e-10);
    assert_values_equal_tol(double(C212[1]), 0., 1.e-10);
    assert_values_equal_rel(double(C212[2]), std::sqrt(0.5), 1.e-10);

    std::vector<float_type> C122 =
        SpecialFunctions::get_clebsch_gordan_coefficients<float_type>(1, 2, 2);
    assert_condition(C213.size() == 3);
    assert_values_equal_rel(double(C122[0]), -std::sqrt(0.5), 1.e-10);
    assert_values_equal_tol(double(C122[1]), 0., 1.e-10);
    assert_values_equal_rel(double(C122[2]), std::sqrt(0.5), 1.e-10);
  }

  return 0;
}
