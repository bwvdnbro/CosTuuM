/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2021 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * CosTuuM is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * CosTuuM is distributed in the hope that it will be useful, but WITOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with CosTuuM. If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

/**
 * @file FullScatteringMatrixTask.hpp
 *
 * @brief Task that computes full ScatteringMatrices.
 *
 * @author Bert Vandenbroucke (vandenbroucke@strw.leidenuniv.nl)
 */
#ifndef FULLSCATTERINGMATRIXTASK_HPP
#define FULLSCATTERINGMATRIXTASK_HPP

#include "Configuration.hpp"
#include "InteractionResource.hpp"
#include "QuickSchedWrapper.hpp"
#include "Result.hpp"
#include "TMatrixResource.hpp"

#include <vector>

/**
 * @brief Result of a full scattering matrix calculation.
 */
class FullScatteringMatrixResult : public Result {

  /*! @brief Give access to the computation task. */
  friend class FullScatteringMatrixTask;

  /*! @brief Give access to the averaging task. */
  friend class FullScatteringMatrixShapeAveragingTask;

private:
  /*! @brief Scattering matrices. */
  std::vector<float_type> _Z;

public:
  /**
   * @brief Constructor.
   *
   * @param composition Composition parameter value for the result.
   * @param size Particle size parameter value for the result (in m).
   * @param wavelength Wavelength value for the result (in m).
   * @param ntheta_in Number of input zenith angles.
   * @param ntheta_out Number of output zenith angles.
   * @param nphi Number of relative azimuth angles between the input and output
   * direction.
   */
  inline FullScatteringMatrixResult(const int_fast32_t composition,
                                    const float_type size,
                                    const float_type wavelength,
                                    const uint_fast32_t ntheta_in,
                                    const uint_fast32_t ntheta_out,
                                    const uint_fast32_t nphi)
      : Result(composition, size, wavelength, RESULTTYPE_FULLSCATTERINGMATRIX),
        _Z(ntheta_in * ntheta_out * nphi * 4) {}

  virtual ~FullScatteringMatrixResult() {}

  /**
   * @brief Get the size in memory of a hypothetical FullScatteringMatrixResult
   * object with the given parameters.
   *
   * @param ntheta_in Number of input zenith angles.
   * @param ntheta_out Number of output zenith angles.
   * @param nphi Number of relative azimuth angles between the input and output
   * direction.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size(const uint_fast32_t ntheta_in,
                                       const uint_fast32_t ntheta_out,
                                       const uint_fast32_t nphi) {
    size_t size = sizeof(FullScatteringMatrixResult);
    size += 4 * ntheta_in * ntheta_out * nphi * sizeof(float_type);
    return size;
  }

  /**
   * @brief Get the scattering matrix components for the given indices.
   *
   * @param itheta_in Input zenith angle index.
   * @param itheta_out Output zenith angle index.
   * @param iphi Relative azimuth angle index.
   * @param Z Array (of at least size 4) to store the non-trivial components
   * of the corresponding scattering matrix in.
   */
  inline void get_scattering_matrix(const uint_fast32_t itheta_in,
                                    const uint_fast32_t itheta_out,
                                    const uint_fast32_t iphi,
                                    float_type *Z) const {
    ctm_assert(4 * itheta_in * itheta_out * iphi + 3 < _Z.size());
    const uint_fast32_t index = 4 * itheta_in * itheta_out * iphi;
    Z[0] = _Z[index];
    Z[1] = _Z[index + 1];
    Z[2] = _Z[index + 2];
    Z[3] = _Z[index + 3];
  }
};

/**
 * @brief Angular grid used to compute full scattering matrices.
 *
 * The grid stores three angles and some useful function values based on
 * these:
 *  - the input and output zenith angles that correspond to the directions for
 *    which we want to compute scattering matrices.
 *  - the relative azimuth angle between the incoming and outgoing direction.
 */
class FullScatteringMatrixGrid : public Resource,
                                 public Task,
                                 public Computable {

  /*! @brief Give access to the computation task. */
  friend class FullScatteringMatrixTask;

  /*! @brief Give access to special Wigner D resources. */
  friend class FullScatteringMatrixSpecialWignerDResources;

private:
  /*! @brief Input zenith angles (in radians). */
  std::vector<float_type> _theta_in;

  /*! @brief Output zenith angles (in radians). */
  std::vector<float_type> _theta_out;

  /*! @brief Relative azimuth angles (in radians). */
  std::vector<float_type> _phi;

  /*! @brief Cosines of the input zenith angles. */
  std::vector<float_type> _cos_theta_in;

  /*! @brief Sines of the input zenith angles. */
  std::vector<float_type> _sin_theta_in;

  /*! @brief Inverse sines of the input zenith angles. */
  std::vector<float_type> _sin_theta_in_inverse;

  /*! @brief Cosines of the output zenith angles. */
  std::vector<float_type> _cos_theta_out;

  /*! @brief Sines of the output zenith angles. */
  std::vector<float_type> _sin_theta_out;

  /*! @brief Inverse sines of the output zenith angles. */
  std::vector<float_type> _sin_theta_out_inverse;

  /*! @brief Cosines of the relative azimuth angles. */
  std::vector<float_type> _cos_phi;

  /*! @brief Sines of the relative azimuth angles. */
  std::vector<float_type> _sin_phi;

public:
  /**
   * @brief Constructor.
   *
   * @param ntheta_in Number of input zenith angles.
   * @param theta_in Grid of input zenith angles (in radians, of size ntheta_in
   * or more).
   * @param ntheta_out Number of output zenith angles.
   * @param theta_out Grid of input zenith angles (in radians, of size
   * ntheta_out or more).
   * @param nphi Number of relative azimuth angles.
   * @param phi Grid of relative azimuth angles (in radians, of size nphi or
   * more).
   */
  inline FullScatteringMatrixGrid(const uint_fast32_t ntheta_in,
                                  const float_type *theta_in,
                                  const uint_fast32_t ntheta_out,
                                  const float_type *theta_out,
                                  const uint_fast32_t nphi,
                                  const float_type *phi)
      : _theta_in(ntheta_in), _theta_out(ntheta_out), _phi(nphi),
        _cos_theta_in(ntheta_in), _sin_theta_in(ntheta_in),
        _sin_theta_in_inverse(ntheta_in), _cos_theta_out(ntheta_out),
        _sin_theta_out(ntheta_out), _sin_theta_out_inverse(ntheta_out),
        _cos_phi(nphi), _sin_phi(nphi) {

    for (uint_fast32_t i = 0; i < ntheta_in; ++i) {
      _theta_in[i] = theta_in[i];
    }
    for (uint_fast32_t i = 0; i < ntheta_out; ++i) {
      _theta_out[i] = theta_out[i];
    }
    for (uint_fast32_t i = 0; i < nphi; ++i) {
      _phi[i] = phi[i];
    }
  }

  virtual ~FullScatteringMatrixGrid() {}

  /**
   * @brief Get the size in memory of a hypothetical FullScatteringMatrixGrid
   * object with the given parameters.
   *
   * @param ntheta_in Number of input zenith angles.
   * @param theta_in Grid of input zenith angles (in radians, of size ntheta_in
   * or more).
   * @param ntheta_out Number of output zenith angles.
   * @param theta_out Grid of input zenith angles (in radians, of size
   * ntheta_out or more).
   * @param nphi Number of relative azimuth angles.
   * @param phi Grid of relative azimuth angles (in radians, of size nphi or
   * more).
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t
  get_memory_size(const uint_fast32_t ntheta_in, const float_type *theta_in,
                  const uint_fast32_t ntheta_out, const float_type *theta_out,
                  const uint_fast32_t nphi, const float_type *phi) {
    size_t size = sizeof(FullScatteringMatrixGrid);
    size += 4 * ntheta_in * sizeof(float_type);
    size += 4 * ntheta_out * sizeof(float_type);
    size += 3 * nphi * sizeof(float_type);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, *this, true);
  }

  /**
   * @brief Get the number of read/write resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readwrite_resources() { return 1; }

  /**
   * @brief Get the number of read only resources for this task.
   *
   * @return 0.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 0; }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t ntheta_in = _cos_theta_in.size();
    for (uint_fast32_t itheta = 0; itheta < ntheta_in; ++itheta) {
      ctm_assert(_theta_in[itheta] >= 0.);
      ctm_assert(_theta_in[itheta] <= M_PI);
      _cos_theta_in[itheta] = cos(_theta_in[itheta]);
      // note that we only use the positive root since the sin(theta) values
      // are guaranteed to be in the range [0, 1]
      _sin_theta_in[itheta] =
          sqrt((1. - _cos_theta_in[itheta]) * (1. + _cos_theta_in[itheta]));
      if (_sin_theta_in[itheta] != 0.) {
        _sin_theta_in_inverse[itheta] = 1. / _sin_theta_in[itheta];
      } else {
        // this value should be ignored, we just set it to a safe and
        // recognisable value
        _sin_theta_in_inverse[itheta] = 9000.;
      }
    }

    const uint_fast32_t ntheta_out = _cos_theta_out.size();
    for (uint_fast32_t itheta = 0; itheta < ntheta_out; ++itheta) {
      ctm_assert(_theta_out[itheta] >= 0.);
      ctm_assert(_theta_out[itheta] <= M_PI);
      _cos_theta_out[itheta] = cos(_theta_out[itheta]);
      // note that we only use the positive root since the sin(theta) values
      // are guaranteed to be in the range [0, 1]
      _sin_theta_out[itheta] =
          sqrt((1. - _cos_theta_out[itheta]) * (1. + _cos_theta_out[itheta]));
      if (_sin_theta_out[itheta] != 0.) {
        _sin_theta_out_inverse[itheta] = 1. / _sin_theta_out[itheta];
      } else {
        // this value should be ignored, we just set it to a safe and
        // recognisable value
        _sin_theta_out_inverse[itheta] = 9000.;
      }
    }

    const uint_fast32_t nphi = _cos_phi.size();
    for (uint_fast32_t iphi = 0; iphi < nphi; ++iphi) {
      const float_type phi = _phi[iphi];
      _cos_phi[iphi] = cos(phi);
      _sin_phi[iphi] = sin(phi);
    }
  }
};

/**
 * @brief Precomputed special Wigner D functions that depend on a specific value
 * of @f$n_{max}@f$ and @f$n_{GL}@f$.
 */
class FullScatteringMatrixSpecialWignerDResources : public Resource,
                                                    public Task,
                                                    public Computable {
private:
  /*! @brief Maximum order, @f$n_{max}@f$. */
  const uint_fast32_t _nmax;

  /*! @brief Full scattering matrix grid to use. */
  const FullScatteringMatrixGrid &_grid;

  /*! @brief Wigner D functions divided by sine. */
  std::vector<Matrix<float_type>> _wigner_d_sinx[2];

  /*! @brief Derivatives of the Wigner D functions. */
  std::vector<Matrix<float_type>> _dwigner_d[2];

public:
  /**
   * @brief Constructor.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param grid Full scattering matrix grid.
   */
  inline FullScatteringMatrixSpecialWignerDResources(
      const uint_fast32_t nmax, const FullScatteringMatrixGrid &grid)
      : _nmax(nmax), _grid(grid) {

    const uint_fast32_t number_of_elements = (2 + nmax + 1) * nmax / 2;

    _wigner_d_sinx[0].reserve(number_of_elements);
    _wigner_d_sinx[1].reserve(number_of_elements);
    _dwigner_d[0].reserve(number_of_elements);
    _dwigner_d[1].reserve(number_of_elements);
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      for (uint_fast32_t m = 0; m < n + 1; ++m) {
        _wigner_d_sinx[0].push_back(
            Matrix<float_type>(grid._cos_theta_out.size(), n));
        _dwigner_d[0].push_back(
            Matrix<float_type>(grid._cos_theta_out.size(), n));
        _wigner_d_sinx[1].push_back(
            Matrix<float_type>(grid._cos_theta_in.size(), n));
        _dwigner_d[1].push_back(
            Matrix<float_type>(grid._cos_theta_in.size(), n));
      }
    }
    ctm_assert(_wigner_d_sinx[0].size() == number_of_elements);
  }

  virtual ~FullScatteringMatrixSpecialWignerDResources() {}

  /**
   * @brief Get the size in memory of a hypothetical SpecialWignerDResources
   * object with the given parameters.
   *
   * @param nmax Maximum order, @f$n_{max}@f$.
   * @param ntheta Number of zenith angles.
   * @return Size in bytes of the hypothetical object.
   */
  static inline size_t get_memory_size(const uint_fast32_t nmax,
                                       const uint_fast32_t ntheta) {
    size_t size = sizeof(FullScatteringMatrixSpecialWignerDResources);
    for (uint_fast32_t n = 1; n < nmax + 1; ++n) {
      for (uint_fast32_t m = 0; m < n + 1; ++m) {
        size += 2 * ntheta * n * sizeof(float_type);
        size += 2 * ntheta * n * sizeof(float_type);
      }
    }
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, *this, true);

    // read access
    quicksched.link_task_and_resource(*this, _grid, false);
  }

  /**
   * @brief Get the number of read/write resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readwrite_resources() { return 1; }

  /**
   * @brief Get the number of read only resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 1; }

  /**
   * @brief Compute the factors.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id = 0) {

    const uint_fast32_t ntheta_in = _grid._cos_theta_out.size();
    const uint_fast32_t ntheta_out = _grid._cos_theta_in.size();
    for (uint_fast32_t n = 1; n < _nmax + 1; ++n) {
      for (uint_fast32_t m = 0; m < n + 1; ++m) {
        for (uint_fast32_t itheta_in = 0; itheta_in < ntheta_in; ++itheta_in) {
          const uint_fast32_t index = (2 + n) * (n - 1) / 2 + m;
          ctm_assert(index < _wigner_d_sinx[0].size());
          ctm_assert_message(
              _wigner_d_sinx[0][index].get_number_of_columns() == n,
              "cols: %" PRIuFAST32 ", n: %" PRIuFAST32 ", m: %" PRIuFAST32
              ", index: %" PRIuFAST32,
              _wigner_d_sinx[0][index].get_number_of_columns(), n, m, index);
          SpecialFunctions::wigner_dn_0m_sinx(
              _grid._cos_theta_out[itheta_in], _grid._sin_theta_out[itheta_in],
              _grid._sin_theta_out_inverse[itheta_in], n, m,
              &_wigner_d_sinx[0][index].get_row(itheta_in)[0],
              &_dwigner_d[0][index].get_row(itheta_in)[0]);
        }
        for (uint_fast32_t itheta_out = 0; itheta_out < ntheta_out;
             ++itheta_out) {
          const uint_fast32_t index = (2 + n) * (n - 1) / 2 + m;
          ctm_assert(index < _wigner_d_sinx[1].size());
          ctm_assert_message(
              _wigner_d_sinx[1][index].get_number_of_columns() == n,
              "cols: %" PRIuFAST32 ", n: %" PRIuFAST32 ", m: %" PRIuFAST32
              ", index: %" PRIuFAST32,
              _wigner_d_sinx[1][index].get_number_of_columns(), n, m, index);
          SpecialFunctions::wigner_dn_0m_sinx(
              _grid._cos_theta_in[itheta_out], _grid._sin_theta_in[itheta_out],
              _grid._sin_theta_in_inverse[itheta_out], n, m,
              &_wigner_d_sinx[1][index].get_row(itheta_out)[0],
              &_dwigner_d[1][index].get_row(itheta_out)[0]);
        }
      }
    }
    make_available();
  }

  /**
   * @brief Get the special Wigner D function for the given input angle.
   *
   * @param igrid Internal grid to sample.
   * @param m @f$m@f$ value.
   * @param itheta_in Index of the input angle.
   * @param n Order, @f$n@f$.
   * @param nmax Maximum order.
   * @return Corresponding special Wigner D function value.
   */
  inline float_type get_wigner_d_sinx(const uint_fast8_t igrid,
                                      const uint_fast32_t m,
                                      const uint_fast32_t itheta_in,
                                      const uint_fast32_t n,
                                      const uint_fast32_t nmax) const {

    ctm_assert(nmax > 0);
    const uint_fast32_t index = (2 + nmax) * (nmax - 1) / 2 + m;
    ctm_assert(index < _wigner_d_sinx[igrid].size());
    ctm_assert(n > 0);
    ctm_assert(m <= n);
    ctm_assert(itheta_in < _wigner_d_sinx[igrid][index].get_number_of_rows());
    ctm_assert_message(n - 1 <
                           _wigner_d_sinx[igrid][index].get_number_of_columns(),
                       "m: %" PRIuFAST32 ", n: %" PRIuFAST32
                       ", nmax: %" PRIuFAST32 ", index: %" PRIuFAST32,
                       m, n, nmax, index);
    // check that the resource was actually computed
    check_use();
    return _wigner_d_sinx[igrid][index](itheta_in, n - 1);
  }

  /**
   * @brief Get the derivative of the Wigner D function for the given input
   * angle.
   *
   * @param igrid Internal grid to sample.
   * @param m @f$m@f$ value.
   * @param itheta_in Index of the input angle.
   * @param n Order, @f$n@f$.
   * @param nmax Maximum order.
   * @return Corresponding derivative value.
   */
  inline float_type get_dwigner_d(const uint_fast8_t igrid,
                                  const uint_fast32_t m,
                                  const uint_fast32_t itheta_in,
                                  const uint_fast32_t n,
                                  const uint_fast32_t nmax) const {

    ctm_assert(nmax > 0);
    const uint_fast32_t index = (2 + nmax) * (nmax - 1) / 2 + m;
    ctm_assert(index < _dwigner_d[igrid].size());
    ctm_assert(n > 0);
    ctm_assert(m <= n);
    ctm_assert(itheta_in < _dwigner_d[igrid][index].get_number_of_rows());
    ctm_assert(n - 1 < _dwigner_d[igrid][index].get_number_of_columns());
    // check that the resource was actually computed
    check_use();
    return _dwigner_d[igrid][index](itheta_in, n - 1);
  }
};

/**
 * @brief Task that computes AbsorptionCoefficients.
 */
class FullScatteringMatrixTask : public Task {
private:
  /*! @brief Angular grid to use. */
  const FullScatteringMatrixGrid &_grid;

  /*! @brief Interaction variables. */
  const InteractionVariables &_interaction_variables;

  /*! @brief T-matrix to use (read only). */
  const TMatrixResource &_Tmatrix;

  /*! @brief N based resources to use (read only). */
  const NBasedResources &_nfactors;

  /*! @brief Special Wigner D resources to use (read only). */
  const FullScatteringMatrixSpecialWignerDResources &_wigner_d;

  /*! @brief Resource in which the result is stored. */
  FullScatteringMatrixResult &_result;

public:
  /**
   * @brief Constructor.
   *
   * @param grid Angular grid to use.
   * @param interaction_variables Interaction variables.
   * @param Tmatrix T-matrix to use.
   * @param nfactors N based resources to use.
   * @param wigner_d Special Wigner D resources to use.
   * @param result Resource in which the result is stored.
   */
  inline FullScatteringMatrixTask(
      const FullScatteringMatrixGrid &grid,
      const InteractionVariables &interaction_variables,
      const TMatrixResource &Tmatrix, const NBasedResources &nfactors,
      const FullScatteringMatrixSpecialWignerDResources &wigner_d,
      FullScatteringMatrixResult &result)
      : _grid(grid), _interaction_variables(interaction_variables),
        _Tmatrix(Tmatrix), _nfactors(nfactors), _wigner_d(wigner_d),
        _result(result) {}

  virtual ~FullScatteringMatrixTask() {}

  /**
   * @brief Get the size in memory of a hypothetical FullScatteringMatrixTask
   * object with the given parameters.
   *
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t get_memory_size() {
    size_t size = sizeof(FullScatteringMatrixTask);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _result, true);

    // read access
    quicksched.link_task_and_resource(*this, _interaction_variables, false);
    quicksched.link_task_and_resource(*this, _Tmatrix, false);
    quicksched.link_task_and_resource(*this, _nfactors, false);
    quicksched.link_task_and_resource(*this, _grid, false);
    quicksched.link_task_and_resource(*this, _wigner_d, false);
  }

  /**
   * @brief Get the number of read/write resources for this task.
   *
   * @return 1.
   */
  inline static uint_fast32_t number_of_readwrite_resources() { return 1; }

  /**
   * @brief Get the number of read only resources for this task.
   *
   * @return 5.
   */
  inline static uint_fast32_t number_of_readonly_resources() { return 5; }

  /**
   * @brief Get the forward scattering matrix @f$S@f$ for a scattering event
   * from the given input angles to the given output angles at a particle with
   * its symmetry axis fixed to the @f$z@f$-axis of the reference frame.
   *
   * @param grid_out Internal output grid to sample.
   * @param itheta_out Index of the output zenith angle.
   * @param grid_in Internal input grid to sample.
   * @param itheta_in Index of the input zenith angle.
   * @param cosphi_in Cosine of the input azimuth angle, @f$\cos(\phi{}_i)@f$.
   * @param sinphi_in Sine of the input azimuth angle, @f$\sin(\phi{}_i)@f$.
   * @param S Scattering matrix for this scattering event (of size at least 4).
   */
  inline void get_forward_scattering_matrix(const uint_fast8_t grid_in,
                                            const uint_fast32_t itheta_in,
                                            const uint_fast8_t grid_out,
                                            const uint_fast32_t itheta_out,
                                            const float_type cosphi_in,
                                            const float_type sinphi_in,
                                            std::complex<float_type> *S) const {

    // initialize the scattering matrix
    S[0] = 0.;
    S[1] = 0.;
    S[2] = 0.;
    S[3] = 0.;

    const uint_fast32_t nmax = _Tmatrix.get_nmax();

    // now compute the matrix S^P
    // we precompute e^{i(phi_out-phi_in)}
    const std::complex<float_type> expiphi_p_out_m_in(cosphi_in, -sinphi_in);
    // e^{im(phi_out-phi_in)} is computed recursively, starting with the value
    // for m=0: 1
    std::complex<float_type> expimphi_p_out_m_in(1., 0.);
    // instead of summing over n and n', we sum over m, since then we can reuse
    // the e^{im(phi_out-phi_in)}, pi and tau factors
    for (uint_fast32_t m = 0; m < nmax + 1; ++m) {
      // only n and n' values larger than or equal to m have non-trivial
      // contributions to the S matrix
      const uint_fast32_t nmin = std::max(m, static_cast<uint_fast32_t>(1));

      // we get the real and imaginary part of e^{im\phi{}} and multiply with
      // 2 to account for both m and -m
      const float_type fcos = 2. * expimphi_p_out_m_in.real();
      const float_type fsin = 2. * expimphi_p_out_m_in.imag();
      // recurse the exponential for the next iteration
      expimphi_p_out_m_in *= expiphi_p_out_m_in;

      // now perform the actual sums over n and n'
      for (uint_fast32_t nn = nmin; nn < nmax + 1; ++nn) {

        // get the specific pi and tau for this n'
        const float_type pi_nn =
            m * _wigner_d.get_wigner_d_sinx(grid_in, m, itheta_in, nn, nmax);
        const float_type tau_nn =
            _wigner_d.get_dwigner_d(grid_in, m, itheta_in, nn, nmax);

        for (uint_fast32_t n = nmin; n < nmax + 1; ++n) {

          // get the specific pi and tau for this n
          const float_type pi_n =
              m * _wigner_d.get_wigner_d_sinx(grid_out, m, itheta_out, n, nmax);
          const float_type tau_n =
              _wigner_d.get_dwigner_d(grid_out, m, itheta_out, n, nmax);

          // get the c factor for these values of n and n'
          const std::complex<float_type> c_nnn = _nfactors.get_cnn(nn, n);

          // get the T11 and T22 elements for this m, n and n' (we need these
          // in all cases)
          const std::complex<float_type> T11nmnnm = _Tmatrix(0, n, m, 0, nn, m);
          const std::complex<float_type> T22nmnnm = _Tmatrix(1, n, m, 1, nn, m);
          // if m=0, the T12 and T21 matrices are trivially zero, and we can
          // simplify the expression for S
          if (m == 0) {
            const std::complex<float_type> factor = c_nnn * tau_n * tau_nn;
            S[0] += factor * T22nmnnm;
            S[3] += factor * T11nmnnm;
          } else {
            // in the general case m=/=0, we also need the T12 and T21 elements
            // for this m, n and n'
            const std::complex<float_type> T12nmnnm =
                _Tmatrix(0, n, m, 1, nn, m);
            const std::complex<float_type> T21nmnnm =
                _Tmatrix(1, n, m, 0, nn, m);

            // due to m symmetry, S11 and S22 only have the cosine factor,
            // while S12 and S21 only have the sine factor
            const std::complex<float_type> real_factor = c_nnn * fcos;
            const std::complex<float_type> imag_factor = c_nnn * fsin;

            // precompute the pi and tau factor combinations
            const float_type pi_pi = pi_n * pi_nn;
            const float_type pi_tau = pi_n * tau_nn;
            const float_type tau_pi = tau_n * pi_nn;
            const float_type tau_tau = tau_n * tau_nn;

            S[0] += real_factor * (T11nmnnm * pi_pi + T21nmnnm * tau_pi +
                                   T12nmnnm * pi_tau + T22nmnnm * tau_tau);
            S[1] += imag_factor * (T11nmnnm * pi_tau + T21nmnnm * tau_tau +
                                   T12nmnnm * pi_pi + T22nmnnm * tau_pi);
            S[2] -= imag_factor * (T11nmnnm * tau_pi + T21nmnnm * pi_pi +
                                   T12nmnnm * tau_tau + T22nmnnm * pi_tau);
            S[3] += real_factor * (T11nmnnm * tau_tau + T21nmnnm * pi_tau +
                                   T12nmnnm * tau_pi + T22nmnnm * pi_pi);
          }
        }
      }
    }
    // now divide all expressions by the wavenumber
    const float_type kinv = 1. / _interaction_variables.get_wavenumber();
    S[0] *= kinv;
    S[1] *= kinv;
    S[2] *= kinv;
    S[3] *= kinv;
  }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t ntheta_in = _grid._theta_in.size();
    const uint_fast32_t ntheta_out = _grid._theta_out.size();
    const uint_fast32_t nphi = _grid._phi.size();

    const std::complex<float_type> icompl(0., 1.);
    const float_type half(0.5);
    for (uint_fast32_t itheta_out = 0; itheta_out < ntheta_out; ++itheta_out) {

      for (uint_fast32_t itheta_in = 0; itheta_in < ntheta_in; ++itheta_in) {
        for (uint_fast32_t iphi = 0; iphi < nphi; ++iphi) {

          const float_type cos_phi = _grid._cos_phi[iphi];
          const float_type sin_phi = _grid._sin_phi[iphi];

          std::complex<float_type> S[4];
          get_forward_scattering_matrix(1, itheta_in, 0, itheta_out, cos_phi,
                                        sin_phi, S);

          float_type *Z = &_result._Z[4 * itheta_out * itheta_in * iphi];
          Z[0] = (half * (S[0] * conj(S[0]) + S[1] * conj(S[1]) +
                          S[2] * conj(S[2]) + S[3] * conj(S[3])))
                     .real();
          Z[1] = (half * (S[0] * conj(S[0]) - S[1] * conj(S[1]) +
                          S[2] * conj(S[2]) - S[3] * conj(S[3])))
                     .real();
          Z[2] = (S[0] * conj(S[3]) + S[1] * conj(S[2])).real();
          Z[3] = (-icompl * (S[0] * conj(S[3]) + S[2] * conj(S[1]))).real();
        }
      }
    }
  }
};

#endif // FULLSCATTERINGMATRIXTASK_HPP
