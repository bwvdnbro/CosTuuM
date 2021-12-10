/*******************************************************************************
 * This file is part of CosTuuM
 * Copyright (C) 2019, 2020 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
 * @file ShapeAveragingTask.hpp
 *
 * @brief Task that computes the shape distribution average of the absorption
 * coefficients.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */
#ifndef SHAPEAVERAGINGTASK_HPP
#define SHAPEAVERAGINGTASK_HPP

#include "AbsorptionCoefficientTask.hpp"
#include "Configuration.hpp"
#include "ExtinctionCoefficientTask.hpp"
#include "FullScatteringMatrixTask.hpp"
#include "QuickSchedWrapper.hpp"
#include "ScatteringMatrixTask.hpp"
#include "ShapeDistribution.hpp"

/**
 * @brief Task that computes the shape distribution average of the absorption
 * coefficients.
 */
class AbsorptionShapeAveragingTask : public Task {
private:
  /*! @brief Shape distribution. */
  const ShapeDistribution &_shape_distribution;

  /*! @brief Input absorption coefficients. */
  std::vector<AbsorptionCoefficientResult *> _input_coefficients;

  /*! @brief Output (averaged) absorption coefficients. */
  AbsorptionCoefficientResult &_output_coefficients;

public:
  /**
   * @brief Constructor.
   *
   * @param shape_distribution Shape distribution.
   * @param output_coefficients Output coefficients.
   */
  inline AbsorptionShapeAveragingTask(
      const ShapeDistribution &shape_distribution,
      AbsorptionCoefficientResult &output_coefficients)
      : _shape_distribution(shape_distribution),
        _input_coefficients(shape_distribution.get_number_of_points(), nullptr),
        _output_coefficients(output_coefficients) {}

  virtual ~AbsorptionShapeAveragingTask() {}

  /**
   * @brief Get the size in memory of a hypothetical
   * AbsorptionShapeAveragingTask object with the given parameters.
   *
   * @param shape_distribution Shape distribution.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t
  get_memory_size(const ShapeDistribution &shape_distribution) {
    size_t size = sizeof(AbsorptionShapeAveragingTask);
    size += shape_distribution.get_number_of_points() *
            sizeof(AbsorptionCoefficientResult *);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _output_coefficients, true);
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
   * @brief Add input coefficients for the calculation.
   *
   * @param quicksched QuickSched library.
   * @param ishape Index of the shape for which the coefficients are computed.
   * @param input_coefficient Input coefficients.
   */
  inline void
  add_input_coefficient(QuickSched &quicksched, const uint_fast32_t ishape,
                        AbsorptionCoefficientResult *input_coefficient) {
    ctm_assert(ishape < _input_coefficients.size());
    _input_coefficients[ishape] = input_coefficient;
    quicksched.link_task_and_resource(*this, *input_coefficient, false);
  }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t ntheta = _output_coefficients._Qabs.size();
    const uint_fast32_t nshape = _shape_distribution.get_number_of_points();

    // make sure the average values are set to 0
    for (uint_fast32_t itheta = 0; itheta < ntheta; ++itheta) {
      _output_coefficients._Qabs[itheta] = 0.;
      _output_coefficients._Qabspol[itheta] = 0.;
    }

    // compute the nominator and denominator in the expression for the
    // average
    float_type norm = 0.;
    for (uint_fast32_t ishape = 0; ishape < nshape; ++ishape) {

      // some sanity checks
      ctm_assert(_input_coefficients[ishape] != nullptr);
      ctm_assert(_input_coefficients[ishape]->_Qabs.size() == ntheta);

      // get the value of the shape distribution at this evaluation point
      const float_type weight = _shape_distribution.get_weight(ishape);
      // add it to the norm (denominator in expression)
      norm += weight;
      // add the contributions to the absorption coefficients from this shape
      for (uint_fast32_t itheta = 0; itheta < ntheta; ++itheta) {
        _output_coefficients._Qabs[itheta] +=
            weight * _input_coefficients[ishape]->_Qabs[itheta];
        _output_coefficients._Qabspol[itheta] +=
            weight * _input_coefficients[ishape]->_Qabspol[itheta];
      }
    }

    // normalise the average quantities
    const float_type norm_inv = 1. / norm;
    for (uint_fast32_t itheta = 0; itheta < ntheta; ++itheta) {
      _output_coefficients._Qabs[itheta] *= norm_inv;
      _output_coefficients._Qabspol[itheta] *= norm_inv;
    }
  }
};

/**
 * @brief Task that computes the shape distribution average of the extinction
 * coefficients.
 */
class ExtinctionShapeAveragingTask : public Task {
private:
  /*! @brief Shape distribution. */
  const ShapeDistribution &_shape_distribution;

  /*! @brief Input extinction coefficients. */
  std::vector<ExtinctionCoefficientResult *> _input_coefficients;

  /*! @brief Output (averaged) extinction coefficients. */
  ExtinctionCoefficientResult &_output_coefficients;

public:
  /**
   * @brief Constructor.
   *
   * @param shape_distribution Shape distribution.
   * @param output_coefficients Output coefficients.
   */
  inline ExtinctionShapeAveragingTask(
      const ShapeDistribution &shape_distribution,
      ExtinctionCoefficientResult &output_coefficients)
      : _shape_distribution(shape_distribution),
        _input_coefficients(shape_distribution.get_number_of_points(), nullptr),
        _output_coefficients(output_coefficients) {}

  virtual ~ExtinctionShapeAveragingTask() {}

  /**
   * @brief Get the size in memory of a hypothetical
   * ExtinctionShapeAveragingTask object with the given parameters.
   *
   * @param shape_distribution Shape distribution.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t
  get_memory_size(const ShapeDistribution &shape_distribution) {
    size_t size = sizeof(ExtinctionShapeAveragingTask);
    size += shape_distribution.get_number_of_points() *
            sizeof(ExtinctionCoefficientResult *);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _output_coefficients, true);
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
   * @brief Add input coefficients for the calculation.
   *
   * @param quicksched QuickSched library.
   * @param ishape Index of the shape for which the coefficients are computed.
   * @param input_coefficient Input coefficients.
   */
  inline void
  add_input_coefficient(QuickSched &quicksched, const uint_fast32_t ishape,
                        ExtinctionCoefficientResult *input_coefficient) {
    ctm_assert(ishape < _input_coefficients.size());
    _input_coefficients[ishape] = input_coefficient;
    quicksched.link_task_and_resource(*this, *input_coefficient, false);
  }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t ntheta = _output_coefficients._Qext.size();
    const uint_fast32_t nshape = _shape_distribution.get_number_of_points();

    // make sure the average values are set to 0
    for (uint_fast32_t itheta = 0; itheta < ntheta; ++itheta) {
      _output_coefficients._Qext[itheta] = 0.;
      _output_coefficients._Qextpol[itheta] = 0.;
      _output_coefficients._Qextcpol[itheta] = 0.;
    }

    // compute the nominator and denominator in the expression for the
    // average
    float_type norm = 0.;
    for (uint_fast32_t ishape = 0; ishape < nshape; ++ishape) {

      // some sanity checks
      ctm_assert(_input_coefficients[ishape] != nullptr);
      ctm_assert(_input_coefficients[ishape]->_Qext.size() == ntheta);

      // get the value of the shape distribution at this evaluation point
      const float_type weight = _shape_distribution.get_weight(ishape);
      // add it to the norm (denominator in expression)
      norm += weight;
      // add the contributions to the absorption coefficients from this shape
      for (uint_fast32_t itheta = 0; itheta < ntheta; ++itheta) {
        _output_coefficients._Qext[itheta] +=
            weight * _input_coefficients[ishape]->_Qext[itheta];
        _output_coefficients._Qextpol[itheta] +=
            weight * _input_coefficients[ishape]->_Qextpol[itheta];
        _output_coefficients._Qextcpol[itheta] +=
            weight * _input_coefficients[ishape]->_Qextcpol[itheta];
      }
    }

    // normalise the average quantities
    const float_type norm_inv = 1. / norm;
    for (uint_fast32_t itheta = 0; itheta < ntheta; ++itheta) {
      _output_coefficients._Qext[itheta] *= norm_inv;
      _output_coefficients._Qextpol[itheta] *= norm_inv;
      _output_coefficients._Qextcpol[itheta] *= norm_inv;
    }
  }
};

/**
 * @brief Task that computes the shape distribution average of the scattering
 * matrix.
 */
class ScatteringMatrixShapeAveragingTask : public Task {
private:
  /*! @brief Shape distribution. */
  const ShapeDistribution &_shape_distribution;

  /*! @brief Input scattering matrix. */
  std::vector<ScatteringMatrixResult *> _input_matrices;

  /*! @brief Output (averaged) extinction coefficients. */
  ScatteringMatrixResult &_output_matrices;

public:
  /**
   * @brief Constructor.
   *
   * @param shape_distribution Shape distribution.
   * @param output_matrices Output matrices.
   */
  inline ScatteringMatrixShapeAveragingTask(
      const ShapeDistribution &shape_distribution,
      ScatteringMatrixResult &output_matrices)
      : _shape_distribution(shape_distribution),
        _input_matrices(shape_distribution.get_number_of_points(), nullptr),
        _output_matrices(output_matrices) {}

  virtual ~ScatteringMatrixShapeAveragingTask() {}

  /**
   * @brief Get the size in memory of a hypothetical
   * ScatteringMatrixShapeAveragingTask object with the given parameters.
   *
   * @param shape_distribution Shape distribution.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t
  get_memory_size(const ShapeDistribution &shape_distribution) {
    size_t size = sizeof(ScatteringMatrixShapeAveragingTask);
    size += shape_distribution.get_number_of_points() *
            sizeof(ScatteringMatrixResult *);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _output_matrices, true);
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
   * @brief Add input matrices for the calculation.
   *
   * @param quicksched QuickSched library.
   * @param ishape Index of the shape for which the matrices are computed.
   * @param input_matrices Input coefficients.
   */
  inline void add_input_matrices(QuickSched &quicksched,
                                 const uint_fast32_t ishape,
                                 ScatteringMatrixResult *input_matrices) {
    ctm_assert(ishape < _input_matrices.size());
    _input_matrices[ishape] = input_matrices;
    quicksched.link_task_and_resource(*this, *input_matrices, false);
  }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t nmat = _output_matrices._scattering_matrix.size();
    const uint_fast32_t nshape = _shape_distribution.get_number_of_points();

    // make sure the average values are set to 0
    for (uint_fast32_t imat = 0; imat < nmat; ++imat) {
      _output_matrices._scattering_matrix[imat].reset();
    }

    // compute the nominator and denominator in the expression for the
    // average
    float_type norm = 0.;
    for (uint_fast32_t ishape = 0; ishape < nshape; ++ishape) {

      // some sanity checks
      ctm_assert(_input_matrices[ishape] != nullptr);
      ctm_assert(_input_matrices[ishape]->_scattering_matrix.size() == nmat);

      // get the value of the shape distribution at this evaluation point
      const float_type weight = _shape_distribution.get_weight(ishape);
      // add it to the norm (denominator in expression)
      norm += weight;
      // add the contributions from this shape
      for (uint_fast32_t imat = 0; imat < nmat; ++imat) {
        for (uint_fast8_t i = 0; i < 4; ++i) {
          for (uint_fast8_t j = 0; j < 4; ++j) {
            _output_matrices._scattering_matrix[imat](i, j) +=
                weight *
                _input_matrices[ishape]->_scattering_matrix[imat](i, j);
          }
        }
      }
    }

    // normalise the average quantities
    const float_type norm_inv = 1. / norm;
    for (uint_fast32_t imat = 0; imat < nmat; ++imat) {
      for (uint_fast8_t i = 0; i < 4; ++i) {
        for (uint_fast8_t j = 0; j < 4; ++j) {
          _output_matrices._scattering_matrix[imat](i, j) *= norm_inv;
        }
      }
    }
  }
};

/**
 * @brief Task that computes the shape distribution average of a full scattering
 * matrix.
 */
class FullScatteringMatrixShapeAveragingTask : public Task {
private:
  /*! @brief Shape distribution. */
  const ShapeDistribution &_shape_distribution;

  /*! @brief Input scattering matrices. */
  std::vector<FullScatteringMatrixResult *> _input_matrices;

  /*! @brief Output (averaged) scattering matrices. */
  FullScatteringMatrixResult &_output_matrices;

public:
  /**
   * @brief Constructor.
   *
   * @param shape_distribution Shape distribution.
   * @param output_matrices Output matrices.
   */
  inline FullScatteringMatrixShapeAveragingTask(
      const ShapeDistribution &shape_distribution,
      FullScatteringMatrixResult &output_matrices)
      : _shape_distribution(shape_distribution),
        _input_matrices(shape_distribution.get_number_of_points(), nullptr),
        _output_matrices(output_matrices) {}

  virtual ~FullScatteringMatrixShapeAveragingTask() {}

  /**
   * @brief Get the size in memory of a hypothetical
   * FullScatteringMatrixShapeAveragingTask object with the given parameters.
   *
   * @param shape_distribution Shape distribution.
   * @return Size in bytes that the object would occupy.
   */
  static inline size_t
  get_memory_size(const ShapeDistribution &shape_distribution) {
    size_t size = sizeof(FullScatteringMatrixShapeAveragingTask);
    size += shape_distribution.get_number_of_points() *
            sizeof(FullScatteringMatrixResult *);
    return size;
  }

  /**
   * @brief Link the resources for this task.
   *
   * @param quicksched QuickSched library.
   */
  inline void link_resources(QuickSched &quicksched) {
    // write access
    quicksched.link_task_and_resource(*this, _output_matrices, true);
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
   * @brief Add input matrices for the calculation.
   *
   * @param quicksched QuickSched library.
   * @param ishape Index of the shape for which the matrices are computed.
   * @param input_matrices Input coefficients.
   */
  inline void add_input_matrices(QuickSched &quicksched,
                                 const uint_fast32_t ishape,
                                 FullScatteringMatrixResult *input_matrices) {
    ctm_assert(ishape < _input_matrices.size());
    _input_matrices[ishape] = input_matrices;
    quicksched.link_task_and_resource(*this, *input_matrices, false);
  }

  /**
   * @brief Execute the task.
   *
   * @param thread_id ID of the thread that executes the task.
   */
  virtual void execute(const int_fast32_t thread_id) {

    const uint_fast32_t nmat = _output_matrices._Z.size() / 4;
    const uint_fast32_t nshape = _shape_distribution.get_number_of_points();

    // make sure the average values are set to 0
    for (uint_fast32_t imat = 0; imat < nmat; ++imat) {
      _output_matrices._Z[4 * imat] = 0.;
      _output_matrices._Z[4 * imat + 1] = 0.;
      _output_matrices._Z[4 * imat + 2] = 0.;
      _output_matrices._Z[4 * imat + 3] = 0.;
    }

    // compute the nominator and denominator in the expression for the
    // average
    float_type norm = 0.;
    for (uint_fast32_t ishape = 0; ishape < nshape; ++ishape) {

      // some sanity checks
      ctm_assert(_input_matrices[ishape] != nullptr);
      ctm_assert(_input_matrices[ishape]->_Z.size() / 4 == nmat);

      // get the value of the shape distribution at this evaluation point
      const float_type weight = _shape_distribution.get_weight(ishape);
      // add it to the norm (denominator in expression)
      norm += weight;
      // add the contributions from this shape
      for (uint_fast32_t imat = 0; imat < nmat; ++imat) {
        for (uint_fast8_t i = 0; i < 4; ++i) {
          _output_matrices._Z[4 * imat + i] +=
              weight * _input_matrices[ishape]->_Z[4 * imat + i];
        }
      }
    }

    // normalise the average quantities
    const float_type norm_inv = 1. / norm;
    for (uint_fast32_t imat = 0; imat < nmat; ++imat) {
      for (uint_fast8_t i = 0; i < 4; ++i) {
        _output_matrices._Z[4 * imat + i] *= norm_inv;
      }
    }
  }
};

#endif // SHAPEAVERAGINGTASK_HPP
