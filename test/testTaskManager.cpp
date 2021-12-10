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
 * @file testTaskManager.cpp
 *
 * @brief Unit test for the TaskManager class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Assert.hpp"
#include "DraineDustProperties.hpp"
#include "SizeBasedAlignmentDistribution.hpp"
#include "TaskManager.hpp"

#include <fstream>
#include <string.h>

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Delete all pointers in the given pointer vector.
 *
 * @param vec Pointer vector.
 */
template <typename T> inline void clear_vector(std::vector<T *> &vec) {
  for (uint_fast32_t i = 0; i < vec.size(); ++i) {
    // make sure we used all elements in the vectors
    ctm_assert(vec[i] != nullptr);
    delete vec[i];
  }
}

/**
 * @brief Print all tasks in the given vector to the given file.
 *
 * @param vec Vector containing task pointers.
 * @param quicksched QuickSched library wrapper.
 * @param ofile Output file.
 */
template <typename T>
inline void print_vector(std::vector<T *> &vec, QuickSched &quicksched,
                         std::ofstream &ofile) {
  for (uint_fast32_t i = 0; i < vec.size(); ++i) {
    quicksched.print_task(*vec[i], ofile);
  }
}

/**
 * @brief Write the given integer to the given binary file.
 *
 * @param ofile File to write to.
 * @param value Value to write.
 * @tparam T Integer type.
 */
template <typename T>
inline void int_to_file(std::ofstream &ofile, const T value) {
  const uint_fast64_t intval = value;
  ofile.write(reinterpret_cast<const char *>(&intval), sizeof(intval));
}

/**
 * @brief Write the given string to the given binary file.
 *
 * @param ofile File to write to.
 * @param value Value to write.
 */
inline void string_to_file(std::ofstream &ofile, const std::string value) {
  char strval[] = "        ";
  memcpy(strval, value.c_str(), value.size());
  ofile.write(strval, 8);
}

/**
 * @brief Write the given array to the given binary file.
 *
 * @param ofile File to write to.
 * @param value Value to write.
 */
inline void array_to_file(std::ofstream &ofile,
                          const std::vector<float_type> &value) {
  for (uint_fast32_t i = 0; i < value.size(); ++i) {
    const double ival = static_cast<double>(value[i]);
    ofile.write(reinterpret_cast<const char *>(&ival), sizeof(ival));
  }
}

/**
 * @brief Write the given float to the given binary file.
 *
 * @param ofile File to write to.
 * @param value Value to write.
 * @tparam T Floating point type.
 */
template <typename T>
inline void float_to_file(std::ofstream &ofile, const T &value) {
  const double fltval = value;
  ofile.write(reinterpret_cast<const char *>(&fltval), sizeof(fltval));
}

/**
 * @brief Unit test for the generate_tasks() function.
 */
inline int test_generate_tasks() {
  uint_fast32_t shape_distribution_type = 2;
  ShapeDistribution *shape_distribution;
  if (shape_distribution_type == 0) {
    shape_distribution = new ShapeDistribution();
    shape_distribution->evaluate(20u);
  } else if (shape_distribution_type == 1) {
    shape_distribution = new DraineHensleyShapeDistribution(20u);
  } else if (shape_distribution_type == 2) {
    shape_distribution = new SingleShapeShapeDistribution(1.00001);
  } else {
    ctm_error("Invalid shape distribution type!");
    shape_distribution = nullptr;
  }
  SizeBasedAlignmentDistribution alignment_distribution(1.e-5, 0, 100);
  const DraineDustProperties dust_properties;
  TaskManager task_manager(10, 100, 2, 1.e-4, 1e10, *shape_distribution,
                           alignment_distribution, dust_properties);

  const float_type log_min_size = -9.;
  const float_type log_max_size = -5.;
  const uint_fast32_t num_sizes = 3;
  std::vector<float_type> sizes(num_sizes);
  for (uint_fast32_t isize = 0; isize < num_sizes; ++isize) {
    sizes[isize] =
        pow(10., log_min_size +
                     isize * (log_max_size - log_min_size) / (num_sizes - 1.));
  }

  const float_type log_min_wavelength = -5.;
  const float_type log_max_wavelength = -3.;
  const uint_fast32_t num_wavelengths = 3;
  std::vector<float_type> wavelengths(num_wavelengths);
  for (uint_fast32_t ilambda = 0; ilambda < num_wavelengths; ++ilambda) {
    wavelengths[ilambda] =
        pow(10., log_min_wavelength +
                     ilambda * (log_max_wavelength - log_min_wavelength) /
                         (num_wavelengths - 1.));
  }

  task_manager.add_composition(DUSTGRAINTYPE_SILICON);
  for (uint_fast32_t i = 0; i < sizes.size(); ++i) {
    task_manager.add_size(sizes[i]);
  }
  for (uint_fast32_t i = 0; i < wavelengths.size(); ++i) {
    task_manager.add_wavelength(wavelengths[i]);
  }

  QuickSched quicksched(4, true, "test_TaskManager.log");

  const uint_fast32_t ntheta = 10;
  std::vector<float_type> thetas(ntheta);
  for (uint_fast32_t i = 0; i < ntheta; ++i) {
    thetas[i] = (i + 0.5) * M_PI / ntheta;
  }

  std::vector<Task *> tasks;
  std::vector<Resource *> resources;
  ResultKey *result_key = nullptr;
  std::vector<Result *> results;
  TMatrixAuxiliarySpaceManager *space_manager = nullptr;
  task_manager.generate_tasks(thetas, 20, quicksched, tasks, resources,
                              result_key, results, space_manager, false, true,
                              false, false);

  quicksched.execute_tasks();

  std::ofstream taskfile("test_taskmanager_tasks.txt");
  taskfile << "# thread\tstart\tend\ttype\ttask id\n";
  print_vector(tasks, quicksched, taskfile);
  std::ofstream typefile("test_taskmanager_types.txt");
  typefile << "# type\tlabel\n";
  quicksched.print_type_dict(typefile);

  // output some example results along each axis
  {
    std::ofstream ofile("test_taskmanager_size.txt");
    ofile << "# size\tQabs\tQabspol\n";
    for (uint_fast32_t isize = 0; isize < sizes.size(); ++isize) {
      const uint_fast32_t result_index = isize * wavelengths.size();
      const AbsorptionCoefficientResult &result =
          *static_cast<AbsorptionCoefficientResult *>(results[result_index]);
      ofile << sizes[isize] << "\t" << result.get_Qabs(0) << "\t"
            << result.get_Qabspol(0) << "\n";
    }
  }
  {
    std::ofstream ofile("test_taskmanager_wavelength.txt");
    ofile << "# wavelength\tQabs\tQabspol\n";
    for (uint_fast32_t ilambda = 0; ilambda < wavelengths.size(); ++ilambda) {
      const uint_fast32_t result_index = ilambda;
      const AbsorptionCoefficientResult &result =
          *static_cast<AbsorptionCoefficientResult *>(results[result_index]);
      ofile << wavelengths[ilambda] << "\t" << result.get_Qabs(0) << "\t"
            << result.get_Qabspol(0) << "\n";
    }
  }
  {
    std::ofstream ofile("test_taskmanager_theta.txt");
    ofile << "# theta\tQabs\tQabspol\n";
    const AbsorptionCoefficientResult &result =
        *static_cast<AbsorptionCoefficientResult *>(results[0]);
    for (uint_fast32_t itheta = 0; itheta < thetas.size(); ++itheta) {
      ofile << thetas[itheta] << "\t" << result.get_Qabs(itheta) << "\t"
            << result.get_Qabspol(itheta) << "\n";
    }
  }

  {
    std::ofstream ofile("test_table.stab", std::ios::binary);
    const std::string skirt_tag("SKIRT X\n");
    ofile.write(skirt_tag.c_str(), skirt_tag.size());
    int_to_file(ofile, 0x010203040A0BFEFF);

    int_to_file(ofile, 3);
    string_to_file(ofile, "a");
    string_to_file(ofile, "lambda");
    string_to_file(ofile, "theta");
    string_to_file(ofile, "m");
    string_to_file(ofile, "m");
    string_to_file(ofile, "rad");
    string_to_file(ofile, "log");
    string_to_file(ofile, "log");
    string_to_file(ofile, "lin");

    int_to_file(ofile, sizes.size());
    array_to_file(ofile, sizes);
    int_to_file(ofile, wavelengths.size());
    array_to_file(ofile, wavelengths);
    int_to_file(ofile, thetas.size());
    array_to_file(ofile, thetas);

    int_to_file(ofile, 2);
    string_to_file(ofile, "Qabs");
    string_to_file(ofile, "Qabspol");
    string_to_file(ofile, "1");
    string_to_file(ofile, "1");
    string_to_file(ofile, "log");
    string_to_file(ofile, "lin");
    for (uint_fast32_t itheta = 0; itheta < thetas.size(); ++itheta) {

      for (uint_fast32_t ilambda = 0; ilambda < result_key->wavelength_size();
           ++ilambda) {

        for (uint_fast32_t isize = 0; isize < result_key->size_size();
             ++isize) {

          const uint_fast32_t result_index =
              result_key->get_result_index(0, isize, ilambda);
          const AbsorptionCoefficientResult &result =
              *static_cast<AbsorptionCoefficientResult *>(
                  results[result_index]);

          ctm_warning("size: %g vs %g", double(result_key->get_size(isize)),
                      double(result.get_size()));
          ctm_warning("lambda: %g vs %g",
                      double(result_key->get_wavelength(isize)),
                      double(result.get_wavelength()));
          assert_condition(result_key->get_size(isize) == result.get_size());
          assert_condition(result_key->get_wavelength(ilambda) ==
                           result.get_wavelength());

          float_to_file(ofile, result.get_Qabs(itheta));
          float_to_file(ofile, result.get_Qabspol(itheta));
        }
      }
    }

    const std::string trail_tag("STABEND\n");
    ofile.write(trail_tag.c_str(), trail_tag.size());
  }

  clear_vector(tasks);
  clear_vector(resources);
  delete result_key;
  clear_vector(results);
  delete space_manager;

  delete shape_distribution;

  return 0;
}

/**
 * @brief Unit test for the generate_scattering_tasks() function.
 */
inline int test_generate_scattering_tasks() {
  uint_fast32_t shape_distribution_type = 2;
  ShapeDistribution *shape_distribution;
  if (shape_distribution_type == 0) {
    shape_distribution = new ShapeDistribution();
    shape_distribution->evaluate(20u);
  } else if (shape_distribution_type == 1) {
    shape_distribution = new DraineHensleyShapeDistribution(20u);
  } else if (shape_distribution_type == 2) {
    shape_distribution = new SingleShapeShapeDistribution(1.00001);
  } else {
    ctm_error("Invalid shape distribution type!");
    shape_distribution = nullptr;
  }
  SizeBasedAlignmentDistribution alignment_distribution(1.e-5, 0, 100);
  const DraineDustProperties dust_properties;
  TaskManager task_manager(10, 100, 2, 1.e-4, 1e10, *shape_distribution,
                           alignment_distribution, dust_properties);

  const float_type log_min_size = -9.;
  const float_type log_max_size = -5.;
  const uint_fast32_t num_sizes = 3;
  std::vector<float_type> sizes(num_sizes);
  for (uint_fast32_t isize = 0; isize < num_sizes; ++isize) {
    sizes[isize] =
        pow(10., log_min_size +
                     isize * (log_max_size - log_min_size) / (num_sizes - 1.));
  }

  const float_type log_min_wavelength = -5.;
  const float_type log_max_wavelength = -3.;
  const uint_fast32_t num_wavelengths = 3;
  std::vector<float_type> wavelengths(num_wavelengths);
  for (uint_fast32_t ilambda = 0; ilambda < num_wavelengths; ++ilambda) {
    wavelengths[ilambda] =
        pow(10., log_min_wavelength +
                     ilambda * (log_max_wavelength - log_min_wavelength) /
                         (num_wavelengths - 1.));
  }

  task_manager.add_composition(DUSTGRAINTYPE_SILICON);
  for (uint_fast32_t i = 0; i < sizes.size(); ++i) {
    task_manager.add_size(sizes[i]);
  }
  for (uint_fast32_t i = 0; i < wavelengths.size(); ++i) {
    task_manager.add_wavelength(wavelengths[i]);
  }

  QuickSched quicksched(4, true, "test_TaskManager.log");

  const uint_fast32_t ntheta = 10;
  std::vector<float_type> theta_in(ntheta);
  std::vector<float_type> theta_out(ntheta);
  std::vector<float_type> phi(ntheta);
  for (uint_fast32_t i = 0; i < ntheta; ++i) {
    theta_in[i] = (i + 0.5) * M_PI / ntheta;
    theta_out[i] = (i + 0.5) * M_PI / ntheta;
    phi[i] = (i + 0.5) * 2. * M_PI / ntheta;
  }

  std::vector<Task *> tasks;
  std::vector<Resource *> resources;
  ResultKey *result_key = nullptr;
  std::vector<Result *> results;
  TMatrixAuxiliarySpaceManager *space_manager = nullptr;
  task_manager.generate_scattering_tasks(theta_in, theta_out, phi, quicksched,
                                         tasks, resources, result_key, results,
                                         space_manager);

  quicksched.execute_tasks();

  std::ofstream taskfile("test_taskmanager_scattering_tasks.txt");
  taskfile << "# thread\tstart\tend\ttype\ttask id\n";
  print_vector(tasks, quicksched, taskfile);
  std::ofstream typefile("test_taskmanager_scattering_types.txt");
  typefile << "# type\tlabel\n";
  quicksched.print_type_dict(typefile);

  std::ofstream ofile("test_full_scattering_matrix.txt");
  ofile << "#lambda\tsize\ttheta_in\ttheta_out\tphi\tZ[0-16]\n";
  for (uint_fast32_t ilambda = 0; ilambda < result_key->wavelength_size();
       ++ilambda) {
    for (uint_fast32_t isize = 0; isize < result_key->size_size(); ++isize) {
      const uint_fast32_t result_index =
          result_key->get_result_index(0, isize, ilambda);
      const FullScatteringMatrixResult &result =
          *static_cast<FullScatteringMatrixResult *>(results[result_index]);
      for (uint_fast32_t itheta_in = 0; itheta_in < ntheta; ++itheta_in) {
        for (uint_fast32_t itheta_out = 0; itheta_out < ntheta; ++itheta_out) {
          for (uint_fast32_t iphi = 0; iphi < ntheta; ++iphi) {
            ofile << wavelengths[ilambda] << "\t" << sizes[isize] << "\t"
                  << theta_in[itheta_in] << "\t" << theta_out[itheta_out]
                  << "\t" << phi[iphi];
            float_type Z[16];
            result.get_scattering_matrix(itheta_in, itheta_out, iphi, Z);
            for (uint_fast8_t i = 0; i < 16; ++i) {
              ofile << "\t" << Z[i];
            }
            ofile << "\n";
          }
        }
      }
    }
  }
  clear_vector(tasks);
  clear_vector(resources);
  delete result_key;
  clear_vector(results);
  delete space_manager;

  delete shape_distribution;

  return 0;
}

/**
 * @brief Unit test for the TaskManager class.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  if (test_generate_tasks() != 0) {
    ctm_error("Error running test_generate_tasks()");
    return 1;
  }

  if (test_generate_scattering_tasks() != 0) {
    ctm_error("Error running test_generate_scattering_tasks()");
    return 1;
  }

  return 0;
}
