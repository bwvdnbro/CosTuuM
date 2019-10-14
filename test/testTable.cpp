/**
 * @file testTable.cpp
 *
 * @brief Unit test for the Table class.
 *
 * @author Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
 */

#include "Configuration.hpp"
#include "Table.hpp"

#if defined(HAVE_MULTIPRECISION) && defined(HAVE_QUAD_PRECISION)
using namespace boost::multiprecision;
#else
using namespace std;
#endif

/**
 * @brief Unit test for the Table class.
 *
 * @param argc Number of command line arguments (ignored).
 * @param argv Command line arguments (ignored).
 * @return Exit code: 0 on success.
 */
int main(int argc, char **argv) {

  LinearInterpolatedTable<3, float_type> table;
  table.from_ascii_file("test_table.txt");

  for (uint_fast32_t i = 0; i < table.size(); ++i) {
    ctm_warning("%g %g %g", double(table[i][0]), double(table[i][1]),
                double(table[i][2]));
  }

  ctm_warning("%g %g %g", double(table.get_row<0>(1.5)[0]),
              double(table.get_row<0>(1.5)[1]),
              double(table.get_row<0>(1.5)[2]));

  return 0;
}
