/*******************************************************************************
 * This file is part of QuickSched.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

/*******************************************************************************
 * CHANGELOG
 *  - 22/10/2019, Bert Vandenbroucke
 *    changed names of macros to solve some conflicts with existing functions
 *    with the same name.
 ******************************************************************************/

/* Error macro. */
#define quicksched_error(s, ...)                                               \
  {                                                                            \
    fprintf(stderr, "%s:%s():%i: " s "\n", __FILE__, __FUNCTION__, __LINE__,   \
            ##__VA_ARGS__);                                                    \
    abort();                                                                   \
  }

/* Message macro. */
#define quicksched_message(s, ...)                                             \
  {                                                                            \
    printf("%s: " s "\n", __FUNCTION__, ##__VA_ARGS__);                        \
    fflush(stdout);                                                            \
  }
