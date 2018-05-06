/**************************************************************************
 *   This file is part of GIJSBERTIUS.                                    *
 *                                                                        *
 *   Author: Ivo Filot <ivo@ivofilot.nl>                                  *
 *                                                                        *
 *   GIJSBERTIUS is free software:                                        *
 *   you can redistribute it and/or modify it under the terms of the      *
 *   GNU General Public License as published by the Free Software         *
 *   Foundation, either version 3 of the License, or (at your option)     *
 *   any later version.                                                   *
 *                                                                        *
 *   GIJSBERTIUS is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty          *
 *   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.              *
 *   See the GNU General Public License for more details.                 *
 *                                                                        *
 *   You should have received a copy of the GNU General Public License    *
 *   along with this program.  If not, see http://www.gnu.org/licenses/.  *
 *                                                                        *
 **************************************************************************/

#ifndef _SCALAR_FIELD_H
#define _SCALAR_FIELD_H

#include <string>
#include <vector>
#include <iostream>
#include <ios>
#include <sstream>
#include <fstream>
#include <math.h>

#include <glm/glm.hpp>

#include "wavefunction.h"

class ScalarField{
private:
    double scalar;

    glm::mat3 mat33;                //!< glm version of the matrix
    glm::mat3 imat33;               //!< glm version of the inverse matrix

    unsigned int grid_dimensions[3];
    std::vector<double> gridptr;    //!< grid to first pos of double array

    unsigned int gridsize;
    double resolution;

public:

    ScalarField(unsigned int _gridsize, double _resolution);

    void load_wavefunction(const WaveFunction &wf);

    /*
     * double get_value_interp(x,y,z)
     *
     * Grabs a value from the 3D scalar field. Calculate the value
     * by using a trilinear interpolation.
     *
     * The trilinear interpolation algorithm has been extracted from:
     * http://paulbourke.net/miscellaneous/interpolation/
     *
     * Future algorithm can make use of a cubic interpolation.
     *
     */
    double get_value_interp(double x, double y, double z) const;

    /**
     * @brief      test whether point is inside unit cell
     *
     * @param[in]  x     x position
     * @param[in]  y     y position
     * @param[in]  z     z position
     *
     * @return     True if inside, False otherwise.
     */
    bool is_inside(double x, double y, double z) const;

    double get_value(unsigned int i, unsigned int j, unsigned int k) const;

    glm::dvec3 grid_to_realspace(double i, double j, double k) const;

    glm::dvec3 realspace_to_grid(double i, double j, double k) const;

    glm::dvec3 realspace_to_direct(double i, double j, double k) const;

    void copy_grid_dimensions(unsigned int _grid_dimensions[]) const;

    inline const glm::mat3& get_mat_unitcell() const {
        return this->mat33;
    }

    inline const glm::mat3& get_mat_unitcell_inverse() const {
        return this->imat33;
    }

    inline auto* get_grid_ptr() const {
        return &this->gridptr[0];
    }

    unsigned int get_size() const {
        return this->gridptr.size();
    }

private:
    void calculate_inverse();

    unsigned int get_idx(unsigned int i, unsigned int j, unsigned int k) const;

    void set_value(unsigned int i, unsigned int j, unsigned int k, double _value);

    void init();
};

#endif //_SCALAR_FIELD_H
