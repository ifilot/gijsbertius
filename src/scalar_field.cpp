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

#include "scalar_field.h"

/**
 * @brief      constructor
 *
 * @param[in]  _filename   url to filename
 * @param[in]  _flag_is_locpot  whether this file is a locpot
 */
ScalarField::ScalarField(unsigned int _gridsize, double _resolution) :
    gridsize(_gridsize),
    resolution(_resolution) {

    this->init();
}

void ScalarField::load_wavefunction(const WaveFunction &wf, bool sgnd) {
    double r = 0;
    double phi = 0;
    double theta = 0;

    unsigned int i=0, j=0, k=0;

    if(sgnd) {
        for(double x = -(double)(gridsize - 1)/2.0 * resolution; x <= (gridsize - 1)/2.0 * resolution; x += resolution) {
            j = 0;
            for(double y = -(double)(gridsize - 1)/2.0 * resolution; y <= (gridsize - 1)/2.0 * resolution; y += resolution) {
                k = 0;
                for(double z = -(double)(gridsize - 1)/2.0 * resolution; z <= (gridsize - 1)/2.0 * resolution; z += resolution) {
                    /* transform to spherical coordinates */
                    r = std::sqrt(x*x + y*y + z*z);
                    phi = std::atan2(y, x);
                    theta = std::acos(z / r);

                    /* obtain wavefunction value from object */
                    this->set_value(i, j, k, wf.value(r, theta, phi));
                    k++;
                }
                j++;
            }
            i++;
        }
    } else {
        for(double x = -(double)(gridsize - 1)/2.0 * resolution; x <= (gridsize - 1)/2.0 * resolution; x += resolution) {
            j = 0;
            for(double y = -(double)(gridsize - 1)/2.0 * resolution; y <= (gridsize - 1)/2.0 * resolution; y += resolution) {
                k = 0;
                for(double z = -(double)(gridsize - 1)/2.0 * resolution; z <= (gridsize - 1)/2.0 * resolution; z += resolution) {
                    /* transform to spherical coordinates */
                    r = std::sqrt(x*x + y*y + z*z);
                    phi = std::atan2(y, x);
                    theta = std::acos(z / r);

                    /* obtain wavefunction value from object */
                    this->set_value(i, j, k, std::abs(wf.value(r, theta, phi)));
                    k++;
                }
                j++;
            }
            i++;
        }
    }
}

/**
 * @brief      Saves to df3 density file
 *
 * @param[in]  filename  The filename
 */
void ScalarField::save_to_df3(const std::string& filename) {
    // get largest and smallest value
    auto minmax = std::minmax_element(this->gridptr.begin(), this->gridptr.end());
    double min = *minmax.first;
    double max = *minmax.second;
    double ival = max - min;

    std::ofstream out(filename, std::ios::binary);

    // introduce dummy variable to store byteswapped value;
    // note that df3 needs to be stored in lower endian!
    uint16_t byteswap = 0;

    // write data size
    uint16_t nx = this->grid_dimensions[0];
    uint16_t ny = this->grid_dimensions[1];
    uint16_t nz = this->grid_dimensions[2];

    byteswap = (nx>>8) | (nx<<8);
    out.write((char*)&byteswap, sizeof(uint16_t));
    byteswap = (ny>>8) | (ny<<8);
    out.write((char*)&byteswap, sizeof(uint16_t));
    byteswap = (nz>>8) | (nz<<8);
    out.write((char*)&byteswap, sizeof(uint16_t));

    for(unsigned int i=0; i<this->gridptr.size(); i++) {
        uint16_t val = (uint16_t)((this->gridptr[i] - min) / ival * 65536);
        byteswap = (val>>8) | (val<<8);
        out.write((char*)&byteswap, sizeof(uint16_t));
    }

    out.close();
}

/* set a value in the grid by specifying its grid location */
void ScalarField::set_value(unsigned int i, unsigned int j, unsigned int k, double _value) {
    this->gridptr[this->get_idx(i, j, k)] = _value;
}

unsigned int ScalarField::get_idx(unsigned int i, unsigned int j, unsigned int k) const {
    return k * this->gridsize * this->gridsize + j * this->gridsize + i;
}

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
double ScalarField::get_value_interp(double x, double y, double z) const {
    if(!this->is_inside(x,y,z)) {
        return 0.0f;
    }

    // cast the input to the nearest grid point
    glm::dvec3 r = this->realspace_to_grid(x,y,z);
    glm::dvec3 d = this->realspace_to_direct(x,y,z);

    // calculate value using trilinear interpolation
    double xd = remainderf(r[0], 1.0);
    double yd = remainderf(r[1], 1.0);
    double zd = remainderf(r[2], 1.0);

    if(xd < 0) xd += 1.0;
    if(yd < 0) yd += 1.0;
    if(zd < 0) zd += 1.0;

    double x0 = floor(r[0]);
    double x1 = ceil(r[0]);
    double y0 = floor(r[1]);
    double y1 = ceil(r[1]);
    double z0 = floor(r[2]);
    double z1 = ceil(r[2]);

    return
    this->get_value(x0, y0, z0) * (1.0 - xd) * (1.0 - yd) * (1.0 - zd) +
    this->get_value(x1, y0, z0) * xd                 * (1.0 - yd) * (1.0 - zd) +
    this->get_value(x0, y1, z0) * (1.0 - xd) * yd                 * (1.0 - zd) +
    this->get_value(x0, y0, z1) * (1.0 - xd) * (1.0 - yd) * zd                 +
    this->get_value(x1, y0, z1) * xd                 * (1.0 - yd) * zd                 +
    this->get_value(x0, y1, z1) * (1.0 - xd) * yd                 * zd                 +
    this->get_value(x1, y1, z0) * xd                 * yd                 * (1.0 - zd) +
    this->get_value(x1, y1, z1) * xd                 * yd                 * zd;
}

/**
 * @brief      test whether point is inside unit cell
 *
 * @param[in]  x     x position
 * @param[in]  y     y position
 * @param[in]  z     z position
 *
 * @return     True if inside, False otherwise.
 */
bool ScalarField::is_inside(double x, double y, double z) const {
    // cast the input to the nearest grid point
    glm::dvec3 d = this->realspace_to_direct(x,y,z);

    if(d[0] < 0 || d[0] > 1.0) {
        return false;
    }
    if(d[1] < 0 || d[1] > 1.0) {
        return false;
    }
    if(d[2] < 0 || d[2] > 1.0) {
        return false;
    }

    return true;
}

/*
 * void calculate_inverse()
 *
 * Calculates the inverse of a 3x3 matrix. This is a convenience
 * function for the read_matrix() function.
 *
 */
void ScalarField::calculate_inverse() {
    this->imat33 = glm::inverse(mat33);
}

/*
 * double get_value(i,j,k)
 *
 * Grabs the value at a particular grid point.
 *
 * This is a convenience function for the get_value_interp() function
 *
 */
double ScalarField::get_value(unsigned int i, unsigned int j, unsigned int k) const {
    unsigned int idx = k * this->grid_dimensions[0] * this->grid_dimensions[1] +
                       j * this->grid_dimensions[0] +
                       i;

    return this->gridptr[idx];
}

/*
 * glm::dvec3 grid_to_realspace(i,j,k)
 *
 * Converts a grid point to a realspace vector. This function
 * is not being used at the moment.
 *
 */
glm::dvec3 ScalarField::grid_to_realspace(double i, double j, double k) const {
    double dx = (double)i / (double)grid_dimensions[0];
    double dy = (double)j / (double)grid_dimensions[1];
    double dz = (double)k / (double)grid_dimensions[2];

    glm::dvec3 r(dx, dy, dz);

    return mat33 * r;
}

/*
 * glm::dvec3 realspace_to_grid(i,j,k)
 *
 * Convert 3d realspace vector to a position on the grid. Non-integer
 * values (i.e. doubleing point) are given as the result.
 *
 * This is a convenience function for the get_value_interp() function
 *
 */
glm::dvec3 ScalarField::realspace_to_grid(double i, double j, double k) const {
    glm::dvec3 r;
    r[0] = imat33[0][0] * i + imat33[0][1] * j + imat33[0][2] * k;
    r[1] = imat33[1][0] * i + imat33[1][1] * j + imat33[1][2] * k;
    r[2] = imat33[2][0] * i + imat33[2][1] * j + imat33[2][2] * k;

    r[0] *= double(this->grid_dimensions[0]-1);
    r[1] *= double(this->grid_dimensions[1]-1);
    r[2] *= double(this->grid_dimensions[2]-1);

    return r;
}

/*
 * glm::dvec3 realspace_to_direct(i,j,k)
 *
 * Convert 3d realspace vector to direct position.
 *
 */
glm::dvec3 ScalarField::realspace_to_direct(double i, double j, double k) const {
    glm::dvec3 r;
    r[0] = imat33[0][0] * i + imat33[0][1] * j + imat33[0][2] * k;
    r[1] = imat33[1][0] * i + imat33[1][1] * j + imat33[1][2] * k;
    r[2] = imat33[2][0] * i + imat33[2][1] * j + imat33[2][2] * k;

    return r;
}

void ScalarField::copy_grid_dimensions(unsigned int _grid_dimensions[]) const {
    for(unsigned int i=0; i<3; i++) {
        _grid_dimensions[i] = this->grid_dimensions[i];
    }
}

void ScalarField::init() {
    this->gridptr.resize(this->gridsize * this->gridsize * this->gridsize, 0.0);
    this->grid_dimensions[0] = this->gridsize;
    this->grid_dimensions[1] = this->gridsize;
    this->grid_dimensions[2] = this->gridsize;

    this->mat33[0][0] = this->gridsize * this->resolution;
    this->mat33[0][1] = 0.0;
    this->mat33[0][2] = 0.0;

    this->mat33[1][0] = 0.0;
    this->mat33[1][1] = this->gridsize * this->resolution;
    this->mat33[1][2] = 0.0;

    this->mat33[2][0] = 0.0;
    this->mat33[2][1] = 0.0;
    this->mat33[2][2] = this->gridsize * this->resolution;

    this->calculate_inverse();
}
