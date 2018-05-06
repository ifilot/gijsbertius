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

#include "integrator.h"

Integrator::Integrator(WaveFunction *_wf) {
    this->wf = _wf;
};

double Integrator::integrate_radial_density(double a, double b, int steps, double* r) {
    double s = 0;
    double h = (b-a)/steps;
    *r = 0;

    for (int i = 0; i < steps; ++i) {
        *r = a + h * (double)i;

        s += (this->wf->radial_distribution_function(*r) +
              4 * this->wf->radial_distribution_function(*r + h / 2.0) +
              this->wf->radial_distribution_function(*r + h)) / 6.0;
    }

    return h * s;
}

double Integrator::find_isosurface_volume(double a, double b, int steps, double* r) {
    double s = 0;
    double h = (b-a)/steps;
    *r = 0;

    for (int i = 0; i < steps; ++i) {
        *r = a + h * (double)i;

        s += (this->wf->radial_distribution_function(*r) +
              4 * this->wf->radial_distribution_function(*r + h / 2.0) +
              this->wf->radial_distribution_function(*r + h)) / 6.0;
        if(h * s >= 0.95) {
            return std::sqrt(this->wf->radial_distribution_function(*r) /
                       (4.0 * M_PI * (*r) * (*r)));
        }
    }

    return -1;
}

double Integrator::integrate_spherical_harmonic(double a, double b, int steps) {
    double s = 0;
    double h = (b-a)/steps;
    double r = 0;
    for (int i = 0; i < steps; ++i) {
        r = a + h * (double)i;

        s += (this->wf->polar_density(r) +
              4 * this->wf->polar_density(r + h / 2.0) +
              this->wf->polar_density(r + h)) / 6.0;
    }
    return h * s / 2.0; // due to azimuthal function
}

double Integrator::integrate_x2(double a, double b, int steps) {
    double s = 0;
    double h = (b-a)/steps;
    double r = 0;

    for (int i = 0; i < steps; ++i) {
        r = a + h * (double)i;

        s += ((std::pow(r, 2.0)) +
              4.0 * (std::pow(r + h / 2.0, 2.0)) +
              (std::pow(r + h, 2.0))) / 6.0;
    }

    return h * s; // due to azimuthal function
}
