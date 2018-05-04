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

#ifndef _INTEGRATOR_H
#define _INTEGRATOR_H

#include "wavefunction.h"

class Integrator{
private:
    //unsigned int steps;
    WaveFunction *wf;

public:
    Integrator(WaveFunction *_wf);
    double integrate_radial_density(double a, double b, int steps);
    double find_isosurface_volume(double a, double b, int steps);
    double integrate_spherical_harmonic(double a, double b, int steps);
    double integrate_x2(double a, double b, int steps); //testroutine
};

#endif
