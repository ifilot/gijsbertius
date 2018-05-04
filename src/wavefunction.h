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

#ifndef _WAVEFUNCTION_H
#define _WAVEFUNCTION_H

#include <iostream>
#include "stdlib.h"
#include "factorial.h"
#include "legendre.h"
#include "laguerre.h"

class WaveFunction{
private:
    // principal, angular and magnetic quantum number
    int n, l, m;

public:
    WaveFunction(int _n, int _l, int _m);
    double value(double r, double theta,
                 double phi) const;

    double radial_distribution_function(double r) const;

    double polar_density(double theta) const;

private:
    double radial_function(double r) const;

    double polar_function(double theta) const;

    double azimuthal_function(double phi) const;
};

#endif
