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

 #include "wavefunction.h"

WaveFunction::WaveFunction(int _n, int _l, int _m) {
    if(_n < 1) {
        throw std::runtime_error("Invalid principal quantum number, exiting...");
    }

    if(_l < 0) {
        throw std::runtime_error("Invalid angular quantum number ( l < 0), exiting...");
    }

    if(_l > _n - 1) {
        throw std::runtime_error("Invalid angular quantum number (l > n - 1), exiting...");
    }

    if(_m < -_l) {
        throw std::runtime_error("Invalid magnetic quantum number, exiting...");
    }

    if(_m > _l) {
        throw std::runtime_error("Invalid magnetic quantum number, exiting...");
    }

    this->n = _n;
    this->l = _l;
    this->m = _m;
}

double WaveFunction::value(double r, double theta, double phi) const {
    return radial_function(r) * polar_function(theta) * azimuthal_function(phi);
}

double WaveFunction::radial_distribution_function(double r) const {
    double val = this->radial_function(r);
    return r * r * val * val;
}

double WaveFunction::polar_density(double theta) const {
    double val = this->polar_function(theta);
    return val * val * sin(theta);
}

/* calculate the radial wave functions using atomic units
 *
 * PLEASE READ THIS CAREFULLY!!
 * I used two literature sources to construct this function,
 * being "Introduction to Quantum Mechanics" from Griffiths and
 * "Intoduction to Quantum Mechanics With Applications to Chemistry" from
 * Linus Pauling and E. Bright Wilson, Jr.
 * Both books discuss the radial wavefunction and use an associated
 * Laguerre polynomial. Neither of the formulas in these books seem to
 * give me the correct results and I believe this has to do with the definition
 * of the Laguerre polynomial and/or the normalization constant. In fact,
 * Griffiths makes a point of this (ref 16, page 152).
 *
 * In the end, I checkedmy Laguerre and Legendre polynomials by comparing
 * these at Wolfram Alpha (www.wolframalpha.com) and
 * use some reverse engineering to find the proper normalisation constants. I
 * am not 100% sure, but at least the radial distribution function integrates
 * to unity for all the wave functions.
 *
 * If you can resolve this issue, please contact me and I will gladly give
 * you a very nice beer (or other drink of your choice) for your troubles!!
 */

double WaveFunction::radial_function(double r) const {
    double pre = sqrt(pow(2.0 / (double)n, 3.0) * (double)factorial(n - l - 1) /
                     ( 2.0 * (double)n * (double)factorial(n + l)));
    double rho = 2.0 * r / (double)n;

    return pre * exp(-rho / 2.0) * pow(rho, (double)l) *
                 laguerre_l(n-l-1, 2*l+1, rho);
}

double WaveFunction::polar_function(double theta) const {
    double pre = sqrt((double)(2 * l + 1) * (double)factorial(l - std::abs(m)) /
        (double)factorial(l + std::abs(m)));

    return pre * legendre_p(l, std::abs(m), cos(theta));
}

double WaveFunction::azimuthal_function(double phi) const {
    const double pre = 1.0 / sqrt(4 * M_PI);

    if(m == 0) return pre;

    if(m > 0) {
        return pre * cos((double)m * phi);
    } else {
        return pre * sin(-(double)m * phi);
    }
}
