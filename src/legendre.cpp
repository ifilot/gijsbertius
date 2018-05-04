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

#include "legendre.h"


/**
 * @brief      legendre function
 *
 * @param[in]  n     { parameter_description }
 * @param[in]  x     { parameter_description }
 *
 * @return     { description_of_the_return_value }
 */
double legendre (int n, double x) {
    int i;

    if(n < 0) {
        return -1;
    }

    double v[n];
    v[0] = 1.0;

    if(n < 1) {
        return 1.0;
    }

    v[1] = x;

    for ( i = 2; i <= n; i++ ) {
            v[i] =         ( ( double ) ( 2 * i - 1 ) * x        * v[i-1]
                                 - ( double ) (         i - 1 ) *                v[i-2] )
                                 / ( double ) (         i         );
    }

    return v[n];
}

/**
 * @brief      associated legendre function
 *
 * @param[in]  n     n
 * @param[in]  m     m
 * @param[in]  x     value x
 *
 * @return     value
 */
double legendre_p (int n, int m, double x) {
    double fact;
    int i;
    int j;
    int k;
    double v[n+1];

    for ( i = 0; i < n + 1; i++ ) {
        v[i] = 0.0;
    }

    if ( m <= n ) {
        v[m] = 1.0;

        fact = 1.0;
        for ( k = 0; k < m; k++ ) {
            v[m] *= - fact * sqrt ( 1.0 - x * x);
            fact += 2.0;
        }
    }

    if ( m + 1 <= n ) {
        v[m+1] = x * ( double ) ( 2 * m + 1 ) * v[m];
    }

    for ( j = m + 2; j <= n; j++ ) {
            v[j] = ( ( double ) ( 2 * j         - 1 ) * x * v[(j-1)]
                    + ( double ) (     - j - m + 1 ) *                v[(j-2)] )
                    / ( double ) (         j - m         );
    }

    return v[n];
}
