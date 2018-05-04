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

#include "laguerre.h"

/* calculate the Laguerre polynomial */

double laguerre(int n, double x) {
    int i;

    if ( n < 0 ) {
        return -1;
    }

    double v[n+1];
    v[0] = 1.0;

    if ( n == 0 ) {
        return 1.0;
    }

    v[1] = 1.0 - x;

    for ( i = 2; i <= n; i++ ) {
            v[i] = ( ( ( double ) ( 2 * i - 1 ) - x ) * v[(i-1)]
                                     + ( double ) (     - i + 1 )                    * v[(i-2)] )
                                     / ( double ) (         i         );
    }

    return v[n];
}

/* calculate the associated Laguerre polynomial */

double laguerre_l(int n, int m, double x) {
    int i;
    double v[n+1];

    if ( n < 0 ) {
        return -1;
    }

    v[0] = 1.0;
    for (i = 1; i <= n; i++) {
        v[i] = 0.0;
    }

    if ( n == 0 ) {
        return v[0];
    }

    v[1] = ( double ) ( m + 1 ) - x;

    for ( i = 2; i <= n; i++ ) {
            v[i]            = ( ( ( double ) (     m + 2 * i - 1 ) - x        ) * v[i-1]
                                        + ( double ) ( - m         - i + 1 )                    * v[i-2] )
                                        / ( double ) (                     i         );
    }

    return v[n];
}
