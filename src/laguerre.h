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

#ifndef _LAGUERRE_H
#define _LAGUERRE_H

/**
 * @brief      Laguerre polynomial
 *
 * @param[in]  n     nth Laguerre polynomial
 * @param[in]  x     value x to calculate
 *
 * @return     value at x of the nth Laguerre polynomial
 */
double laguerre(int n, double x);

/**
 * @brief      Associated Laguerre polynomial
 *
 * @param[in]  n     nth Laguerre polynomial
 * @param[in]  m     value m
 * @param[in]  x     value x to calculate
 *
 * @return     value at x of the nth/mth Associated Laguerre polynomial
 */
double laguerre_l(int n, int m, double x);

#endif
