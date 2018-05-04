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

#include <chrono>
#include <string>
#include <tclap/CmdLine.h>

#define BOUNDING_BOX 75.0

#include "config.h"
#include "scalar_field.h"
#include "isosurface_mesh.h"
#include "wavefunction.h"
#include "integrator.h"

int main(int argc, char **argv) {
    try {
        TCLAP::CmdLine cmd("Construct hydrogen-like orbital.", ' ', PROGRAM_VERSION);

        // input filename
        TCLAP::ValueArg<int> arg_n("n","n","Primary quantum number", true, 1, "1");
        TCLAP::ValueArg<int> arg_l("l","l","Angular quantum number", true, 0, "0");
        TCLAP::ValueArg<int> arg_m("m","m","Magnetic quantum number", true, 0, "0");
        TCLAP::ValueArg<std::string> arg_o("o","o","Output file", true, "test.obj", "test.obj");

        cmd.add(arg_n);
        cmd.add(arg_l);
        cmd.add(arg_m);
        cmd.add(arg_o);

        cmd.parse(argc, argv);

        int n = arg_n.getValue();
        int l = arg_l.getValue();
        int m = arg_m.getValue();
        const std::string output_filename = arg_o.getValue();

        WaveFunction wf(n,l,m);
        Integrator in(&wf);
        double isovalue = in.find_isosurface_volume(0, 100, 1e5);

        unsigned int gridsize = 201;
        double resolution = BOUNDING_BOX * 2.0 / (double)gridsize;
        ScalarField sf(gridsize, resolution);
        sf.load_wavefunction(wf);

        // analyze the grid and generate the isosurface using the isovalue
        IsoSurface is(&sf);
        is.marching_cubes(isovalue);

        IsoSurfaceMesh ism(&sf, &is);
        ism.construct_mesh(true);
        ism.write_obj(output_filename, "test", "test");

        return 0;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
