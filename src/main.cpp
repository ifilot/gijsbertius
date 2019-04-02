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
#include <boost/format.hpp>

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
        TCLAP::ValueArg<std::string> arg_d("d","d","Density file", false, "", "density.d3d");

        cmd.add(arg_n);
        cmd.add(arg_l);
        cmd.add(arg_m);
        cmd.add(arg_o);
        cmd.add(arg_d);

        cmd.parse(argc, argv);

        std::cout << "--------------------------------------------------------------" << std::endl;
        std::cout << "Executing Gijsbertius v." << PROGRAM_VERSION << std::endl;
        std::cout << "Author: Ivo Filot <ivo@ivofilot.nl>" << std::endl;
        std::cout << "--------------------------------------------------------------" << std::endl;

        int n = arg_n.getValue();
        int l = arg_l.getValue();
        int m = arg_m.getValue();
        const std::string output_filename = arg_o.getValue();

        auto start = std::chrono::system_clock::now();

        WaveFunction wf(n,l,m);
        Integrator in(&wf);

        double r = 0;
        double isovalue = in.find_isosurface_volume(0, 100, 1e5, &r);

        unsigned int gridsize = 301;
        double resolution = std::ceil(r * 2.0) * 2.0 / (double)gridsize;
        ScalarField sf(gridsize, resolution);
        sf.load_wavefunction(wf);

        // analyze the grid and generate the isosurface using the isovalue
        IsoSurface is(&sf);
        is.marching_cubes(isovalue);

        IsoSurfaceMesh ism(&sf, &is);
        ism.construct_mesh(true);
        const std::string name = (boost::format("%i%i%i") % n % l % m).str();
        ism.write_obj(output_filename, name, name);

        std::string densityfile = arg_d.getValue();
        if(!densityfile.empty()) {
            sf.save_to_df3(densityfile);
        }

        auto end = std::chrono::system_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << boost::format("Total elapsed time: %f ms\n") % elapsed.count();
        std::cout << "--------------------------------------------------------------" << std::endl;

        return 0;

    } catch (TCLAP::ArgException &e) {
        std::cerr << "error: " << e.error() <<
                     " for arg " << e.argId() << std::endl;
        return -1;
    }
}
