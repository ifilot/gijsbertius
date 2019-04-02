//
// Example POVRAY rendering script
//
// To create the input for this file, run in the build folder
// ./gijsbertius -n 5 -l 1 -m 0 -o 510.obj -d 510.df3
//

#declare NX = 301;
#declare NY = 301;
#declare NZ = 301;
#declare DIAG = <NX, NY, NZ>;

global_settings {
    ambient_light <1,1,1>
    assumed_gamma 1
}

camera {
    location <0, -5/4 * NY, 0>
    up z
    right x
    sky <0,0,1>
    look_at <0,0,0>
}

light_source {
    <2 * NX, -NY, 2*NZ>
    color rgb <1,1,1>
    media_interaction on
    media_attenuation on
    shadowless
}

#declare DENS = interior {
    media {
        intervals 100
        ratio 0.5
        samples 3,3
        method 3
        emission 3*<1,1,1>/100
        absorption <1,1,1>/1000
        scattering {1, <0,0,0>}
        confidence 0.99
        variance 1/1000
        density {
            density_file df3 "../build/510.df3"
            interpolate 3
            color_map {
                #include "YlGn.inc"
            }
        }
    }
}

box {
    <0,0,0>, <1,1,1>
    pigment {rgbt <0,0,0,1>}
    hollow
    interior {DENS}
    scale DIAG
    translate -DIAG / 2
    rotate <360*clock, 0, 0>
}
