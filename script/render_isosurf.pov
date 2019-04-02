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
    location <0, -5, 0>
    up z
    right x
    sky <0,0,1>
    look_at <0,0,0>
}

light_source {
    <0, -NY, 2*NZ>
    color rgb <1,1,1>
    media_interaction on
    media_attenuation on
    shadowless
}

#declare DENS = function {
    pigment {
        density_file df3 "../build/510.df3"
        interpolate 3
    }
}

isosurface {
    function { DENS(x,y,z).gray }
    contained_by {
        box {0.0,1.0}
    }
    threshold 0.01
    open
    accuracy 0.001
    max_gradient 50
    scale 5
    pigment {rgb <1,1,1>}
    finish {
        phong 0.8
        specular 1.0
    }
    translate -2.5
}
