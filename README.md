# gijsbertius

## Purpose
Generate hydrogen-like orbitals that can be imported in a 3D rendering program such as Blender

## Dependencies
Gijsbertius depends on the following packages being installed
* GLM
* Boost (format)
* TCLAP

```
sudo apt-get install build-essential cmake libglm-dev libboost-all-dev libtclap-dev
```

## Compilation

```
mkdir build
cd build
cmake ../src
make -j9
```

## Usage
```
./gijsbertius -n 5 -l 1 -m 0 -o 510.ply [-s]
```
The `-s` tag splits positive and negative contributions and generates two files.

## Converting to stl
Install `openctm-tools`
```
sudo apt install openctm-tools
```

and run (e.g.)

```
ctmconv 100.ply 100.stl

```
