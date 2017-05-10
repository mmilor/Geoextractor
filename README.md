# Geoextractor
**An interface to connect MaGe simulations to pulse shape simulations of ADL and SIGGEN.**

Achieved by read-in of MaGe .ROOT files into Geoextractor, resulting in the generation of output files for each detector in the native formats of ADL and SIGGEN. JSON files with all detector and hit information are generated as well. (This is also helpful for programs in the future that need the MaGe geometry extracted.)

## Description

MaGe simulations : .ROOT output into Geoextractor
* Geoextractor generates all necessary files for ADL and SIGGEN in their native format. 

Geoextractor **extracts** automatically:
* Detector absolute position inside of MaGe. [.GDML]
* Detector parameters (radius, length, dead layer size, etc.). [.GDML]
* Hit positions transformed into the individual detector frames. [.ROOT]

**Input** for Geoextractor:
* .GDML file extracted from MaGe geometry through MaGe macro.
* detsettings file that contains additional information (e.g. HV).
* .ROOT files of a Monte Carlo simulation with the same MaGe geometry.

**Output** automatically generated by Geoextractor:
* .JSON file for each detector (transformed hit positions/energy/particleID, detector parameters, detector absolute positions)
* .config files for SIGGEN for each detector (detector parameters)
* DET & DET_SETUP files for ADL (detector parameters)
* (hit parameters in simple .txt file for test purposes, but will be replaced fully by JSON usage for all programs)

## Installation

Installation from github in directory of choice:
```
git clone git@github.com:mmilor/Geoextractor.git
```
## Usage
The working directory after installation already contains:
* GDML file 'gdml.gdml' suitable for GERDA Phase II created with the detector- and matrixfiles:
  * geometry_PhaseII_FCCD_growth.dat
  * matrix_phase_ii_StatusDec2015_pos_vis.dat
  * both found here: https://github.com/mppmu/MaGe/tree/GERDAphaseII/gerdageometry
* 'detsettings.txt' file that contains additional information about the GERDA Phase II detectors. Any upgrade with addition channels, change of HV etc. means that this file needs to be edited, but Geoextractor does not need to be recompiled.

These default input files are usually enough for the GERDA Phase II geometry. 
If a different arrangement or detector- and matrixfile combination is desired, a new gdml file can be produced as described later on in the section "Producing New .GDML Files".

Geoextractor can be run directly as CINT ROOT program or compiled standalone. 

### Compile and run
Geoextractor is compiled with the following command:
```
cd /path/to/Geoextractor
make
```

Options:
*	**-m, --mcdir PATH**:	Specify input path for ROOT files to be merged and used for extraction. Default: MC_DIR shell variable.
*	**-a, --adl PATH**:	Specify output path for ADL files. Default: ./Output/adl
*	**-j, --json PATH**:	Specify output path for JSON files. Default: ./Output/json
*	**-s, --siggen PATH**:	Specify output path for SIGGEN files. Default: ./Output/siggen
*	**-g, --gdml PATH**:	Specify path to GDML file. Default: ./defaultgdml.gdml
*	**-d, --detsettings PATH**:	Specify path to Detector Settings file. Default: ./detsettings.txt
*	**-x, --maxhits**: INT	Limit amount of hits per detector to be extracted and written to file to a fixed number. Default: No limit.
*	**-n, --nohits**: 	Skip all hit extraction, only extract geometry.
* **-c, --checkgdml**: 	Visualise .GDML file with OGL drivers. Best to use: root Geoextractor.cxx+(1)
*	**-h, --help**:		Show help message

To extract only the geometries without any hits, use the command: 
```
./Geoextractor --nohits
```

In order to extract hits out of .ROOT files, the directory of the files needs to be specified (The files should be simulated in the same detector- and matrixfile configuration as the .GDML file.). 

There are two options (order of priority):
* Via option: **-m,--mcdir PATH**
* Via shell variable $MC_DIR pointing to the directory

To specify $MC_DIR, use the following shell commands:  
```
cd /path/to/MCSimulations
export MC_DIR=$(pwd)
```
Simply evoking the program is then enough to extract into the default folders:
```
cd /path/to/Geoextractor
./Geoextractor
```

### Run with CINT ROOT 
Geoextractor can also be run directly by CINT ROOT.
First the shell variable $MC_DIR should be set to the directory of the .ROOT MC files from where the hits should be extracted. (The files should be simulated in the same detector- and matrixfile configuration as the .GDML file.) 
```
cd /path/to/MCSimulations
export MC_DIR=$(pwd)
cd /path/to/Geoextractor
```

When in the Geoextractor directory, the program will be executed by simply evoking:
```
root Geoextractor.cxx
```
In this configuration, Geoextractor will use default options and use as Input:
* 'defaultgdml.gdml' in the Geoextractor directory
* 'detsettings.txt' in the Geoextractor directory
* All .ROOT files contained in the directory pointed at by shell variable $MC_DIR 

The Output will be saved in the directories: 
* /Output/adl
* /Output/json
* /Output/siggen

This is an easy method to get results.

### Producing New .GDML Files
If a different arrangement or detector- and matrixfile combination is desired, a new gdml file can be produced with a GDML macro provided in the working directory: 'gdmlextractor.mac'. It needs to be edited by changing only three lines to the appropriate new file names:
```
In gdmlextractor.mac:
(...)
/MG/geometry/detector/geometryfile geometry_PhaseII_FCCD_growth.dat
/MG/geometry/detector/matrixfile matrix_phase_ii_StatusDec2015_pos_vis.dat
(...)
/MG/geometry/GDML/outputName extractedgdml
```
Make sure to use the _pos_vis.dat matrixfile and check the gdml later on visually.
Then MaGe needs to be run with the macro. gdmlextractor.mac can for example be copied to the MaGe path and executed from there.
```
cd /path/to/MaGe
./MaGe gdmlextractor.mac
```
This results in a file 'extractedgdml.gdml' (depending on how the macro was edited) in the MaGe folder. Geoextractor can link to the .gdml file with the command **-g, --gdml PATH**. Another possibility is to rename the defaultfile in Geoextractor folder and copy and rename your new GDML to 'defaultgdml.gdml' to automatically use this file without needing any options.

To visually check the .GDML file, use ROOT CINT with the command (explained in the corresponding section): 
```
root Geoextractor.cxx+(1)
```
This sets the checkGDML option to true and enables you to look at a OGL rendering of the GDML file.

## File formatting
The ADL and SIGGEN file structure for the detector geometries/parameters are given by the individual programs. The JSON formatting is created for Geoextractor and looks as follows (here for GD00A and two hits):
```
{ 
    "detector": "GD00A",
    "magechannel": "33",
    "type": "BEGe",
    "mageposition": 
        [
            { "detcenter_xpos": "55", "detcenter_ypos": "-95.2628", "detcenter_zpos": "67.7497" }
        ],
    "parameters": 
        [
            { "detdimx": "34.25", "detdimy": "34.25", "detdimz": "12.7497", "xtal_radius": "34.25", "xtal_length": "25.4994", "pc_length": "3e-4", "pc_radius": "7.5", "wrap_around_radius": "10.5", "ditch_depth": "2", "ditch_thickness": "3", "outer_taper_length": "12.06", "outer_taper_width": "11.92", "Li_thickness": "0.91", "impurity_z0": "-1", "impurity_gradient": "0.0454", "xtal_HV": "2500" }
        ],
    "hits": 
        [
            { "id": "0", "eventnumber": "2276", "xpos": "-29.9293", "ypos": "-2.55369", "zpos": "13.0999", "edep": "11.067", "trackpdg": "22" },
            { "id": "1", "eventnumber": "27604", "xpos": "19.1338", "ypos": "-25.4563", "zpos": "11.7349", "edep": "0.12938", "trackpdg": "22" }
        ]
}
``` 
* **detcenter_xpos/ypos/zpos** describe the absolute position of the detector center in the MaGe coordinate system.
* **detdimx/y/z** describe the dimensions of the detector in each direction from the center point. Point contact position thus is detcenter_xpos, detcenter_ypos, detcenter_zpos-detdimz. Inverted detectors: " " detcenter_zpos+detdimz.
* xtal_radius, length etc are described here: radware.phy.ornl.gov/MJ/mjd_siggen/fieldgen_geometry.pdf
* **xpos, ypos, zpos** are the coordinates of each hit already transformed into the detector frame with the point contact being (0,0,0). (Inverted detectors are already taken into consideration here).

To be added:
* inverted (already used/determined/known, but not yet written into JSON file)
* coax Li_thickness faulty atm, needs to be fixed
* ...

## Plans
* Unify JSON usage for hits read-in for SIGGEN and ADL.
* Work out optimal file formatting.
* SSE/MSE identification and categorisation inside JSON.
* Proper geometry extraction for Coax as well. (Currently only verified with BEGes)
* Bugtesting.

### Prerequisites
Geoextractor:
* ROOT

GDML macro:
* MaGe

MaGe for MC simulations:
* GEANT4


## Author

Michael Miloradovic 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
