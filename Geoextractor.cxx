// Created: Michael Miloradovic March 2017
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "TGeoManager.h"
#include "TApplication.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TFile.h"
#include "TChain.h"

#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TSystem.h"
#include "TRegexp.h"

#include <time.h>


using namespace std;


//Function to search through GDML node levels to find nodes matching germanium detectors. Then extract all information about them. This is run iteratively.
void print_node(TGeoNode *node, string path, unsigned int& MaGeChannelNumIterator, vector<string> convDetNameSortedMC, vector<int> convIsGrooveOnTopSortedMC, vector<int> convIsBEGeSortedMC, vector<int> convIsPSSSortedMC, vector<double> convTaperLengthSortedMC, vector<double> convTaperWidthSortedMC, vector<double> convImpurityZSortedMC, vector<double> convImpurityGradSortedMC, vector<int> convHVSortedMC, vector<double> convFCCD, string jsonOutputPath, string adlOutputPath, string siggenOutputPath, ofstream & jsonOutputFile, ofstream & siggenOutputFile, ofstream & adlOutputFile, ofstream & hitsOutputFile, ofstream & hitsADLOutputFile, TChain * fTree, TChain * aTree, int maxHits, int useHits, int isVisHitsEnabled, string visHitsCond)
{
    string str = "";
    string convString = "debug";

    if (node != NULL)
    {
        const Char_t *node_name = node->GetName();
        path += Form("/%s", node_name);
        TGeoVolume *volume = node->GetVolume();
        const char * material = volume->GetMaterial()->GetName();
        string strMaterial(material);
        if (volume != NULL)
        {
            //Find all germanium detectors.
            if (strMaterial == "NaturalGe" || strMaterial == "EnrichedGe")
            {
                convString = node_name;
                //Ignore passivation layer nodes: Only deadlayer and active nodes/volumes remain.
                if ((convString.find("Passivation") != string::npos) == 0)
                {
                    //Determine if groove/ditch/point contact is on top of the detector.
                    double grooveOnTop = 0.;
                    if(convIsGrooveOnTopSortedMC[MaGeChannelNumIterator] == 1) {
                        grooveOnTop = 1.;
                    };
                    int isCurrentQActiveV = 0;
                    if (convString.find("Active") != string::npos)
                    {
                        cout << "Detector / MaGechannel: " << convDetNameSortedMC[MaGeChannelNumIterator] << " / " << MaGeChannelNumIterator << endl;
                        isCurrentQActiveV = 1;
                    }
                    cout << "Node: " << node_name << endl;
                    //cout << "PATH: " << path << endl;
                    //cout << "VOLUME: " << volume->GetName() << endl;
                    //cout << "CAPACITY: " << volume->Capacity() << endl;
                    //cout << "MATERIAL: " << strMaterial << endl;

                    //Shape properties to get absolute coordinate centers of detectors & detector dimensions
                    TGeoBBox *shape = (TGeoBBox *)volume->GetShape();
                    cout << "Dimensions dx/dy/dz from center: " << shape->GetDX() << " / " << shape->GetDY() << " / " << shape->GetDZ() << endl;
                    cout << "Shape: " << shape->GetName() << endl;
                    if (shape != NULL)
                    {

                        Double_t master[3];
                        const Double_t *local = shape->GetOrigin();
                        gGeoManager->cd(path.c_str());
                        gGeoManager->LocalToMaster(local, master);
                        str += Form("{%.3f, %.3f, %.3f}", master[0], master[1], master[2]);
                        //Terminal output of coordinate center, point contact, of germanium detectors and their extensions (dx/dy/dz)
                        cout << "Coordinate center x/y/z (abs): = " << master[0] << "/"<< master[1] << "/" << master[2] << endl;
                        if(grooveOnTop == 0) {
                            cout << "Pointcontact coordinates x/y/z (abs): = " << master[0] << "/"<< master[1] << "/" << master[2] - shape->GetDZ() << endl;
                        }
                        if(grooveOnTop == 1) {
                            cout << "Inverted pointcontact coordinates (abs): = " << master[0] << "/"<< master[1] << "/" << master[2] + shape->GetDZ() << endl;
                        }

                        const int MaxArrayLength(5000);
                        // Declaration of leaf types
                        Int_t           eventnumber;
                        /*Int_t           vertex_totnum;
                        Float_t         vertex_xpos[MaxArrayLength];   //[vertex_totnum]
                        Float_t         vertex_ypos[MaxArrayLength];   //[vertex_totnum]
                        Float_t         vertex_zpos[MaxArrayLength];   //[vertex_totnum]
                        Float_t         vertex_time[MaxArrayLength];   //[vertex_totnum]
                        Int_t           vertex_numparticle[MaxArrayLength];   //[vertex_totnum]
                        Int_t           mc_totnumparticles;
                        Int_t           mc_iparticle[MaxArrayLength];   //[mc_totnumparticles]
                        Float_t         mc_ivertex[MaxArrayLength];   //[mc_totnumparticles]
                        Float_t         mc_px[MaxArrayLength];   //[mc_totnumparticles]
                        Float_t         mc_py[MaxArrayLength];   //[mc_totnumparticles]
                        Float_t         mc_pz[MaxArrayLength];   //[mc_totnumparticles]
                        Float_t         mc_pe[MaxArrayLength];   //[mc_totnumparticles]
                        Int_t           mc_id[MaxArrayLength];   //[mc_totnumparticles]
                        Float_t         nuclei_baryonnumber;
                        Float_t         nuclei_charge;*/
                        Int_t           hits_totnum[MaxArrayLength];
                        //Float_t         hits_tote; 
                        Float_t         hits_edep[MaxArrayLength];   //[hits_totnum]
                        Float_t         hits_xpos[MaxArrayLength];   //[MaxArrayLength]
                        Float_t         hits_ypos[MaxArrayLength];   //[hits_totnum]
                        Float_t         hits_zpos[MaxArrayLength];   //[hits_totnum]
                        Int_t           hits_iddet[MaxArrayLength];   //[hits_totnum]
                        //Int_t           hits_idseg[MaxArrayLength];   //[hits_totnum]
                        //Float_t         hits_time[MaxArrayLength];   //[hits_totnum]
                        //Int_t           hits_trackid[MaxArrayLength];   //[hits_totnum]
                        Int_t           hits_trackpdg[MaxArrayLength];   //[hits_totnum]
                        /*Int_t           hits_passivation_totnum;
                        Float_t         hits_passivation_tote;
                        Float_t         hits_passivation_edep[MaxArrayLength];   //[hits_passivation_totnum]
                        Float_t         hits_passivation_xpos[MaxArrayLength];   //[hits_passivation_totnum]
                        Float_t         hits_passivation_ypos[MaxArrayLength];   //[hits_passivation_totnum]
                        Float_t         hits_passivation_zpos[MaxArrayLength];   //[hits_passivation_totnum]
                        Int_t           hits_deadlayer_totnum;
                        Float_t         hits_deadlayer_tote;
                        Float_t         hits_deadlayer_edep[MaxArrayLength];   //[hits_deadlayer_totnum]
                        Float_t         hits_deadlayer_xpos[MaxArrayLength];   //[hits_deadlayer_totnum]
                        Float_t         hits_deadlayer_ypos[MaxArrayLength];   //[hits_deadlayer_totnum]
                        Float_t         hits_deadlayer_zpos[MaxArrayLength];   //[hits_deadlayer_totnum]
                        Int_t           seg_totnum;
                        Int_t           det_totnum;
                        Int_t           seg_id[3];   //[seg_totnum]
                        Int_t           seg_numhits[3];   //[seg_totnum]
                        Float_t         seg_edep[3];   //[seg_totnum]
                        Int_t           det_id[3];   //[det_totnum]
                        Int_t           det_numhits[3];   //[det_totnum]
                        Float_t         det_edep[3];   //[det_totnum]
                        Float_t         ene_in_water;
                        Float_t         ene_in_nitrogen;
                        Float_t         ene_in_scint;
                        Int_t           neutronflag; */

                        // Address of branches
                        fTree->SetBranchAddress("eventnumber",&eventnumber);
                        /*fTree->SetBranchAddress("vertex_totnum",&vertex_totnum);
                        fTree->SetBranchAddress("vertex_xpos",&vertex_xpos);
                        fTree->SetBranchAddress("vertex_ypos",&vertex_ypos);
                        fTree->SetBranchAddress("vertex_zpos",&vertex_zpos); */
                        fTree->SetBranchAddress("hits_totnum",&hits_totnum);
                        fTree->SetBranchAddress("hits_edep",&hits_edep);
                        fTree->SetBranchAddress("hits_xpos",&hits_xpos);
                        fTree->SetBranchAddress("hits_ypos",&hits_ypos);
                        fTree->SetBranchAddress("hits_zpos",&hits_zpos);
                        fTree->SetBranchAddress("hits_trackpdg",&hits_trackpdg);
                        fTree->SetBranchAddress("hits_iddet",&hits_iddet);

                        //fTree->SetBranchAddress("hits_deadlayer_totnum",&hits_deadlayer_totnum);
                        //fTree->SetBranchAddress("hits_deadlayer_edep",&hits_deadlayer_edep);


                        ////-----------WRITE OUTPUT FILES:------------
                        ostringstream condTempStream;
                        ostringstream jsonContentStream;
                        ostringstream siggenContentStream;
                        ostringstream adlContentStream;
                        ostringstream hitsContentStream;
                        ostringstream hitsADLContentStream;

                        ofstream adlSetupOutputFile;

                        string fileContentString = "";
                        string jsonContentString = "";
                        string siggenContentString = "";
                        string adlContentString = "";
                        string hitsContentString = "";
                        string hitsADLContentString = "";

                        //CHECK DETECTOR SUBSET & if textfiles are even wanted:
                        //if(convIsPSSSortedMC[MaGeChannelNumIterator] == 1) {
                        if(true) {
                            // ORDER: First write all dead layer parts, then all active volume parts:
                            //CHECK DEADLAYER is CURRENT VOLUME:
                            if(isCurrentQActiveV == 0) {
                                //BEGEs:
                                if(convIsBEGeSortedMC[MaGeChannelNumIterator] == 1 ) {

                                    //JSON: Not used for now.
                                    //Just console output
                                    /*
                                    if(MaGeChannelNumIterator<10){cout << "For JSON: Writing BEGe Geometry into DET_0" << MaGeChannelNumIterator << ".json" << endl;}
                                    else{cout << "For JSON: Writing BEGe Geometry into DET_" << MaGeChannelNumIterator << ".json" << endl;}
                                    //Flush TempStream
                                    condTempStream.str(string());
                                    //Add jsonOutputPath
                                    condTempStream << jsonOutputPath;
                                    condTempStream << "/DET_";
                                    //Create Output string:
                                    if(MaGeChannelNumIterator<10){condTempStream << "0" << MaGeChannelNumIterator;}else{condTempStream << MaGeChannelNumIterator;}
                                    condTempStream << ".json";
                                    fileContentString = condTempStream.str();
                                    //JSON FILE CREATION:
                                    jsonOutputFile.open(fileContentString);
                                    //Flush TempStream
                                    condTempStream.str(string());
                                    //PREPARE JSON
                                    jsonContentStream.str(string());
                                    jsonContentStream <<
                                    "TEST JSON" << endl;
                                    jsonContentString = jsonContentStream.str();
                                    jsonOutputFile << jsonContentString << endl;
                                    //FILE LEFT OPEN TO BE CLOSED DURING ACTIVE VOLUME READOUT LATER ON
                                    */

                                    //SIGGEN: Config files for siggen
                                    //Just console output
                                    if(MaGeChannelNumIterator<10) {
                                        cout << "For SigGen: Writing BEGe Geometry into DET_0" << MaGeChannelNumIterator << ".config" << endl;
                                    }
                                    else {
                                        cout << "For SigGen: Writing BEGe Geometry into DET_" << MaGeChannelNumIterator << ".config" << endl;
                                    }
                                    condTempStream.str(string());
                                    condTempStream << siggenOutputPath;
                                    condTempStream << "/DET_";
                                    //Create Output string:
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".config";
                                    fileContentString = condTempStream.str();
                                    //SIGGEN CONFIG FILE CREATION: ------------------------- OPEN FOR DEAD LAYERS
                                    siggenOutputFile.open(fileContentString);
                                    condTempStream.str(string());
                                    //Prepare string for siggen config output file creation:
                                    siggenContentStream.str(string());
                                    siggenContentStream <<
                                                        "#Configuration file generated from MaGe geometry by Geoextractor for fieldgen/siggen by David Radford" << endl <<
                                                        "#Format: <key_word> <value>  #comment; all lengths are in mm" << endl << endl <<
                                                        "#" << convDetNameSortedMC[MaGeChannelNumIterator]  << endl <<
                                                        "#Deadlayer volume lengths from center (DX/DY/DZ): " << shape->GetDX()*10. << " / " << shape->GetDY()*10. << " / " << shape->GetDZ()*10. << endl <<
                                                        "#Deadlayer volume center coordinates (X/Y/Z): " << master[0]*10. << " / "<< master[1]*10. << " / " << master[2]*10. << endl <<
                                                        "# general" << endl <<
                                                        "verbosity_level 0        #  0 = terse, 1 = normal, 2 = chatty/verbose" << endl <<
                                                        "# detector geometry" << endl <<
                                                        //2x the height from center of detector and convert to mm from cm: ASSUMING ACTIVE VOLUME SIZE?
                                                        "xtal_length " << shape->GetDZ()*2.*10.       <<  "    # z length" << endl <<
                                                        "xtal_radius " << shape->GetDX()*10.          <<  "    # radius" << endl;
                                    convFCCD[MaGeChannelNumIterator] = shape->GetDX()*10.;
                                    siggenContentString = siggenContentStream.str();
                                    siggenOutputFile << siggenContentString << endl;
                                    //FILE LEFT OPEN TO BE CLOSED DURING ACTIVE VOLUME READOUT LATER ON


                                    //ADL:
                                    //------------WRITE ADL DET_XX.txt for BEGes:
                                    //Just console output
                                    if(MaGeChannelNumIterator<10) {
                                        cout << "For ADL: Writing BEGe Geometry into DET_0" << MaGeChannelNumIterator << ".txt" << endl;
                                    }
                                    else {
                                        cout << "For ADL: Writing BEGe Geometry into DET_" << MaGeChannelNumIterator << ".txt" << endl;
                                    }
                                    condTempStream.str(string());
                                    condTempStream << adlOutputPath << "ConfigFiles/";
                                    condTempStream << "/DET_";
                                    //Create Output string:
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".txt";
                                    fileContentString = condTempStream.str();
                                    adlOutputFile.open(fileContentString);
                                    condTempStream.str(string());
                                    adlContentStream.str(string());

                                    //------------Prepare SETUP_DET_XX.txt:
                                    //Contains no ACTIVE/DEADLAYER info, thus fully written in deadlayer queue.
                                    //Just console output
                                    if(MaGeChannelNumIterator<10) {
                                        cout << "Writing BEGe Geometry into SETUP_DET_0" << MaGeChannelNumIterator << ".txt" << endl;
                                    }
                                    else {
                                        cout << "Writing BEGe Geometry into SETUP_DET_" << MaGeChannelNumIterator << ".txt" << endl;
                                    }
                                    condTempStream.str(string());
                                    condTempStream << adlOutputPath << "ConfigFiles/";
                                    condTempStream << "/SETUP_DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".txt";
                                    fileContentString = condTempStream.str();
                                    adlSetupOutputFile.open(fileContentString);
                                    condTempStream.str(string());
                                    //Prepare string for SETUP_DET_XX.txt
                                    condTempStream <<
                                                   "#ADL SETUP & DET files generated from MaGe geometry by Geoextractor for ADL 3.0" << endl <<
                                                   "#All lengths are in cm" << endl << endl <<
                                                   "#" << convDetNameSortedMC[MaGeChannelNumIterator]  << endl <<
                                                   "#Deadlayer volume lengths from center (DX/DY/DZ): " << shape->GetDX()*10. << " / " << shape->GetDY()*10. << " / " << shape->GetDZ()*10. << endl <<
                                                   "#Deadlayer volume center coordinates (X/Y/Z): " << master[0]*10. << " / "<< master[1]*10. << " / " << master[2]*10. << endl <<
                                                   "ADL_G_VERSION 3.0" << endl <<
                                                   "ADL_G_DEBUG 0" << endl <<
                                                   "#BEGe Geometry File:" << endl <<
                                                   "SIMION_GEOMETRY_BEGE ConfigFiles/DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".txt" << endl <<
                                                   "#BEGe Geometry File:" << endl <<
                                                   "SIMION_SOLVER_INHOMOG ConfigFiles/DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".txt" << endl <<
                                                   "ADL_FIELDS_SIMION	ConfigFiles/DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".txt" << endl <<
                                                   "# How many interactions maximum in an event, how many samples in the traces, how many segments..." << endl <<
                                                   "ADL_EVENT		ConfigFiles/EVENT.txt" << endl <<
                                                   "# Response function for preamp can be convolution with dummy function:" << endl <<
                                                   "# This setup file will invert channel 0 (core)" << endl <<
                                                   "#ADL_CONVL_SAMPLES		ConfigFiles/CONV_SAMPLES.txt" << endl <<
                                                   "#Inserted AGATA Convl Dummy:" << endl <<
                                                   "ADL_CONVL_DUMMY ConfigFiles/Template_CONVL_DUMMY.txt" << endl <<
                                                   "# Drift velocity parameters, see:" << endl <<
                                                   "# Bruyneel et al NIM A 569 (2006) 764-773" << endl <<
                                                   "# Mihailescu et al NIM A 447 (2000) 350-360" << endl <<
                                                   "ADL_DRIFT_GE		ConfigFiles/DRIFT_GE.txt" << endl <<
                                                   "# Choices are binary BIN, text TXT, dino's tkt TKT (write only)" << endl <<
                                                   "ADL_READWRITE		ConfigFiles/READWRITE_TXT.txt" << endl <<
                                                   "ADL_TRACES_NUMRES	ConfigFiles/TRACES_NUMRES.txt"
                                                   "ADL_TRAPPING		ConfigFiles/TRAPPING.txt" << endl;
                                    fileContentString = condTempStream.str();
                                    adlSetupOutputFile << fileContentString << endl;
                                    adlSetupOutputFile.close();
                                }



                                //------------WRITE DET_XX.txt for Coax:
                                if(convIsBEGeSortedMC[MaGeChannelNumIterator] == 0) {
                                    if(MaGeChannelNumIterator<10) {
                                        cout << "Writing Coax Geometry into DET_0" << MaGeChannelNumIterator << ".txt" << endl;
                                    }
                                    else {
                                        cout << "Writing Coax Geometry into DET_" << MaGeChannelNumIterator << ".txt" << endl;
                                    }
                                    //------------WRITE SETUP_DET_XX.txt:
                                    if(MaGeChannelNumIterator<10) {
                                        cout << "Writing Coax Geometry into SETUP_DET_0" << MaGeChannelNumIterator << ".txt" << endl;
                                    }
                                    else {
                                        cout << "Writing Coax Geometry into SETUP_DET_" << MaGeChannelNumIterator << ".txt" << endl;
                                    }
                                }

                              //End of dead layer loop:
                            }
                            //Active volume parts & closing of file
                            //Check active is current volume
                            if( isCurrentQActiveV == 1) {
                                //BEGES:
                                if(convIsBEGeSortedMC[MaGeChannelNumIterator] == 1 ) {

                                    //JSON:
                                    //Just console output
                                    /*
                                    if(MaGeChannelNumIterator<10){cout << "Active Volume For JSON: Writing BEGe Geometry into DET_0" << MaGeChannelNumIterator << ".json" << endl;}
                                    else{cout << "Active Volume For JSON: Writing BEGe Geometry into DET_" << MaGeChannelNumIterator << ".json" << endl;}
                                    //Flush TempStream
                                    condTempStream.str(string());
                                    jsonContentStream <<
                                    "ACTIVE VOLUME STUFF: TEST JSON" << endl <<
                                    "etc pp" << endl;
                                    jsonContentString = jsonContentStream.str();
                                    jsonOutputFile << jsonContentString << endl;
                                    //JSON FILE CLOSE
                                    jsonOutputFile.close();
                                    */

                                    //SIGGEN:
                                    //Just console output
                                    if(MaGeChannelNumIterator<10) {
                                        cout << "Active Volume For SigGen: Writing BEGe Geometry into DET_0" << MaGeChannelNumIterator << ".config" << endl;
                                    }
                                    else {
                                        cout << "Active Volume For SigGen: Writing BEGe Geometry into DET_" << MaGeChannelNumIterator << ".config" << endl;
                                    }
                                    condTempStream.str(string());
                                    //PREPARE STRING FOR SIGGEN CONFIG FILE CONTENT:
                                    siggenContentStream <<
                                                        "#Active volume lengths from center (DX/DY/DZ): " << shape->GetDX()*10. << " / " << shape->GetDY()*10. << " / " << shape->GetDZ()*10. << endl <<
                                                        "#Active volume center coordinates (X/Y/Z): " << master[0]*10. << " / "<< master[1]*10. << " / " << master[2]*10. << endl << endl <<
                                                        "pc_length    3e-4         # point contact length " << endl <<
                                                        "pc_radius    7.5         # point contact radius" << endl <<
                                                        "wrap_around_radius 10.5     # wrap-around radius for BEGes. Set to zero for ORTEC" << endl <<
                                                        "ditch_depth        2     # depth of ditch next to wrap-around for BEGes. Set to zero for ORTEC" << endl <<
                                                        "ditch_thickness    3     # width of ditch next to wrap-around for BEGes. Set to zero for ORTEC" << endl <<
                                                        //Taper from detSettings
                                                        "outer_taper_length " << convTaperLengthSortedMC[MaGeChannelNumIterator] << endl <<
                                                        "outer_taper_width  " << convTaperWidthSortedMC[MaGeChannelNumIterator] << endl <<
                                                        //Additional options useful for Coax:
                                                        //"bottom_taper_length 0    #size of 45-degree taper at bottom of ORTEC-type crystal (for r=z)" << endl <<
                                                        //"bulletize_PC    0        # set to 1 for point contact hemispherical, 0 for cylindrical" << endl <<
                                                        //"hole_length 55 #length of hole for inverted-coax style" << endl <<
                                                        //"hole_radius #radius of hole, for inverted-coax style" << endl <<
                                                        //"taper_angle  10 #taper angle in degrees, for inner or outer taper" << endl <<
                                                        //Calculate Deadlayer thickness:
                                                        //convFCCD[MaGeChannelNumIterator] = convFCCD[MaGeChannelNumIterator] - shape->GetDX()*10.;
                                                        "Li_thickness " << convFCCD[MaGeChannelNumIterator] - shape->GetDX()*10. << " # FCCD boundary for Li contact (not currently used)" << endl <<
                                                        "# configuration for mjd_fieldgen (calculates electric fields & weighing potentials)" << endl <<
                                                        "xtal_grid         0.1    # grid size in mm for field files (usually 0.5 or 0.1 mm)" << endl <<
                                                        //impurity from detSettings
                                                        "impurity_z0       " << convImpurityZSortedMC[MaGeChannelNumIterator] << " # net impurity concentration at Z=0, in 1e10 e/cm3" << endl <<
                                                        "impurity_gradient  " << convImpurityGradSortedMC[MaGeChannelNumIterator]  << " # net impurity gardient, in 1e10 e/cm4" << endl <<
                                                        //High Voltage from detSettings
                                                        "xtal_HV " << convHVSortedMC[MaGeChannelNumIterator] << " # detector bias for fieldgen, in Volts" << endl <<
                                                        "# options for mjd_fieldgen:" << endl <<
                                                        "max_iterations    30000  # maximum number of iterations to use in mjd_fieldgen" << endl <<
                                                        "write_field       1      # 0/1: do_not/do write the standard field output file" << endl <<
                                                        "write_WP          1      # 0/1: do_not/do calculate the weighting potential and write it to the file" << endl <<
                                                        "# file names" << endl <<
                                                        "drift_name ./drift_vel_tcorr.tab    # drift velocity lookup table" << endl;
                                    condTempStream.str(string());
                                    condTempStream << siggenOutputPath;
                                    condTempStream << "ev_DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".dat";
                                    fileContentString = condTempStream.str();
                                    siggenContentStream <<
                                                        "field_name " << fileContentString << " # potential/efield file name; no included spaces allowed" << endl;
                                    condTempStream.str(string());
                                    condTempStream << siggenOutputPath;
                                    condTempStream << "wp_DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".dat";
                                    fileContentString = condTempStream.str();
                                    siggenContentStream <<
                                                        "wp_name " << fileContentString << " # weighting potential file name; no included spaces allowed" << endl <<
                                                        "# configuration for signal calculation" << endl <<
                                                        //Temperature from GSTR
                                                        "xtal_temp         88     # crystal temperature in Kelvin, temperature in GERDA cryostat 185°C = 88K " << endl <<
                                                        "preamp_tau        0     # integration time constant for preamplifier, in ns" << endl <<
                                                        "time_steps_calc   2000   # number of time steps used in calculations" << endl <<
                                                        "step_time_calc    1    # length of time step used for calculation, in ns" << endl <<
                                                        "step_time_out     1   # length of time step for output signal, in ns" << endl <<
                                                        "# nonzero values in the next few lines significantly slows down the code" << endl <<
                                                        "charge_cloud_size 0      # initial FWHM of charge cloud, in mm" << endl <<
                                                        "use_diffusion     0      # set to 0/1 for ignore/add diffusion as the charges drift" << endl;
                                    siggenContentString = siggenContentStream.str();
                                    siggenOutputFile << siggenContentString << endl;
                                    //FINISH SIGGEN CONFIG FILE CREATION: ----------------------- CLOSE AFTER ACTIVE VOLUME
                                    siggenOutputFile.close();

                                    //ADL:
                                    //------------WRITE ADL DET_XX.txt for BEGes:
                                    //Just console output
                                    if(MaGeChannelNumIterator<10) {
                                        cout << "Active Volume For ADL: Writing BEGe Geometry into DET_0" << MaGeChannelNumIterator << ".txt" << endl;
                                    }
                                    else {
                                        cout << "Active Volume For ADL: Writing BEGe Geometry into DET_" << MaGeChannelNumIterator << ".txt" << endl;
                                    }
                                    //Flush TempStream
                                    condTempStream.str(string());
                                    adlContentStream.str(string());
                                    //PREPARE STRING FOR ADL DET FILE CONTENT:
                                    adlContentStream <<
                                                     "BEGE_G_GrooveDepth               0.20 ! Standard Canberra size" << endl <<
                                                     "BEGE_G_GrooveWidth               0.30 ! Standard Canberra size" << endl <<
                                                     "BEGE_G_PointContactRadius        0.75 ! Standard Canberra size" << endl <<
                                                     "BEGE_G_GrooveInnerRadius         0.75 ! usually = POINT_CONTACT_RADIUS" << endl <<
                                                     "BEGE_G_PointContactDepth        -0.01 ! smallest value possible (negative value = one gridsize)" << endl <<
                                                     "BEGE_G_PasLayThickness           0.00 ! " << endl <<
                                                     "BEGE_G_Radius                   " << shape->GetDX()   <<  " !" << endl <<
                                                     "BEGE_G_Height                   " << shape->GetDZ()*2.    <<  " !" << endl <<
                                                     "BEGE_G_SurfaceContactDepth       0.10 ! " << endl <<
                                                     "BEGE_G_ImpTop                   -1.00 ! " << endl <<
                                                     "BEGE_G_ImpBot                   -1.00 ! " << endl <<
                                                     "BEGE_G_Spacing                   0.00 ! No spacing " << endl <<
                                                     "BEGE_G_ExtGroundWidth            0.00 ! External ground potential (if in a grounded container) " << endl <<
                                                     "SIMION_G_GridSize                0.01 ! grid size in cm " << endl <<
                                                     "SIMION_G_EpsScale                16.0 ! epsilon scale relative to epsilon_0 " << endl <<
                                                     "SIMION_G_EpsExtScale              1.0 ! external permittivity (usually 1=vacuum) " << endl <<
                                                     "SIMION_G_Description             Bege ! " << endl <<
                                                     "SIMION_G_Dimension                  2 ! 2 for 2D 3 for 3D " << endl <<
                                                     "SIMION_G_Voltage                " << convHVSortedMC[MaGeChannelNumIterator] << " ! From DetSettings File" << endl <<
                                                     "SIMION_G_RhoScale                 1.0 ! space charge density scale, relative to 10^10/cm^3 " << endl <<
                                                     "SIMION_G_Tol                    1e-10 ! required tolerance for solution " << endl <<
                                                     "SIMION_G_DEBUG                      0 ! print extra information if != 0 " << endl <<
                                                     "ADL_G_SIMION_SmallPot            1e-6 ![V] Small potential, Defines e.g. ADL_In_Detector_SIMION " << endl <<
                                                     //Field Locations
                                                     "ADL_G_Wpot              " << "Fields/DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        adlContentStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        adlContentStream << MaGeChannelNumIterator;
                                    }
                                    adlContentStream << "_Wpot.pa ! Bege_Wpot.pa  Location where the weighting potential is saved " << endl <<
                                                     "ADL_G_Epot              " << "Fields/DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        adlContentStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        adlContentStream << MaGeChannelNumIterator;
                                    }
                                    adlContentStream << "_Epot.pa ! Bege_Epot.pa  Location where the electric potential is saved " << endl <<
                                                     "ADL_G_Stru              " << "Fields/DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        adlContentStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        adlContentStream << MaGeChannelNumIterator;
                                    }
                                    adlContentStream << "_Stru.pa ! Bege_Epot.pa  Location where the structural potential is saved " << endl <<
                                                     "ADL_SCALE_0                         1 !Epot.pa0 (Electrical field space charge) " << endl;
                                    adlContentString = adlContentStream.str();
                                    adlOutputFile << adlContentString << endl;
                                    //FINISH SIGGEN CONFIG FILE CREATION: ----------------------- CLOSE AFTER ACTIVE VOLUME
                                    adlOutputFile.close();
                                }
                            }
                            //////////////////////////////////////////////////////
                            //---------------- HITS ---------------------------- 
                            /////////////////////////////////////////////////////
                            if(useHits == 1 && convIsPSSSortedMC[MaGeChannelNumIterator] == 1) {
                        

                                //------------CUT TO ACTIVE ITERATION:
                                condTempStream.str(string());
                                condTempStream <<
                                               //"det_id==" << MaGeChannelNumIterator << " && " <<
                                               "hits_zpos <= " << master[2] + shape->GetDZ() << " && " <<
                                               "hits_zpos >= " << master[2] - shape->GetDZ() << " && " <<
                                               "((hits_xpos-" << master[0] << ")^2+(hits_ypos-" << master[1] << ")^2) <= " << (shape->GetDX()) * (shape->GetDX());

                                string condTempCut = condTempStream.str();
                                //char const *condCut = condTempCut.c_str();
                                cout << "Active volume fiducial cut: " << condTempCut << endl;
                                //--------------------------------------------------------------------------------------------------------------------
                                ///// READOUT POSITION DATA FROM EACH POINT OF INTEREST (CUT & TRANSFORM)

                                //------------------ACTIVE VOLUME ONLY ATM!!
                                if(isCurrentQActiveV == 1) {
                                    Long64_t nentries = fTree->GetEntries();
                                    //Long64_t nentries = fTree->GetEntries("hits_iddet == 33");

                                    //----------------CREATE FILES FOR HITS OUTPUT
                                    //Commandline output
                                    if(MaGeChannelNumIterator<10) {
                                        cout << "Writing hit output of " << convDetNameSortedMC[MaGeChannelNumIterator] << " into HITS_DET_0" << MaGeChannelNumIterator << ".json" << endl;
                                    }
                                    else {
                                        cout << "Writing hit output of " << convDetNameSortedMC[MaGeChannelNumIterator] << " into HITS_DET_" << MaGeChannelNumIterator << ".json" << endl;
                                    }
                                    //ROOT HITS OUTPUT filepath and name generation:
                                    condTempStream.str(string());
                                    condTempStream << jsonOutputPath;
                                    condTempStream << "/HITS_DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".root";
                                    fileContentString = condTempStream.str();
                                    //ROOT HITS OPEN FILE
                        
                                 //   TFile hfile(fileContentString.c_str(),"recreate");
                                    TFile hfile(fileContentString.c_str(),"RECREATE","Hit positions in cm (0,0,0) at point contact");
                                    // Create a TTree
                                            TTree *hTree = new TTree("fTree","Hit positions in cm (0,0,0) at point contact");


                                            Float_t         cHits_edep;   //[hits_totnum]
                                            Float_t         cHits_xpos;   //[MaxArrayLength]
                                            Float_t         cHits_ypos;   //[hits_totnum]
                                            Float_t         cHits_zpos;   //[hits_totnum]
                                            Int_t           cHits_iddet;   //[hits_totnum]
                                           // Float_t         cHits_time;   //[hits_totnum]
                                           // Int_t           cHits_trackid;   //[hits_totnum]
                                            Int_t           cHits_trackpdg;   //[hits_totnum]    
                                     
                                            hTree->Branch("hits_edep", &cHits_edep, "edep");
                                            hTree->Branch("hits_xpos", &cHits_xpos, "xpos");
                                            hTree->Branch("hits_ypos", &cHits_ypos, "ypos");
                                            hTree->Branch("hits_zpos", &cHits_zpos, "zpos");
                                            hTree->Branch("hits_iddet", &cHits_iddet);
                                            hTree->Branch("hits_trackpdg", &cHits_trackpdg);





                                    //HITS Filepath stream writing:
                                    condTempStream.str(string());
                                    condTempStream << jsonOutputPath;
                                    condTempStream << "HITS_DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".json";
                                    fileContentString = condTempStream.str();
                                    //HITS JSON FILE CREATION with the stream:
                                    hitsOutputFile.open(fileContentString);

                                    //ADLOutputPath generation for file name
                                    condTempStream.str(string());
                                    condTempStream << adlOutputPath;
                                    condTempStream << "HITS_DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".txt";
                                    fileContentString = condTempStream.str();
                                    //ADL HITS FILE PATH:
                                    hitsADLOutputFile.open(fileContentString);
                                    condTempStream.str(string());
                                    hitsContentStream.str(string());
                                    hitsADLContentStream.str(string());
                                    //Write Headers for JSON and ADL HIT files:
                                    hitsContentStream << "{ " << endl
                                                      << "    \"detector\": " << "\"" << convDetNameSortedMC[MaGeChannelNumIterator] << "\"," << endl
                                                      << "    \"magechannel\": " << "\"" << MaGeChannelNumIterator << "\"," << endl;
                                    if(convIsBEGeSortedMC[MaGeChannelNumIterator] == 1) {
                                        hitsContentStream << "    \"type\": " << "\"" << "BEGe" << "\"," << endl;
                                    }
                                    if(convIsBEGeSortedMC[MaGeChannelNumIterator] == 0) {
                                        hitsContentStream << "    \"type\": " << "\"" << "COAX" << "\"," << endl;
                                    }
                                    hitsContentStream << "    \"mageposition\": " << endl
                                                      << "        [" << endl
                                                      << "            { \"detcenter_xpos\": \"" << master[0]*10.
                                                      << "\", \"detcenter_ypos\": \"" << master[1]*10.
                                                      << "\", \"detcenter_zpos\": \"" << master[2]*10.
                                                      << "\" }" << endl << "        ]," << endl
                                                      << "    \"parameters\": " << endl
                                                      << "        [" << endl
                                                      << "            { \"detdimx\": \"" << shape->GetDX()*10.
                                                      << "\", \"detdimy\": \"" << shape->GetDY()*10.
                                                      << "\", \"detdimz\": \"" << shape->GetDZ()*10.
                                                      << "\", \"xtal_radius\": \"" << shape->GetDX()*10.
                                                      << "\", \"xtal_length\": \"" << shape->GetDZ()*2.*10.
                                                      << "\", \"pc_length\": \"" << "3e-4"
                                                      << "\", \"pc_radius\": \"" << "7.5"
                                                      << "\", \"wrap_around_radius\": \"" << "10.5"
                                                      << "\", \"ditch_depth\": \"" << "2"
                                                      << "\", \"ditch_thickness\": \"" << "3"
                                                      << "\", \"outer_taper_length\": \"" << convTaperLengthSortedMC[MaGeChannelNumIterator]
                                                      << "\", \"outer_taper_width\": \"" << convTaperWidthSortedMC[MaGeChannelNumIterator]
                                                      << "\", \"Li_thickness\": \"" << convFCCD[MaGeChannelNumIterator] - shape->GetDX()*10.
                                                      << "\", \"impurity_z0\": \"" << convImpurityZSortedMC[MaGeChannelNumIterator]
                                                      << "\", \"impurity_gradient\": \"" << convImpurityGradSortedMC[MaGeChannelNumIterator]
                                                      << "\", \"xtal_HV\": \"" << convHVSortedMC[MaGeChannelNumIterator]
                                                      << "\" }" << endl << "        ]," << endl
                                                      << "    \"hits\": " << endl
                                                      << "        [" << endl;
                                    //"#Active volume lengths from center (DX/DY/DZ): " << endl << shape->GetDX()*10. << "  " << shape->GetDY()*10. << "  " << shape->GetDZ()*10. << endl;

                                    hitsADLContentStream <<
                                                         "#Transformed HITS for ADL into detector frame (point contact as (0,0,0)) (commented line absolute coordinates: x,y,z, dx,dy,dz) of " 
                                                         << convDetNameSortedMC[MaGeChannelNumIterator] << " by Geoextractor: (X,Y,Z [cm], Energy [keV])" << endl 
                                                         << "#" << "  " << master[0]*10. << "  " << master[1]*10. << "  " << master[2]*10. << "  " << shape->GetDX()*10. << "  " << shape->GetDY()*10. << "  " << shape->GetDZ()*10. << endl;
                                                         //"#Active volume center coordinates (X/Y/Z): " << endl << master[0]*10. << "  "<< master[1]*10. << "  " << master[2]*10. << endl <<
                                                         //"#Active volume lengths from center (DX/DY/DZ): " << endl << shape->GetDX()*10. << "  " << shape->GetDY()*10. << "  " << shape->GetDZ()*10. << endl;

                                    double condPointCZPos = 0.;
                                    double inverseFrame = 0.;
                                    if(grooveOnTop == 0) {
                                        condPointCZPos = (master[2] - shape->GetDZ());
                                        inverseFrame=1.;
                                    }
                                    if(grooveOnTop == 1) {
                                        condPointCZPos = (master[2] + shape->GetDZ());
                                        inverseFrame=-1.;
                                    }
                                    ////---- LOOP TO GET EVERY EVENT THAT MEETS CONDITION
                                    int idcounter = 0;
                                    for(Long64_t i=0; i<nentries && idcounter<maxHits; i++) {
                                        fTree->GetEntry(i);
                                        //Spatial cut conditions:
                                        if(*hits_zpos <= (master[2] + shape->GetDZ())  &&
                                           *hits_zpos >= (master[2] - shape->GetDZ())  &&
                                           (*hits_xpos - master[0])*(*hits_xpos - master[0]) + (*hits_ypos - master[1])*(*hits_ypos - master[1]) <= (shape->GetDX())*(shape->GetDX())  ) {
                                            //ROOT FILE OUTPUT FOR HIT POINTS &HITS& Point contact as (0,0,0):
                                             

                                            cHits_edep = *hits_edep*1000.  ; 
                                            cHits_xpos = inverseFrame*(*hits_xpos - master[0])*10.;
                                            cHits_ypos = inverseFrame*(*hits_ypos - master[1])*10.;
                                            cHits_zpos = inverseFrame*(*hits_zpos - condPointCZPos)*10.;
                                            cHits_iddet = *hits_iddet;  
                                           
                                            //cout << "cHits_xpos: " << cHits_xpos << " &cHits_xpos: "<< &cHits_xpos << " hits_xpos: " << hits_xpos << " *hits_xpos: " << *hits_xpos << endl;
                                            hTree->Fill();
                                           



                                            //JSON Hit output into files: Point contact as (0,0,0):
                                            hitsContentStream <<
                                                              "            { \"id\": \"" << idcounter << "\", \"eventnumber\": \"" << i
                                                              << "\", \"xpos\": \"" << inverseFrame*(*hits_xpos - master[0])*10.
                                                              << "\", \"ypos\": \"" << inverseFrame*(*hits_ypos - master[1])*10.
                                                              << "\", \"zpos\": \"" << inverseFrame*(*hits_zpos - condPointCZPos)*10.
                                                              << "\", \"edep\": \"" << *hits_edep*1000.
                                                              << "\", \"trackpdg\": \"" << *hits_trackpdg
                                                              << "\" }," << endl;
                                            //Terminal output
                                            cout <<
                                                 inverseFrame*(*hits_xpos - master[0]) << "  " <<
                                                 inverseFrame*(*hits_ypos - master[1]) << "  " <<
                                                 inverseFrame*(*hits_zpos - condPointCZPos) << "  " <<
                                                 *hits_edep*1000. <<
                                                 //"hits_tote: " << *hits_tote << ") " <<
                                                 //"hits_totnum: " << *hits_totnum << "} " <<
                                                 //"hits_time: " << *hits_time << "} " <<
                                                 //"hits_trackpdg: " << *hits_trackpdg << "} " <<
                                                 //"hits_trackid: " << *hits_trackid << "} " <<
                                                 " hits_iddet: " << *hits_iddet << endl;
                                            //HITS Files for ADL output files:   
                                            hitsADLContentStream <<
                                                                 //*hits_xpos*10. << "  " << *hits_ypos*10. << "  " << *hits_zpos*10. << "  " << *hits_edep*1000. << endl;   //absolute frame
                                                                 inverseFrame*(*hits_xpos - master[0]) << "  " << inverseFrame*(*hits_ypos - master[1]) << "  " << inverseFrame*(*hits_zpos - condPointCZPos) << "  " << *hits_edep*1000. << endl; //point contact as (0,0,0)
                                            idcounter++;

                                        }



                                        //END OF ALL EVENTS FOR CURRENT DETECTOR
                                    }   
                                    hitsContentStream.seekp(-2,std::ios_base::end);
                                    hitsContentStream << "" << endl << "        ]" << endl << "}" << endl;
                                    hitsContentString = hitsContentStream.str();
                                    hitsADLContentString = hitsADLContentStream.str();
                                    hitsOutputFile << hitsContentString << endl;
                                    hitsADLOutputFile << hitsADLContentString << endl;
                                    //Hit Files close
                                    hfile.Write();
                                    hfile.Close();

                                    //ROOT HITS TCHAIN ADD
                                    condTempStream.str(string());
                                    condTempStream << jsonOutputPath;
                                    condTempStream << "/HITS_DET_";
                                    if(MaGeChannelNumIterator<10) {
                                        condTempStream << "0" << MaGeChannelNumIterator;
                                    }
                                    else {
                                        condTempStream << MaGeChannelNumIterator;
                                    }
                                    condTempStream << ".root";
                                    fileContentString = condTempStream.str();
                                    aTree->Add(fileContentString.c_str());
                                    if(MaGeChannelNumIterator+1 == convDetNameSortedMC.size()){
                                    //ALL FILE CREATION:
                                    condTempStream.str(string());
                                    condTempStream << jsonOutputPath;
                                    condTempStream << "HITS_DET_ALL.root";
                                    fileContentString = condTempStream.str();
                                    TFile afile(fileContentString.c_str(),"RECREATE","Hit positions in cm (0,0,0) at point contact");
                                    aTree->CloneTree(-1,"fast");                                    
                                    afile.Write();
                                    afile.Close();                                    
                                    }
                                                                        

                                    hitsOutputFile.close();
                                    hitsADLOutputFile.close();

                                }



                                //-------------------- 3D SCATTER PLOTS TO CHECK GEOMETRY ------------- Could be used for detectorwise Visualisation
/*
                                if(isVisHitsEnabled == 1){
                                    //TApplication *app = new TApplication("App", &argc, argv )
                                    TApplication *app = new TApplication("App", 0, 0 );  
                                    //TApplication *app = new TApplication("App", appargc, appargv );
                                    //cout << "3D SCATTER PLOT OF DET/CHAN: " << convDetNameSortedMC[MaGeChannelNumIterator] << "/" << MaGeChannelNumIterator << endl;
                                    //char const *condCanvasName = condTempCanvasName.c_str();
                                    TCanvas *cX = new TCanvas("Visualised Hits");
                                    fTree->SetMarkerStyle(6);

                                    fTree->Draw("hits_zpos:hits_ypos:hits_xpos",visHitsCond.c_str());
                                    cX->Update();
                                    app->Run();
                                    
                                    if(isCurrentQActiveV == 0){
                                    //----------------------------------------- Hit position 3D plot: ------------------------------------
                                    //CANVAS NAMING ITERATION for 3D PLOT:
                                    condTempStream.str(string());
                                    condTempStream << "DET/MaGeChan: " << convDetNameSortedMC[MaGeChannelNumIterator] << "/" <<MaGeChannelNumIterator;
                                    string condCanvasString = condTempStream.str();
                                    string condTempCanvasName = condTempStream.str();
                                    char const *condCanvasName = condTempCanvasName.c_str();
                                    TCanvas *cX = new TCanvas(condCanvasName);
                                    }

                                    //accurate, no hits outside -- DOESNT WORK FOR HIGHER ORDERS
                                    //fTree->SetMarkerColor(kGreen);
                                    //fTree->Draw("hits_zpos:hits_ypos:hits_xpos","hits_idseg == 1");

                                    //NOT accurate
                                    //fTree->Draw("hits_zpos:hits_ypos:hits_xpos","hits_iddet == 1");


                                    //fTree->Draw("hits_zpos:hits_ypos:hits_xpos",condDet_id);
                                    //DEADLAYER
                                    //fTree->SetMarkerColor(kBlue);
                                    //fTree->Draw("hits_deadlayer_zpos:hits_deadlayer_ypos:hits_deadlayer_xpos","hits_iddet == 1","same");
                                    //fTree->Draw("hits_deadlayer_zpos:hits_deadlayer_ypos:hits_deadlayer_xpos",condDet_id,"same");

                                    //DEAD LAYER -RED
                                    if(isCurrentQActiveV == 0){

                                    //fTree->SetMarkerColor(kRed);
                                    //fTree->SetMarkerSize(6);
                                    //fTree->Draw("hits_zpos:hits_ypos:hits_xpos",condCut);
                                    //fTree->Draw("hits_zpos:hits_ypos:hits_xpos",condDet_id);


                                    condTree->SetMarkerColor(kRed);
                                    condTree->SetMarkerStyle(6);
                                    //Draw Deadlayer
                                    condTree->Draw("hits_zpos:hits_ypos:hits_xpos");
                                    condTree->SetMarkerColor(kBlue);
                                    condTree->Draw("hits_zpos:hits_ypos:hits_xpos",condCut,"same");

                                    cout << "**********************************************************************DEADLAYER**************** " << endl;
                                    }

                                    //ACTIVE VOLUME -GREEN
                                    if(isCurrentQActiveV == 1){
                                    //fTree->SetMarkerColor(kGreen);
                                    //fTree->SetMarkerSize(6);
                                    //fTree->Draw("hits_zpos:hits_ypos:hits_xpos",condCut,"same");

                                    condTree->SetMarkerColor(kRed);
                                    condTree->SetMarkerStyle(6);
                                    condTree->SetMarkerColor(kGreen);
                                    condTree->Draw("hits_zpos:hits_ypos:hits_xpos","","same");
                                    condTree->SetMarkerColor(kPink);
                                    condTree->Draw("hits_zpos:hits_ypos:hits_xpos",condCut,"same");


                                    cout << "**********************************************************************ACTIVE VOLUME**************** " << endl;
                                    } 
                                } */


                                //-------------------------------- 3D SCATTER PLOTS FINISHED ---------------

                                


                            }
                        }
                        /////////////////////////////////////////////////////////////////////////////////////
                        //--------------------------------------------------END HITS USAGE------------------------------------------------------
                        /////////////////////////////////////////////////////////////////////////////////////


                        //Check when to move to next Det Channel
                        if(isCurrentQActiveV == 1) {
                            MaGeChannelNumIterator++;
                            cout << "----------------------------------------------------";
                        };
                        cout << "" << endl;
                    }
                }
            }
        }
        //Iteration to the next Node by calling print_node function again:
        TIter next(node->GetNodes());

        while (TGeoNode *n = (TGeoNode *)next())
        {
            print_node(n, path, MaGeChannelNumIterator, convDetNameSortedMC, convIsGrooveOnTopSortedMC, convIsBEGeSortedMC, convIsPSSSortedMC, convTaperLengthSortedMC, convTaperWidthSortedMC, convImpurityZSortedMC, convImpurityGradSortedMC, convHVSortedMC, convFCCD, jsonOutputPath, adlOutputPath, siggenOutputPath, jsonOutputFile, siggenOutputFile, adlOutputFile, hitsOutputFile, hitsADLOutputFile, fTree, aTree, maxHits, useHits, isVisHitsEnabled, visHitsCond);
        }

    }
}

//Function to be used directly with "root Geoextractor.cxx" terminal command without need to compile. Main function also basically runs this function.
void Geoextractor(int checkGDML = 0, string InputPath = "", string jsonOutputPath = "./Output/json/", string adlOutputPath = "./Output/adl/", string siggenOutputPath = "./Output/siggen/", string gdmlFilePath = "defaultgdml.gdml", string convFilePath = "detsettings.txt", int maxHits = 1000000000, int useHits = 1, int isVisHitsEnabled = 0, string visHitsCond = "") {

    cout << "Geoextractor running..." << endl;

    /////////////////////////////////////////////////////////////////////////
    ////// JOIN ROOT FILES to TChain named fTree IN CHOSEN DIRECTORY ($MC_DIR) TO BE ANALYSED
    /////////////////////////////////////////////////////////////////////////
    TChain* aTree = new TChain("fTree");    
    TChain* fTree = new TChain("fTree");
    string MC_DIR = "";
    if(useHits == 1)
    {
        if(InputPath == "")
        {
            cout << "No -m,-mcdir path given, looking into directory defined by shell variable MC_DIR... " << endl;
            MC_DIR = getenv("MC_DIR");
            if( MC_DIR == "" )
            {
                cout << "Please insert the MC directory where there root files of your simulation are located. " << endl;
                cin >> MC_DIR;
            }
        } else {
            MC_DIR = InputPath; 
        }
        // read in the output of simulation and chain the files
        
        string pattern = "*.root";
        vector<string> flist;
        void *dir = gSystem->OpenDirectory( MC_DIR.c_str() );
        if (!dir)
        {
            cout << "Couldn't open directory '" << MC_DIR << "'" << std::endl;
            vector<string> MC_FileList = flist;
        }

        const Char_t * c_fileForList;
        TRegexp regexp(pattern.c_str(), true);

        while ( (c_fileForList = gSystem->GetDirEntry(dir)) )
        {
            string fileForList = c_fileForList;
            if ( (TString (fileForList)).Index(regexp) != kNPOS ) flist.push_back( fileForList );
        }
        gSystem->FreeDirectory(dir);
        sort( flist.begin(),flist.end() );
        vector<string> MC_FileList = flist;


        cout << "Adding all Monte Carlo .ROOT files from MC_DIR directory: " << endl;

        for(unsigned int ifile = 0; ifile < MC_FileList.size(); ifile++)
        {
            fTree->Add( Form("%s/%s", MC_DIR.c_str(), MC_FileList.at( ifile ).c_str()) );
            cout << "\t " << MC_DIR << "/" << MC_FileList.at( ifile ) << endl;
        }

        //-------------------- 3D SCATTER PLOTS TO CHECK GEOMETRY ------------- Used in VisHits Option
        if(isVisHitsEnabled == 1){  
          TApplication app("appKey",0,0);                             
          fTree->SetMarkerStyle(6);
    
          fTree->Draw("hits_zpos:hits_ypos:hits_xpos",visHitsCond.c_str());
          app.Run();

       }










    }





    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////
    /// TGEOMANAGER START:
    //////////////////////////////////////////////////////
    // Input: GDML
    char const *gdmlFilePathcc = gdmlFilePath.c_str();
    cout << endl << "Reading in .GDML file from path: " << gdmlFilePath << endl;
    //gSystem->Load("libGeom");
    //TGeoManager *geom = new TGeoManager("simple1", "Simple geometry");
    //gGeoManager = new TGeoManager("Geometry", "default geometry");
    TGeoManager::Import(gdmlFilePathcc);

    //---------// DRAWING options; Can be used to check GDML geometry
    if ( checkGDML == 1 ) 
    {
        gGeoManager->SetVisLevel(10);
        gGeoManager->GetTopVolume()->SetVisContainers(kTRUE);
        gGeoManager->GetTopVolume()->Draw("ogl");
    }
  
    //////////////////////////////////////////////////////////////////////////////
    //------------------------------------------DetSettings & Det Channels CONVERSION:
    ///////////////////////////////////////////////////////////////////////////////
    // Input: DetSettings
    //Organised to detector name: DetName (string), DataChannel (int), MaGeChannel (int), GrooveOnTop? (int), COAX=0/BEGe=1 (int), Should be analysed for PSS?
    cout << endl <<"Reading in DetectorSettings from path: " << convFilePath << endl;
    cout << "DetectorSettings structure: Extract Hits=1, DetName, DataChannel, MaGeChannel, GrooveOnTop=1, COAX=0/BEGe=1, TaperLength, TaperWidth, ImpurityZ0, ImpurityZGradient, HV" << endl;

    vector<int> convIsPSSSortedData;
    vector<int> convIsPSSSortedMC;
    vector<string> convDetNameSortedData;
    vector<string> convDetNameSortedMC;
    vector<int> convDataChan;
    vector<unsigned int> convMCChan;
    vector<int> convIsGrooveOnTopSortedData;
    vector<int> convIsGrooveOnTopSortedMC;
    vector<int> convIsBEGeSortedData;
    vector<int> convIsBEGeSortedMC;
    vector<double> convTaperLengthSortedData;
    vector<double> convTaperLengthSortedMC;
    vector<double> convTaperWidthSortedData;
    vector<double> convTaperWidthSortedMC;
    vector<double> convImpurityZSortedData;
    vector<double> convImpurityZSortedMC;
    vector<double> convImpurityGradSortedData;
    vector<double> convImpurityGradSortedMC;
    vector<int> convHVSortedData;
    vector<int> convHVSortedMC;

    vector<double> convFCCD;


    int debugInt = -1;

    convDataChan.push_back(debugInt);           //DataChannel
    convMCChan.push_back(debugInt);              //MaGeChannel

    convIsPSSSortedData.push_back(debugInt); //ExtractHits?
    convDetNameSortedData.push_back("DEBUG");   //DetName
    convIsGrooveOnTopSortedData.push_back(debugInt); //GrooveOnTop?
    convIsBEGeSortedData.push_back(debugInt); //COAX=0/BEGe=1

    convTaperLengthSortedData.push_back(debugInt); //TaperLength
    convTaperWidthSortedData.push_back(debugInt);      //TaperWidth
    convImpurityZSortedData.push_back(debugInt);       //Impurity Z0
    convImpurityGradSortedData.push_back(debugInt);    //Impurity Gradient
    convHVSortedData.push_back(debugInt); //HV

    convIsPSSSortedMC.push_back(debugInt);
    convDetNameSortedMC.push_back("DEBUG");
    convIsGrooveOnTopSortedMC.push_back(debugInt);
    convIsBEGeSortedMC.push_back(debugInt);   
    convTaperLengthSortedMC.push_back(debugInt);
    convTaperWidthSortedMC.push_back(debugInt);
    convImpurityZSortedMC.push_back(debugInt);
    convImpurityGradSortedMC.push_back(debugInt);
    convHVSortedMC.push_back(debugInt);

    convFCCD.push_back(debugInt);

    unsigned int convIterator = 0;
    ifstream file(convFilePath, std::ios::in);
    //Generate Vectors
    if( !file )
        cerr << "Can't open Detector Channel File!" << endl;

    string commentLine;
    getline(file, commentLine);
    while( file >> convIsPSSSortedData[convIterator] //ExtractHits?
            >> convDetNameSortedData[convIterator]   //DetName
            >> convDataChan[convIterator]           //DataChannel
            >> convMCChan[convIterator]              //MaGeChannel
            >> convIsGrooveOnTopSortedData[convIterator] //GrooveOnTop?
            >> convIsBEGeSortedData[convIterator] //COAX=0/BEGe=1            
            >> convTaperLengthSortedData[convIterator]
            >> convTaperWidthSortedData[convIterator]
            >> convImpurityZSortedData[convIterator]
            >> convImpurityGradSortedData[convIterator]
            >> convHVSortedData[convIterator] //HV
         )
    {
        if(convMCChan.size() <= convMCChan[convIterator]) {     //Check if vector is at least big enough to accomodate entries, otherwise resize
            convDataChan.resize(convMCChan[convIterator]+1);
            convMCChan.resize(convMCChan[convIterator]+1);

            convDetNameSortedData.resize(convMCChan[convIterator]+1);
            convIsGrooveOnTopSortedData.resize(convMCChan[convIterator]+1);
            convIsBEGeSortedData.resize(convMCChan[convIterator]+1);
            convIsPSSSortedData.resize(convMCChan[convIterator]+1);
            convTaperLengthSortedData.resize(convMCChan[convIterator]+1);
            convTaperWidthSortedData.resize(convMCChan[convIterator]+1);
            convImpurityZSortedData.resize(convMCChan[convIterator]+1);
            convImpurityGradSortedData.resize(convMCChan[convIterator]+1);
            convHVSortedData.resize(convMCChan[convIterator]+1);

            convDetNameSortedMC.resize(convMCChan[convIterator]+1);
            convIsGrooveOnTopSortedMC.resize(convMCChan[convIterator]+1);
            convIsBEGeSortedMC.resize(convMCChan[convIterator]+1);
            convIsPSSSortedMC.resize(convMCChan[convIterator]+1);
            convTaperLengthSortedMC.resize(convMCChan[convIterator]+1);
            convTaperWidthSortedMC.resize(convMCChan[convIterator]+1);
            convImpurityZSortedMC.resize(convMCChan[convIterator]+1);
            convImpurityGradSortedMC.resize(convMCChan[convIterator]+1);
            convHVSortedMC.resize(convMCChan[convIterator]+1);

            convFCCD.resize(convMCChan[convIterator]+1);

        }
        if(convMCChan.size() <= convIterator) {            //Pushback increase vector size to maximum number of detectors
            convDataChan.push_back(debugInt);           //DataChannel
            convMCChan.push_back(debugInt);              //MaGeChannel

            convDetNameSortedData.push_back("DEBUG");   //DetName
            convIsGrooveOnTopSortedData.push_back(debugInt); //GrooveOnTop?
            convIsBEGeSortedData.push_back(debugInt); //COAX=0/BEGe=1
            convIsPSSSortedData.push_back(debugInt); //PSS?
            convTaperLengthSortedData.push_back(debugInt);
            convTaperWidthSortedData.push_back(debugInt);
            convImpurityZSortedData.push_back(debugInt);
            convImpurityGradSortedData.push_back(debugInt);
            convHVSortedData.push_back(debugInt); //HV

            convDetNameSortedMC.push_back("DEBUG");
            convIsGrooveOnTopSortedMC.push_back(debugInt);
            convIsBEGeSortedMC.push_back(debugInt);
            convIsPSSSortedMC.push_back(debugInt);
            convTaperLengthSortedMC.push_back(debugInt);
            convTaperWidthSortedMC.push_back(debugInt);
            convImpurityZSortedMC.push_back(debugInt);
            convImpurityGradSortedMC.push_back(debugInt);
            convHVSortedMC.push_back(debugInt);

            convFCCD.push_back(debugInt);

        }

        convDetNameSortedMC[convMCChan[convIterator]] = convDetNameSortedData[convIterator];    //Fill up MaGe sorted vectors
        convIsGrooveOnTopSortedMC[convMCChan[convIterator]] = convIsGrooveOnTopSortedData[convIterator];
        convIsBEGeSortedMC[convMCChan[convIterator]] = convIsBEGeSortedData[convIterator];
        convIsPSSSortedMC[convMCChan[convIterator]] = convIsPSSSortedData[convIterator];
        convTaperLengthSortedMC[convMCChan[convIterator]] = convTaperLengthSortedData[convIterator];
        convTaperWidthSortedMC[convMCChan[convIterator]] = convTaperWidthSortedData[convIterator];
        convImpurityZSortedMC[convMCChan[convIterator]] = convImpurityZSortedData[convIterator];
        convImpurityGradSortedMC[convMCChan[convIterator]] = convImpurityGradSortedData[convIterator];
        convHVSortedMC[convMCChan[convIterator]] = convHVSortedData[convIterator];
        cout << convIsPSSSortedData[convIterator]
             << " " << convDetNameSortedData[convIterator]
             << " " << convDataChan[convIterator]
             << " " << convMCChan[convIterator]
             << " " << convIsGrooveOnTopSortedData[convIterator]
             << " " << convIsBEGeSortedData[convIterator]
             << " " << convTaperLengthSortedData[convIterator]
             << " " << convTaperWidthSortedData[convIterator]
             << " " << convImpurityZSortedData[convIterator]
             << " " << convImpurityGradSortedData[convIterator]
             << " " << convHVSortedData[convIterator]
             //             << " " << convDetNameSortedMC[convIterator]
             //             << " " << convIsGrooveOnTopSortedMC[convIterator]
             << endl;


        convIterator++;
    }

    file.close();
    if(InputPath != "") 
    {
        cout << "Input path for ROOT files to be merged for extraction: " << InputPath << endl;
    }
    ////////// Output PATHS: JSON Output, ADL Output, SigGen Output
    /////JSON:
    cout << endl << "JSON output path: " << jsonOutputPath << endl;

    /////ADL:
    cout << "ADL output path: " << adlOutputPath << endl;

    /////SIGGEN:
    cout << "SIGGEN output path: " << siggenOutputPath << endl << endl;
    cout << "---------------------------------------------------------------------" <<  endl; 
    cout << "Extracted Detector geometries: " << endl;
    cout << "---------------------------------------------------------------------" << endl << endl;

    //Textfiles to be written inside Output documents by the loops
    ofstream jsonOutputFile;
    ofstream siggenOutputFile;
    ofstream adlOutputFile;
    ofstream hitsOutputFile;
    ofstream hitsADLOutputFile;



    clock_t start, end;
    start = clock();
    
    //Set counter of channel numbers outside of loop to increase with every iteration
    unsigned int MaGeChannelNumIterator = 0;
    print_node(gGeoManager->GetTopNode(), "", MaGeChannelNumIterator, convDetNameSortedMC, convIsGrooveOnTopSortedMC, convIsBEGeSortedMC, convIsPSSSortedMC,
               convTaperLengthSortedMC, convTaperWidthSortedMC, convImpurityZSortedMC, convImpurityGradSortedMC, convHVSortedMC, convFCCD, jsonOutputPath, adlOutputPath, siggenOutputPath, jsonOutputFile, siggenOutputFile, adlOutputFile, hitsOutputFile, hitsADLOutputFile, fTree, aTree, maxHits, useHits, isVisHitsEnabled, visHitsCond);

    end=clock();
    cout << "Time required for execution: "
         << (double)(end-start)/CLOCKS_PER_SEC
         << " seconds." << "\n\n";
}


int main(int argc, char* argv[])
{
    //Define default values to be passed to Geoextractor function if unchanged by arguments.
    string InputPath = "";
    string jsonOutputPath = "./Output/json/"; 
    string adlOutputPath = "./Output/adl/";
    string siggenOutputPath = "./Output/siggen/";
    string gdmlFilePath = "defaultgdml.gdml";
    string convFilePath = "detsettings.txt";
    string visHitsCond = "";
    //Maximum of number of hits default value: 100 billion
    int maxHits = 1000000000;
    int useHits = 1;
    int checkGDML = 0;
    int isVisHitsEnabled = 0;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            cerr << "Options:\n"
                 << "\t-h,--help\t\tShow help message\n"
                 << "\t-m,--mcdir PATH\t\tSpecify input path for ROOT files to be merged and used for extraction. Default: MC_DIR shell variable." << endl
                 << "\t-a,--adl PATH\t\tSpecify output path for ADL files. Default: ./Output/adl" << endl
                 << "\t-j,--json PATH\t\tSpecify output path for JSON files. Default: ./Output/json" << endl
                 << "\t-s,--siggen PATH\tSpecify output path for SIGGEN files. Default: ./Output/siggen" << endl
                 << "\t-g,--gdml PATH\t\tSpecify path to GDML file. Default: ./defaultgdml.gdml" << endl
                 << "\t-d,--detsettings PATH\tSpecify path to Detector Settings file. Default: ./detsettings.txt" << endl
                 << "\t-x,--maxhits INT\tLimit amount of hits per detector to be extracted and written to file to a fixed number. Default: No limit." << endl
                 << "\t-v,--vishits STRING\tVisualise hit extraction with condition set by STRING." << endl
                 << "\t-n,--nohits \t\tSkip all hit extraction, only extract geometry." << endl
                 << "\t-c,--checkgdml \t\tVisualise .GDML file with OGL drivers. Best to use: root Geoextractor.cxx+(1)" << endl;
            return 0;
        } else if ((arg == "-a") || (arg == "--adl")) {
            if (i + 1 < argc) {
                i++;
                adlOutputPath = argv[i];
            } 
        } else if ((arg == "-j") || (arg == "--json")) {
            if (i + 1 < argc) {
                i++;
                jsonOutputPath = argv[i];
            } 
        } else if ((arg == "-s") || (arg == "--siggen")) {
            if (i + 1 < argc) { 
                i++;
                siggenOutputPath = argv[i];
            } 
        }  else if ((arg == "-g") || (arg == "--gdml")) {
            if (i + 1 < argc) { 
                i++;
                cout << "GDML path defined as: " << argv[i] << endl;
                gdmlFilePath = argv[i];
            }
        }  else if ((arg == "-d") || (arg == "--detsettings")) {
            if (i + 1 < argc) { 
                i++;
                cout << "Detector settings file path defined as: " << argv[i] << endl;
                convFilePath = argv[i];
            } 
        }  else if ((arg == "-m") || (arg == "--mcdir")) {
            if (i + 1 < argc) { 
                i++;
                InputPath = argv[i];
                cout << "MC Input path defined by -m,-mc command as: " << argv[i] << endl;
            } 
        }  else if ((arg == "-x") || (arg == "--maxhits")) {
            if (i + 1 < argc) { 
                i++;
                maxHits = stoi( argv[i] );
                cout << "Maximum number of hits per detector chosen to be: " << argv[i] << endl;
            } 
        }  else if ((arg == "-v") || (arg == "--vishits")) {
            isVisHitsEnabled = 1;
            cout << "Visualisation of hits enabled. " << endl;
            if (i + 1 < argc) { 
                i++;
                visHitsCond = argv[i];
                cout << "Selection cut of hits: " << argv[i] << endl;
            }            
        }  else if ((arg == "-n") || (arg == "--nohits")) {
             useHits = 0;
             cout << "--nohits enabled: Hits are not extracted or written to file. Geometries/Parameters are still extracted & written." << endl << endl;
            
        } else if ((arg == "-c") || (arg == "--checkgdml")) {
             checkGDML = 1;
             cout << "--checkgdml enabled: Visualising .GDML file with OGL drivers. Best use: root Geoextractor.cxx+(1)" << endl << endl;
            
        } 
    }
    
    cout << "Running of compiled Geoextractor:" << endl;
    Geoextractor(checkGDML, InputPath, jsonOutputPath, adlOutputPath, siggenOutputPath, gdmlFilePath, convFilePath, maxHits, useHits, isVisHitsEnabled, visHitsCond);
    cout << "End of Geoextractor reached." << endl;
}




