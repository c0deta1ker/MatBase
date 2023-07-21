# ADRESSTools: A suite of data analysis tools for photoelectron spectroscopy at ADRESS

A compilation of MATLAB scripts used for the analysis of soft X-ray angle-resolved photoemission spectroscopy (SX-ARPES) experiments that give direct access to the electronic band-structure of a material. Designed to be directly compatible with the data format of SX-ARPES experiments at the ADRESS beamline, at the Swiss Light Source (SLS) in the Paul Scherrer Institute (PSI), but can be generalised to other data formats if required. The scripts can also be used for perform x-ray photoelectron spectroscopy (PES) curve fitting, using a variety of background shapes.

## Installation  
1. Download the *ADRESSTools* repository.
2. Open MATLAB and use *Set Path* in the *Home* tab to add the *ADRESSTools* repository and all its sub-folders into its saved search paths.
3. Make sure you also use *Set Path* to add the repository / folder that contains all of your data to be loaded in.
4. In the *PESTools* folder, there is an *Examples* folder, which contains many ARPES / XPS data processing and curve fitting templates / examples that can be used. You should use this as a starting point.

## Materials Database Tools
**(1) Material Properties Database (MPD)**: This is a local MATLAB database that compiles the most useful physical, electronic, optical and magnetic material properties of elements / compounds. The data is taken from a range of sources, where the 'average' values are used for parameters that had more than 1 unique value. Accessible in MATLAB via 'get_mpd_props()'.
![000_MPD](ADRESSTools_vEos/PESTools_PCC/0_ReadMeImages/000_MPD.png)

**(2) Electron Inelastic Mean Free Path Database (eIMFPD)**: This is a local MATLAB database that compiles the optical data from the NIST Electron Inelastic-Mean-Free-Path Database (http://dx.doi.org/10.18434/T48C78) so that the results can be easily called and accessed within MATLAB. This is accessible in MATLAB via 'get_eimfpd_props()'. Furthermore, predictive eIMFP formulas's are also available using the (1) Universal, (2) TPP-2M and (3) S1 & S2 formalisms. 
![000_eIMFPD](ADRESSTools_vEos/PESTools_PCC/0_ReadMeImages/000_eIMFPD.png)

**(3) Photoionisation Cross-Section and Asymmetry Database (PIXSAD)**:
This is a local MATLAB database that compiles the photoionisation cross-section and asymmetry parameter data of the elements. This is useful when modelling the total photoelectron intensity that originates from a given layer of a sample when performing ARPES / XPS experiments. Accessible in MATLAB via 'get_pixsad_props()'.
![000_PIXSAD](ADRESSTools_vEos/PESTools_PCC/0_ReadMeImages/000_PIXSAD.png)

**(4) Photoionisation Energy and Fluorescence Database (PIEFD)**:
This is a local MATLAB database that compiles the photoionisation energy and fluorescence yield data of the elements. Accessible in MATLAB via 'get_piefd_props()'.
![000_PIEFD](ADRESSTools_vEos/PESTools_PCC/0_ReadMeImages/000_PIEFD.png)

**(5) Physics Constants Database**:
This is a MATLAB function that loads in many physics constants that can be used for data processing, or modelling of physical systems. Accessible in MATLAB via 'physics_constants()'.


## MATLAB Apps
**(1) IMFP Calculator App**: A simple app that allows the user to calculate the electron inelastic mean free path (IMFP) by selecting a predefined element or compound, or by manually entering the known material parameters. The IMFP is calculated using the (i) Universal, (ii) TPP-2M and (iii) S1 and (iv) S2 formalisms. Accessible in MATLAB via 'IMFP_UI'.
![IMFP_UI_Snapshot](ADRESSTools_vEos/MaterialsDatabase_PCC/0_ReadMeImages/IMFP_UI_Snapshot.png)

**(2) Photoionisation Binding Energy, Cross-Section and Asymmetry App**: A simple app that allows the user to plot the (i) binding energy spectrum, (ii) photoionisation cross-sections and (iii) photoionisation asymmetry parameters of a selected element. This is useful for on-the-fly experiments. Accessible in MATLAB via 'PES_UI()'.
![PES_UI_Snapshot](ADRESSTools_vEos/MaterialsDatabase_PCC/0_ReadMeImages/PES_UI_Snapshot.png)



**ADRESS Job File Generators**:
Allows the user to generate custom job-files in the form of a text-file that is compatible to use in the 'Restore Tab' button on 'SmartGUI' at the ADRESS beamline. This allows the user to very quickly create a list of ARPES / PES scans to be taken as a function of photon energy, time or position, where each scan can be referenced to a known Fermi-edge, if desired. List of available job file generators:  
- 'ADRESS_JobFileGenerator_CLvsTIME()'  
- 'ADRESS_JobFileGenerator_CLvsHV()'  
- 'ADRESS_JobFileGenerator_CLvsPOS()'  

**Data Loaders**:
A MATLAB data loader for the ARPES / PES data acquired at the ADRESS, SIS and PEARL beamlines in the Swiss Light Source (SLS) are all available. 
- 'load_adress_data()'  
- 'load_sis_data()'  
- 'load_pearl_data()'  


## MATLAB Version control  
MATLAB version:   2022a  
MATLAB add-ons (recommended): Database Toolbox, Image Processing Toolbox, Global Optimization Toolbox, Optimization Toolbox, Curve Fitting Toolbox, Parallel Processing Toolbox, Statistics and Machine Learning Toolbox, Signal Processing Toolbox.


## Authors
MatBase Scripts:  
**Dr. Procopios Constantinou**,  
Swiss Light Source (SLS),  
Paul Scherrer Institute (PSI),  
email: procopios.constantinou@psi.ch

## Acknowledgments

## License  
This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details.

--PCC, July 2023