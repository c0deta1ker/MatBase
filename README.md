# MatBase: A comprehensive MATerials DataBASE for photoelectron spectroscopy

[MatBase](https://github.com/c0deta1ker/MatBase) is an app that allows you to explore, analyze and extract the properties and photoelectron spectroscopy behaviors of various materials using MATLAB. You can use this app to:  
- Access a comprehensive database of material properties and photoelectron spectroscopy parameters, and edit it as you wish.  
- Visualize and manipulate crystal structures in both real and reciprocal space, and calculate the Brillouin zone and 2D slices. 
- Compute the electron inelastic mean free path for different materials using various formalisms, such as the Universal, TPP-2M, S1, and S2 methods.  
- Extract the binding energy spectrum, photoionization cross-sections, and photoionization asymmetry parameters of all elements from 1 to 118.  .  
- Model an N-layer sample stack to determine the PES intensity contribution from each layer.  


This app is designed for researchers, students, and enthusiasts who are interested in learning more about the physical and chemical properties of materials. Whether you are working on a project, studying for an exam, or just curious about the world of materials science and photoelectron spectroscopy, this app can help you find the information and tools you need. Download MatBase today and discover the fascinating world of materials science!  


If you are interested in photoelectron spectroscopy data analysis tools, we recommend downloading [ARPEScape](https://github.com/c0deta1ker/ARPEScape), a powerful and convenient MATLAB-based analysis codebase that allows you to load, process, simulate, analyse and fit data from several PES techniques, such as X-ray photoelectron spectroscopy (XPS) and angle-resolved photoelectron spectroscopy (ARPES).  


## Installation  
1. Download the *MatBase* repository.
2. Open MATLAB and use *Set Path* in the *Home* tab to add the *MatBase* repository and all its sub-folders into its saved search paths.
3. Make sure you also use *Set Path* to add the repository / folder that contains all of your data to be loaded in.
4. Type 'App_MatBase' in the MATLAB Command Prompt to boot up the Main Menu App.

## MatBase Apps
**(1) MatBase Main Menu**: The main MATLAB App that provides seamless navigation to all other MATLAB apps using the comprehensive Materials Properties Database. Accessible in MATLAB via 'App_MatBase'.
![App_MainMenu](MatBase-v1.0.1/0_ReadMeImages/App_MainMenu.png)  

**(2) Materials Database Editor App**: Effortlessly manage the Materials Properties Database with this user-friendly app. Add, delete, and edit entries with ease, streamlining your workflow and enhancing the efficiency of other apps that utilize the database for calculations. Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_DatabaseEditor'.
![App_MaterialsDatabaseEditor](MatBase-v1.0.1/0_ReadMeImages/App_MaterialsDatabaseEditor.png)  

**(3) Crystallography App**: Explore the intricacies of crystal structures with our powerful visualization app. View unit cells in both real and reciprocal space, calculate the Brillouin zone, and extract 2D slices for in-depth analysis. Conveniently access material parameters from the Materials Properties Database or manually enter data for materials not included in the database. Unlock new insights and enhance your understanding of crystal structures with ease! Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_Crystallography'.
![App_Crystallography](MatBase-v1.0.1/0_ReadMeImages/App_Crystallography.png)  

**(4) IMFP Calculator App**: Achieve unparalleled precision in your calculations of electron inelastic mean free path with our advanced app. Utilize a range of formalisms, including the Universal, TPP-2M, S1, and S2 methods, to achieve the most accurate results. Conveniently access material parameters from the Materials Properties Database or manually enter data for materials not included in the database. Enhance your workflow and achieve new levels of accuracy with ease! Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_IMFP'.
![App_IMFP](MatBase-v1.0.1/0_ReadMeImages/App_IMFP.png)  

**(5) Photoionization Database App**: Extract the binding energy spectrum, photoionization cross-sections, and photoionization asymmetry parameters of all elements between 1-118 with our powerful app. Easily identify the best core-levels to probe by finding the maximum cross-sections and gain a comprehensive overview of the core-level spectra of all elements! Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_Photoionization'.
![App_Photoionization](MatBase-v1.0.1/0_ReadMeImages/App_Photoionization.png)  

**(6) N-Layer PES Intensity Model App**: Create a sample stack of N layers with different materials and thicknesses from the Materials Properties Database. Then, use the Photoionisation Cross-Section and Asymmetry Database to get the core-levels, cross-sections, and asymmetry parameters for each layer. The PES intensities for each layer are then determined using the Beer-Lambert law, modulated by the electron inelastic mean free path calculated using the TPP-2M formalism whilst considering the experimental geometry. This will help you compare XPS depth profile data (spectral intensity versus photon energy or emission angle), and determine the layer thicknesses in your samples from XPS! Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_NLayerPESIntensities'.
![App_NLayerPESIntensities](MatBase-v1.0.1/0_ReadMeImages/App_NLayerPESIntensities.png)  


## List of Available Databases  
- **Material Properties Database (MPD)**: This is a local MATLAB database that compiles the most useful physical, electronic, optical and magnetic material properties of elements / compounds. The data is taken from a range of sources, where the 'average' values are used for parameters that had more than 1 unique value. Accessible in MATLAB via 'get_mpd_props()'.  

- **Physics Constants Database**: This is a MATLAB function that loads in many physics constants that can be used for data processing, or modelling of physical systems. Accessible in MATLAB via 'physics_constants()'.  

- **Crystallography Database**: Contains some MATLAB functions for extracting / viewing the unit cells in both real and reciprocal space, calculate the Brillouin zone, and extract 2D slices through the Brilluoin Zone for data analysis. You can translate and rotate the 2D slices as you wish. For some examples of how to use these functions, see the ‘CrystallographyExamples.mlx’ file.  

- **Electron Inelastic Mean Free Path Database (eIMFPD)**: This is a local MATLAB database that compiles the optical data from the NIST Electron Inelastic-Mean-Free-Path Database (http://dx.doi.org/10.18434/T48C78) so that the results can be easily called and accessed within MATLAB. This is accessible in MATLAB via 'get_eimfpd_props()'. Furthermore, predictive eIMFP formulas's are also available using the (1) Universal, (2) TPP-2M and (3) S1 & S2 formalisms.  

- **Photoionisation Cross-Section and Asymmetry Database (PIXSAD)**: This is a local MATLAB database that compiles the photoionisation cross-section and asymmetry parameter data of the elements. This is useful when modelling the total photoelectron intensity that originates from a given layer of a sample when performing ARPES / XPS experiments. Accessible in MATLAB via 'get_pixsad_props()'.   

- **Photoionisation Energy and Fluorescence Database (PIEFD)**: This is a local MATLAB database that compiles the photoionisation energy and fluorescence yield data of the elements. Accessible in MATLAB via 'get_piefd_props()'.  


## MATLAB Version control  
MATLAB version:   2022a  
MATLAB add-ons (recommended): Database Toolbox, Image Processing Toolbox.

## Authors
**Dr. Procopios Constantinou**,  
Swiss Light Source (SLS),  
Paul Scherrer Institute (PSI),  
email: procopios.constantinou@psi.ch

## Acknowledgments

## License  
This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details.

--PCC, July 2023