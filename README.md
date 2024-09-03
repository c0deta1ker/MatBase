# MatBase: A comprehensive MATerials DataBASE for photoelectron spectroscopy analysis

[MatBase](https://github.com/c0deta1ker/MatBase) is a MATLAB-integrated application that offers a comprehensive database and versatile toolbox for photoelectron spectroscopy (PES) analysis. Users can access and expand a material properties database, visualize crystal structures, compute electron inelastic mean free paths (IMFP) using up to 8 different formalisms, and extract binding energy spectra, photoionization cross-sections, and photoionization asymmetry parameters for all elements from 1 to 118 from 4 different database sources. Additionally, [MatBase](https://github.com/c0deta1ker/MatBase) allows for the simulation of layered sample stacks and provides an intuitive interface for curve fitting in XPS and PES data, supporting up to 10 curves, various curve shapes (Gaussian, Lorentzian, and Voigt), background subtraction methods (Polynomial, Shirley, and Tougaard), and parameter constraints.   


This app is designed for researchers, students, and enthusiasts who are interested in learning more about the physical and chemical properties of materials. Whether you are working on a project, studying for an exam, or just curious about the world of materials science and photoelectron spectroscopy, this app can help you find the information and tools you need. Download [MatBase](https://github.com/c0deta1ker/MatBase) today and discover the fascinating world of materials science!  

## Installation  
1. Download the *MatBase* repository.
2. Open MATLAB and use *Set Path* in the *Home* tab to add the *MatBase* repository and all its sub-folders into its saved search paths.
3. Make sure you also use *Set Path* to add the repository / folder that contains all of your data to be loaded in.
4. Type 'App_MatBase' in the MATLAB Command Prompt to boot up the Main Menu App.


## List of Available Databases  
- **Material Properties Database (MPD)**: This is a local MATLAB database that compiles the most useful physical, electronic, optical and magnetic material properties of elements / compounds. The data is taken from a range of sources, where the 'average' values are used for parameters that had more than 1 unique value. Accessible in MATLAB via the function 'get_mpd_props()' or a dedicated UI App by calling 'App_MatBase'.    

- **Physics Constants Database**: This is a MATLAB function that loads in many physics constants that can be used for data processing, or modelling of physical systems. Accessible in MATLAB via the function 'physics_constants()'.    

- **Crystallography Database**: This tool offers a suite of functions for extracting and viewing crystal unit cells in both real and reciprocal space. The materials database includes the lattice types and vectors for most elemental materials and some compounds. Additionally, you can extract 2D slices through the Brillouin Zone and obtain a projected view of the tessellated extended Brillouin Zone structure; you can also translate and rotate these 2D slices to match data overlays as needed.      

- **Electron Inelastic Mean Free Path Database (IMFPD)**: This is a local MATLAB database that allows you to calculate the electron inelastic mean free path (IMFP) easily. A range of IMFP formalisms are included, such as the Optical (NIST), Universal, TPP-2M, S1, S2, S3, S4 and JTP methods to achieve the most accurate results. This is interfaced with the materials database includes the lattice types and vectors for most elemental materials and some compounds, meaning that many parameters for inputs are already known. Accessible in MATLAB via the function 'calc_imfp()' or a dedicated UI App by calling 'App_MatBase'.    

- **Photoionization Cross-Section and Asymmetry Database (PIXSAD)**: This is a local MATLAB database that compiles the photoionization cross-sections and asymmetry parameters of all elements between 1 - 98. A range of formalisms are available, including the Scofield (1973), Yeh & Lindau (1985), Trzhaskovskaya (2018) and Cant (2022) methods. Accessible in MATLAB via the function 'calc_xsect()' or a dedicated UI App by calling 'App_MatBase'.      

- **Photoionization Energy and Fluorescence Database (PIEFD)**: This local MATLAB database compiles photoionization energy and fluorescence yield data for all elements from 1 to 98. It includes data from various sources, such as the Moulder (1999), Trzhaskovskaya (2018), and Cant (2022) databases. To find the binding energy of a specific element and core level, you can use the 'calc_be()' function or the dedicated UI app by calling 'App_MatBase'.      

- **PES Quantification, Modelling & Fitting**: Discover a powerful toolbox designed for handling photoelectron spectroscopy (PES) data. This comprehensive suite includes functions for data processing, compositional analysis, overlayer thickness determination, sensitivity factor determination, and curve fitting using various models. You’ll find a range of curve fitting functions, including Peak-Like Models (Gaussian, Lorentzian, Voigt, Doniach-Sunjic, and Pseudo-Voigt), as well as TopHat- and Step-Like functions for plotting atomic concentration curves. Additionally, various background subtraction methods (Polynomial, Shirley, and Tougaard) are available. For an enhanced user experience, a dedicated UI App for curve fitting can be accessed by calling ‘App_MatBase’.  


## MatBase Apps
**(1) MatBase Main Menu**: The main MATLAB App that provides seamless navigation to all other MATLAB apps using the comprehensive Materials Properties Database. Accessible in MATLAB via 'App_MatBase'.
![App_MatBase](MatBase-v1.2/0_ReadMeImages/App_MatBase.png)  

**(2) Materials Database Editor App**: Effortlessly manage the Materials Properties Database with this user-friendly app. Add, delete, and edit entries with ease, streamlining your workflow and enhancing the efficiency of other apps that utilize the database for calculations. Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_DatabaseEditor'.
![App_MatBase_DatabaseEditor](MatBase-v1.2/0_ReadMeImages/App_MatBase_DatabaseEditor.png)  

**(3) Crystallography App**: Explore the intricacies of crystal structures with our powerful visualization app. View unit cells in both real and reciprocal space, calculate the Brillouin zone, and extract 2D slices for in-depth analysis. Conveniently access material parameters from the Materials Properties Database or manually enter data for materials not included in the database. Unlock new insights and enhance your understanding of crystal structures with ease! Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_Crystallography'.
![App_MatBase_Crystallography](MatBase-v1.2/0_ReadMeImages/App_MatBase_Crystallography.png)  

**(4) IMFP Calculator App**: Achieve unparalleled precision in your calculations of electron inelastic mean free path with our advanced app. Utilize a range of formalisms, including the Universal, TPP-2M, S1, S2 and JTP methods, to achieve the most accurate results. Conveniently access material parameters from the Materials Properties Database or manually enter data for materials not included in the database. Enhance your workflow and achieve new levels of accuracy with ease! Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_IMFP'.
![App_MatBase_IMFP](MatBase-v1.2/0_ReadMeImages/App_MatBase_IMFP.png)  

**(5) Photoionization Database App**: Extract the binding energy spectrum, photoionization cross-sections, and photoionization asymmetry parameters of all elements between 1-118 with our powerful app. Easily identify the best core-levels to probe by finding the maximum cross-sections and gain a comprehensive overview of the core-level spectra of all elements! Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_Photoionization'.
![App_MatBase_Photoionization](MatBase-v1.2/0_ReadMeImages/App_MatBase_Photoionization.png)  

**(6) PES Curve Fitting App**: Easily load your XPS data from a standard .txt file and construct a series of curves, including principal and spin-orbit split components. Begin by developing an initial model (using Gaussian, Lorentzian, Voigt and Doniach-Sunjic Peak Shapes), applying background subtraction (Polynomial, Shirley, and Tougaard) and setting constraints, an optimization algorithm is then run to achieve the best fit. Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_PESCurveFitter'.  
![App_MatBase_PESCurveFitter](MatBase-v1.2/0_ReadMeImages/App_MatBase_PESCurveFitter.png)     


**(7) N-Layer PES Intensity Model App**: Create a sample stack of N layers with different materials and thicknesses from the Materials Properties Database. Then, use the Photoionization Cross-Section and Asymmetry Database to get the core-levels, cross-sections, and asymmetry parameters for each layer. The PES intensities for each layer are then determined using the Beer-Lambert law, modulated by the electron inelastic mean free path calculated using the JTP formalism whilst considering the experimental geometry. This will help you compare XPS depth profile data (spectral intensity versus photon energy or emission angle), and determine the layer thicknesses in your samples from XPS! Accessible in MATLAB via the MatBase Main Menu, or explicitely by typing 'App_MatBase_NLayerPESIntensities'.  
![App_MatBase_NLayerModeller](MatBase-v1.2/0_ReadMeImages/App_MatBase_NLayerModeller.png)       


## MATLAB Version control  
MATLAB version:   2024a  
MATLAB add-ons (recommended): Database Toolbox, Image Processing Toolbox.


## Authors
**Dr. Procopios Constantinou**,  
Swiss Light Source (SLS),  
Paul Scherrer Institute (PSI),  
email: procopios.constantinou@psi.ch


## Acknowledgments  


## References

**Inelastic Mean-Free Path (IMFP) Sources**  
Universal (1979):    
[[1](http://dx.doi.org/10.18434/T48C78)] _Seah, M. Pl, and W. A. Dench. "Quantitative electron spectroscopy of surfaces: A standard data base for electron inelastic mean free paths in solids." Surface and interface analysis 1.1 (1979): 2-11_  

TPP-2M Formalism (1994):      
[[2](https://doi.org/10.1002/sia.740010103)] _Tanuma, Shigeo, Cedric J. Powell, and David R. Penn. "Calculations of electron inelastic mean free paths. V. Data for 14 organic compounds over the 50–2000 eV range." Surface and interface analysis 21.3 (1994): 165-176_  
[[3](https://doi.org/10.1002/sia.1526)] _Tanuma, Shigeo, Cedric J. Powell, and David R. Penn. "Calculation of electron inelastic mean free paths (IMFPs) VII. Reliability of the TPP‐2M IMFP predictive equation." Surface and interface analysis 35.3 (2003): 268-275_    
[[4](https://doi.org/10.1002/sia.4816)] _Seah, M. P. "An accurate and simple universal curve for the energy‐dependent electron inelastic mean free path." Surface and interface analysis 44.4 (2012): 497-503_    
[[5](https://doi.org/10.1002/sia.3522)] _Tanuma, Shigeo, C. J. Powell, and D. R. Penn. "Calculations of electron inelastic mean free paths. IX. Data for 41 elemental solids over the 50 eV to 30 keV range." Surface and interface analysis 43.3 (2011): 689-713_    

NIST Electron IMFP Database (1999):     
[[6](http://dx.doi.org/10.18434/T48C78)] NIST Optical Experiment Data  

S1 & S2 Formalism (2011):     
[[4](https://doi.org/10.1002/sia.4816)] _Seah, M. P. "An accurate and simple universal curve for the energy‐dependent electron inelastic mean free path." Surface and interface analysis 44.4 (2012): 497-503_    

S3 & S4 Formalism (2012):    
[[7](https://doi.org/10.1002/sia.5033)] _Seah, M. P. "Simple universal curve for the energy‐dependent electron attenuation length for all materials." Surface and interface analysis 44.10 (2012): 1353-1359_   

JTP Formalism (2023):    
[[8](https://doi.org/10.1002/sia.7217)] _Jablonski, Aleksander, Shigeo Tanuma, and Cedric J. Powell. "Calculations of electron inelastic mean free paths (IMFPs). XIV. Calculated IMFPs for LiF and Si3N4 and development of an improved predictive IMFP formula." Surface and Interface Analysis 55.8 (2023): 609-637_   

**Photoionization Energies, Cross-Sections & Asymmetry Sources**  
Scofield Database (1973):  
[[9](https://doi.org/10.2172/4545040)] _Scofield, J. H. "Theoretical photoionization cross sections from 1 to 1500 keV"_   

Yeh & Lindau Database (1985):   
[[10](https://doi.org/10.1016/0092-640X(85)90016-6)] _Yeh, J. J., and I. Lindau. "Atomic subshell photoionization cross sections and asymmetry parameters: 1⩽ Z⩽ 103." Atomic data and nuclear data tables 32.1 (1985): 1-155_   

Moulder Binding Energies Database (1993):   
[[11](https://scholar.google.com/scholar_url?url=https://www.researchgate.net/profile/Akif-Zeb/post/How_can_I_evaluate_at_of_N_in_TiO2_using_XPS_technique/attachment/5f3ebac4ce377e00016e3bc5/AS%253A926669195993088%25401597946561045/download/MANXPS.pdf&hl=en&sa=T&oi=gsb-ggp&ct=res&cd=0&d=11053645406167494942&ei=Aw3XZr_mFv-Xy9YPgcCeoQo&scisig=AFWwaeajp-vE3wtFLu1NvP33L_uI)] _Moulder, John F. et al. “Handbook of X-Ray Photoelectron Spectroscopy.” (1992)._   

Trzhaskovskaya Database (2018-2019):    
[[12](https://doi.org/10.1016/j.adt.2017.04.003)] _Trzhaskovskaya, M. B., and V. G. Yarzhemsky. "Dirac–Fock photoionization parameters for HAXPES applications." Atomic Data and Nuclear Data Tables 119 (2018): 99-174._   
[[13](https://doi.org/10.1016/j.adt.2019.05.001)] _Trzhaskovskaya, M. B., and V. G. Yarzhemsky. "Dirac–Fock photoionization parameters for HAXPES applications, Part II: Inner atomic shells." Atomic Data and Nuclear Data Tables 129 (2019): 101280_   

Cant Database (2018-2019):    
[[14](https://doi.org/10.1002/sia.7059)] _Cant, David JH, et al. "Quantification of hard X‐ray photoelectron spectroscopy: Calculating relative sensitivity factors for 1.5‐to 10‐keV photons in any instrument geometry." Surface and Interface Analysis 54.4 (2022): 442-454_   

**Materials Properties Database Sources**  
[[15](https://www.wolframalpha.com/)] _Electronegativity, electron affinity, ionisation energies, temperatures and crystal structures_   
[[16](https://www.schoolmykids.com/learn/periodic-table-of-elements/)] _Temperatures, crystal structures, unit cell parameters, electronic and magnetic properties_   
[[17](https://doi.org/10.1063/1.3253115/)] _Strehlow, W. H., and Earl L. Cook. "Compilation of energy band gaps in elemental and binary compound semiconductors and insulators." (1973): 163-200._   
[[18](https://en.wikipedia.org/)] _Band-gap estimates,  temperatures,  crystal structures,  unit cell parameters,  electronic and magnetic properties_   
[[19](https://pubchem.ncbi.nlm.nih.gov/)] _Atomic weights of elements and compounds_   

**Useful PES Quantification Literature**  
[[20](https://doi.org/10.1016/0009-2614(76)80496-4)] _Hill, J. M., et al. "Properties of oxidized silicon as determined by angular-dependent X-ray photoelectron spectroscopy." Chemical Physics Letters 44.2 (1976): 225-231_   
[[21](https://doi.org/10.1002/sia.5934)] _Walton, J., et al. "Film thickness measurement and contamination layer correction for quantitative XPS." Surface and Interface Analysis 48.3 (2016): 164-172_   
[[22](https://doi.org/10.1116/1.5141395)] _Shard, Alexander G. "Practical guides for x-ray photoelectron spectroscopy: Quantitative XPS." Journal of Vacuum Science & Technology A 38.4 (2020)_   

**Useful Websites for XPS Information & Application Notes**    
[[23](https://srdata.nist.gov/xps/)] NIST X-ray Photoelectron Spectroscopy Database  
[[24](https://www.xpsfitting.com/search/label/About%20This%20Site)] Surface Science Western laboratories (XPS Reference Pages)  
[[25](https://www.xpsdata.com/xpsdata.htm)] XPS Information and Application Notes   
[[26](https://xpsdatabase.net/)] B. Vincent Crist: International XPS Database  
[[27](https://xpslibrary.com/)] B. Vincent Crist: XPS Information and Application Notes  
[[28](https://a-x-s.org/research/cross-sections/)] A. Regoutz: Source of the Digitized Photoionization Parameters  


## License  
This project is licensed under the GNU General Public License v3.0 - see the [LICENSE.md](LICENSE.md) file for details.

--PCC, September 2024