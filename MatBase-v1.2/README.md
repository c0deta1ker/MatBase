## List of known bugs


## Log of updates
**v1.3 - //:** 

**v1.2 - 2024.09.03:** 
-> Updated all MatBase Apps to make them consistent and more modular so they run more efficiently.
--> App_MatBase_IMFP: 
----> Updated and added a new field to show the standard deviation of the IMFP for a given electron kientic energy.
----> Added a new IMFP formalism (JTP method).
--> App_MatBase_Photoionization: 
----> Large update which now allows the user to select the cross-section formalism (Scofield 1973, Yeh & Lindau 1985, Trzhaskovskaya 2018, Cant 2022) and export tabulated cross-section data as a *.txt file over a user-defined photo energy range. 
----> Added new buttons to make plotting of individual or all core-levels more efficient. 
----> Added new buttons to plot the binding energy spectrum of each element, as well as the option to export this as a *.txt file.
--> App_MatBase_PESModelViewer: 
----> A new app that lets the user view all the available PES Models.
--> App_MatBase_PESCurveFitter: 
----> A new app that lets the user curve-fit XPS / PES data.

**v1.0.2 - 2023.11.11:** 
-> Updated the Materials Properties Database (MPD)
---> Added new column called stoichiometry. This is an integer scalar of the stoiciometry of the material (e.g. for elements, stoic = 1, for molecules of the form G_gH_h, the stoiciometry is g + h). Used in the IMFP calculatios.
-> Updated the Photoionisation Cross-Section and Asymmetry Database (PIXSAD)
---> Added 1973 Scofield Database
---> Added 1985 Yeh & Lindau Database
---> Added 2018 Trzhaskovskaya HAXPES Database
---> Added 2022 Cant Database
-> Updated the Photoionisation Energy and Fluorescence Database (PIEFD)
---> Added 1993 Moulder Database
---> Added 2022 Cant Database
-> Updated the IMFP Database
---> Added literature sources for each calculation method.
---> Added the more recent S3 and S4 formalism to determine the effective attenuation length (EAL).

**v1.0.1 - 2023.09.06:** 
--> New App: PES Intensity Modelling of an N-Layer system.

**v1.0.0 - 2023.07.21:** MatBase published.