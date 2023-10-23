# simSAS
Simulation code for a self-adaptive system

This is simSAS, a Matlab tool to simulate a certain self-adaptive system.


REMARKS

For copyright and licensing, see header of m-files.
If you use simSAS, please refer to the corresponding article(s).
For details on the theory, please see the articles referenced in the headers of the m-files.
All technical drawings (including material specifications) and CAD files are provided here: https://doi.org/10.18419/darus-2462.
Feel free to use those files to reproduce the experiments and/or carry out further reserarch.


INSTRUCTIONS

To test simSAS, it is recommended to use publicly avialable experimental data, and compare
with simulation data:
1. Download experimental data from http://dx.doi.org/10.18419/darus-2883. 
    The file names are 'Time_VelBeam_RelSliderPos__fex104_ExcLevel14_RelSliderPosStart03.txt' and so on.
2. Put it into a subfolder 'DAT' within the folder containing the main m-files and the 'SRC' folder.
3. Run the 'simulate_FreeSliderModel.m' script.
    Here you have to select a scenario, labeled after the Figures in reference [1] 
    (see header of m-files).
    This will save simulation data in the 'DAT' folder.
4. Run the 'simulate_sSIM.m' script.
    Here you have to select a scenario, labeled after the Figures in reference [1] 
    (see header of m-files).
    This will save simulation data in the 'DAT' folder.
5. Run the 'illustrate_results.m' script.
    Here you have to select a scenario, labeled after the Figures in reference [1] 
    (see header of m-files).
    This requires that experimental and simulation data are available in the 'DAT' folder.
6. Have fun with the tool. Adjust the parameters or the model or the simulation to your 
    needs. 


Do not hesitate to contact me if you have any questions.

Malte Krack

(malte.krack@ila.uni-stuttgart.de)
