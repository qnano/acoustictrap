General information (2DGAUSS-CHARACTERIZATION)

This software belongs to the article: Low-cost acoustic force trap in a
microfluidic channel. This software is implemented in MATLAB R2022a
requiring the Image processing, curve fitting, optimization and Parallel
computing toolboxes.

Data availability:

The experimental data underlying the figures in the article can be found
at: https://figshare.com/s/dd6a9df122444dfbc496

Operation instruction calibration:
1. Download the calibration data from https://figshare.com/s/dd6a9df122444dfbc496 (data are in the zip folder "1ms dial led 30 2.zip")
2. Give the path to the files in the folder 1ms dial led 30 2.zip in the code bead_calibration_fwhm.m
3. Run the code in bead_calibration_fwhm.m  and the calibration curve will be generated


Operation instruction trapping experiments (doing your own trap analysis): 
1. Download the experimental trapping videos from
https://figshare.com/s/dd6a9df122444dfbc496 called last_mearumentdschip...vpp where ... stand for a number between 3 to 9. 
2. Give the path to the files named "last_measurementdschip...vpp" in bead_experiments_fwhm.m, where ... needs to be replaced by number between 3 to 9. 
3. Download the folder COMSOL_data.zip and Calibration data.zip and put the data in the same directory as the MATLAB file bead_experiments_fwhm.m
4. Give the path to the right trap video called last_mearumentdschip...vpp where ... stand for a number between 3 to 9 in bead_experiments_fwhm.m. 
5. Run this code
6. Save the data called output after running bead_experiments_fwhm.m 
7. Repeat step 4 to 6 for all the trapping videos files 
8. Use all the saved output data in trapping_stiffness.m to calculate the trap stiffnes and create the trap
stiffness curve vs. peak-to-peak-voltage

Operation instruction trapping experiments (using avaialble dataset):
1. Open the MATLAB file
trapping_stiffness.m
2. Run this file with the recorded experimental dataset called "dschip_...vpp.mat" in the folder experimental data.zip,
where ... needs to be replace by a number between 3 to 9.


Reproducing data:

The trap stiffness from experiments is calculated in
trapping_stiffness.m in this code figure 7a and 7b are generated. The
code in trapping_stiffness.m uses the output calculated (sigma_x) in
bead_experiments_fwhm.m and the calibration curve calculated in
bead_calibration_fwhm.m, the tracking software in the tracking
software folder and 2D Gaussian function coded in D2GaussFunction.m. The
trap stiffness from simulation (COMSOL) is derived in
comsol_stiffness_dschip.m in the folder COMSOL_stiffness_dschip
(figure 5a and 5b) The calibration curve figure 6d is generated in the
MATLAB file bead_calibration.m. This is summarized below:

1.  COMSOL_stiffness_dschip: figure 5a and 5b
2.  bead_calibration.m: figure 6d
3.  trapping_stiffness.m: figure 7a and 7b are

Reference Software The tracking software is developed by Daniel Blair
and Eric Dufresne and can be found using the following link:

https://site.physics.georgetown.edu/matlab/ 

The code for publication quality plots can be found using the following link:
https://www.mathworks.com/matlabcentral/fileexchange/47921-plotpub-publication-quality-graphs-in-matlab
