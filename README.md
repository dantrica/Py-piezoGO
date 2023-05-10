# Cement-based-composites

This repository includes scripts to analyze physicochemical, electrical, and mechanical data, such as UV-Vis, DLS, Z potential, impedance spectroscopy, OCP determination, transients, compressive strength, and others.

## The repository contains:

- The Python module called [analytics_GO.ipynb](analytics_GO.ipynb) to compute the electromechanical characterization.
- All [sample data](/data) acquired during the experimental campaign. Herein is located the main DataFrame "Medidas_GO.xlsx", which contains the information about the measurements performed and specimens fabrication procedures.
- The [output files](/outputs), where the new dataframes (Excel files), figures (.png files), or new series (.txt files) are saved.
- Some python [scripts](scripts) that can be called by [analytics_GO.ipynb](analytics_GO.ipynb) to improve the mechanical or electrical characterization.
- The [Simulink](/Simulink) folder conotains the Simulink scripts (files with extension .slx) and the parameter estimator files for each rGO-cement composites frabricated.

## Project Description: 

$\texttt{Pypiezo-GO}$ comprises a sequence of methods developed in Python language and implemented on a Google Colaboratory notebook apt for efficient cloud computing. The experimental database is organised with a DataFrame called $\texttt{Medidas GO.xlsx}$ including the paths and metadata of the experimental measurement files. In particular, the experimental database includes the measurement records extracted from cyclic voltammetry (CV), open circuit potential (OCP), and compressive testing (CT), as reported in reference [Triana-Camacho2023](https://doi.org/10.1016/j.cemconcomp.2023.105063). On this basis, $\texttt{Pypiezo-GO}$ synchronizes the mechanical and electrical records and calculates the effective capacitance and the piezoelectric parameters of voltage $g_{33}$ and charge $d_{33}$, which represent the key features of these materials for their use as strain sensors. The code also includes the possibility of calibrating a lumped equivalent circuiy developed in MATLAB/Simulink for signal processing applications.

## Usage instructions:

The module must be [analytics_GO.ipynb](analytics_GO.ipynb) uploaded to Google Drive and opened with Google Collaboratory. Then, the folders [data](/data), [scripts](/scripts), and [outputs](/outputs) must be uploaded to the Google Colab folder. 

## License information:

The module [analytics_GO.ipynb](analytics_GO.ipynb) is open source and published under the GPL option. Nevertheless, the Matlab files require the user get has a LICENCE agreement with Matlab.

## Troubleshooting:

The users must be take into account that the pandas' library requires the openpyxl module. Google collaboratory uses to update the version of this library periodically. Then, the users must modify the version of openpyxl on the section <b>0.3. Importing the libraries: lmfit, pro_data, ... from the Colab directory in Google Drive.</b> into the module [analytics_GO.ipynb](analytics_GO.ipynb)

### Matlab:
The authors have tested the Matlab (Simulink model) scripts with Matlab 9.12 (R2022a), but it should be compatible also with older versions.


## Contact information:

For support, comments or sugestions feel free to contact dantrica@saber.uis.edu.co
