# Cell registration across multiple sessions in large-scale calcium imaging data
This package is an implementation of a probabilistic approach for tracking the same neurons (cell registration) across multiple sessions 
in Ca2+ imaging data, developed by Sheintuch et al., 2016. The package includes a GUI that supports the entire registration procedure. 

## Usage and documentation
To run the cell registration procedure run CellReg.m from CellReg directory.


The inputs for the cell registration method are the spatial footprints of cellular activity (ROIs) of the cells that were detected in the different sessions. 


The main output for the cell registration method is the obtained mapping of cell identity across all registered sessions.
Other outputs include a log file with all the relevant information regarding the data, registration
configurations, and a summary of the registration results and quality and a figures directory with important figures that are saved automatically. 

For more information refer to the user manual.

## Main stages of the cell registration procedure

1. Loading the spatial footprints of cellular activity from the different sessions.

2. Transforming all sessions to a reference coordinate system using rigid-body transformation.

3. Computing a probabilistic model of the spatial footprints similarities
of neighboring cell-pairs from different sessions using the centroid
distances and spatial correlations.

4. Obtaining an initial cell registration according to an optimized registration threshold.

5. Obtaining the final cell registration based on a clustering algorithm.