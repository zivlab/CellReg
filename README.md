# Cell registration across multiple sessions in large-scale calcium imaging data
This package is an implementation of a probabilistic approach for tracking the same neurons (cell registration) across multiple sessions 
in Ca2+ imaging data, developed by Sheintuch et al., 2017. The package includes a GUI that supports the entire registration procedure. 

## Setting up the repository
1. Cloning:
`git clone https://github.com/zivlab/CellReg.git`
2. Checkout version:
`git checkout v<major>.<minor>.<bugfix> (e.g., v1.1.1)`
3. Run `setup.m`

## Usage and documentation
To run the cell registration procedure you can either use the full GUI version or access the API directly.
To use the GUI run *GUI\CellReg.m*.
An example of how to use the API is provided in the file *CellReg\demo.m*.


The inputs for the cell registration method are the spatial footprints of cellular activity (ROIs) of the cells that were detected in the different sessions. 


The main output for the cell registration method is the obtained mapping of cell identity across all registered sessions.
Other outputs include a log file with all the relevant information regarding the data, registration
configurations, and a summary of the registration results and quality. In addition, important figures are saved automatically in a figures directory. 


An example data set and cell registration results are provided in the *SampleData* directory.


For more information refer to the user manual found in the *Docs* directory.

## Main stages of the cell registration procedure

1. Loading the spatial footprints of cellular activity from the different sessions.

2. Aligning all sessions to a reference coordinate system using rigid-body transformation.

3. Computing a probabilistic model of the spatial footprints similarities
of neighboring cell-pairs from different sessions using the centroid
distances and spatial correlations.

4. Obtaining an initial cell registration according to an optimized registration threshold.

5. Obtaining the final cell registration based on a clustering algorithm.