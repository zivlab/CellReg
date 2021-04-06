# Cell registration across multiple sessions in large-scale calcium imaging data
This package is an implementation of a probabilistic approach for tracking the same neurons (cell registration) across multiple sessions 
in Ca2+ imaging data. The package includes a GUI that supports the entire registration procedure. 

For more information contact lironsheintuch@gmail.com or join our [slack channel](https://cellreg.slack.com).

## Setting up the repository
We encourage the use of official versions (e.g., v1.4.2) for easier debugging processes. Switch to the releases tab on GitHub and checkout the latest version.

1. Cloning:
`git clone https://github.com/zivlab/CellReg.git`
2. Checkout version:
`git checkout v<major>.<minor>.<bugfix> (e.g., v1.4.2)`
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

## References
Sheintuch, L., Rubin, A., Brande-Eilat, N., Geva, N., Sadeh, N., Pinchasof, O., Ziv, Y. (2017). Tracking the Same Neurons across Multiple Days in Ca2+ Imaging Data. *Cell Reports*, 21(4), pp. 1102–1115. doi: 10.1016/j.celrep.2017.10.013.
