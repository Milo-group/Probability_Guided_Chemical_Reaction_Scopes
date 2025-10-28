In this folder there is one subfolder and one ```*.R``` file:

**DeuterationData**
Contains the Training_Data.csv used for the deuteration analysis.

**Extreme_values.R**
R script performing five analyses on the dataset:
- Z-score: identifies outliers based on standard deviation from the mean.  
- Electronic z-score: highlights extreme values in electronic features.  
- Steric z-score: highlights extreme values in steric features.  
- Convex hull: checks the boundaries of the dataset in feature space.  
- Similarity-based extremes: detected data points that are extreme in each class according to our similarity-base under sampling.