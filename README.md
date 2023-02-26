# OrdinalScale_Statistics_PowerAnalyses

Main functionalities are finalized. Follow-ups: R Markdown documents and Shiny app for ease of demonstration and more practicality for the users with limited background in coding and R scripts.

The pipeline allows statistical analyses of ordinal scale data, compares distributions of measurements in control groups to test samples. 'MASS' (polr()) and 'emmeans' packages are used for these purposes. Furthermore, power analyses were employed by utilizing simulation based approaches. Hence, the functionalities are:
- The users can import experimental results that are based on ordinal scale measurements, 
- Produce spine plots from their data and obtain statistical descriptions (effect sizes, significance values etc)
- Collect curve graphs for power analyses
- Obtain minimum amount of samples needed to achieve user specified power levels
