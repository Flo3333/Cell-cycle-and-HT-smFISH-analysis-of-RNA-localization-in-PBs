# Cell cycle and HT-smFISH analysis of RNA localization in PBs
Source code from pipeline used for Cell cycle and HT-smFISH analysis of RNA localization in PBs - python 3.8.16  
*Author : Slimani Floric Institut de Génétique Humaine - CNRS - Montpellier*  
*2023-11-06*



**Dependencies**

- numpy  
- pandas  
- matplotlib  
- scipy  
- skimage  
- cellpose 2.0 (https://github.com/MouseLand/cellpose/blob/main/README.md)  
Pachitariu, M., Stringer, C. Cellpose 2.0: how to train your own model. Nat Methods 19, 1634–1641 (2022).  
- bigFish(https://github.com/fish-quant/big-fish)  
Imbert A, Ouyang W, Safieddine A, Coleno E, Zimmer C, Bertrand E, Walter T, Mueller F. FISH-quant v2: a scalable and modular tool for smFISH image analysis. RNA. 2022 Jun;28(6):786-795. doi: 10.1261/rna.079073.121. Epub 2022 Mar 28. PMID: 35347070; PMCID: PMC9074904.


 **Repository content**

 - **Pbody_analysis_byGenes.py**  
  Pipeline processing field of views and saves raw data in feather dataframes.
   
 - **Pbody_results_plots.py**  
   Read raw data in feather dataframes and compute a set of quantification graphs.

 - **Pbody_table_extracts.py**  
   Read raw data in feather dataframes and extracts set of data into excels files.
   
 - **CustomPandasFramework**  
   Custom python package used to handle data management during et after analysis pipeline.
   
 - **pbwrap**  
   Custom functions used to analyse images such as wrappers around cellpose and bigFish.
   
 - **cellpose model**  
   Retrained cellpose2 model.
