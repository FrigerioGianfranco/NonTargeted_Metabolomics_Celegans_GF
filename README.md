
This repository contains the R-code prepared by Gianfranco Frigerio for the elaboration of data presented in the following paper:

**Leveraging open cheminformatics tools for non-targeted metabolomics analysis of C. elegans: a workflow comparison and application to strains related to xenobiotic metabolism and neurodegeneration**

Gianfranco Frigerio<sup>1,2</sup>, Yunjia Lai<sup>3</sup>, Emma L. Schymanski<sup>1</sup>, Gary W. Miller<sup>3</sup>.


<sup>1</sup> Luxembourg Centre for Systems Biomedicine (LCSB), University of Luxembourg, 6, Avenue du Swing, L-4367 Belvaux, Luxembourg.

<sup>2</sup> Center for Omics Sciences (COSR), IRCCS San Raffaele Scientific Institute, Milan, Italy.

<sup>3</sup> Department of Environmental Health Sciences, Mailman School of Public Health at Columbia University, New York, NY, USA.


Corresponding authors: GF & ELS.



The expanded WormJam chemical list and related MSP libraries, the raw files of LC-MS/MS analyses, and the main input and output tables for running the R-script are reported in the Zenodo repository: https://doi.org/10.5281/zenodo.14975586.


The codes are reported in the folders following the chronological order of the data elaboration. In particular:
- The folder “01_patRoon”, containing the subfolders “HILIC POS” and “RPLC NEG”, contains the codes used to run the feature extractions and compound annotation with patRoon. The scripts include the separate use of PubChemLite and the expanded WormJam as chemical databases. The .mzmL files are reported in the repository on Zenodo in the folder "00_mzmL_files_raw_data_HILIC_RPLC".
- The folder “02_creation_of_WJ_MSP_libraries” contains the code used for the preparation of the MSP libraries by in-silico fragmenting the compounds from the expanded WormJam list with CFM-ID. The in-silico fragmentation itself was performed outside R (with the Docker Destkp software) but the codes used also for that are reported within the R-script as comments. The expanded WormJam database and the output MSP libraries are reported in the repository on Zenodo in the folder "01_WJ_and_MSPlibraries".
- The folder “03_following_elaboration1” contains a first step of following data elaboration, performed after the retrieval of the feature dataset and annotations with both patRoon (with the PubChemLite, PCL, and expanded WormJam, WJ, dataset) and MS-DIAL (with the publicMSP libraries and with the WJ MSP libraires prepared as reported in the previous folder). The elaborations include: filtering the features considering the separate pooled quality control samples, transformations of feature intensities, ANOVA analyses, visualisations with Eulero-Venn graphs (then reported in the paper in the Fig. 2 and supplementary figures S2, S5, S8, S11), assignment of levels of annotations. These elaborations were separately performed for each of the subfolders: “MS_DIAL_HILIC_POS_publicMSP”, “MS_DIAL_HILIC_POS_WJ”, “MS_DIAL_RPLC_NEG_publicMSP”, “MS_DIAL_RPLC_NEG_WJ”, “patRoon_HILIC_POS_PCL_WJ”, and “patRoon_RPLC_NEG_PCL_WJ”. The input data can be taken from the folder "02_feature_data" of the repository on Zenodo.
- Lastly, the folder “04_following_elaboration_all”, contains in a single script all the following elaboration together, including: the retrieval of chemical information (such as PubChem chemical identifier, KEGG codes), the creation of the summary table with annotated compounds at a level 3 or above (Table S1 and S2), the compound classification with Classyfire, the data visualisation building upset plots (Fig. S14), Eulero-Venn graphs (then reported in the Fig.3 of the paper), and Sankey Diagram with categories of molecules (Fig. S15, S16, S17), the graph reported in the Fig. 4 of the paper, the box-plots (Supplementary data 2), and finally the pathway enrichment analysis with the FELLA package (Fig. S18 and tables S3, S4, S5, S6, and S7). The input data can be found inside the folders "03_feature_data_elaborated_featINFO" and "04_feature_data_elaborated_feat_transf" in the repository on Zenodo.

In order to run the code with the same versions of dependencies used the following code just before running the code. 

```r
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::restore()
```

