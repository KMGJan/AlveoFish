# AlveoFish

The project contains code and data.

# **code/**
This folder contains all the code needed to reproduce this study.

-   `analyses.R` contains the analyses

## **data/raw/**
This folder contains the asv table assigned to the PR2 v.5.0.0 database and the associated metadata

#### ASV table

- ASV table (**asv_18S.csv**) from [Novotny et al. (2022)](https://doi.org/10.1038/s41598-022-15116-7)[^1] and [Jan et al. (2025)](https://doi.org/10.1093/icesjms/fsaf122)[^2].
  - Each **row** represents a unique ASV.  
  - Column **1** contains a unique ID for each ASV
  - Columns **2–10** contain taxonomic annotations.  
  - Column **11** contains the ASV sequence.  
  - Column **12** contains the reference database used for taxonomy assignment
  - All subsequent columns represent unique samples, with cell values corresponding to the number of reads for each ASV in that sample.


#### Sample Metadata

The corresponding metadata (**metadata_18S.csv**) file contains one row per sample.  

| Field Name  | Description |
|---|-------------|
| `library_ID` | Unique identifier for each sequencing library. |
| `title` | Short descriptive title for the sample. |
| `organism`   | Taxonomic name or environmental sample descriptor (e.g., "Clupea harengus"). |
| `collection_date` | Date when the sample was collected (YYYY-MM-DD). |
| `geo_loc_name` | Geographic location, correspond to the ICES statistical rectangle or to the monitoring station name. |
| `depth` | Sampling depth interval in meters. |
| `samp_size` | Size or volume of the sample collected (e.g., 1 L, 1 individual). |
| `size_frac` | Size fraction (length for fish species, mesh size for WP2 samples, filter size for water samples) |
| `lat_lon` | Latitude and longitude in decimal degrees. |
| `design_description` | Short description of the sampling design or context. |
| `env_broad_scale` | Broad environmental context (e.g., "Pelagic Baltic Sea"). |
| `env_local_scale` | More specific local context (e.g., "ICES statistical rectangle 44G7"). |
| `env_medium` | Type of environmental material sampled (e.g., "Seawater"). |
| `collection_method` | Description or name of the sampling protocol or method. |
| `samp_collection_device` | Equipment or device used to collect the sample (e.g., "Niskin bottle", "Pelagic trawl", "WP2 , 90um"). |
| `samp_mat_process` |  Description of any material processing steps (e.g., "Bulk DNA extraction, QIAmp Micro Kit"). |
| `source_material_ID` | Identifier linking to a parent or source material, here matches ASV table column. |
| `sample_name` | Unique sample name used in the study. |
| `sample_accession` | ENA BioSample accession number. |
| `study_accession` | Associated ENA BioProject or study accession number. |


[^1]: Novotny A, Jan KMG, Dierking J, Winder M. 2022. *Niche partitioning between planktivorous fish in the pelagic Baltic Sea assessed by DNA metabarcoding, qPCR and microscopy*. Scientific Reports. https://doi.org/10.1038/s41598-022-15116-7
[^2]: Jan KMG, Hentati-Sundberg J, Larson N, Winder M. 2025. *Limited resource use overlaps among small pelagic fish species in the central Baltic Sea*. ICES Journal of Marine Science. https://doi.org/10.1093/icesjms/fsaf122
