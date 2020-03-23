# Side-Effect-Target-statistic-search
This repository contain all the file necessary to download, compute and statistically validate the pairwise interaction between Drug side effects and Drug Targets

In order to make it works a couple of pre-processing steps are necessary.
1) Set up a postgresql database where to run the FAERS cleaning procedure. It's necessary just create the server, all other procedure, like the database creation are done automatically in the bash script
2) Download an updated version of ATHENA vocabularies. These vocabularies are used to map and retrieve unique DRUG namses in the FAERS procedure. The total number of used vocabularies are 23, included the SNOWMED and MEDDRA dictionaries.
  - To download that you need to be registered in the ATHENA website http://athena.ohdsi.org and have a license for the access the MEDDRA dictionary, the latter can be requested registering in the UMLS website.
  - Once the registering procedures are completed, to download the dictionaries  just go again to the ATHENA website at http://athena.ohdsi.org/vocabulary/list, add to the already selected databases the MEDRRA and click on "download vocabulaires".
  - Save and extract the files in a folder called "athena" in the same parent folder in which are saved all the other repository files.
3) launch the run.sh file
