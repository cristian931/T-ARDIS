# Side-Effect-Target-statistic-search
This repository contain all the file necessary to download, compute and statistically validate the pairwise interaction between Drug side effects and Drug Targets

In order to make it works a couple of pre-processing steps are necessary.
1) Set up a postgresql database where to run the FAERS cleaning procedure. It's necessary just create the server, all other procedure, like the database creation are done automatically in the bash script
2) Download an updated version of ATHENA vocabularies. These vocabularies are used to map and retrieve unique DRUG namses in the FAERS procedure. The total number of used vocabularies are 23, included the SNOWMED and MEDDRA dictionaries.
  - To download that you need to be registered in the ATHENA website http://athena.ohdsi.org and have a license for the access the MEDDRA dictionary, the latter can be requested registering in the UMLS website.
  - Once the registering procedures are completed, to download the dictionaries  just go again to the ATHENA website at http://athena.ohdsi.org/vocabulary/list, add to the already selected databases the MEDDRA and click on "download vocabulaires".
  - Save and extract the files in a folder called "athena" in the same parent folder in which are saved all the other repository files.
3) Launch the run_file.sh script
  - The script will start downlaoding the updated files of the different DRUG - SIDE EFFECT databases:
  
      - FAERS (https://www.fda.gov/drugs/questions-and-answers-fdas-adverse-event-reporting-system-faers/fda-adverse-event-reporting-system-faers-latest-quarterly-data-files)
      
      - MEDEFFECT (https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database/medeffect-canada-caveat-privacy-statement-interpretation-data-search-canada-vigilance-adverse-reaction-online-database.html)
      
      - OFFSIDE (http://tatonettilab.org/offsides/)
      
      - SIDER_4.1 (http://sideeffects.embl.de/)
      
  - And the updated versions of the DRUG - TARGETs databases:
  
      - DRUGBANK (https://www.drugbank.ca/)
      
      - DRUGCENTRAL (http://drugcentral.org/)
      
      - DGidb (http://www.dgidb.org/)

4) The procedure will standardize, clean and relate the information in the databases to obtain a list of pairwise relationships between the Drug's targets and side effects expliting the drug name as bridge.
In order to confirm the relationship the results have been statistically validated using the Fisher exact test as explained in the paper of Kuhn et al. (http://europepmc.org/article/MED/23632385)

# Output
Different files of output are obtained from the procedure:
  - The databases cleaned versions, in the folder relationship_analysis_input_files. The latter are tsv files containing the access code of the different entries (if available) and the information of the drugs and the associated side effects 
  - The pairwise relationship file containing the computed p-value and q-avalue correction of the association target - side effect, called "qvalues_all_interactions"

# Usage
1) Set up the postgresql server
2) clone or download all the file in the repository and save it to a folder
3) Download the Athena dictionaries and save them into a folder called "athena" inside the repository folder
4) run the run_file.sh script (you need to have the permission to create/remove folder)
