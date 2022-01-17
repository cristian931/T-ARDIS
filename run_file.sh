#!/usr/bin/env bash



#####################################################
#
# This Script is used to Automatize the FAERS cleaning procedure (That rely on a Postgresql and different vocabularies), 
# the download and cleaning procedure of three different drug-se databases (MEDEFFECT
# OFFSIDE and SIDER (using a python script), and a compendium of different DRUG-TARGET databases called DRUG-TARGETS COMMON and STITCH
# cleaned with a python script too.
# The script manages the downloading and unpacking of file, the folder creation and file / database management.
#
# The script is divided in sections based on the work done
# The entire procedure will require at least two days and AT LEAST 45 GB OF FREE MEMORY, 
# since the files are huge and unfortunately need to be processed at the same time
#
# USAGE bash download_data.sh
#
#####################################################

# First thing: Set up and install a postgresql database on the local machine or server, 
# optional download and activate the pgadmin utility to explore with a graphic interface the server.
# The procedure won't require it since each database command is run by script or line.
#
# service httpd start
# service postgresql-11 start
# service pgadmin4 start
#
#

# create a python virtual environment
echo
echo setting up Python Virtual Enviroment and installing required packages

conda update conda --all
conda update anaconda

conda create --name conda_env python=3.7

echo Installing packages
conda install -q -y -n conda_env requirements.txt
conda install -q -y -n conda_env -c rdkit rdkit
conda install -q -y -n conda_env -c bjrn pandarallel

conda activate conda_env

git clone https://github.com/puolival/multipy.git
ipython multipy/setup.py install


#create a database called DRUG_ADR_polishing_procedure

echo
echo Setting up Database
echo

# if already exist the database disconnect everyone from it
psql -U postgres \
     -c "SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE datname = 'DRUG_ADR_polishing_procedure';"

# drop the entire database for a clean start
psql -U postgres \
     -c "drop database if exists DRUG_ADR_polishing_procedure;"  > /dev/null 2>&1

createdb -U postgres \
         -h localhost DRUG_ADR_polishing_procedure \
         > /dev/null 2>&1 # create the database

psql -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -c "CREATE SCHEMA faers"  \
     > /dev/null 2>&1 # create the faers schema


# FAERS and LAERS download and preparation

echo
echo FAERS and LAERS download and preparation
echo
# use a python script to get the links for downloading the FAERS data
python3.7 all_scripts/get_urls.py | grep ascii > urls

foldvar=faers_files
mkdir $foldvar
cd $foldvar || { echo "Error ${foldvar} not found"; exit 1; }

for url in $(cat ../urls)
do
	name=$(echo ${url} | sed 's|https://fis.fda.gov/content/Exports/||g')
	echo Downloading "${name}"
	wget "${url}" > /dev/null 2>&1
	unzip "${name}" > /dev/null 2>&1

	# remove the non necessary file for each data package
	mv ascii/* . 2> /dev/null
	mv ASCII/* . 2> /dev/null
	rm -rf ASCII/ 2> /dev/null
	rm -rf deleted/ 2> /dev/null
	rm -rf DELETED/ 2> /dev/null
	rmdir ascii/ 2> /dev/null
	rm -rf *.pdf 2> /dev/null
	rm -rf *.PDF 2> /dev/null
	rm -rf FAQs* 2> /dev/null
	rm -rf Readme* 2> /dev/null
	rm -rf ASC_NTS* 2> /dev/null
	rm -rf *.zip 2> /dev/null
	rm -rf *.doc 2> /dev/null
done

mv DEMO18Q1_new.txt ./DEMO18Q1.txt 2> /dev/null

# convert all file names in uppercase (even extension)
for i in *
do
    mv -i "$i" "$(echo "$i" | tr '[:lower:]' '[:upper:]')"
done 2> /dev/null


# launch the github preprocessing on the files
echo
echo wrapping up files
echo

bash ../all_scripts/Legacy_data_wrapper
bash ../all_scripts/Current_data_wrapper # This scripts take the data until 2018Q4

ls | grep "with_filename" | grep -v all | xargs rm # remove intermediate non necessary files


# update the current files with new data if present
echo
echo Updating with files later than 2018
echo

filetype="DEMO  DRUG  INDI  OUTC  REAC  RPSR  THER" # define the different files class

for type in ${filetype}
do
	ls | grep ${type} > tmp # grep the file of interest

	# transform the type to lowercase to match the name of the new files
	type_lower=$(echo "${type}" | tr '[:upper:]' '[:lower:]')
	
	for file in $(cat tmp) # check the file if are later than 2018q4
	do
		value=$(echo "${file}" | sed "s/${type}//g" | sed 's/Q//g' | sed 's/.TXT//g')
		if [[ ${value#0} -gt 184 ]] # 184 is the suffix of file from year 18q4 without the q, 
					    # so if the suffix on the file is bigger than 184 it means that the file is later
					    # 2018 and have to be processed and appended.
					    # The #0 means that the variable value has to be interpreted as decimal
    			then
    			sed 's/\r$//' "${file}"| sed '1,1d' | sed "1,$ s/$/\$${file}/" \
    			>> all_${type_lower}_data_with_filename.TXT
		fi
	done
done


# append to the correct file the new data and erase the secondary file
cat all_demo_data_with_filename.TXT >> all_version_B_demo_data_with_filename.TXT
rm all_demo_data_with_filename.TXT

cat all_drug_data_with_filename.TXT >> all_version_B_drug_data_with_filename.TXT
rm all_drug_data_with_filename.TXT

cat all_reac_data_with_filename.TXT >> all_version_B_reac_data_with_filename.TXT
rm all_reac_data_with_filename.TXT

rm -rf tmp
cd ..

rm -rf urls 2> /dev/null


#DOWNLOAD OF ORANGE BOOK
echo
echo Orange book download and Load
echo
foldvar=orange_book
mkdir $foldvar
cd $foldvar || { echo "Error ${foldvar} not found"; exit 1; }
wget https://www.fda.gov/media/76860/download > /dev/null 2>&1
mv download ./orange_book.zip
unzip orange_book.zip > /dev/null 2>&1
psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f ../all_scripts/load_nda_table.sql \
     > /dev/null 2>&1
cd ..


#load the alternative ISO codes
echo
echo Loading ISO Codes
echo
psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/load_country_code_table.sql \
     > /dev/null 2>&1


#load the european active drug references
echo
echo Loading EU references
echo
psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/load_eu_drug_name_active_ingredient_table.sql \
     > /dev/null 2>&1


#Update Athena vocabularies
echo
echo Updating Athena Vocabularies
echo
cd athena || { echo "Error folder not found"; exit 1; }
read -p "Enter UMLS apikey \(retrivable at https://uts.nlm.nih.gov/uts/profile\): " varname
java -Dumls-apikey="$varname" -jar cpt4.jar 5
cd ..


# Load Athena files in a new schema called cdmv5
echo
echo Loading Athena files in Database
echo

psql -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -c "CREATE SCHEMA cdmv5" \
     > /dev/null

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/data_table_creation_athena_step_1.sql \
      > /dev/null 2>&1

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/data_table_creation_athena_step_2.sql \
     > /dev/null 2>&1

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/data_table_creation_athena_step_3.sql \
     > /dev/null 2>&1

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/data_table_creation_athena_step_4.sql \
     > /dev/null 2>&1

psql -h localhost  \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/data_table_creation_athena_step_5.sql \
     > /dev/null 2>&1


# Create the meddra / snomed mapping
echo
echo Creating Meddra-Snowmed reference Map
echo
psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/create_meddra_snomed_mapping_table.sql \
     > /dev/null 2>&1


# LOAD THE FAERS FILE INTO THE DATABASE
echo
echo Loading Faers file into Database
echo

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/load_legacy_faers.sql \
     > /dev/null 2>&1

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/load_current_faers.sql \
     > /dev/null 2>&1


# De-duplicate cases
echo
echo De-duplicating cases
echo

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/derive_unique_all_case.sql \
     > /dev/null 2>&1


# Map Drugs from Rx-norm and other databases (without USAGI)
echo
echo Mapping drug with Rx-norm db \(without Usagi procedure\)
echo

psql -h localhost -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/map_all_drugname_to_rxnorm_without_usagi.sql \
     > /dev/null 2>&1


#Download cleaned tables and create final table with python
mkdir FAERS_almost_clean

psql -h localhost -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/standardize_combined_drug_mapping.sql \
     > /dev/null 2>&1

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -c "\copy faers.standard_case_drug TO 'FAERS_almost_clean/cleaned_faers_drugs.csv' DELIMITER ',' CSV HEADER;" \
     > /dev/null 2>&1

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -c "\copy faers.reac_pt_legacy_list TO 'FAERS_almost_clean/legacy_side_effects.csv' DELIMITER ',' CSV HEADER;" \
     > /dev/null 2>&1

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -c "\copy faers.reac_pt_list TO 'FAERS_almost_clean/current_side_effects.csv' DELIMITER ',' CSV HEADER;" \
     > /dev/null 2>&1

echo
echo Faers final cleaning
echo

python3.7 all_scripts/faers_final_polishing.py \
          FAERS_almost_clean/cleaned_faers_drugs.csv \
          FAERS_almost_clean/legacy_side_effects.csv \
          FAERS_almost_clean/current_side_effects.csv > /dev/null 2>&1

# This concludes the FAERS clenaing procedure

########################################################################################################################

#OTHER DATABASES

echo
echo Starting Downloading Other Drug-Side Effects Databases
echo

# medeffect Download

echo
echo Downloading Medeffect
echo

foldvar=MEDEFFECT
mkdir $foldvar
cd $foldvar || { echo "Error ${foldvar} not found"; exit 1; }

wget \
https://www.canada.ca/content/dam/hc-sc/migration/hc-sc/dhp-mps/alt_formats/zip/medeff/databasdon/extract_extrait.zip \
> /dev/null 2>&1

unzip extract_extrait.zip > /dev/null 2>&1
rm *zip
mv cvponline*/* .
rmdir cvponline*
cd ..

# MEDEFFECT Cleaning Procedure (Applying the same steps of FAERS)
echo
echo "MEDEFFECT cleaning and statistical procedure"
echo

psql -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -c "CREATE SCHEMA medeffect"  \
     > /dev/null 2>&1 # create the medeffect schema

cd orange_book || { echo "Error orange_book folder not found"; exit 1; }
psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f ../all_scripts/load_nda_table_med.sql \
     > /dev/null 2>&1
cd ..

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/load_eu_drug_name_active_ingredient_table_med.sql \
     > /dev/null 2>&1

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/load_medeffect.sql \
     > /dev/null 2>&1

psql -h localhost -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/map_drugname_to_rxnorm_medeffect.sql \
     > /dev/null 2>&1

psql -h localhost -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -f all_scripts/standardize_combined_drug_mapping_med.sql \
     > /dev/null 2>&1

psql -h localhost \
     -U postgres \
     -d DRUG_ADR_polishing_procedure \
     -c "\copy medeffect.standard_case_drug TO 'MEDEFFECT/MEDEFFECT_DRUG_CLEANED.csv' DELIMITER ',' CSV HEADER;" \
     > /dev/null 2>&1

sed -i 's/"//g' MEDEFFECT/MEDEFFECT_DRUG_CLEANED.csv



# SIDER Download

echo
echo Dowloading SIDER 4.1
echo

foldvar=SIDER_4.1

mkdir $foldvar
cd $foldvar || { echo "Error ${foldvar} not found"; exit 1; }
wget sideeffects.embl.de/media/download/drug_names.tsv > /dev/null 2>&1
wget sideeffects.embl.de/media/download/meddra_all_se.tsv.gz > /dev/null 2>&1
gunzip meddra_all_se.tsv.gz > /dev/null 2>&1
cd ..


# OFFSIDE DOWNLOAD
echo
echo Downloading OFFSIDE
echo

foldvar=OFFSIDE
mkdir $foldvar
cd $foldvar || { echo "Error ${foldvar} not found"; exit 1; }
wget tatonettilab.org/resources/nsides/OFFSIDES.csv.gz > /dev/null 2>&1
gunzip OFFSIDES.csv.gz > /dev/null 2>&1
cd ..


# Cleaning procedure for MEDEFFECT, SIDER, OFFSIDE
echo
echo Cleaning Other Drug-Side effects Databases
echo
python3.7 all_scripts/Cleaning_procedure.py > /dev/null 2>&1


# Statistical Validation of FAERS and MEDEFFECT data

python3.7 all_scripts/stat_validation_Community_DRUG_ADR.py FAERS_DRUG_SE.input FAERS
python3.7 all_scripts/stat_validation_Community_DRUG_ADR.py MEDEFFECT_DRUG_SE.input MEDEFFECT
#######################################################################


# Download drug-targets file

echo
echo Starting Downloading Drug-Targets Databases
echo

#Download DRUG_TARGETS_COMMONS
echo
echo Downloading DRUG_TARGETS_COMMONS
echo

foldvar=DRUG_TARGETS_COMMONS
mkdir $foldvar
cd $foldvar || { echo "Error ${foldvar} not found"; exit 1; }
wget --no-check-certificate https://drugtargetcommons.fimm.fi/static/Excell_files/DTC_data.csv > /dev/null 2>&1
cd ..

echo
echo Downloading STITCH
echo

foldvar=STITCH
mkdir $foldvar
cd $foldvar || { echo "Error ${foldvar} not found"; exit 1; }
wget http://stitch.embl.de/download/protein_chemical.links.v5.0/9606.protein_chemical.links.v5.0.tsv.gz
wget http://stitch.embl.de/download/chemicals.v5.0.tsv.gz
wget http://stitch.embl.de/download/chemicals.inchikeys.v5.0.tsv.gz

gunzip *tsv.gz

cd ..


# cleaning targets
echo
echo Cleaning Drug Targets Databases
echo
python3.7 all_scripts/DTC_cleaning > /dev/null 2>&1
python3.7 all_scripts/STITCH_cleaning > /dev/null 2>&1

##################################################################################


# Merging, cleaning, relate the different databases.
# The output file will be a tab delimited file containing the relationship between target and side effect with
# the respective p-values and q-values
# A supplementary file with only the interactions <= 0.05 will be created as well.


mkdir relationship_analysis_input_files
mv *.input relationship_analysis_input_files/

echo
echo Final computation in progress
echo

python3.7 all_scripts/drug_target_se_computation.py


################################################################################
