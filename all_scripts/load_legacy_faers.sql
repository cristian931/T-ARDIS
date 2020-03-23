--######################################################
--# Create legacy demo staging tables DDL and
--# Load legacy FAERS data files into the demo_legacy table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists demo_legacy_staging_version_A;
create table demo_legacy_staging_version_A
(
ISR varchar,
"CASE" varchar,
I_F_COD varchar,
FOLL_SEQ varchar,
IMAGE varchar,
EVENT_DT varchar,
MFR_DT varchar,
FDA_DT varchar,
REPT_COD varchar,
MFR_NUM varchar,
MFR_SNDR varchar,
AGE varchar,
AGE_COD varchar,
GNDR_COD varchar,
E_SUB varchar,
WT varchar,
WT_COD varchar,
REPT_DT varchar,
OCCP_COD varchar,
DEATH_DT varchar,
TO_MFR varchar,
CONFID varchar,
FILENAME varchar
);
truncate demo_legacy_staging_version_A;

\copy demo_legacy_staging_version_A FROM 'faers_files/all_version_A_demo_legacy_data_with_filename.txt' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from demo_legacy_staging_version_A order by 1 limit 10;

drop table if exists demo_legacy_staging_version_B;
create table demo_legacy_staging_version_B
(
ISR varchar,
"CASE" varchar,
I_F_COD varchar,
FOLL_SEQ varchar,
IMAGE varchar,
EVENT_DT varchar,
MFR_DT varchar,
FDA_DT varchar,
REPT_COD varchar,
MFR_NUM varchar,
MFR_SNDR varchar,
AGE varchar,
AGE_COD varchar,
GNDR_COD varchar,
E_SUB varchar,
WT varchar,
WT_COD varchar,
REPT_DT varchar,
OCCP_COD varchar,
DEATH_DT varchar,
TO_MFR varchar,
CONFID varchar,
REPORTER_COUNTRY varchar,
FILENAME varchar
);
truncate demo_legacy_staging_version_B;

\copy demo_legacy_staging_version_B FROM 'faers_files/all_version_B_demo_legacy_data_with_filename.txt' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from demo_legacy_staging_version_B order by 1 ;

drop table if exists demo_legacy ;
create table demo_legacy as
select
ISR,
"CASE",
I_F_COD,
FOLL_SEQ,
IMAGE,
EVENT_DT,
MFR_DT,
FDA_DT,
REPT_COD,
MFR_NUM,
MFR_SNDR,
AGE,
AGE_COD,
GNDR_COD,
E_SUB,
WT,
WT_COD,
REPT_DT,
OCCP_COD,
DEATH_DT,
TO_MFR,
CONFID,
null as REPORTER_COUNTRY,
FILENAME
from demo_legacy_staging_version_A
union all
select * from demo_legacy_staging_version_B;
--######################################################
--# Create legacy drug staging tables DDL and
--# Load legacy FAERS data files into the drug_legacy table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists drug_legacy;
create table drug_legacy
(
ISR varchar,
DRUG_SEQ varchar,
ROLE_COD varchar,
DRUGNAME varchar,
VAL_VBM varchar,
ROUTE varchar,
DOSE_VBM varchar,
DECHAL varchar,
RECHAL varchar,
LOT_NUM varchar,
EXP_DT varchar,
NDA_NUM varchar,
FILENAME varchar
);
truncate drug_legacy;

\copy drug_legacy FROM 'faers_files/all_drug_legacy_data_with_filename.txt' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from drug_legacy order by 1 ;

--######################################################
--# Create legacy indication staging tables DDL and
--# Load legacy FAERS data files into the indi_legacy table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists indi_legacy;
create table indi_legacy
(
ISR varchar,
DRUG_SEQ varchar,
INDI_PT varchar,
FILENAME varchar
);
truncate indi_legacy;

\copy indi_legacy FROM 'faers_files/all_indi_legacy_data_with_filename.txt' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from indi_legacy order by 1 ;

--######################################################
--# Create legacy outc staging tables DDL and
--# Load legacy FAERS data files into the outc_legacy table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists outc_legacy;
create table outc_legacy
(
ISR varchar,
OUTC_COD varchar,
FILENAME varchar
);
truncate outc_legacy;

\copy outc_legacy FROM 'faers_files/all_outc_legacy_data_with_filename.txt' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from outc_legacy order by 1 ;

--######################################################
--# Create legacy reaction staging tables DDL and
--# Load legacy FAERS data files into the reac_legacy table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists reac_legacy;
create table reac_legacy
(
ISR varchar,
PT varchar,
FILENAME varchar
);
truncate reac_legacy;

\copy reac_legacy FROM 'faers_files/all_reac_legacy_data_with_filename.txt' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from reac_legacy order by 1 ;

--######################################################
--# Create legacy rpsrcation staging tables DDL and
--# Load legacy FAERS data files into the rpsr_legacy table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists rpsr_legacy;
create table rpsr_legacy
(
ISR varchar,
RPSR_COD varchar,
FILENAME varchar
);
truncate rpsr_legacy;

\copy rpsr_legacy FROM 'faers_files/all_rpsr_legacy_data_with_filename.txt' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from rpsr_legacy order by 1 ;
--######################################################
--# Create legacy therapy staging tables DDL and
--# Load legacy FAERS data files into the ther_legacy table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists ther_legacy;
create table ther_legacy
(
ISR varchar,
DRUG_SEQ varchar,
START_DT varchar,
END_DT varchar,
DUR varchar,
DUR_COD varchar,
FILENAME varchar
);
truncate ther_legacy;

\copy ther_legacy FROM 'faers_files/all_ther_legacy_data_with_filename.txt' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from ther_legacy order by 1 ;

