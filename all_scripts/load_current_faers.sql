--######################################################
--# Create demo staging tables DDL and
--# Load current FAERS data files into the demo table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists demo_staging_version_A;
create table demo_staging_version_A
(
primaryid varchar,
caseid varchar,
caseversion varchar,
i_f_code varchar,
event_dt varchar,
mfr_dt varchar,
init_fda_dt varchar,
fda_dt varchar,
rept_cod varchar,
mfr_num varchar,
mfr_sndr varchar,
age varchar,
age_cod varchar,
gndr_cod varchar,
e_sub varchar,
wt varchar,
wt_cod varchar,
rept_dt varchar,
to_mfr varchar,
occp_cod varchar,
reporter_country varchar,
occr_country varchar,
filename varchar
);
truncate demo_staging_version_A;

\copy demo_staging_version_A FROM 'faers_files/all_version_A_demo_data_with_filename.TXT' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from demo_staging_version_A order by 1;

drop table if exists demo_staging_version_B;
create table demo_staging_version_B
(
primaryid varchar,
caseid varchar,
caseversion varchar,
i_f_code varchar,
event_dt varchar,
mfr_dt varchar,
init_fda_dt varchar,
fda_dt varchar,
rept_cod varchar,
auth_num varchar,
mfr_num varchar,
mfr_sndr varchar,
lit_ref varchar,
age varchar,
age_cod varchar,
age_grp varchar,
sex varchar,
e_sub varchar,
wt varchar,
wt_cod varchar,
rept_dt varchar,
to_mfr varchar,
occp_cod varchar,
reporter_country varchar,
occr_country varchar,
filename varchar
);
truncate demo_staging_version_B;

\copy demo_staging_version_B FROM 'faers_files/all_version_B_demo_data_with_filename.TXT' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from demo_staging_version_B order by 1;

drop table if exists demo;
create table demo as
select
primaryid,
caseid,
caseversion,
i_f_code,
event_dt,
mfr_dt,
init_fda_dt,
fda_dt,
rept_cod,
null as auth_num,
mfr_num,
mfr_sndr,
null as lit_ref,
age,
age_cod,
null as age_grp,
gndr_cod as sex,
e_sub,
wt,
wt_cod,
rept_dt,
to_mfr,
occp_cod,
reporter_country,
occr_country,
filename
from demo_staging_version_A
union all
select * from demo_staging_version_B;

select distinct filename from demo order by 1;
--######################################################
--# Create drug staging tables DDL and
--# Load current FAERS data files into the drug table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists drug_staging_version_A ;
create table drug_staging_version_A
(
primaryid varchar,
caseid varchar,
drug_seq varchar,
role_cod varchar,
drugname varchar,
val_vbm varchar,
route varchar,
dose_vbm varchar,
cum_dose_chr varchar,
cum_dose_unit varchar,
dechal varchar,
rechal varchar,
lot_nbr varchar,
exp_dt varchar,
nda_num varchar,
dose_amt varchar,
dose_unit varchar,
dose_form varchar,
dose_freq varchar,
filename varchar
);
truncate drug_staging_version_A;

\copy drug_staging_version_A FROM 'faers_files/all_version_A_drug_data_with_filename.TXT' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from drug_staging_version_A order by 1;

drop table if exists drug_staging_version_B;
create table drug_staging_version_B
(
primaryid varchar,
caseid varchar,
drug_seq varchar,
role_cod varchar,
drugname varchar,
prod_ai varchar,
val_vbm varchar,
route varchar,
dose_vbm varchar,
cum_dose_chr varchar,
cum_dose_unit varchar,
dechal varchar,
rechal varchar,
lot_num varchar,
exp_dt varchar,
nda_num varchar,
dose_amt varchar,
dose_unit varchar,
dose_form varchar,
dose_freq varchar,
filename varchar
);
truncate drug_staging_version_B;

\copy drug_staging_version_B FROM 'faers_files/all_version_B_drug_data_with_filename.TXT' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from drug_staging_version_B order by 1;

drop table if exists drug;
create table drug as
select
primaryid,
caseid,
drug_seq,
role_cod,
drugname,
null as prod_ai,
val_vbm,
route,
dose_vbm,
cum_dose_chr,
cum_dose_unit,
dechal,
rechal,
lot_nbr as lot_num,
exp_dt,
nda_num,
dose_amt,
dose_unit,
dose_form,
dose_freq,
filename
from drug_staging_version_A
union all
select * from drug_staging_version_B;

select filename, count(*) from drug group by filename order by filename;
--######################################################
--# Create indi staging tables DDL and
--# Load current FAERS data files into the indi table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists indi ;
create table indi
(
primaryid varchar,
caseid varchar,
indi_drug_seq varchar,
indi_pt varchar,
filename varchar
);
truncate indi;

\copy indi FROM 'faers_files/all_indi_data_with_filename.TXT' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select filename, count(*) from indi group by filename order by 1;

--######################################################
--# Create outc staging tables DDL and
--# Load current FAERS data files into the outc table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists outc;
create table outc
(
primaryid varchar,
caseid varchar,
outc_code varchar,
filename varchar
);
truncate outc;

\copy outc FROM 'faers_files/all_outc_data_with_filename.TXT' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select filename, count(*) from outc group by filename order by 1

--######################################################
--# Create reac staging tables DDL and
--# Load current FAERS data files into the reac table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists reac_staging_version_A ;
create table reac_staging_version_A
(
primaryid varchar,
caseid varchar,
pt varchar,
filename varchar
);
truncate reac_staging_version_A;

\copy reac_staging_version_A FROM 'faers_files/all_version_A_reac_data_with_filename.TXT' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from reac_staging_version_A order by 1;

drop table if exists reac_staging_version_B;
create table reac_staging_version_B
(
primaryid varchar,
caseid varchar,
pt varchar,
drug_rec_act varchar,
filename varchar
);
truncate reac_staging_version_B;

\copy reac_staging_version_B FROM 'faers_files/all_version_B_reac_data_with_filename.TXT' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select distinct filename from reac_staging_version_B order by 1;

drop table if exists reac;
create table reac as
select
primaryid,
caseid,
pt,
null as drug_rec_act,
filename
from reac_staging_version_A
union all
select * from reac_staging_version_B;

select filename, count(*) from reac group by filename order by filename;
--######################################################
--# Create rpsr staging tables DDL and
--# Load current FAERS data files into the rpsr table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists rpsr;
create table rpsr
(
primaryid varchar,
caseid varchar,
rpsr_cod varchar,
filename varchar
);
truncate rpsr;

\copy rpsr FROM 'faers_files/all_rpsr_data_with_filename.TXT' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select filename, count(*) from rpsr group by filename order by 1;

--######################################################
--# Create therapy staging tables DDL and
--# Load current FAERS data files into the ther table
--#
--# LTS Computing LLC
--######################################################

set search_path = faers;

drop table if exists ther;
create table ther
(
primaryid varchar,
caseid varchar,
dsg_drug_seq varchar,
start_dt varchar,
end_dt varchar,
dur varchar,
dur_cod varchar,
filename varchar
);
truncate ther;

\copy ther FROM 'faers_files/all_ther_data_with_filename.TXT' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b' ;
select filename, count(*) from ther group by filename order by 1 ;
