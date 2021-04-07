set search_path = medeffect;

drop table if exists medeffect_drugs;
create table medeffect_drugs
(
DRUG_PRODUCT_INGREDIENT_ID varchar,
DRUG_PRODUCT_ID varchar,
DRUGNAME varchar,
ACTIVE_INGREDIENT_ID varchar,
PROD_AI varchar
);
truncate medeffect_drugs;

\copy medeffect_drugs FROM 'MEDEFFECT/drug_product_ingredients.txt' WITH DELIMITER E'$' CSV HEADER QUOTE E'\b';