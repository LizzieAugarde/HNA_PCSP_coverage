################ HNA/PCSP coverage analysis - data extraction ###########################

#SQL extracts for patient cohort and linked HNA/PCSP records from RCRD

#Created September 2024 by Lizzie Augarde 
###################################################################################


library(NDRSAfunctions)

casref01 <- createConnection()
cas2407 <- createConnection(port = 1525, sid = "cas2407")

###### EXTRACTING PATIENT COHORT #######
#patient cohort (all 2021 diagnoses exc C44, all ages, non DCO and not diagnosed same day as death), including demogs----
query <- "select * from (
  select pat.patientid,
  tum.tumourid,
  tum.diagnosisdatebest,
  tum.site_icd10r4_o2_3char_from2013,
  tum.stage_pi_detail,
  tum.stage_best_system,
  tum.stage_best,
  tum.gender,
  tum.age,
  tum.ethnicity,
  tum.deathdatebest,
  tum.diag_trust,
  imd.imd19_decile_lsoas,
  geo.icb_2022_name,
  geo.icb_2022_code,
  site.ndrs_main,
  rank () over (partition by pat.patientid order by tum.diagnosisdatebest, tum.tumourid asc) as rank
  from av2021.at_patient_england@casref01 pat 
  left join av2021.at_tumour_england@casref01 tum on pat.patientid = tum.patientid
  left join analysispollyjeffrey.at_site_england@casref01 site on tum.tumourid = site.tumourid --added 25/07/2024
  left join av2021.at_geography_england@casref01 geo on tum.tumourid = geo.tumourid
  left join imd.imd2019_equal_lsoas@casref01 imd on geo.lsoa11_code = imd.lsoa11_code
  where tum.diagnosisyear = 2021
  and tum.cascade_inci_flag = 1 
  and tum.dco = 'N'
  and tum.dedup_flag = 1
  --and tum.age >17 --removed 17/07/2024 (all ages due to lack of specificity in guidance)
  and ((tum.deathdatebest is null) or (tum.diagnosisdatebest != tum.deathdatebest)) 
  and ((tum.deathdatebest is null) or (tum.deathdatebest - tum.diagnosisdatebest > 0)) 
  and substr(tum.site_icd10_o2_3char, 1, 1) = 'C' 
  and tum.site_icd10_o2_3char <> 'C44'
  and ((tum.gender = '1' and tum.site_icd10_3char not in ('C51','C52','C53','C54','C55','C56','C57','C58') --added 16/09/2024 following QA, limiting to tumours matching gender
        and tum.site_icdo3rev2011_3char not in ('C51','C52','C53','C54','C55','C56','C57','C58'))
      or 
        (tum.gender = '2' and tum.site_icd10_3char not in ('C60','C61','C62','C63')
        and tum.site_icdo3rev2011_3char not in ('C60','C61','C62','C63'))
      or tum.site_icd10 is null or tum.site_icdo3rev2011_3char is null))
      subquery
where rank = 1"

patient_cohort <- dbGetQueryOracle(casref01, query, rowlimit = NA)


###### EXTRACTING PATIENT COHORT LINKED TO RCRD TO FIND THOSE MISSING FROM RCRD #######
query <- "select * from (
  select pat.patientid, 
  rp.nhsnumber
  from av2021.at_patient_england@casref01 pat 
  left join av2021.at_tumour_england@casref01 tum on pat.patientid = tum.patientid
  left join analysisncr.at_rapid_pathway@cas2407 rp on pat.nhsnumber = rp.nhsnumber
  where tum.diagnosisyear = 2021
  and tum.cascade_inci_flag = 1 
  and tum.dco = 'N'
  and tum.dedup_flag = 1
  and ((tum.deathdatebest is null) or (tum.diagnosisdatebest != tum.deathdatebest)) 
  and ((tum.deathdatebest is null) or (tum.deathdatebest - tum.diagnosisdatebest > 0)) 
  and substr(tum.site_icd10_o2_3char, 1, 1) = 'C' 
  and tum.site_icd10_o2_3char <> 'C44'
  and ((tum.gender = '1' and tum.site_icd10_3char not in ('C51','C52','C53','C54','C55','C56','C57','C58') --added 16/09/2024 following QA, limiting to tumours matching gender
        and tum.site_icdo3rev2011_3char not in ('C51','C52','C53','C54','C55','C56','C57','C58'))
      or 
        (tum.gender = '2' and tum.site_icd10_3char not in ('C60','C61','C62','C63')
        and tum.site_icdo3rev2011_3char not in ('C60','C61','C62','C63'))
      or tum.site_icd10 is null or tum.site_icdo3rev2011_3char is null))"

rcrd_missing_diags <- dbGetQueryOracle(casref01, query, rowlimit = NA)


###### EXTRACTING LINKED HNA AND PCSP RECORDS #######
query <- "select * from (
  select pat.patientid,
  rp.nhsnumber as rp_nhsnumber,
  rp.patientid as rp_patientid,       
  tum.tumourid,
  tum.diagnosisdatebest,
  rp.event_type,
  rp.event_property_1,
  rp.event_date,
  rp.trust_code,
  rank () over (partition by pat.patientid order by tum.diagnosisdatebest, tum.tumourid asc) as rank
  from av2021.at_patient_england@casref01 pat 
  left join av2021.at_tumour_england@casref01 tum on pat.patientid = tum.patientid
  left join analysisncr.at_rapid_pathway@cas2407 rp on pat.nhsnumber = rp.nhsnumber 
  where tum.diagnosisyear = 2021
  and tum.cascade_inci_flag = 1 
  and tum.dco = 'N'
  and tum.dedup_flag = 1
  --and tum.age >17 ---removed 17/07/2024
  and ((tum.deathdatebest is null) or (tum.diagnosisdatebest != tum.deathdatebest)) 
  and ((tum.deathdatebest is null) or (tum.deathdatebest - tum.diagnosisdatebest > 0)) 
  and substr(tum.site_icd10_o2_3char, 1, 1) = 'C' 
  and tum.site_icd10_o2_3char <> 'C44'
  and ((tum.gender = '1' and tum.site_icd10_3char not in ('C51','C52','C53','C54','C55','C56','C57','C58') --added 16/09/2024 following QA, limiting to tumours matching gender
        and tum.site_icdo3rev2011_3char not in ('C51','C52','C53','C54','C55','C56','C57','C58'))
      or 
        (tum.gender = '2' and tum.site_icd10_3char not in ('C60','C61','C62','C63')
        and tum.site_icdo3rev2011_3char not in ('C60','C61','C62','C63'))
      or tum.site_icd10 is null or tum.site_icdo3rev2011_3char is null)
  and rp.event_type in (20, 24))
where rank = 1"

raw_hna_pcsp_data <- dbGetQueryOracle(casref01, query, rowlimit = NA)
