/****************** HNA/PCSP coverage analysis ******************

* Extracts 2021 COSD Level 3 registrations and link to HNA/PCSP
* data stored in RCRD. Testing linkage coverage and creating 
* table for inequalities analysis
* Created November 2023 by Lizzie Augarde 

* Change log:
16/01/2024 added test lines for looking at end of treatment dates. 
Reran main table to include demographic information
*****************************************************************/

create table hna_pcsp_2021diags compress basic nologging as 
select * from (
select pat.patientid,
       rp.nhsnumber as rp_nhsnumber,
       rp.patientid as rp_patientid,       
       tum.tumourid,
       tum.diagnosisdatebest,
       tum.diagnosisyear,
       tum.site_icd10_o2_3char,
       tum.stage_best,
       tum.gender,
       tum.age,
       tum.ethnicity,
       tum.deathdatebest,
       tum.diag_trust,
       rp.event_type,
       rp.event_property_1,
       rp.event_property_2,
       rp.event_property_3,
       rp.event_date,
       rp.event_end,
        --need to rank by tumourid too as could have 2 tumours diagnosed on same day so both of these would have rank of 1 and be pulled through, but only want 1 tumour per patient so that tumour will be randomly picked based on which of the 2 tumour ids is smaller
        rank () over (partition by pat.patientid order by tum.diagnosisdatebest, tum.tumourid asc) as rank
from av2021.at_patient_england@casref01 pat 
        left join av2021.at_tumour_england@casref01 tum on pat.patientid = tum.patientid
        left join analysisncr.at_rapid_pathway@cas2312 rp on pat.nhsnumber = rp.nhsnumber
and rp.event_type in (20, 24)
where tum.diagnosisyear = 2021
and tum.cascade_inci_flag = 1 --standard CAS exclusions
and tum.dco = 'N' --excluding those diagnosed on death certificate only 
and tum.dedup_flag = 1 --excluding duplicate records
and ((tum.deathdatebest is null) or (tum.diagnosisdatebest != tum.deathdatebest)) --excluding those diagnosed on the same day as death (not all captured by DCO = N)
and ((tum.deathdatebest is null) or (tum.deathdatebest - tum.diagnosisdatebest > 0)) --excluding those diagnosed after death (not all captured by DCO = N)
and substr(tum.site_icd10_o2_3char, 1, 1) = 'C' --include C codes only 
and tum.site_icd10_o2_3char <> 'C44') --exclude skin cancer
where rank = 1;


create table hna_pcsp_2021diags_eot_test compress basic nologging as 
select * from (
select a.*,
       c.nhsnumber as rp_nhsnumber,
       c.patientid as rp_patientid,       
       b.tumourid,
       b.diagnosisdatebest,
       b.diagnosisyear,
       b.site_icd10_o2_3char,
       b.stage_best,
       c.event_type,
       c.event_property_1,
       c.event_property_2,
       c.event_property_3,
       c.event_date,
       c.event_end,
        --need to rank by tumourid too as could have 2 tumours diagnosed on same day so both of these would have rank of 1 and be pulled through, but only want 1 tumour per patient so that tumour will be randomly picked based on which of the 2 tumour ids is smaller
        rank () over (partition by a.patientid order by b.diagnosisdatebest, b.tumourid asc) as rank
from av2021.at_patient_england@casref01 a 
        left join av2021.at_tumour_england@casref01 b on a.patientid = b.patientid
        left join analysisncr.at_rapid_pathway@cas2312 c on a.nhsnumber = c.nhsnumber
        left join (
            select t.patientid,
                   t.event_type,
                   t.event_date,
                   max(t.event_end) as max_event_end ---finding the latest treatment date 
            from analysisncr.at_rapid_pathway@cas2312 t 
            where t.event_type in (17,18,19,20,21,22,23,24,10,37,38,39,50,105)
            group by t.patientid,
                     t.event_type,
                     t.event_date
        ) t on b.patientid = t.patientid
where b.diagnosisyear = 2021
and c.event_type in (20, 24)
and b.cascade_inci_flag = 1 --standard CAS exclusions
and b.dco = 'N' --excluding those diagnosed on death certificate only 
and b.dedup_flag = 1 --excluding duplicate records
and ((b.deathdatebest is null) or (b.diagnosisdatebest != b.deathdatebest)) --excluding those diagnosed on the same day as death (not all captured by DCO = N)
and ((b.deathdatebest is null) or (b.deathdatebest - b.diagnosisdatebest > 0)) --excluding those diagnosed after death (not all captured by DCO = N)
and substr(b.site_icd10_o2_3char, 1, 1) = 'C' --include C codes only 
and b.site_icd10_o2_3char <> 'C44') --exclude skin cancer
where rank = 1;
