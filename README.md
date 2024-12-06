# HNA_PCSP_coverage
NCRAS-Macmillan project on the coverage of holistic needs assessments (HNAs) and personalised care and support plans (PCSPs)

This project explores data on HNAs and PCSPs from the Cancer Outcomes and Services Dataset (COSD). These data flow into the Rapid Cancer Registrations Datset (RCRD), which is linked to cancer registration data, allowing for patients to be linked to their HNA and PCSP records. The project linked patients diagnosed in 2021 to their offers of HNAs and PCSPs in the 2 years post-diagnosis. It cleans and prepares the data, to create a patient-level dataset containing the date of the earliest offer of each event per patient, and the total number of each recorded for the patient in the 2 years post-diagnosis. It then conducts various analyses of this dataset, exploring data quality and completeness, delivery of HNAs and PCSPs across the cancer pathway and in relation to diagnosis and death dates, and by demographic, clinical, and geographic variables. The project produces HTML files using Quarto, containing the results and relevant background and narrative. The project code can also be used for other similar analysis, or to rerun in future years. 

This project is part of the partnership between the National Cancer Registration and Analysis Service (NCRAS) in NDRS, and Macmillan Cancer Support. 

Project status: 
This project is closed. The code may be further adapted and used for similar projects and/or follow on analyses, but there are currently no specific plans for this. 


Point of contact: 
Lizzie Augarde lizzie994@googlemail.com as the original analyst, Carolynn Gildea carolynn.gildea@nhs.net for further NDRS work using this project


Data requirements: 
- cancer registrations in AV202 and related tables for IMD, cancer site 
- Rapid Cancer Registrations Dataset pathway table in CAS snapshot
- extract of Macmillan eHNA data, stored on NDRS drives


Outputs: 
The project contains Quarto files to output HTML slide decks of the results. Quarto files for 3 slide decks are included, containing different parts of the analysis. 


Prerequisites:
None beyond standard Level 2 CAS access


How to install and use:
Clone the Github repo and run the scripts as required for analysis purposes.


License:
MIT (see license file)


Other legal or regulatory requirements:
None


Other notes:
There are some additional scripts in this project containing exploratory and supplementary analysis. These are not QA'd and were used for project development and results interpretation. 



