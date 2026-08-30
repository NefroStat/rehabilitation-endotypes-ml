**1. Data Availability Statement**

The raw clinical and biochemical dataset analyzed in this repository is restricted from open public access to ensure compliance with medical confidentiality and ethical standards regarding military personnel data. This restriction is mandated by the Bioethics Committee of Danylo Halytsky Lviv National Medical University (Protocol No. 9, dated September 17, 2025).  Researchers seeking to replicate the analysis or conduct secondary studies may request access to the de-identified dataset by contacting the corresponding author (Vladyslav Bardash, v.bardash@1tmolviv.com). Access will be granted subject to the execution of a formal data-sharing agreement and institutional ethical approval.



**2. Script Execution Framework**

To verify the computational pipeline (scripts/01\_data\_preprocessing.R through 06\_data\_visualization.R) without the original clinical dataset, researchers must construct a synthetic CSV file named markers\_results\_comma.csv. The file must strictly adhere to the column nomenclature and data types specified in the dictionary below to ensure algorithmic compatibility.



**3. Data Dictionary**

The input dataset structure relies on the following parameters. Continuous biochemical variables represent serum concentrations quantified via enzyme-linked immunosorbent assay (ELISA).

Variable Name       Data Type           Description

\------------------------------------------------------------------------------------------------------------------

ID                  Integer             Unique anonymized patient identifier.

Вік                 Numeric             Patient age in years (cohort inclusion range: 20-60 years).

ISS                 Numeric             Injury Severity Score (control variable).

NISS                Numeric             New Injury Severity Score (primary baseline anatomical severity metric).

Дата.травми         Date (DD.MM.YYYY)   Date of the initial combat trauma.

Дата                Date (DD.MM.YYYY)   Date of peripheral venous blood sampling (minimum 3 months post-injury).

HSPG2               Numeric             Heparan sulfate proteoglycan 2 serum concentration.

GPC1                Numeric             Glypican-1 serum concentration.

SDC1                Numeric             Syndecan-1 serum concentration.

TIMP4               Numeric             Tissue inhibitor of metalloproteinases 4 serum concentration.

MMP9                Numeric             Matrix metalloproteinase-9 serum concentration.

MMP2                Numeric             Matrix metalloproteinase-2 serum concentration.

eNOS                Numeric             Endothelial nitric oxide synthase serum concentration.

Fibronectin         Numeric             Fibronectin serum concentration.

