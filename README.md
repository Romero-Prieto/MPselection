ata Access: The data for this analysis was obtained from the Demographic and Health Surveys (DHS) program. The DHS data is publicly available and can be requested directly from the DHS website.
Conditions for Access: To access the data, users must create an account, submit a data request, and specify the purpose of the research. Upon approval, datasets are made available for download in multiple formats (e.g., Stata, CSV). DHS data is subject to a license agreement, which includes restrictions on sharing and mandates appropriate citation of the data source.

Number of Records: The number of records (observations) used in the analysis varies based on the datasets selected (e.g., individual or household datasets). The sample sizes for each country/survey are derived from the dhs_datasets and vary by the survey ID and year. In the mobile ownership analysis, the sample sizes are specific to the IR (individual recode) and PR (household recode) datasets.
Variables Used in the Analysis:
	1	Demographic Information:
	◦	v021: Primary sampling unit (PSU) used for survey stratification.
	◦	v022: Stratification variable for survey design.
	◦	v005: Sample weight for respondents, divided by 1,000,000 to create normalized weights.
	◦	v169a: Mobile phone ownership indicator for women.
	◦	v190: Wealth quintile of the household.
	◦	v106: Educational level of the respondent (none, primary, secondary, higher).
	◦	v501: Marital status (single, married, cohabiting, divorced, widowed).
	◦	v011: Date of birth in CMC (Century Month Code).
	◦	v008: Date of interview in CMC, used to calculate current age.
	◦	v140: Type of residence (urban/rural).
	◦	v139: Region of residence (geographical area).
	2	Mortality Analysis:
	◦	mm1: Sibling sex (used for sibling mortality estimates).
	◦	mm2: Sibling survival status.
	◦	mm4: Month of sibling birth.
	◦	mm8: Month of sibling death (if applicable).
	◦	mmidx: Sibling ID within the household, used to index siblings in mortality analysis.
 
o A basic meta-description of these variables
  Demographic variables: Used to characterize the population (age, sex, education, marital status) and to apply sample weights, stratification, and clustering corrections.
  Mobile ownership (v169a): Binary variable (1 for owning a mobile phone, 0 for not) used to assess socio-economic and demographic factors associated with mobile phone access.
  Mortality variables (mm1, mm2, mm4, mm8): Used to estimate adult mortality rates based on sibling survival histories, by age group and gender.
