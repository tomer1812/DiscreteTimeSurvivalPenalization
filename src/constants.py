ADMISSION_TIME_COL = 'admittime'
DISCHARGE_TIME_COL = 'dischtime'
DEATH_TIME_COL = 'deathtime'
ED_REG_TIME = 'edregtime'
ED_OUT_TIME = 'edouttime'
AGE_COL = 'anchor_age' 
GENDER_COL = 'gender'

AGE_BINS = list(range(0, 125, 5))
AGE_LABELS = [f'{AGE_BINS[a]}' for a in range(len(AGE_BINS)-1)]

font_sz = 14
title_sz = 18

YEAR_GROUP_COL = 'anchor_year_group'
SUBSET_YEAR_GROUP = '2017 - 2019'
SUBJECT_ID_COL = 'subject_id'
ADMISSION_ID_COL = 'hadm_id'
ADMISSION_TYPE_COL = 'admission_type'
CHART_TIME_COL = 'charttime'
STORE_TIME_COL = 'storetime'
LOS_EXACT_COL = 'LOS exact'
LOS_DAYS_COL = 'LOS days'
ADMISSION_LOCATION_COL = 'admission_location'
DISCHARGE_LOCATION_COL = 'discharge_location'
RACE_COL = 'race'
INSURANCE_COL = 'insurance'
ADMISSION_TO_RESULT_COL = 'admission_to_result_time'
ADMISSION_AGE_COL = 'admission_age'
ADMISSION_YEAR_COL = 'admission_year'
ADMISSION_COUNT_COL = 'admissions_count'
ITEM_ID_COL = 'itemid'
NIGHT_ADMISSION_FLAG = 'night_admission' 
MARITAL_STATUS_COL = 'marital_status'
STANDARDIZED_AGE_COL = 'standardized_age'
COEF_COL = '   coef   '
STDERR_COL = ' std err '
DIRECT_IND_COL = 'direct_emrgency_flag'
PREV_ADMISSION_IND_COL = 'last_less_than_diff'
ADMISSION_COUNT_GROUP_COL = ADMISSION_COUNT_COL + '_group'


DISCHARGE_REGROUPING_DICT = {
    'HOME': 'HOME',
    'HOME HEALTH CARE': 'HOME',
    'SKILLED NURSING FACILITY': 'FURTHER TREATMENT',
    'DIED': 'DIED',
    'REHAB': 'HOME',
    'CHRONIC/LONG TERM ACUTE CARE': 'FURTHER TREATMENT',
    'HOSPICE': 'FURTHER TREATMENT',
    'AGAINST ADVICE': 'CENSORED',
    'ACUTE HOSPITAL': 'FURTHER TREATMENT',
    'PSYCH FACILITY': 'FURTHER TREATMENT',
    'OTHER FACILITY': 'FURTHER TREATMENT',
    'ASSISTED LIVING': 'HOME',
    'HEALTHCARE FACILITY': 'FURTHER TREATMENT',
}

RACE_REGROUPING_DICT = {
    'WHITE': 'WHITE',
    'UNKNOWN': 'OTHER',
    'BLACK/AFRICAN AMERICAN': 'BLACK',
    'OTHER': 'OTHER',
    'ASIAN': 'ASIAN',
    'WHITE - OTHER EUROPEAN': 'WHITE',
    'HISPANIC/LATINO - PUERTO RICAN': 'HISPANIC',
    'HISPANIC/LATINO - DOMINICAN': 'HISPANIC',
    'ASIAN - CHINESE': 'ASIAN',
    'BLACK/CARIBBEAN ISLAND': 'BLACK',
    'BLACK/AFRICAN': 'BLACK',
    'BLACK/CAPE VERDEAN': 'BLACK',
    'PATIENT DECLINED TO ANSWER': 'OTHER',
    'WHITE - BRAZILIAN': 'WHITE',
    'PORTUGUESE': 'HISPANIC', 
    'ASIAN - SOUTH EAST ASIAN': 'ASIAN',
    'WHITE - RUSSIAN': 'WHITE',
    'ASIAN - ASIAN INDIAN': 'ASIAN',
    'WHITE - EASTERN EUROPEAN': 'WHITE',
    'AMERICAN INDIAN/ALASKA NATIVE': 'OTHER',
    'HISPANIC/LATINO - GUATEMALAN': 'HISPANIC',
    'HISPANIC/LATINO - MEXICAN': 'HISPANIC',
    'HISPANIC/LATINO - SALVADORAN': 'HISPANIC',
    'SOUTH AMERICAN': 'HISPANIC',
    'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER': 'OTHER',
    'HISPANIC/LATINO - COLUMBIAN': 'HISPANIC',
    'HISPANIC/LATINO - CUBAN': 'HISPANIC',
    'ASIAN - KOREAN': 'ASIAN',
    'HISPANIC/LATINO - HONDURAN': 'HISPANIC',
    'HISPANIC/LATINO - CENTRAL AMERICAN': 'HISPANIC',
    'UNABLE TO OBTAIN': 'OTHER',
    'HISPANIC OR LATINO': 'HISPANIC'
}

# 'MCH': 'Mean Cell Hemoglobin', 
# 'MCHC': 'Mean Cell Hemoglobin Concentration',
table1_rename_columns = {
    'AnionGap': 'Anion Gap', 
    'Bicarbonate': 'Bicarbonate', 
    'CalciumTotal': 'Calcium Total', 
    'Chloride': 'Chloride', 
    'Creatinine': 'Creatinine',
    'Glucose': 'Glucose', 
    'Magnesium': 'Magnesium', 
    'Phosphate': 'Phosphate', 
    'Potassium': 'Potassium', 
    'Sodium': 'Sodium',
    'UreaNitrogen': 'Urea Nitrogen', 
    'Hematocrit': 'Hematocrit', 
    'Hemoglobin': 'Hemoglobin', 
    'MCH': 'MCH', 
    'MCHC': 'MCHC', 
    'MCV': 'MCV',
    'PlateletCount': 'Platelet Count', 
    'RDW': 'RDW', 
    'RedBloodCells': 'Red Blood Cells', 
    'WhiteBloodCells': 'White Blood Cells',
    NIGHT_ADMISSION_FLAG: 'Night Admission', 
    GENDER_COL: 'Sex', 
    DIRECT_IND_COL: 'Direct Emergency',
    PREV_ADMISSION_IND_COL: 'Previous Admission This Month', 
    ADMISSION_AGE_COL: 'Admission Age', 
    INSURANCE_COL: 'Insurance', 
    MARITAL_STATUS_COL: 'Marital Status',
    RACE_COL: 'Race', 
    ADMISSION_COUNT_GROUP_COL: 'Admissions Number', 
    LOS_DAYS_COL: 'LOS (days)', 
    DISCHARGE_LOCATION_COL: 'Discharge Location'
}

table1_rename_sex = {0: 'Male', 1: 'Female'}
table1_rename_race = {'ASIAN': 'Asian', 'BLACK': 'Black', 'HISPANIC': 'Hispanic', 'OTHER': 'Other', 'WHITE': 'White'}
table1_rename_marital = {'SINGLE': 'Single', 'MARRIED': 'Married', 'DIVORCED': 'Divorced', 'WIDOWED': 'Widowed'}
table1_rename_yes_no = {0: 'No', 1: 'Yes'}
table1_rename_normal_abnormal = {0: 'Normal', 1: 'Abnormal'}
table1_rename_discharge = {1: 'Home', 2: 'Further Treatment', 3: 'Died', 0: 'Censored'} 


rename_beta_index = {
 'AdmsCount 2': 'Admissions Number 2',
 'AdmsCount 3up': 'Admissions Number 3+',
 'AnionGap': 'Anion Gap',
 'Bicarbonate': 'Bicarbonate',
 'CalciumTotal': 'Calcium Total',
 'Chloride': 'Chloride',
 'Creatinine': 'Creatinine',
 'Ethnicity BLACK': 'Ethnicity Black',
 'Ethnicity HISPANIC': 'Ethnicity Hispanic',
 'Ethnicity OTHER': 'Ethnicity Other',
 'Ethnicity WHITE': 'Ethnicity White',
 'Glucose': 'Glucose',
 'Hematocrit': 'Hematocrit',
 'Hemoglobin': 'Hemoglobin',
 'Insurance Medicare': 'Insurance Medicare',
 'Insurance Other': 'Insurance Other',
 'MCH': 'MCH',
 'MCHC': 'MCHC',
 'MCV': 'MCV',
 'Magnesium': 'Magnesium',
 'Marital MARRIED': 'Marital Married',
 'Marital SINGLE': 'Marital Single',
 'Marital WIDOWED': 'Marital Widowed',
 'Phosphate': 'Phosphate',
 'PlateletCount': 'Platelet Count',
 'Potassium': 'Potassium',
 'RDW': 'RDW',
 'RedBloodCells': 'Red Blood Cells',
 'Sodium': 'Sodium',
 'UreaNitrogen': 'Urea Nitrogen',
 'WhiteBloodCells': 'White Blood Cells',
 'direct emrgency flag': 'Direct Emergency',
 'gender': 'Sex',
 'last less than diff': 'Recent Admission',
 'night admission': 'Night Admission',
 'standardized age': 'Standardized Age',
}

beta_units = {
    'Admissions Number 2': '2',
    'Admissions Number 3+': '3+',
    'Anion Gap': 'Abnormal',
    'Bicarbonate': 'Abnormal',
    'Calcium Total': 'Abnormal',
    'Chloride': 'Abnormal',
    'Creatinine': 'Abnormal',
    'Ethnicity Black': 'Black',
    'Ethnicity Hispanic': 'Hispanic',
    'Ethnicity Other': 'Other',
    'Ethnicity White': 'White',
    'Glucose': 'Abnormal',
    'Hematocrit': 'Abnormal',
    'Hemoglobin': 'Abnormal',
    'Insurance Medicare': 'Medicare',
    'Insurance Other': 'Other',
    'MCH': 'Abnormal',
    'MCHC': 'Abnormal',
    'MCV': 'Abnormal',
    'Magnesium': 'Abnormal',
    'Marital Married': 'Married',
    'Marital Single': 'Single',
    'Marital Widowed': 'Widowed',
    'Phosphate': 'Abnormal',
    'Platelet Count': 'Abnormal',
    'Potassium': 'Abnormal',
    'RDW': 'Abnormal',
    'Red Blood Cells': 'Abnormal',
    'Sodium': 'Abnormal',
    'Urea Nitrogen': 'Abnormal',
    'White Blood Cells': 'Abnormal',
    'Direct Emergency': 'Yes',
    'Sex': 'Female',
    'Recent Admission': 'Yes',
    'Night Admission': 'Yes',
    'Standardized Age': '',
}