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