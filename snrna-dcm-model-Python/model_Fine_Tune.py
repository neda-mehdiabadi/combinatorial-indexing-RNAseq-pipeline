import datetime
import os
from geneformer import Classifier


current_date = datetime.datetime.now()
datestamp = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}{current_date.hour:02d}{current_date.minute:02d}{current_date.second:02d}"
datestamp_min = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}"

output_prefix = "cm_classifier_test"
output_dir = f"/data/mcri_heartv2/mcri_nedar/DCM-Foundation-Model-data/output_dir/{datestamp}"
os.makedirs(output_dir)

filter_data_dict={"cell_type":"Cardiomyocyte"}

training_args = {
    "num_train_epochs": 0.9,
    "learning_rate": 0.000804,
    "lr_scheduler_type": "polynomial",
    "warmup_steps": 1812,
    "weight_decay":0.258828,
    "per_device_train_batch_size": 12,
    "seed": 73,
}


cc = Classifier(classifier="cell",
                cell_state_dict = {"state_key": "disease", "states": "all"},
                filter_data=filter_data_dict,
                training_args=training_args,
                max_ncells=None,
                freeze_layers = 10,
                num_crossval_splits = 1,
                forward_batch_size=200,
                nproc=16)

train_ids = ["1472_1", "1678_1", "1539_2", "1558_1", "1582_1", "1300_2", "1504_2", "1558_2", "1371_2", "1437_1", "1430_2", "1702_1", "1678_2", "1516_2", "1702_2", "1437_2", "1617_1", "1561_2", "1539_1", "1304_1", "1515_1", "1430_1", "1718_1", "1472_2", "1561_1", "1515_2", "1516_1", "1617_2", "1600_1", "1600_2", "1371_1", "1304_2", "1300_1", "1718_2", "1582_2", "GSM4923404", "GSM4923403", "GSM4923401", "TWCM-11-78", "TWCM-13-104", "TWCM-11-74", "TWCM-13-152", "TWCM-13-181", "TWCM-13-102", "TWCM-14-173", "TWCM-13-208", "TWCM-10-68", "TWCM-11-42", "TWCM-LVAD3", "TWCM-13-198", "TWCM-13-1", "TWCM-13-285", "TWCM-11-103", "TWCM-11-264", "TWCM-11-104", "TWCM-13-101", "H_ZC-LVAD-1", "TWCM-11-41", "BO_H20_LVN", "BS_DS3_LV", "BO_H01_LV0", "HCAHeart7664652", "BS_DP2_LV", "HCAHeart7757636", "HCAHeart7664654", "H0015_LV", "BO_H25_LV0", "BO_H12_LV0", "HCAHeart7985086", "HCAHeart8287124", "BO_H39_LVW", "H0025_LV", "BS_DP1_LV", "BO_H80_LV0", "H0035_LV", "BO_H06_LV0", "BS_DT4_LV", "BO_H40_LV2", "BO_H24_LV0", "BO_H40_LVW", "BO_H27_LV0", "BO_H13_LV0", "HCAHeart7702873", "H0037_LV", "HCAHeart7664653", "BS_H26_LV", "BS_DL3_LV", "HCAHeart7698015", "BS_DO1_LV", "BO_H31_LV0", "BO_H02_LV0", "BO_H03_LV0", "BO_H37_LV0", "HCAHeart7835148", "BO_H40_LV1_rep", "BO_H14_LV0", "BO_H34_LVW", "BO_H26_LV0", "P75_snrna_rep1_GSM6165805", "P36_snrna_rep1_GSM6165799", "UK1_snrna_rep2_GSM6165813", "UK2_snrna_rep1_GSM6165809", "UK2_snrna_rep3_GSM6165811", "UK1_snrna_rep1_GSM6165812", "UK2_snrna_rep2_GSM6165810", "WU198LV_snrna_rep1_GSM6165815", "P36_snrna_rep2_GSM6165800", "WU198LV_snrna_rep2_GSM6165816", "P75_snrna_rep2_GSM6165806", "UK1_snrna_rep3_GSM6165814", "DCM4_GSM5606409", "Adult1_GSM4742854", "Adult2_GSM4742855", "DCM1_GSM5606406", "DCM3_GSM5606408"]
eval_ids = ["1547_1", "1549_2", "1547_2", "1540_2", "1540_1", "1358_1", "1549_1", "1358_2", "GSM4923402", "TWCM-LVAD2", "TWCM-11-93", "TWCM-13-235", "TWCM-11-256", "H_ZC-11-292", "TWCM-13-192", "TWCM-13-132", "TWCM-13-17", "BO_H33_LV0", "BO_H22_LV0", "BO_H19_LV0", "BS_DL2_LV", "BO_H79_LV0", "BO_H41_LV0", "DCM2_GSM5606407", "Young1_GSM5606413"]
test_ids = ["1290_2", "1290_1", "1606_1", "1610_1", "1610_2", "1622_2", "1603_2", "TWCM-10-5", "TWCM-13-36", "TWCM-11-3", "TWCM-13-96", "TWCM-13-80", "TWCM-13-168", "TWCM-11-192", "BO_H28_LV0", "HCAHeart7880860", "HCAHeart7829978", "IC_H02_LV0", "BO_H35_LV0", "BO_H36_LV0", "H0020_LV", "IC_H01_LV0", "Young3_GSM5606415", "Adult3_GSM4742856", "Young2_GSM5606414"]

train_test_id_split_dict = {"attr_key": "individual",
                            "train": train_ids+eval_ids,
                            "test": test_ids}

cc.prepare_data(input_data_file="/data/mcri_heartv2/mcri_nedar/DCM-Foundation-Model-data/output_dir_dataset/20250410_165samples_human_dcm_control_pediatric_adult_CMonly_4096_w_length.dataset",
                output_directory=output_dir,
                output_prefix=output_prefix,
                split_id_dict=train_test_id_split_dict)


train_valid_id_split_dict = {"attr_key": "individual",
                            "train": train_ids,
                            "eval": eval_ids}

all_metrics = cc.validate(model_directory="/data/mcri_heartv2/mcri_nedar/Geneformer/gf-12L-95M-i4096",
                          prepared_input_data_file=f"{output_dir}/{output_prefix}_labeled_train.dataset",
                          id_class_dict_file=f"{output_dir}/{output_prefix}_id_class_dict.pkl",
                          output_directory=output_dir,
                          output_prefix=output_prefix,
                          split_id_dict=train_valid_id_split_dict,
                          n_hyperopt_trials=8)