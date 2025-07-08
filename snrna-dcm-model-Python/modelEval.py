from geneformer import Classifier

output_prefix = "cm_classifier_test"
output_dir = "/data/mcri_heartv2/mcri_nedar/DCM-Foundation-Model-data/output_dir/run-4d746651_modelEval"

filter_data_dict={"cell_type":"Cardiomyocyte"}

training_args = {
    "num_train_epochs": 1,
    "learning_rate": 0.00022,
    "lr_scheduler_type": "polynomial",
    "warmup_steps": 303,
    "weight_decay":0.28786,
    "per_device_train_batch_size": 12,
    "seed": 64,
}


cc = Classifier(classifier="cell",
                cell_state_dict = {"state_key": "disease", "states": "all"},
                filter_data=filter_data_dict,
                training_args=training_args,
                forward_batch_size=64,
                nproc=16)


all_metrics_test = cc.evaluate_saved_model(
        model_directory="/data/mcri_heartv2/mcri_nedar/DCM-Foundation-Model-data/output_dir/250218173556/250218_geneformer_cellClassifier_cm_classifier_test/ksplit1/run-4d746651/checkpoint-21688",
        id_class_dict_file="/data/mcri_heartv2/mcri_nedar/DCM-Foundation-Model-data/output_dir/250218173556/cm_classifier_test_id_class_dict.pkl",
        test_data_file="/data/mcri_heartv2/mcri_nedar/DCM-Foundation-Model-data/output_dir/250218173556/cm_classifier_test_labeled_test.dataset",
        output_directory=output_dir,
        output_prefix=output_prefix,
    )


cc.plot_conf_mat(
        conf_mat_dict={"Geneformer": all_metrics_test["conf_matrix"]},
        output_directory=output_dir,
        output_prefix=output_prefix,
        custom_class_order=["Control","DCM"],
)

cc.plot_predictions(
    predictions_file=f"{output_dir}/{output_prefix}_pred_dict.pkl",
    id_class_dict_file="/data/mcri_heartv2/mcri_nedar/DCM-Foundation-Model-data/output_dir/250218173556/cm_classifier_test_id_class_dict.pkl",
    title="disease",
    output_directory=output_dir,
    output_prefix=output_prefix,
    custom_class_order=["Control","DCM"],
)