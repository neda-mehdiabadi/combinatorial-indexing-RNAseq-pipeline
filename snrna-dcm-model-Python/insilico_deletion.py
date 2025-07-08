from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
from geneformer import EmbExtractor
from multiprocess import set_start_method

if __name__ == '__main__':

    set_start_method('spawn')

    cell_states_to_model={"state_key": "disease", 
                          "start_state": "DCM", 
                          "goal_state": "Control"}

    embex = EmbExtractor(model_type= "CellClassifier", 
                         num_classes=2,
                         filter_data=None,
                         max_ncells=10000,
                         emb_layer=0,
                         summary_stat="exact_mean",
                         forward_batch_size=10,
                         nproc=60)

    state_embs_dict = embex.get_state_embs(cell_states_to_model,
                            "output/fine_tune_CM/250412_geneformer_cellClassifier_cm_classifier_test/ksplit1/run-9ab57793/checkpoint-21688", 
                            "data/datasetG1/20250412_25samples_human_dcm_control_pediatric_adult_TEST_CMonly_4096_w_length.dataset",
                            "output/in_silico_deletion_CM/state_embs",
                                           "insilico_embs")


    isp = InSilicoPerturber(perturb_type="delete",
                            perturb_rank_shift=None,
                            genes_to_perturb= "all",
                            combos=0,
                            anchor_gene= None,
                            model_type="CellClassifier", 
                            num_classes=2,
                            emb_mode="cls",
                            cell_emb_style="mean_pool",
                            filter_data=None,
                            cell_states_to_model=cell_states_to_model,
                            state_embs_dict=state_embs_dict,
                            max_ncells=10000,
                            emb_layer=0,
                            forward_batch_size=200,
                            nproc=60)

    isp.perturb_data("output/fine_tune_CM/250412_geneformer_cellClassifier_cm_classifier_test/ksplit1/run-9ab57793/checkpoint-21688", 
                    "data/datasetG1/20250412_25samples_human_dcm_control_pediatric_adult_TEST_CMonly_4096_w_length.dataset",
                    "output/in_silico_deletion_CM/isp",
                    "insilico_isp")

    ispstats = InSilicoPerturberStats(mode="goal_state_shift",
                                      genes_perturbed= "all",
                                      combos=0,
                                      anchor_gene= None,
                                      cell_states_to_model=cell_states_to_model)
 
    ispstats.get_stats("output/in_silico_deletion_CM/isp", 
                       None,
                       "output/in_silico_deletion_CM/isp_stats",
                       "insilico_isp_stats")

#for aggregate data use the following: 

#isp = InSilicoPerturber(perturb_type="delete",
#                        perturb_rank_shift=None,
#                        genes_to_perturb= ["ENSG00000152061"],
#                        combos=0,
#                        anchor_gene=None,
#                        model_type="CellClassifier", 
#                       num_classes=2,
#                        emb_mode="cls_and_gene",
#                        cell_emb_style="mean_pool",
#                        filter_data=None,
#                        cell_states_to_model=None,
#                        state_embs_dict=state_embs_dict,
#                        max_ncells=40,
#                        emb_layer=0,
#                        forward_batch_size=200,
#                        nproc=1)

#ispstats = InSilicoPerturberStats(mode="aggregate_data",
#                                  genes_perturbed= ["ENSG00000152061"],
#                                  combos=0,
#                                  anchor_gene= None,
#                                  cell_states_to_model=None)



