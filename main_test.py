import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import decoupler as dc
import os
#from importlib.resources import files



#Will not work until package directory is created, please enter regulon csv manually.
#def load_dorothea_regulon (organism):
#    if (organism == "human"):
#        dorothea_regulon_human = files("LR2TF_py.data").joinpath("human_dorothea_reg.csv").pd.read_csv(index_col=0)
#        regulon =  dorothea_regulon_human.loc[dorothea_regulon_human["confidence"].isin(["A","B", "C", "D"])]
#        regulon = pd.DataFrame.rename(regulon, columns={"source" : "tf"})
#
#    elif (organism == "mouse"):
#        dorothea_regulon_mouse = files("LR2TF_py.data").joinpath("mouse_dorothea_reg.csv").pd.read_csv(index_col=0)
#        regulon =  dorothea_regulon_mouse.loc[dorothea_regulon_mouse["confidence"].isin(["A","B", "C", "D"])]
#        regulon = pd.DataFrame.rename(regulon, columns={"source" : "tf"})
#        
#    else:
#        print("Only human and mouse regulons can be loaded by default!")
#    return regulon


def validate_input_arguments (arguments_list):
    if arguments_list["out_path"] is None:
        print("Please provide an output path")
    elif arguments_list["out_path"][-1] != "/":
        arguments_list["out_path"] = arguments_list["out_path"] + "/"

    if arguments_list["celltype"] is None:
        print("Please provide the name of the metadata field containing cell type annotations")

    if arguments_list["condition"] is None:
        print("Please provide the name of the metadata field containing condition annotations")

    if arguments_list["organism"] is None:
        arguments_list["organism"] = "human"

    if arguments_list["comparison_list"] is None:
        arguments_list["comparison_list"] = np.nan

    if arguments_list["logfc"] is None:
        arguments_list["logfc"] = 0.0

    if arguments_list ["pval"] is None:
        arguments_list["pval"] = 0.05

#    if arguments_list["reg"] is None:
#        arguments_list["reg"] = load_dorothea_regulon(arguments_list["organism"])

    elif isinstance(arguments_list["reg"], str):
        arguments_list["reg"] = pd.read_csv(arguments_list["reg"], index_col=0)
        arguments_list["reg"] = pd.DataFrame.rename(arguments_list["reg"], columns={"source" : "tf"})

    if not "tf" in arguments_list["reg"] and "target" in arguments_list["reg"] and "weight"in arguments_list["reg"]:
        raise Exception("Not all necessary columns found in regulon table! Please make sure that the regulon has the columns source, target and weight!")
    
    return(arguments_list)


def AverageExpression(anndataobject, celltype = None, outpath = None):
    gene_ids = anndataobject.var.index.values
    obs = anndataobject[:,gene_ids].X.toarray()
    sub_object = pd.DataFrame(obs,columns=gene_ids,index= anndataobject.obs[celltype])
    sub_object = sub_object.groupby(level=0, observed=False).mean()
    sub_object.T.to_csv(outpath + "average_gene_expression_by_cluster.csv")

    return sub_object.T


def tf_activity_analysis (anndataobject, tf_activities = None, arguments_list = None):
    
    if (isinstance(anndataobject, str)):
        anndataobject = ad.read_h5ad(anndataobject)

    arguments_list = validate_input_arguments(arguments_list)

    if not os.path.isdir(arguments_list["out_path"]):
        os.mkdir(arguments_list["out_path"])
        tf_path = arguments_list["out_path"] + "TF_results/"
        os.mkdir(tf_path)
    else:
        tf_path = arguments_list["out_path"] + "TF_results/"

    #skipped tf activities part. ignore extra tf data from decoupler

    anndataobject.obs["condition"] = arguments_list["condition"] 
    anndataobject.obs["cell_type"] = arguments_list["celltype"]
    anndataobject.obs["comparison_list"] = arguments_list["comparison_list"]

    if not np.isnan(arguments_list["comparison_list"]):
        if len(arguments_list["comparison_list"]) > 0 & len(anndataobject.obs["comparison_list"]) < 2:
            arguments_list["comparison_list"] <- np.nan
            print("Only one condition was found in the data, although a list of comparisons was provided. The analyses are performed only for the present condition!")

    #code for single condition  analysis

    if np.isnan(arguments_list["comparison_list"]):
        anndataobject.uns["tf_annotation"] = pd.DataFrame({"result_list" : [],
        "gene_expression_list" : [],
        "CTR_cluster_list" : [],
        "intranet_cluster_list" : []})

    #anndataobject_list = split by condition, skipped for now
    sub_object = anndataobject

    sub_object.uns["Average_Expression"] = AverageExpression(sub_object, celltype = arguments_list["celltype"], outpath= arguments_list["out_path"])

    return(sub_object)


#Please edit the following argument list to suit your own data.

#Enter your regulon csv path for reg

sub_object = tf_activity_analysis(anndataobject= "LR2TF_test_run/anndata_object.h5ad", arguments_list= {"out_path" : "folder2", "celltype" : "new_annotation", "condition" : "control", "organism" : None, "comparison_list" : None, "logfc" : None, "pval" : None, "reg" : "human_dorothea_reg.csv"})
