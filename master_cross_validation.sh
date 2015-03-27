#!/bin/bash

# written by Georg Zeller EMBL Heidelberg 2012-2015
version="0.1.0"

### R version

# TODO using Pauls R 3.0
LD_LIBRARY_PATH=/g/bork3/home/costea/gcc-4.8.2/lib64/:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
pathtoR="/g/bork3/home/costea/R-3.0.2/bin/"

### directories and data sets

# directory where source scripts are located
sourcedir="/g/bork4/zeller/dev/siamcat"
# directory where input data is located
datadir="/g/bork2/zeller/project_data/siamcat/original_data"
# working directory where all (intermediate) results will be stored
# (this is a prefix that will be modfied to a more specific subdirectory later)
workdir="/g/bork2/zeller/project_data/siamcat/results"

# data set-specific identifier tag
tag="baboon-groups_N48_tax-profile_metaphlan2-smg"
#tag="cancer-vs-healthy_study-pop-I_N141_tax_profile_mocat_bn_specI_clusters"
#tag="UC_vs_healthy163Spain_tax_profile_mocat_bn_curated_species"
#tag="c-diff_N182_tax-profile_16s"
#tag="Buccal_mucosa-vs-rest_N1431_tax-profile-QIIME-genus_16s"
#tag="Hard_palate-vs-rest_N1431_tax-profile-QIIME-genus_16s"
#tag="Keratinized_gingiva-vs-rest_N1431_tax-profile-QIIME-genus_16s"
#tag="Saliva-vs-rest_N1431_tax-profile-QIIME-genus_16s"
#tag="Subgingival_plaque-vs-rest_N1431_tax-profile-QIIME-genus_16s"
#tag="Supragingival_plaque-vs-rest_N1431_tax-profile-QIIME-genus_16s"
tag="Tongue_dorsum-vs-rest_N1431_tax-profile-QIIME-genus_16s"

# a date tag to identify / trace results
dt=`date +"%F_%Hh%Mm%S"`


### parameters
# TODO copy parameters here
COL="RdYlBu" # other options: "matlab", "RdYlBu", "RdBu", and other diverging schemes from ColorBrewer


### logic for switching optional pipeline steps on or off
PREPARE=true
VALIDATE=true
SELECT=false
CHECKMETA=true
FILTER=true
ADDMETAPRED=false
CHECKASSOC=true
NORM=true
SPLIT=true
TRAIN=true
APPLY=true
EVAL=true
INTERPRET=true

# the method used to build the linear model
LMMETHOD="lasso" # "lasso" # "lasso" for binary classification


### instance names of the modular components of the pipeline
validator="data_validator.r"
selector="sample_selector.r"
metachecker="confounder_check.r"
filter="generic_filter.r"
metapredadder="meta_predictor_adder.r"
assocchecker="association_check.r"
normalizer="generic_normalizer.r"
splitter="data_splitter.r"
builder="${LMMETHOD}_trainer.r"
predictor="${LMMETHOD}_predictor.r"
assessor="model_evaler.r"
interpretor="model_interpretor.r"
frozennormalizer="frozen_normalizer.r"


##############################
### START RUNNING PIPELINE ###
##############################

### PREPARE initial data
# prepare MASTER LOG file
masterlog="${workdir}/${tag}_${dt}/log_master_${tag}.txt"
echo -e "Master log will be written to $masterlog"

if $PREPARE
then
    if mkdir ${workdir}/${tag}_${dt}
    then
        workdir="${workdir}/${tag}_${dt}"
        echo -e "#master log\n#version ${version}\n#run on ${dt}\n#run in ${workdir}" > $masterlog
        echo -e "Successfully created working directory: ${workdir}"
    else
        echo -e "Failed to create working directory: ${workdir}"
        exit 1
    fi
    if cp ${datadir}/feat_${tag}.tsv ${workdir} \
        && cp ${datadir}/label_${tag}.tsv ${workdir}/
    then
        curr_feat="${workdir}/feat_${tag}"
        curr_label="${workdir}/label_${tag}"
        echo -e "Successfully copied initial data from ${datadir} to ${workdir}"
    else
        echo -e "Could not get initial data from ${datadir}"
        exit 1
    fi
    if cp ${datadir}/num_metadata_${tag}.tsv ${workdir}/
    then
        curr_meta="${workdir}/num_metadata_${tag}"
        echo -e "Successfully copied meta-data from ${datadir} to ${workdir}"
    else
        curr_meta="NULL"
        echo -e "Could not get meta-data from ${datadir}, but will continue anyway"
    fi
fi


### VALIDATE data
if $VALIDATE
then
    # output: validated feature file
    vld_feat="${curr_feat}_vld"
    # output: validated label file
    vld_label="${curr_label}_vld"
    # output: validated meta-data file
    vld_meta="${curr_meta}_vld"

    if [[ $curr_meta == "NULL" ]]
    then
        cmd="${validator} --srcdir=${sourcedir} --feat_in=${curr_feat}.tsv --feat_out=${vld_feat}.tsv --label_in=${curr_label}.tsv --label_out=${vld_label}.tsv"
    else
        cmd="${validator} --srcdir=${sourcedir} --feat_in=${curr_feat}.tsv --feat_out=${vld_feat}.tsv --label_in=${curr_label}.tsv --label_out=${vld_label}.tsv --metadata_in=${curr_meta}.tsv --metadata_out=${vld_meta}.tsv"
    fi
    echo "${pathtoR}Rscript $cmd" >> $masterlog
    if "${pathtoR}Rscript" $cmd > ${workdir}/log_${validator}.txt 2>&1
    then
        curr_feat=${vld_feat}
        curr_label=${vld_label}
        if [[ $curr_meta != "NULL" ]]
        then
            curr_meta=${vld_meta}
        fi
        echo "#succeeded" >> $masterlog
        echo "Successfully ran ${validator} in ${workdir}"
    else
        echo "#failed" >> $masterlog
        echo -e "Data validation failed!\nSee ${workdir}/log_${validator}.txt"
        exit 1
    fi
fi


#### SELECT examples / subjects based on meta-data
# (e.g. to better match cases and controls for reduced confounding or to stratify data based on meta-data such as gender)
if $SELECT
then
    # output: feature file with selected samples
    sel_feat="${curr_feat}_sel"
    # output: label file with selected samples
    sel_label="${curr_label}_sel"
    # output: label file with selected samples
    sel_meta="${curr_meta}_sel"
    # meta-variable to be used for sample selection
    filter_var="age"
    # samples can either be selected by specifying a range or a set of allowed values
    toggle_range_set="range"
    allowed_range="[60,120]"
    allowed_set="{55,56,57,58,59,60}"

# TODO doesn't work if $curr_meta == "NULL"

    if [ ${toggle_range_set}="range" ];
    then
        cmd="${selector} --srcdir=${sourcedir} --feat_in=${curr_feat}.tsv --feat_out=${sel_feat}.tsv --label_in=${curr_label}.tsv --label_out=${sel_label}.tsv --metadata_in=${curr_meta}.tsv --metadata_out=${sel_meta}.tsv --filter_var=${filter_var} --allowed_range=${allowed_range}"
    fi
    if [ ${toggle_range_set}="set" ];
    then
        cmd="${selector} --srcdir=${sourcedir} --feat_in=${curr_feat}.tsv --feat_out=${sel_feat}.tsv --label_in=${curr_label}.tsv --label_out=${sel_label}.tsv --metadata_in=${curr_meta}.tsv --metadata_out=${sel_meta}.tsv --filter_var=${filter_var} --allowed_set=${allowed_set}"
    fi

    echo "${pathtoR}Rscript $cmd" >> $masterlog
    if "${pathtoR}Rscript" $cmd > ${workdir}/log_${selector}.txt 2>&1
    then
        curr_feat=${sel_feat}
        curr_label=${sel_label}
        curr_meta=${sel_meta}
        echo "#succeeded" >> $masterlog
        echo "Successfully ran ${selector} in ${workdir}"
    else
        echo "#failed" >> $masterlog
        echo -e "Example selection failed!\nSee ${workdir}/log_${selector}.txt"
        exit 1
    fi
else
    sel_feat="${curr_feat}"
    sel_label="${curr_label}"
    sel_meta="${curr_meta}"
fi


### call script to CHECK METAdata for potential counfounders and biases
if $CHECKMETA
then
    if [[ $curr_meta == "NULL" ]]
    then
        echo "No meta-data found, thus skipping confounder check!"
    else
        # output: pdf file showing the results of confounder analysis
        conf_plot="${workdir}/confounder_plots_${tag}.pdf"

        cmd="${metachecker} --srcdir=${sourcedir} --label_in=${curr_label}.tsv --metadata_in=${curr_meta}.tsv --plot=${conf_plot}"
        echo "${pathtoR}Rscript $cmd" >> $masterlog
        if "${pathtoR}Rscript" $cmd > ${workdir}/log_${metachecker}.txt 2>&1
        then 
            echo "#succeeded" >> $masterlog
            echo -e "Successfully ran ${metachecker} in ${workdir}"
        else
            echo "#failed" >> $masterlog
            echo -e "Confounder check failed!\nSee ${workdir}/log_${metachecker}.txt"
            exit 1
        fi
    fi
fi


### call script to FILTER features in an unsupervised manner (i.e. without looking at the label information)
if $FILTER
then
    # output: filtered features
    filt_feat="${curr_feat}_filt"
    # filtering method
    filt_method="abundance"
    # cutoff on abundance or prevalence to remove features with very low abundance/prevalence (depends on filtering method)
    cutoff="1e-3"
    # option to remove the fraction of unmapped reads from the feature table
    rm_unmapped="TRUE"
    # convert features to relative abundances (i.e. proportions)
    recomp_prop="FALSE"

    cmd="${filter} --feat_in=${curr_feat}.tsv --feat_out=${filt_feat}.tsv --method=${filt_method} --cutoff=${cutoff} --rm_unmapped=${rm_unmapped} --recomp_prop=${recomp_prop}"
    echo "${pathtoR}Rscript $cmd" >> $masterlog
    if "${pathtoR}Rscript" $cmd > ${workdir}/log_${filter}.txt 2>&1
    then 
        curr_feat=$filt_feat
        echo "#succeeded" >> $masterlog
        echo -e "Successfully ran ${filter} in ${workdir}"
    else
        echo "#failed" >> $masterlog
        echo -e "Filtering failed!\nSee ${workdir}/log_${filter}.txt"
        exit 1
    fi
fi


### ADD METAdata as PREDictors into the model (the filtered feature table)
if $ADDMETAPRED
then
    if [[ $curr_meta == "NULL" ]]
    then
        echo "No meta-data found, thus cannot add any as features!"
    else
        # output: combined features
        comb_feat="${curr_feat}_comb"
        pred_names="age,bmi,gender"

        cmd="${metapredadder} --feat_in=${curr_feat}.tsv --metadata_in=${curr_meta}.tsv --feat_out=${comb_feat}.tsv --pred_names=${pred_names}"
        echo "${pathtoR}Rscript $cmd" >> $masterlog
        if "${pathtoR}Rscript" $cmd > ${workdir}/log_${metapredadder}.txt 2>&1
        then 
            curr_feat=$comb_feat
            echo "#succeeded" >> $masterlog
            echo -e "Successfully ran ${metapredadder} in ${workdir}"
        else
            echo "#failed" >> $masterlog
            echo -e "Adding predictors from metadata failed!\nSee ${workdir}/log_${metapredadder}.txt"
            exit 1
	fi
    fi
fi


### call script to CHECK for significant ASSOCiations between features (abundances) and labels
if $CHECKASSOC
then
    # output: pdf file showing the results of confounder analysis
    assoc_plot="${workdir}/associations_plots_${tag}.pdf"
    # correction for multiple testing
    mult_corr="FDR"
    # significance level (p-value cutoff)
    alpha="0.05"
    # fold-change cutoff (on absolute log10(FC))
    min_fc="0.0"
    # lower detection limit on rel. abundance
    detect_lim="1E-6"

    cmd="${assocchecker} --srcdir=${sourcedir} --label_in=${curr_label}.tsv --feat_in=${curr_feat}.tsv --plot=${assoc_plot} --mult_test=${mult_corr} --alpha=${alpha} --min_fc=${min_fc} --detect_limit=${detect_lim} --col_scheme=${COL}"
    echo "${pathtoR}Rscript $cmd" >> $masterlog
    if "${pathtoR}Rscript" $cmd > ${workdir}/log_${assocchecker}.txt 2>&1
    then 
        echo "#succeeded" >> $masterlog
        echo -e "Successfully ran ${assocchecker} in ${workdir}"
    else
        echo "#failed" >> $masterlog
        echo -e "Associations check failed!\nSee ${workdir}/log_${assocchecker}.txt"
        exit 1
    fi
fi


### call script to NORMalize features (again without looking at the labels!)
if $NORM
then
    # output: normalized features
    norm_feat="${curr_feat}_norm"
    # output: paramters of the normalization transformation
    norm_param_file="${workdir}/normalization_parameters.txt"
    # normalization method
    norm_method="log.std"
    # small positive value added to every feature before log-transformation
    log_n0="1E-6" # for species features this should be 10^-8
    # the type of vector norm (1-norm normalizes by sum of absolute entries, 2-norm by sum of squared entries)
    norm="2"
    # shall per-feature normalization be applied?
    n_feat="TRUE"
    # shall per-sample normalization be applied?
    n_sample="TRUE"
    # shall global normalization be applied?
    n_global="FALSE"
    # a quantity that is added to the denominator during standardization
    # (estimated as a quantile of the distribution of s.d. across all features)
    sd_min_q="0.1"

    cmd="${normalizer} --feat_in=${curr_feat}.tsv --feat_out=${norm_feat}.tsv --param_out=${norm_param_file} --method=${norm_method} --log_n0=${log_n0} --sd_min_quantile=${sd_min_q} --norm_feature=${n_feat} --norm_sample=${n_sample} --norm_global=${n_global} --vector_norm=${norm}"
    echo "${pathtoR}Rscript $cmd" >> $masterlog
    if "${pathtoR}Rscript" $cmd > ${workdir}/log_${normalizer}.txt 2>&1
    then 
        curr_feat=${norm_feat}
        echo "#succeeded" >> $masterlog
        echo -e "Successfully ran ${normalizer} in ${workdir}"
    else
        echo "#failed" >> $masterlog
        echo -e "Normalization failed!\nSee ${workdir}/log_${normalizer}.txt"
        exit 1
    fi
fi


### call script to SPLIT data for cross validation
# output of data split: training sets
train_sets="${workdir}/train_sets.tsv"
# output of data split: test sets
test_sets="${workdir}/test_sets.tsv"

if $SPLIT
then
    # number of cross validation folds (subsets)
    n_folds="10"
    # number of resampling rounds
    n_repeats="1"
    # shall cross validation be stratified?
    # (i.e. shall the relative occurrence of positive examples be balanced across subset?)
    stratify="TRUE"

    cmd="${splitter} --srcdir=${sourcedir} --label_in=${curr_label}.tsv --train_sets=${train_sets} --test_sets=${test_sets} --num_folds=${n_folds} --resample=${n_repeats} --stratify=${stratify}"
    echo "${pathtoR}Rscript $cmd" >> $masterlog
    if "${pathtoR}Rscript" $cmd > ${workdir}/log_${splitter}.txt 2>&1
    then 
        echo "Successfully ran ${splitter} in ${workdir}"
    else
        echo -e "Splitting data failed!\nSee ${workdir}/log_${splitter}.txt"
        exit 1
    fi
else
    if [ -f "$train_sets" ] && [ -f "$test_sets" ]
    then
        echo -e "Reusing training and test sets (${train_sets}, ${test_sets})"
    else
        echo -e "No training and test set found! Consider enabling data split module (set \$SPLIT=true)"
        exit 1
    fi
fi


### call script to TRAIN a model for each data split
if $TRAIN
then
    # output: trained model
    trained_model="${workdir}/trained_${LMMETHOD}_models_${tag}.tsv"
    # number of cross validation folds (subsets) for internal model selection
    n_folds="5"
    # shall cross validation be stratified?
    # (i.e. shall the relative occurrence of positive examples be balanced across subset?)
    stratify="TRUE"
    # Criterion used for evaluation in model selection in order to determine the best model
    # (options: 'acc', 'auroc', 'auprc', 'auroc.2')
    sel_criterion="auroc"
    # required minimum number of nonzero coefficients for a model to be considered in model selection
    min_nonzero_coeff=1

    cmd="${builder} --srcdir=${sourcedir} --label_in=${curr_label}.tsv --feat_in=${curr_feat}.tsv --model=${trained_model} --train_sets=${train_sets} --num_folds=${n_folds} --stratify=${stratify} --sel_criterion=${sel_criterion} --min_nonzero_coeff=${min_nonzero_coeff}"
    echo "${pathtoR}Rscript $cmd" >> $masterlog
    if "${pathtoR}Rscript" $cmd > ${workdir}/log_${builder}.txt 2>&1
    then
        echo "#succeeded" >> $masterlog
        echo -e "Successfully ran ${builder} in ${workdir}"
    else
        echo "#failed" >> $masterlog
        echo -e "Training ${LMMETHOD} models failed!\nSee ${workdir}/log_${builder}.txt"
        exit 1
    fi
fi

### call script to APPLY the model to make predictions
if $APPLY
then
    # output: predictions
    pred="${workdir}/${LMMETHOD}_predictions_${tag}.tsv"

    cmd="${predictor} --srcdir=${sourcedir} --label_in=${curr_label}.tsv --feat_in=${curr_feat}.tsv --test_sets=${test_sets} --model=${trained_model} --pred=${pred}"
    echo "${pathtoR}Rscript $cmd" >> $masterlog
    if "${pathtoR}Rscript" $cmd > ${workdir}/log_${predictor}.txt 2>&1
    then
        echo "#succeeded" >> $masterlog
        echo -e "Successfully ran ${predictor} in ${workdir}"
    else
        echo "#failed" >> $masterlog
        echo -e "Applying model failed!\nSee ${workdir}/log_${predictor}.txt"
        exit 1
    fi
fi


### call script to EVALuate model (test) performance
if $EVAL
then
    # output: pdf file with evaluation graphs (ROC & precision-recall curves)
    eval_plot="${workdir}/evaluation_plots_${tag}.pdf"

    cmd="${assessor} --srcdir=${sourcedir} --label=${curr_label}.tsv --pred=${pred} --plot=${eval_plot}"
    echo "${pathtoR}Rscript $cmd" >> $masterlog
    if "${pathtoR}Rscript" $cmd > ${workdir}/log_${assessor}.txt 2>&1
    then
        echo "#succeeded" >> $masterlog
        echo -e "Successfully ran ${assessor} in ${workdir}"
    else
        echo "#failed" >> $masterlog
        echo -e "Assessing model performance failed!\nSee ${workdir}/log_${assessor}.txt"
        exit 1
    fi
fi


### call script that extracts model properties for INTERPRETATION
if $INTERPRET
then
    # output: pdf file showing model and marker properties
    model_plot="${workdir}/model_plots_${tag}.pdf"
    # optionally include heatmap with metadata into the model graphs
    plot_metadata=true
    if $plot_metadata
    then
        metadata=${curr_meta}.tsv
    else
        metadata="NULL"
    fi

    cmd="${interpretor} --srcdir=${sourcedir} --feat=${curr_feat}.tsv --label=${curr_label}.tsv --meta=${metadata} --model=${trained_model} --pred=${pred} --plot=${model_plot} --col_scheme=${COL}"
    echo "${pathtoR}Rscript $cmd" >> $masterlog
    if "${pathtoR}Rscript" $cmd > ${workdir}/log_${interpretor}.txt 2>&1
    then
        echo "#succeeded" >> $masterlog
        echo "Successfully ran ${interpretor} in ${workdir}"
    else
        echo "#failed" >> $masterlog
        echo -e "Interpreting model failed!\nSee ${workdir}/log_${interpretor}.txt"
        exit 1
    fi
fi

echo "#master finished successfully" >> $masterlog
