from goatools.base import get_godag
import numpy as np
from collections import namedtuple
from sklearn import metrics
import sys
import random
from ALL_EC_GOTERMS_IEA_MF import *
from ALL_EC_GOTERMS_NonIEA_MF import *
from ALL_EC_GOTERMS_IEA_BP import *
from ALL_EC_GOTERMS_NonIEA_BP import *
from GOntoSim import (
    Semantic_Value,
    Agglomerative_Clustering,
    AgglomerativeClusteringIC,
    purity_score,
    ARI_Score,
    AMI_Score,
    semantic_value_resnik_lin_IEA_MF,
    semantic_value_resnik_lin_nonIEA_MF,
)

go = get_godag("go-basic.obo", optional_attrs={"relationship"})

# Relationship Semantic Contribution
is_a = 0.8
part_of = 0.6
regulates = 0.7
negatively_regulates = 0.7
positively_regulates = 0.7
reg = "reg0.7"
Gene = namedtuple("Gene", "GeneName GOAnnotations")
c = 0.67

### New Functions
similarities_down = {}


def main():

    # if sys.argv[1] == "wang":
    #     method = "wang"
    # elif sys.argv[1] == "baseline":
    #     method = "Baseline"
    # elif sys.argv[1] == "lca":
    #     method = "Baseline_LCA_avg"
    # elif sys.argv[1] == "baselineDesc":
    #     method = "Baseline_Desc"
    # elif sys.argv[1] == "gontosim":
    #     method = "GOntoSim"
    # elif sys.argv[1] == "GOGO" or sys.argv[1] == "gogo":
    #     method = "GOGO"
    # elif sys.argv[1] == "resnik":
    #     method = "Resnik"
    # elif sys.argv[1] == "lin":
    #     method = "Lin"
    # else:
    #     method = "GOntoSim"
    method_map = {
        "wang": "wang",
        "baseline": "Baseline",
        "lca": "Baseline_LCA_avg",
        "baselineDesc": "Baseline_Desc",
        "gontosim": "GOntoSim",
        "GOGO": "GOGO",
        "gogo": "GOGO",
        "resnik": "Resnik",
        "lin": "Lin",
    }

    method = method_map.get(sys.argv[1].lower(), "GOntoSim")

    if sys.argv[2] == "MF":
        if sys.argv[3] == "IEA":
            # 10890
            # ALL_MF_Classes = [

            # ]
            # ALL_MF_true_Lables = [

            # ]
            # 1366 Unique terms
            unique_goterms_MF = []
            # Class1: 864_MF_I
            class1_mf_IEA = []

            class2_mf_IEA = []
            class3_mf_IEA = []
            class4_mf_IEA = []
            class5_mf_IEA = []
            class6_mf_IEA = []

        if sys.argv[3] == "NONIEA":
            # ALL_MF_Classes = [

            # ]
            unique_goterms_MF = []
            # ALL_MF_true_Lables = [

            # ]
            # 230
            class1_mf = []
            # 1435_MF
            class2_mf = []
            # 1168_MF
            class3_mf = []
            # 171
            class4_mf = []
            # 155
            class5_mf = []
            # 218
            class6_mf = []

        # Classes = ALL_MF_Classes
        if method != "Baseline_Desc":
            if method == "GOntoSim":
                S_values = [
                    (x, Semantic_Value(x, go, "Baseline_LCA_avg"))
                    for x in unique_goterms_MF
                ]
            else:
                print(method)
                S_values = [
                    (x, Semantic_Value(x, go, method)) for x in unique_goterms_MF
                ]

        else:
            S_values = [
                (x, Semantic_Value(x, go, "Baseline")) for x in unique_goterms_MF
            ]

    if sys.argv[2] == "BP":
        if sys.argv[3] == "IEA":
            # 10614
            # ALL_BP_Classes = [

            # ]

            # 4099 Unique terms out 41896
            unique_goterms_BP = []
            # 864_BP_C1_IEA
            class1_bp_IEA = []
            # 4401_BP_C2_IEA
            class2_bp_IEA = []
            # 2887_BP_C3_IEA
            class3_bp_IEA = []
            # 587_BP_C4_IEA
            class4_bp_IEA = []
            # 645_BP_C5_IEA
            class5_bp_IEA = []
            # 1133_BP_C6_IEA
            class6_bp_IEA = []

        if sys.argv[4] == "NONIEA":
            # 3116
            # ALL_BP_Classes = [

            # ]
            # ALL_BP_true_Lables = [

            # ]

            # 3504 out of 10915
            unique_goterms_BP = []
            # 197_BP_C1
            class1_bp = []
            # 1405_BP_C2
            class2_bp = []
            # 978_BP_C3
            class3_bp = []
            # 151_BP_C4
            class4_bp = []
            # 149_BP_C5
            class5_bp = []
            # 212_BP_C6
            class6_bp = []

        # Classes = ALL_BP_Classes
        if method != "Baseline_Desc":
            if method == "GOntoSim":
                S_values = [
                    (x, Semantic_Value(x, go, "Baseline_LCA_avg"))
                    for x in unique_goterms_BP
                ]
            else:
                print(method)
                S_values = [
                    (x, Semantic_Value(x, go, method)) for x in unique_goterms_BP
                ]

        else:
            S_values = [
                (x, Semantic_Value(x, go, "Baseline")) for x in unique_goterms_BP
            ]

    samples = int(sys.argv[4])
    purity = []
    ari = []
    ami = []
    hs = []
    cs = []
    vms = []
    fms = []
    for i in range(1, 11):
        if sys.argv[2] == "MF":
            if sys.argv[3] == "IEA":
                sample_class1 = random.sample(class1_mf_IEA, samples)
                sample_class2 = random.sample(class2_mf_IEA, samples)
                sample_class3 = random.sample(class3_mf_IEA, samples)
                sample_class4 = random.sample(class4_mf_IEA, samples)
                sample_class5 = random.sample(class5_mf_IEA, samples)
                sample_class6 = random.sample(class6_mf_IEA, samples)
            if sys.argv[3] == "NONIEA":
                sample_class1 = random.sample(class1_mf, samples)
                sample_class2 = random.sample(class2_mf, samples)
                sample_class3 = random.sample(class3_mf, samples)
                sample_class4 = random.sample(class4_mf, samples)
                sample_class5 = random.sample(class5_mf, samples)
                sample_class6 = random.sample(class6_mf, samples)
        elif sys.argv[2] == "BP":
            if sys.argv[3] == "IEA":
                sample_class1 = random.sample(class1_bp_IEA, samples)
                sample_class2 = random.sample(class2_bp_IEA, samples)
                sample_class3 = random.sample(class3_bp_IEA, samples)
                sample_class4 = random.sample(class4_bp_IEA, samples)
                sample_class5 = random.sample(class5_bp_IEA, samples)
                sample_class6 = random.sample(class6_bp_IEA, samples)
            if sys.argv[3] == "NONIEA":
                sample_class1 = random.sample(class1_bp, samples)
                sample_class2 = random.sample(class2_bp, samples)
                sample_class3 = random.sample(class3_bp, samples)
                sample_class4 = random.sample(class4_bp, samples)
                sample_class5 = random.sample(class5_bp, samples)
                sample_class6 = random.sample(class6_bp, samples)

        Classes = (
            [x[0] for x in sample_class1]
            + [x[0] for x in sample_class2]
            + [x[0] for x in sample_class3]
            + [x[0] for x in sample_class4]
            + [x[0] for x in sample_class5]
            + [x[0] for x in sample_class6]
        )
        #   labels_true = [x[1] for x in sample_class1] + [x[1] for x in sample_class2] + [x[1] for x in sample_class3] + [x[1] for x in sample_class4] + [x[1] for x in sample_class5] + [x[1] for x in sample_class6]

        labels_true = (
            ([0] * samples)
            + ([1] * samples)
            + ([2] * samples)
            + ([3] * samples)
            + ([4] * samples)
            + ([5] * samples)
        )

        n_clusters = 6
        label = "IEA_MF_" + method + str(samples)
        print(label)
        if method == "Resnik" or method == "Lin":
            if sys.argv[2] == "MF":
                if sys.argv[3] == "IEA":
                    termcounts = semantic_value_resnik_lin_IEA_MF()
                    labels_pred = AgglomerativeClusteringIC(
                        label, Classes, n_clusters, method, termcounts
                    )
                if sys.argv[3] == "NONIEA":
                    termcounts = semantic_value_resnik_lin_nonIEA_MF()
                    labels_pred = AgglomerativeClusteringIC(
                        label, Classes, n_clusters, method, termcounts
                    )

        else:
            _dict = dict(S_values)
            labels_pred = Agglomerative_Clustering(
                label, Classes, n_clusters, method, _dict
            )
        purityscore = purity_score(labels_true, labels_pred)
        ARIScore = ARI_Score(labels_true, labels_pred)
        AMIScore = AMI_Score(labels_true, labels_pred)
        homogeneityscore = metrics.homogeneity_score(labels_true, labels_pred)
        completenessscore = metrics.completeness_score(labels_true, labels_pred)
        v_measurescore = metrics.v_measure_score(labels_true, labels_pred)
        fowlkes_mallowsscore = metrics.fowlkes_mallows_score(labels_true, labels_pred)

        print("purity_score = ", purityscore)
        print("ARI_Score = ", ARIScore)
        print("AMI_Score = ", AMIScore)
        print("homogeneity_score = ", homogeneityscore)
        print("completeness_score = ", completenessscore)
        print("v_measure_score = ", v_measurescore)
        print("fowlkes_mallows_score = ", fowlkes_mallowsscore)

        purity.append(purityscore)
        ari.append(ARIScore)
        ami.append(AMIScore)
        hs.append(homogeneityscore)
        cs.append(completenessscore)
        vms.append(v_measurescore)
        fms.append(fowlkes_mallowsscore)

    print(label)

    print(np.mean(purity))
    print(np.std(purity))

    print(np.mean(ari))
    print(np.std(ari))

    print(np.mean(ami))
    print(np.std(ami))

    print(np.mean(hs))
    print(np.std(hs))

    print(np.mean(cs))
    print(np.std(cs))

    print(np.mean(vms))
    print(np.std(vms))

    print(np.mean(fms))
    print(np.std(fms))


if __name__ == "__main__":
    main()
