import json


base_profiles = {
                "datatype_control":     {
                                         "id_filter"                    : "control",
                                         "mutation_stats_fname"         : "mutation_stats.json",
                                         "variant_filter_fnames"        : ["control_filtered_appended.txt"],
                                         "variant_fnames"               : ["control_filtered_appended.txt"]
                                        },
                "datatype_msk":         {
                                         "id_filter"                    : "MSK",
                                         "mutation_stats_fname"         : "mutation_stats_msk.json",
                                         "variant_filter_fnames"        : ["msk_chord_appended.txt"],
                                         "variant_fnames"               : ["msk_chord_appended.txt"]
                                        },
                "datatype_tcga":        {
                                         "id_filter"                    : "TCGA",
                                         "mutation_stats_fname"         : "mutation_stats.json",
                                         "variant_filter_fnames"        : ["tcga.txt"],
                                         "variant_fnames"               : ["tcga.txt"]
                                        },
                "mutation_msk_full":    {
                                         "masks"                        : ["frameshift", "nonsense"],
                                         "variant_filter_value_filter"  : {"FEATURE:prediction": 0, "ID:multiple indels": [False]},
                                         "variant_value_filter"         : {"FEATURE:prediction": 0, "ID:multiple indels": [False]},
                                        },
                "mutation_msk_colorectal_low_ptc":{
                                         "masks"                        : ["frameshift", "nonsense"],
                                         "use_variant_filter"           : True,
                                         "variant_filter_value_filter"  : {"FEATURE:prediction": 0, "FEATURE:prediction": 20, "ID:multiple indels": [False], "ID:cancer type": ["Colorectal Cancer"]},
                                         "variant_value_filter"         : {"FEATURE:prediction": 0, "FEATURE:prediction_LESS": 20, "ID:multiple indels": [False], "ID:cancer type": ["Colorectal Cancer"]},
                                        },
                "mutation_msk_colorectal_high_ptc":{
                                         "masks"                        : ["frameshift", "nonsense"],
                                         "use_variant_filter"           : True,
                                         "variant_filter_value_filter"  : {"FEATURE:prediction": 0, "FEATURE:prediction_LESS": 20, "ID:multiple indels": [False], "ID:cancer type": ["Colorectal Cancer"]},
                                         "variant_value_filter"         : {"FEATURE:prediction": 0, "FEATURE:prediction": 20, "ID:multiple indels": [False], "ID:cancer type": ["Colorectal Cancer"]},
                                        },
                "mutation_msk_nonsense": {
                                         "masks"                        : ["nonsense"],
                                         "variant_filter_value_filter"  : {"FEATURE:prediction": 0, "ID:multiple indels": [False], "ID:variant classification": ["Nonsense_Mutation"]},
                                         "variant_value_filter"         : {"FEATURE:prediction": 0, "ID:multiple indels": [False], "ID:variant classification": ["Nonsense_Mutation"]},
                                        },
                "mutation_msk_frameshift": {
                                         "masks"                        : ["nonsense"],
                                         "variant_filter_value_filter"  : {"FEATURE:prediction": 0, "ID:multiple indels": [False], "ID:variant classification": ["Frame_Shift_Del", "Frame_Shift_Ins"]},
                                         "variant_value_filter"         : {"FEATURE:prediction": 0, "ID:multiple indels": [False], "ID:variant classification": ["Frame_Shift_Del", "Frame_Shift_Ins"]},
                                        },
                "mutation_tcga_full":   {
                                         "masks"                        : ["frameshift", "nonsense"],
                                         "variant_filter_value_filter"  : {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0], "FEATURE:prediction": 0, "ID:multiple indels": [False]},
                                         "variant_value_filter"         : {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0], "FEATURE:prediction": 0, "ID:multiple indels": [False]},
                                        },
                "mutation_tcga_nonsense": {
                                         "masks"                        : ["nonsense"],
                                         "variant_filter_value_filter"  : {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0], "FEATURE:prediction": 0, "ID:multiple indels": [False], "ID:variant classification": ["Nonsense_Mutation"]},
                                         "variant_value_filter"         : {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0], "FEATURE:prediction": 0, "ID:multiple indels": [False], "ID:variant classification": ["Nonsense_Mutation"]},
                                        },
                "mutation_tcga_frameshift": {
                                         "masks"                        : ["frameshift"],
                                         "variant_filter_value_filter"  : {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0], "FEATURE:prediction": 0, "ID:multiple indels": [False], "ID:variant classification": ["Frame_Shift_Del", "Frame_Shift_Ins"]},
                                         "variant_value_filter"         : {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0], "FEATURE:prediction": 0, "ID:multiple indels": [False], "ID:variant classification": ["Frame_Shift_Del", "Frame_Shift_Ins"]},
                                        },
                "patients_on_msk":      {
                                         "block"                        : "ID:patient id",
                                         "mutation_stats_gene_target"   : None,
                                         "mutation_stats_pair_target"   : "total",
                                         "variants_only"                : False,
                                        },
                "projects_on_msk":      {
                                         "block"                        : "ID:cancer type",
                                         "mutation_stats_gene_target"   : None,
                                         "mutation_stats_pair_target"   : "total",
                                         "variants_only"                : False,
                                        },
                "projects_off_msk":     {
                                         "block"                        : None,
                                         "mutation_stats_gene_target"   : "total",
                                         "mutation_stats_pair_target"   : "total",
                                         "variants_only"                : False,
                                        },
                "patients_on_tcga":     {
                                         "block"                        : "ID:case id",
                                         "mutation_stats_gene_target"   : None,
                                         "mutation_stats_pair_target"   : None,
                                         "variants_only"                : True,
                                        },
                "projects_on_control":  {
                                         "block"                        : "ID:sample id",
                                         "mutation_stats_gene_target"   : "total",
                                         "mutation_stats_pair_target"   : "total",
                                         "variants_only"                : True,
                                        },
                "projects_on_tcga":     {
                                         "block"                        : "ID:project",
                                         "mutation_stats_gene_target"   : None,
                                         "mutation_stats_pair_target"   : None,
                                         "variants_only"                : True,
                                        },
                "projects_off_tcga":    {
                                         "block"                        : None,
                                         "mutation_stats_gene_target"   : "total",
                                         "mutation_stats_pair_target"   : "total",
                                         "variants_only"                : False,
                                        },
                "weights_on":           {
                                         "apply_3mers"                  : True,
                                         "apply_mutation_stats"         : True,
                                         "average_indel_probabilities"  : True,
                                        },
                "weights_off":          {
                                         "apply_3mers"                  : False,
                                         "apply_mutation_stats"         : False,
                                         "average_indel_probabilities"  : False,
                                        },
                }


profiles = {
            "control_projectwise":      {
                                         **{key: base_profiles["datatype_control"][key] for key in base_profiles["datatype_control"]},
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_on_control"][key] for key in base_profiles["projects_on_control"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                         "variant_filter_value_filter": {"ID:base mean": 100},
                                         "variant_value_filter": {"ID:base mean": 100}
                                        },
            "msk":                      {
                                         **{key: base_profiles["datatype_msk"][key] for key in base_profiles["datatype_msk"]},
                                         **{key: base_profiles["mutation_msk_full"][key] for key in base_profiles["mutation_msk_full"]},
                                         **{key: base_profiles["projects_off_msk"][key] for key in base_profiles["projects_off_msk"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "msk_mo":                   {
                                         **{key: base_profiles["datatype_msk"][key] for key in base_profiles["datatype_msk"]},
                                         **{key: base_profiles["mutation_msk_full"][key] for key in base_profiles["mutation_msk_full"]},
                                         **{key: base_profiles["projects_off_msk"][key] for key in base_profiles["projects_off_msk"]},
                                         **{key: base_profiles["weights_off"][key] for key in base_profiles["weights_off"]},
                                        },
            "msk_colorectal_low_ptc":   {
                                         **{key: base_profiles["datatype_msk"][key] for key in base_profiles["datatype_msk"]},
                                         **{key: base_profiles["mutation_msk_colorectal_low_ptc"][key] for key in base_profiles["mutation_msk_colorectal_low_ptc"]},
                                         **{key: base_profiles["projects_off_msk"][key] for key in base_profiles["projects_off_msk"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                         "mutation_stats_gene_target": "Colorectal Cancer", # <- added on 251022, must be re-calculated!
                                        },
            "msk_colorectal_high_ptc":  {
                                         **{key: base_profiles["datatype_msk"][key] for key in base_profiles["datatype_msk"]},
                                         **{key: base_profiles["mutation_msk_colorectal_high_ptc"][key] for key in base_profiles["mutation_msk_colorectal_high_ptc"]},
                                         **{key: base_profiles["projects_off_msk"][key] for key in base_profiles["projects_off_msk"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                         "mutation_stats_gene_target": "Colorectal Cancer", # <- added on 251022, must be re-calculated!
                                        },
            "msk_patientwise":          {
                                         **{key: base_profiles["datatype_msk"][key] for key in base_profiles["datatype_msk"]},
                                         **{key: base_profiles["mutation_msk_full"][key] for key in base_profiles["mutation_msk_full"]},
                                         **{key: base_profiles["patients_on_msk"][key] for key in base_profiles["patients_on_msk"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "msk_projectwise":          {
                                         **{key: base_profiles["datatype_msk"][key] for key in base_profiles["datatype_msk"]},
                                         **{key: base_profiles["mutation_msk_full"][key] for key in base_profiles["mutation_msk_full"]},
                                         **{key: base_profiles["projects_on_msk"][key] for key in base_profiles["projects_on_msk"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "msk_mo_projectwise":       {
                                         **{key: base_profiles["datatype_msk"][key] for key in base_profiles["datatype_msk"]},
                                         **{key: base_profiles["mutation_msk_full"][key] for key in base_profiles["mutation_msk_full"]},
                                         **{key: base_profiles["projects_on_msk"][key] for key in base_profiles["projects_on_msk"]},
                                         **{key: base_profiles["weights_off"][key] for key in base_profiles["weights_off"]},
                                        },
            "tcga":                     {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "tcga_mo":                  {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_off"][key] for key in base_profiles["weights_off"]},
                                        },
            "tcga_nonsense":            {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_nonsense"][key] for key in base_profiles["mutation_tcga_nonsense"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "tcga_nonsense_oncogenes":  {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_nonsense"][key] for key in base_profiles["mutation_tcga_nonsense"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                         "variant_fnames": ["tcga_oncogenes.txt"]
                                        },
            "tcga_nonsense_tsgs":       {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_nonsense"][key] for key in base_profiles["mutation_tcga_nonsense"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                         "variant_fnames": ["tcga_tsgs.txt"]
                                        },
            "tcga_mo_nonsense":         {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_nonsense"][key] for key in base_profiles["mutation_tcga_nonsense"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_off"][key] for key in base_profiles["weights_off"]},
                                        },
            "tcga_frameshift":          {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_frameshift"][key] for key in base_profiles["mutation_tcga_frameshift"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "tcga_frameshift_oncogenes":{
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_frameshift"][key] for key in base_profiles["mutation_tcga_frameshift"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                         "variant_fnames": ["tcga_oncogenes.txt"]
                                        },
            "tcga_frameshift_tsgs":     {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_frameshift"][key] for key in base_profiles["mutation_tcga_frameshift"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                         "variant_fnames": ["tcga_tsgs.txt"]
                                        },
            "tcga_mo_frameshift":       {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_frameshift"][key] for key in base_profiles["mutation_tcga_frameshift"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_off"][key] for key in base_profiles["weights_off"]},
                                        },
            "tcga_patientwise":         {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["patients_on_tcga"][key] for key in base_profiles["patients_on_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "tcga_patientwise_test":    {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                         "mutation_stats_gene_target": "TCGA-UCEC",
                                         "mutation_stats_pair_target": "TCGA-UCEC",
                                         "variant_value_filter": {"ID:ptc reads": 1, "ID:cnv total_EXCLUDED": [0], "FEATURE:prediction": 0, "FEATURE:ptc_mutations": 200, "ID:multiple indels": [False], "ID:project": ["TCGA-UCEC"]},
                                        },
            "tcga_projectwise":         {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_on_tcga"][key] for key in base_profiles["projects_on_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "tcga_mo_projectwise":      {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_on_tcga"][key] for key in base_profiles["projects_on_tcga"]},
                                         **{key: base_profiles["weights_off"][key] for key in base_profiles["weights_off"]},
                                        },
            "tcga_nonsense_projectwise":{
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_nonsense"][key] for key in base_profiles["mutation_tcga_nonsense"]},
                                         **{key: base_profiles["projects_on_tcga"][key] for key in base_profiles["projects_on_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "tcga_mo_nonsense_projectwise":{
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_nonsense"][key] for key in base_profiles["mutation_tcga_nonsense"]},
                                         **{key: base_profiles["projects_on_tcga"][key] for key in base_profiles["projects_on_tcga"]},
                                         **{key: base_profiles["weights_off"][key] for key in base_profiles["weights_off"]},
                                        },
            "tcga_frameshift_projectwise":{
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_frameshift"][key] for key in base_profiles["mutation_tcga_frameshift"]},
                                         **{key: base_profiles["projects_on_tcga"][key] for key in base_profiles["projects_on_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "tcga_mo_frameshift_projectwise":{
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_frameshift"][key] for key in base_profiles["mutation_tcga_frameshift"]},
                                         **{key: base_profiles["projects_on_tcga"][key] for key in base_profiles["projects_on_tcga"]},
                                         **{key: base_profiles["weights_off"][key] for key in base_profiles["weights_off"]},
                                        },
            "tcga_oncogenes":           {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                         "variant_fnames": ["tcga_oncogenes.txt"]
                                        },
            "tcga_tsgs":                {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                         "variant_fnames": ["tcga_tsgs.txt"],
                                        },
            "tcga_test":                {
                                         **{
                                            "id_filter"                    : "TCGA",
                                            "mutation_stats_fname"         : "mutation_stats.json",
                                            "variant_filter_fnames"        : ["tcga.txt"],
                                            "variant_fnames"               : ["tcga_combined.txt"]
                                           },
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
            "tcga_mo":                  {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_off"][key] for key in base_profiles["weights_off"]},
                                        },
           }


def apply_profile(params):
    profile = params["profile"]

    # for CPTAC3, create profile according to TCGA and finally replace target file
    if "cptac3" in params["profile"]:
        profile = params["profile"]
        profile = profile.replace("cptac3", "tcga")

    for key in params:
        if key in profiles[profile]:
            params[key] = profiles[profile][key]
    
    if "cptac3" in params["profile"]:
        if "projectwise" in params["profile"]:
            profiles[profile]["mutation_stats_gene_target"] = params["mutation_stats_gene_target"] = "total"
            profiles[profile]["mutation_stats_pair_target"] = params["mutation_stats_pair_target"] = "total"


        profiles[profile]["mutation_stats_fname"]  = params["mutation_stats_fname"]  = "mutation_stats_cptac3.json"
        profiles[profile]["variant_filter_fnames"] = params["variant_filter_fnames"] = ["cptac3.txt"]
        profiles[profile]["variant_fnames"]        = params["variant_fnames"]        = ["cptac3.txt"]
                  
    print(json.dumps(profiles[profile], indent=4))
    return params