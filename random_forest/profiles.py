import json


base_profiles = {
                "datatype_msk":         {
                                         "id_filter"                    : "MSK",
                                         "mutation_stats_fname"         : "msk_mutation_stats.json",
                                         "variant_filter_fnames"        : ["msk_chord_variants.txt"],
                                         "variant_fnames"               : ["msk_chord_variants.txt"]
                                        },
                "datatype_tcga":        {
                                         "id_filter"                    : "TCGA",
                                         "mutation_stats_fname"         : "tcga_mutation_stats.json",
                                         "variant_filter_fnames"        : ["tcga_variants.txt"],
                                         "variant_fnames"               : ["tcga_variants.txt"]
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
                                        },
            "msk_colorectal_high_ptc":  {
                                         **{key: base_profiles["datatype_msk"][key] for key in base_profiles["datatype_msk"]},
                                         **{key: base_profiles["mutation_msk_colorectal_high_ptc"][key] for key in base_profiles["mutation_msk_colorectal_high_ptc"]},
                                         **{key: base_profiles["projects_off_msk"][key] for key in base_profiles["projects_off_msk"]},
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
            "tcga_mo_frameshift":       {
                                         **{key: base_profiles["datatype_tcga"][key] for key in base_profiles["datatype_tcga"]},
                                         **{key: base_profiles["mutation_tcga_frameshift"][key] for key in base_profiles["mutation_tcga_frameshift"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_off"][key] for key in base_profiles["weights_off"]},
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
                                         "variant_fnames": ["tcga_oncogenes.txt"],
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
                                            "mutation_stats_fname"         : "tcga_mutation_stats.json",
                                            "variant_filter_fnames"        : ["tcga_variants.txt"],
                                            "variant_fnames"               : ["tcga_combined.txt"]
                                           },
                                         **{key: base_profiles["mutation_tcga_full"][key] for key in base_profiles["mutation_tcga_full"]},
                                         **{key: base_profiles["projects_off_tcga"][key] for key in base_profiles["projects_off_tcga"]},
                                         **{key: base_profiles["weights_on"][key] for key in base_profiles["weights_on"]},
                                        },
           }


def apply_profile(params):
    for key in params:
        if key in profiles[params["profile"]]:
            params[key] = profiles[params["profile"]][key]
    
    print(json.dumps(profiles[params["profile"]], indent=4))
    return params