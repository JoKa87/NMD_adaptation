import gzip
import os
from progress.bar import IncrementalBar
import requests
import time

from extract_mutations_utils import *


class Data_retriever_utils():
    def __init__(self, params, log):
        self.params = params
        self.log    = log


    # <- function added on 250807 (currently only submitter_id defined)
    def _append_info(self, case_id, file_ids):
        info = {}
        base_url = f"https://api.gdc.cancer.gov/files/"

        # Build the filter for multiple UUIDs
        filters = {
            "op": "in",
            "content": {
                "field": "file_id",
                "value": file_ids
            }
        }

        if "submitter_id" in self.params["appending_targets"]:
            fields = ["file_id,cases.samples.submitter_id,cases.samples.sample_type"]
    
        # Fields to retrieve
        params = {
                  "filters": str(filters).replace("'", '"'),
                  "fields": fields,
                  "format": "JSON",
                  "size": len(file_ids)
                 }

        response = requests.get(base_url, params=params)

        if response.status_code == 200:
            response = response.json()

            if len(response["data"]["hits"]) != len(file_ids):
                print("< exception 1 occurred @_append_info for case_id", case_id)

            else:
                if "submitter_id" in self.params["appending_targets"]:
                    info = {file_id: {"sample_type": [], "submitter_id": []} for file_id in file_ids}

                    for i in range(len(response["data"]["hits"])):
                        if len(response["data"]["hits"][i]["cases"]) != 1:
                            print("< exception 2 occurred @_append_info for case_id", case_id)

                        else:
                            file_id = response["data"]["hits"][i]["file_id"]
                            for j in range(len(response["data"]["hits"][i]["cases"][0]["samples"])):
                                info[file_id]["sample_type"].append(response["data"]["hits"][i]["cases"][0]["samples"][j]["sample_type"])
                                info[file_id]["submitter_id"].append(response["data"]["hits"][i]["cases"][0]["samples"][j]["submitter_id"])

        return info


    # <- function added on 250807 (currently only submitter_id defined)
    def append_info(self, status):
        bar = IncrementalBar(set_bar("appending info"), max=np.sum([len(status["file_ids"][project]) for project in status["file_ids"]]))
        for project in status["file_ids"]:
            for case_id in status["file_ids"][project]:
                # retrieve all file_ids for a given case
                file_ids = []
                for appending_category in self.params["appending_categories"]:
                    file_ids.extend(status["file_ids"][project][case_id][appending_category])

                full_info = self._append_info(case_id, file_ids)

                for appending_category in self.params["appending_categories"]:
                    updated_targets = {appending_target: [] for appending_target in self.params["appending_targets"]}

                    for i, file_id in enumerate(status["file_ids"][project][case_id][appending_category]):
                        #print("i", i, file_id, appending_category, status["file_ids"][project][case_id][appending_category+"_sample_type"], info)
                        info = full_info[file_id]

                        for appending_target in self.params["appending_targets"]:
                            if (len(status["file_ids"][project][case_id][appending_category+"_sample_type"][i]) != len(info["sample_type"])
                                or len(status["file_ids"][project][case_id][appending_category+"_sample_type"][i]) > 2):
                                print("< inconsistent length @append_info for file_id", file_id)
                                print(info)
                                print(status["file_ids"][project][case_id][appending_category+"_sample_type"][i])
                                input("x")

                            # sort appending_target according to sample type in the request file
                            else:
                                if (len(info["sample_type"]) == 1 or (status["file_ids"][project][case_id][appending_category+"_sample_type"][i][0] == info["sample_type"][0]
                                    and status["file_ids"][project][case_id][appending_category+"_sample_type"][i][1] == info["sample_type"][1])):
                                    updated_targets[appending_target].append(info[appending_target])
                                
                                elif (status["file_ids"][project][case_id][appending_category+"_sample_type"][i][0] == info["sample_type"][1]
                                    and status["file_ids"][project][case_id][appending_category+"_sample_type"][i][1] == info["sample_type"][0]):
                                    updated_targets[appending_target].append([info[appending_target][len(info[appending_target])-1-j] for j in range(len(info[appending_target]))])

                                else:
                                    print("< inconsistent sample types @append_info for file_id", file_id)
                                    print(info)
                                    print(status["file_ids"][project][case_id][appending_category+"_sample_type"][i])
                                    input("x")

                    for appending_target in self.params["appending_targets"]:
                        status["file_ids"][project][case_id][appending_category+"_"+appending_target] = updated_targets[appending_target]

                #print(json.dumps(status["file_ids"][project][case_id], indent=4))
                bar.next()
        bar.finish()

        return status


    def case_request(self, start, target):
        response_json = {}

        fields = [
        "case_id",
        "files.data_category",
        "files.experimental_strategy",
        "samples.tumor_code",
        "samples.tumor_code_id",
        "project.primary_site",
        "project.project_id",
        "submitter_id"
        ]

        fields = ",".join(fields)

        cases_endpt = "https://api.gdc.cancer.gov/cases"

        filters = {
            "op": "and",
            "content":[
                {
                "op": "in",
                "content":{
                    "field": "project.project_id",
                    "value": [target]
                    }
                },

                {
                "op": "or",
                "content":[
                    {"op": "in",
                    "content":{
                        "field": "files.data_category",
                        "value": "Copy Number Variation"
                        }
                    },
                    {"op": "in",
                    "content":{
                    "field": "files.experimental_strategy",
                    "value": "RNA-Seq"
                        }
                    },
                    {"op": "in",
                    "content":{
                    "field": "files.experimental_strategy",
                    "value": "WXS"
                        }
                    }
                    ]
                }
            ]
        }

        # A POST is used, so the filter parameters can be passed directly as a Dict object.
        params = {
            "filters": filters,
            "fields":  fields,
            "format":  "JSON",
            "from":    start,
            "size":    "2000"
            }
        
        try:
            response      = requests.post(cases_endpt, headers={"Content-Type": "application/json"}, json=params)
            response_json = response.json()
        
        except:
            print("< error occurred @case_request.")

        return response_json


    def file_request(self, start, target):
        response_json = {}

        fields = [
        "analysis.input_files.experimental_strategy",
        #"cases.case_id",
        "downstream_analyses.output_files.access",
        "downstream_analyses.output_files.data_type",
        "downstream_analyses.output_files.experimental_strategy",
        "downstream_analyses.output_files.file_id",
        "downstream_analyses.output_files.md5sum",
        "downstream_analyses.workflow_type",
        "cases.samples.sample_type",
        #"cases.samples.tumor_code",
        #"cases.samples.tumor_code_id",
        #"cases.project.primary_site",
        "file_id",
        #"submitter_id"
        ]

        fields = ",".join(fields)

        files_endpt = "https://api.gdc.cancer.gov/files"

        filters = {
            "op": "and",
            "content":[
                {
                "op": "in",
                "content":{
                    "field": "cases.case_id",
                    "value": [target]
                    }
                },
            ]
        }

        # A POST is used, so the filter parameters can be passed directly as a Dict object.
        params = {
            "filters": filters,
            "fields":  fields,
            "format":  "JSON",
            "from":    start,
            "size":    "2000"
            }
        
        try:
            response      = requests.post(files_endpt, headers={"Content-Type": "application/json"}, json=params)
            response_json = response.json()
        
        except:
            print("< error occurred @file_request.")
        
        return response_json


    def download_files(self, targets, thread_index, sub_dir, downloads, bar):
        keys = [key for key in targets]
        for i in range(thread_index[0], thread_index[1]+1, 1):
            case = keys[i]

            # check whether files are present for at least RNA and WXS
            if len(targets[case]["RNA"]) > 0 and len(targets[case]["WXS"]) > 0:
                for key in targets[case]:
                    for j in range(len(targets[case][key])):
                        #print(sub_dir+"\\"+key+"\\"+str(targets[case][key][j]), os.path.isfile(sub_dir+"\\"+key+"\\"+str(targets[case][key][j])))
                        if os.path.isfile(sub_dir+"\\"+key+"\\"+str(targets[case][key][j])) == False and key != "CNV_genes" and "sample_type" not in key:
                            print(sub_dir+"\\"+key+"\\"+str(targets[case][key][j]))
                            data_endpt       = "https://api.gdc.cancer.gov/data/{}".format(targets[case][key][j])
                            #os.system("curl " + data_endpt + " -v")

                            response         = requests.get(data_endpt, headers={"Content-Type": "application/json"}).content
                    
                            success          = True

                            if key == "WXS":
                                try:
                                    response = gzip.decompress(response)

                                except:
                                    success = False

                            if success == True:
                                with open(sub_dir+"\\"+key+"\\"+str(targets[case][key][j]), "wb") as outfile:
                                    outfile.write(response)

                            downloads += 1
            else:
                self.log.write("< not all file types found for case " + json.dumps(targets[case]) + "\n")

            #bar.next()

        return downloads
    

    def get_sample_types(self, response, it):
        sample_types = []
        for i in range(len(response["data"]["hits"][it]["cases"])):
            for j in range(len(response["data"]["hits"][it]["cases"][i]["samples"])):
                sample_types.append(response["data"]["hits"][it]["cases"][i]["samples"][j]["sample_type"])

        return sample_types


    def request_defined(self, targets, mode):
        if mode == "case_ids":    output = {target: {"case_id": [], "submitter_id": [], "tumor_code": [], "tumor_code_id": []} for target in targets}
        if mode == "file_ids":    output = {target: {case: {"CNV_genes": [], "CNV_ranges": [], "RNA": [], "WXS": [],
                                                            "CNV_genes_md5sum": [], "CNV_ranges_md5sum": [], "RNA_md5sum": [], "WXS_md5sum": [],
                                                            "CNV_genes_sample_type": [], "CNV_ranges_sample_type": [], "RNA_sample_type": [], "WXS_sample_type": []} 
                                                             for case in targets[target]["case_id"]} for target in targets}
        if mode == "project_ids": output = []

        bar = IncrementalBar("retrieve " + mode.replace("_", ""), max=len(targets))
        cases = 0; files = 0
        for target in targets:
            if mode == "case_ids":
                start     = 0
                no_update = False
                while no_update == False:
                    response = self.case_request(start, target)                 
                    #print(json.dumps(response, indent=4))

                    if "data" in response and "hits" in response["data"]:
                        for i in range(len(response["data"]["hits"])):   
                            if "files" in response["data"]["hits"][i]:
                                cnv_found = False; rna_found = False; wxs_found = False
                                for j in range(len(response["data"]["hits"][i]["files"])):
                                    #print(j, response["data"]["hits"][i]["files"][j]["data_category"])
                                    if ("data_category" in response["data"]["hits"][i]["files"][j] 
                                        and response["data"]["hits"][i]["files"][j]["data_category"] == "Copy Number Variation"):
                                        cnv_found = True
                                    
                                    elif ("experimental_strategy" in response["data"]["hits"][i]["files"][j]
                                        and response["data"]["hits"][i]["files"][j]["experimental_strategy"] == "RNA-Seq"):
                                        rna_found = True

                                    elif ("experimental_strategy" in response["data"]["hits"][i]["files"][j]
                                        and response["data"]["hits"][i]["files"][j]["experimental_strategy"] == "WXS"):
                                        wxs_found = True

                                # changed 241113
                                filter_passed = True
                                if "cnv" in self.params["filter_criterion"] and cnv_found == False: filter_passed = False
                                if "rna" in self.params["filter_criterion"] and rna_found == False: filter_passed = False
                                if "wxs" in self.params["filter_criterion"] and wxs_found == False: filter_passed = False

                                #print("i", i, response["data"]["hits"][i]["submitter_id"], cnv_found, rna_found, wxs_found, len(output[target]["case_id"]), filter_passed)

                                if filter_passed == True and response["data"]["hits"][i]["case_id"] not in output[target]["case_id"]:
                                #if (cnv_found == True and rna_found == True and wxs_found == True
                                #    and response["data"]["hits"][i]["case_id"] not in output[target]["case_id"]):
                                    cases += 1
                                    output[target]["case_id"].append(response["data"]["hits"][i]["case_id"])
                                    output[target]["submitter_id"].append(response["data"]["hits"][i]["submitter_id"])
                                    #for j in range(len(response["data"]["hits"][i]["samples"])):
                                    #    output[target]["tumor_code"].append(response["data"]["hits"][i]["samples"][j]["tumor_code"])                    
                                    #    output[target]["tumor_code_id"].append(response["data"]["hits"][i]["samples"][j]["tumor_code_id"])

                    else:
                        self.log.write("< error. 'data' or 'hits' keywords not found in" + json.dumps(response))
                        print("< error. 'data' or 'hits' keywords not found in", response)

                    #print("start", start, "count", response["data"]["pagination"]["count"], "total", response["data"]["pagination"]["total"], "l", cases)
                    start += response["data"]["pagination"]["count"]
                    if start >= response["data"]["pagination"]["total"] or response["data"]["pagination"]["count"] == 0: no_update = True
                    #input("x")

            elif mode == "file_ids":
                #print(targets[target]["case_id"])
                for case in targets[target]["case_id"]:
                    start          = 0
                    no_update      = False
                    while no_update == False:
                        response = self.file_request(start, case)
                        #print(json.dumps(response, indent=4))
                        print("case", case, len(response["data"]["hits"]))
                        cnv_sample_type = {"CNV_ranges_sample_type": {}, "CNV_genes_sample_type": {}}

                        if "data" in response and "hits" in response["data"]:
                            for i in range(len(response["data"]["hits"])):
                                if "cases" in response["data"]["hits"][i] and "downstream_analyses" in response["data"]["hits"][i]:
                                    for j in range(len(response["data"]["hits"][i]["downstream_analyses"])):
                                        if ("output_files" in response["data"]["hits"][i]["downstream_analyses"][j]
                                            and "workflow_type" in response["data"]["hits"][i]["downstream_analyses"][j]):

                                            for k in range(len(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"])):
                                                if "data_type" in response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]:

                                                    # flexibility required to cover keywords "Copy Number Segment" and "Allele-specific Copy Number Segment"
                                                    # changed on 250313 to specific "Allele-specific Copy Number Segment" as "Copy Number Segment" gives different file type
                                                    # not used in the routine
                                                    #if "Copy Number Segment" in response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["data_type"]:
                                                    if "Allele-specific Copy Number Segment" in response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["data_type"]:
                                                        if response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"] not in output[target][case]["CNV_ranges"]:
                                                            files += 1
                                                            output[target][case]["CNV_ranges"].append(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"])
                                                            output[target][case]["CNV_ranges_md5sum"].append(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["md5sum"])

                                                        # needs to be appended multiple times because sample type information is distributed over two (cancer and normal) entries (different to WXS data)
                                                        if self.get_sample_types(response, i)[0] not in output[target][case]["CNV_ranges_sample_type"]:
                                                            if response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"] in cnv_sample_type["CNV_ranges_sample_type"]:
                                                                cnv_sample_type["CNV_ranges_sample_type"][response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"]].append(self.get_sample_types(response, i)[0])

                                                            else:
                                                                cnv_sample_type["CNV_ranges_sample_type"][response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"]] = [self.get_sample_types(response, i)[0]]
                                                                
                                                    if "Gene Level Copy Number" in response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["data_type"]:
                                                        if response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"] not in output[target][case]["CNV_genes"]:
                                                            files += 1
                                                            output[target][case]["CNV_genes"].append(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"])
                                                            output[target][case]["CNV_genes_md5sum"].append(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["md5sum"])

                                                        # needs to be appended multiple times because sample type information is distributed over two (cancer and normal) entries (different to WXS data)
                                                        if self.get_sample_types(response, i)[0] not in output[target][case]["CNV_genes_sample_type"]:
                                                            if response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"] in cnv_sample_type["CNV_genes_sample_type"]:
                                                                cnv_sample_type["CNV_genes_sample_type"][response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"]].append(self.get_sample_types(response, i)[0])

                                                            else:
                                                                cnv_sample_type["CNV_genes_sample_type"][response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"]] = [self.get_sample_types(response, i)[0]]
                                                                

                                            if response["data"]["hits"][i]["downstream_analyses"][j]["workflow_type"] == "STAR - Counts":
                                                for k in range(len(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"])):
                                                    if "data_type" in response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]:
                                                        if response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["data_type"] == "Gene Expression Quantification":
                                                            if response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"] not in output[target][case]["RNA"]:
                                                                files    += 1
                                                                output[target][case]["RNA"].append(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"])
                                                                output[target][case]["RNA_md5sum"].append(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["md5sum"])
                                                                output[target][case]["RNA_sample_type"].append(self.get_sample_types(response, i))


                                            if response["data"]["hits"][i]["downstream_analyses"][j]["workflow_type"] == "Aliquot Ensemble Somatic Variant Merging and Masking":
                                                for k in range(len(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"])):
                                                    if "data_type" in response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]:
                                                        if response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["data_type"] == "Masked Somatic Mutation":
                                                            if response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"] not in output[target][case]["WXS"]:
                                                                files    += 1
                                                                output[target][case]["WXS"].append(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["file_id"])
                                                                output[target][case]["WXS_md5sum"].append(response["data"]["hits"][i]["downstream_analyses"][j]["output_files"][k]["md5sum"])
                                                                output[target][case]["WXS_sample_type"].append(self.get_sample_types(response, i))


                        else:
                            self.log.write("< error. 'data' or 'hits' keywords not found in" + json.dumps(response))
                            print("< error. 'data' or 'hits' keywords not found in", response)
                        
                        if (len(output[target][case]["CNV_genes"]) == 0 or len(output[target][case]["CNV_ranges"]) == 0) and (len(output[target][case]["RNA"]) > 0 or len(output[target][case]["WXS"]) > 0):
                            self.log.write("< no CNV files found for case id " + case + "\n")
                            self.log.write(json.dumps(response, indent=4))

                        if len(output[target][case]["CNV_genes"]) == 0 and len(output[target][case]["CNV_ranges"]) == 0 and len(output[target][case]["RNA"]) == 0 and len(output[target][case]["WXS"]) == 0:
                            self.log.write("< no files found for case id " + case + "\n")
                            self.log.write(json.dumps(response, indent=4))

                        start += response["data"]["pagination"]["count"]
                        print("start", start, "count", response["data"]["pagination"]["count"], "total", response["data"]["pagination"]["total"])
                        if start >= response["data"]["pagination"]["total"] or response["data"]["pagination"]["count"] == 0: no_update = True
                    
                    for key in cnv_sample_type:
                        for entry in cnv_sample_type[key]:
                            output[target][case][key].append(cnv_sample_type[key][entry])

                    if len(targets[target]) == 0: no_update = True


            if mode == "project_ids":
                response = requests.get(self.params["url"]+"projects/"+target+"?expand=summary,summary.experimental_strategies,summary.data_categories&pretty=true").json()
                matches = 0
                if "data" in response and "summary" in response["data"] and "experimental_strategies" in response["data"]["summary"]:
                    for k in range(len(response["data"]["summary"]["experimental_strategies"])):
                        if (response["data"]["summary"]["experimental_strategies"][k]["experimental_strategy"] == "RNA-Seq"
                            and response["data"]["summary"]["experimental_strategies"][k]["file_count"] > 0):
                            matches += 1

                        if (response["data"]["summary"]["experimental_strategies"][k]["experimental_strategy"] == "WXS"
                            and response["data"]["summary"]["experimental_strategies"][k]["file_count"] > 0):
                            matches += 1
                    #print(target, k, response["data"]["summary"]["experimental_strategies"][k]["experimental_strategy"], matches)
                else:
                    self.log.write("< error. 'data' or 'summary' or 'experimental_strategies' keywords not found in" + str(response))
                    print("< error. 'data' or 'summary' or 'experimental_strategies' keywords not found in", response)

                if matches == 2:
                    output.append(target)
                    
            bar.next()

        bar.finish()
        if mode == "case_ids": print("<", cases, "detected.")
        if mode == "file_ids": print("<", files, "detected.")
        return output


    def request_undefined(self, mode):
        output = []
        start          = 0

        failed_request = False
        no_update      = False
        while failed_request == False and no_update == False:
            try:
                if mode == "projects":
                    response = requests.get(self.params["url"]+"projects?from="+ str(start) + "&pretty=true").json()

                    for k in range(len(response["data"]["hits"])):
                        #print(k, response["data"]["hits"][k]["id"])
                        if (self.params["project_filter"] == None or (self.params["project_filter"] in response["data"]["hits"][k]["id"]
                        and (len(self.params["projects"]) == 0 or (len(self.params["projects"]) > 0 and response["data"]["hits"][k]["id"] in self.params["projects"])))): # <- modified on 251025 
                            output.append(response["data"]["hits"][k]["id"])

                start += response["data"]["pagination"]["count"]
                if len(response["data"]["hits"]) == 0: no_update = True

            except:
                failed_request = True

        return output