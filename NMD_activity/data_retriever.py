import os
from progress.bar import IncrementalBar
import threading

from data_retriever_utils import *
from extract_mutations_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


# parameters
params = {
         "appending_categories"     : ["CNV_genes", "CNV_ranges", "RNA", "WXS"], # <- added on 250807, applies to append_info and allows selection of categories to be appended
         "appending_targets"        : ["submitter_id"], # <- added on 250807, applies to append_info and allows selection of targets to be append to specified categories
         "data_dir"                 : parent_dir+r"\data",
         "filter_criterion"         : ["wxs"],
         "max_threads"              : 1,
         "mode"                     : ["create_status"], # append_info, create_status, download_files # <- append_info added on 250807
         "os_sep"                   : "\\",
         "project_filter"           : "CPTAC-3", # "TCGA", # None means no filter <- modified on 2501020
         "projects"                 : [],
         "status_fname"             : "filtered_status_w_md5sum.json",
         "url"                      : "https://api.gdc.cancer.gov/"
         }


def main():
    log = open(params["data_dir"]+params["os_sep"]+"data_retriever.log", "w")
    dru = Data_retriever_utils(params, log)

    # <- if-statement added on 250708
    if "append_info" in params["mode"]:
        with open(params["data_dir"]+params["os_sep"]+params["status_fname"], "r") as file:
            status = json.load(file)

        status = dru.append_info(status)

        with open(params["data_dir"]+params["os_sep"]+params["status_fname"].split(".")[0]+"_appended.json", "w") as file:
            file.write(status)


    if "create_status" in params["mode"]:
        # initialize file dictionary
        status = {
                 "case_ids":             [],
                 "file_ids":             [],
                 "filtered_project_ids": [],
                 "project_ids":          [],
                 }

        # load status file (if present)
        status = load_status(params, status)
        print("< status loaded.")
        status_update = False

        # find projects
        if len(status["project_ids"]) == 0:          status["project_ids"]          = dru.request_undefined("projects"); status_update = True
        if status_update == True:                    save_status(status, params)

        # find metadata for all project
        if len(status["filtered_project_ids"]) == 0: status["filtered_project_ids"] = dru.request_defined(status["project_ids"], "project_ids"); status_update = True
        if status_update == True:                    save_status(status, params)

        # find case ids for selected project ids
        if len(status["case_ids"]) == 0:             status["case_ids"]             = dru.request_defined(status["filtered_project_ids"], "case_ids"); status_update = True
        if status_update == True:                    save_status(status, params)

        # find file ids for selected case ids
        if len(status["file_ids"]) == 0:             status["file_ids"]             = dru.request_defined(status["case_ids"], "file_ids"); status_update = True
        if status_update == True:                    save_status(status, params)


    # download files project-wise
    if "download_files" in params["mode"]:
        status = load_status(params, None, "250312_status.json")
        print("< status loaded.")

        if len(status["file_ids"]) > 0:
            for project in status["file_ids"]:
                if len(status["file_ids"][project]) > 0:
                    print(project, len(status["file_ids"][project]))
                    sub_dir = params["data_dir"] + "\\" + project
                    if not os.path.isdir(sub_dir): os.mkdir(sub_dir)
                    
                    for case in status["file_ids"][project]:
                        for key in status["file_ids"][project][case]:
                            if os.path.isdir(sub_dir+"\\"+key) == False: os.mkdir(sub_dir+"\\"+key)

                    threads = []
                    thread_index = split_index(len(status["file_ids"][project]), params["max_threads"])
                    bar          = IncrementalBar("download of "+project+" files.", max=len(status["file_ids"][project]))
                    downloads    = 0
                    
                    if len(thread_index) > 1:
                        for i in range(len(thread_index)):
                            thread = threading.Thread(target=dru.download_files, args=(status["file_ids"][project], thread_index[i], sub_dir, downloads, bar))
                            threads.append(thread)
                            thread.start()

                        for i, thread in enumerate(threads):
                            thread.join()

                    if len(thread_index) == 1:
                        dru.download_files(status["file_ids"][project], thread_index[0], sub_dir, downloads, bar)

                    bar.finish()
                    print("<", downloads, "for", len(status["file_ids"][project]), "cases.")
    
    log.close()


if __name__ == "__main__":
    main()