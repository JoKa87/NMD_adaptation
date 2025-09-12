import pandas as pd

from bravado.client import SwaggerClient


data_dir = r"C:\Programming\Translational_genomics\NMD_analysis\data"


cbioportal = SwaggerClient.from_url('https://www.cbioportal.org/api/v2/api-docs',
                                    config={"validate_requests":False,"validate_responses":False,"validate_swagger_spec": False})

for a in dir(cbioportal):
    cbioportal.__setattr__(a.replace(' ', '_').lower(), cbioportal.__getattr__(a))


muts = cbioportal.mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
    molecularProfileId="msk_impact_2017_mutations", # {study_id}_mutations gives default mutations profile for study 
    sampleListId="msk_impact_2017_all", # {study_id}_all includes all samples
    projection="DETAILED" # include gene info
).result()

mutations = {}
for i in range(len(muts)):
    if len(mutations) > 0:
        for key in muts[i]:
            mutations[key].append(muts[i][key])

    if len(mutations) == 0:
        for key in muts[i]:
            mutations[key] = [muts[i][key]]

mutations = pd.DataFrame(mutations)
mutations.to_csv(data_dir+"\\"+"msk_impact_2017_mutations.txt", sep=",")