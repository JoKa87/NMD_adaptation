import os

from plot_load import *
from plot_utils import *

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def main():
    # params #
    data = [
            {
            "data"          : [],
            "datatype"      : ["pandas"],
            "extensions"    : [None],
            "features"      : {
                              "x"           : "average",
                              "xcol"        : [[None]],
                              "ylabel"      : "feature importance",# "mean decrease in impurity",
                              "z"           : "sem",
                              "zcutoff"     : (0, 0.01),
                              "zlabel"      : "sem"
                              },
            "paths"         : [parent_dir+r"\data\nmdelphi\forest_importances.txt"], # insert model folder name here
            "separators"    : [","],
            "type"          : "random forest coefficients"
            }
        ]
    

    dims       = [3, 4]
    resolution = 600
    run_dir    = parent_dir+r"\data"

    # loading section
    pu = Plot_utils()

    data = [data[i] for i in range(len(data)) if "off" not in data[i] or data[i]["off"] == False]
    data = load(data)


    # define figure
    fig = plt.figure(figsize=(180/25.4, 180/25.4))
    plt.subplots_adjust(hspace=0.25)

    gs = fig.add_gridspec(dims[0], 3*dims[1])

    subplots = []
    subplots.append(fig.add_subplot(gs[0:2, 0:6]))

    for i in range(len(data)):
        if data[i]["type"] == "random forest coefficients":
            data[i]["data"][0].index = data[i]["data"][0][data[i]["data"][0].columns[0]]
            subplots[i]              = pu.plot_coefficients(subplots[i], data[i]["data"][0], data[i]["features"])
        
        subplots[i].text(-0.1, 1.1, string.ascii_lowercase[i], transform=subplots[i].transAxes, size=9, weight='bold')

    plt.show()
    fig.savefig(run_dir+"\\figures\\FigS1.svg", dpi=resolution, transparent=True)
    return


if __name__ == '__main__':
    main()