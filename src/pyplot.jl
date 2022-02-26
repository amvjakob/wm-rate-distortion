using PyPlot, LaTeXStrings
using PyCall
@pyimport seaborn as sns

# define custom color palette
colors = [
    "#1F77B4",
    "#FF7F0E",
    "#2CA02C",
    "#D62728",
    "#9467BD",
    "#8C564B",
    "#E377C2",
    "#7F7F7F",
    "#BCBD22",
    "#17BECF",
]

customRcParams = Dict(
    "text.usetex" => false,
    # font
    "font.family" => "sans-serif",
    "font.sans-serif" => ["Arial"],
    "font.size" => "10",
    "font.weight" => "normal",
    # lines
    "lines.linewidth" => 1.5,
    # axes
    "axes.prop_cycle" => plt.cycler(color = colors),
    "axes.labelsize" => 12,
    "axes.titlelocation" => "left",
    "axes.titlesize" => 14,
    "axes.titleweight" => "bold",
    "axes.spines.top" => false,
    "axes.spines.right" => false,
    # error bars
    "errorbar.capsize" => 3,
    # legend
    "legend.handletextpad" => 0.2,
    "legend.frameon" => false,
)

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")

# apply custom rc params
merge!(rcParams, customRcParams)

# output svg in IJulia
PyPlot.svg(true)

"Get color cycle."
get_palette() = rcParams["axes.prop_cycle"].by_key()["color"]

"Get color `i` in color cycle."
get_palette(index) = rcParams["axes.prop_cycle"].by_key()["color"][index]

"Shorten legend lines and make them thicker."
function legend_shorten_lines!(legend)
    for (i, line) in enumerate(legend.get_lines())
        line.set_xdata([5, 15])
        line.set_linewidth(2)
    end
    return legend
end

"Set legend text color."
function legend_set_text_color!(legend, colors)
    for (i, text) in enumerate(legend.get_texts())
        text.set(color = colors[i]) # usetex=true
    end
    return legend
end

"Add 'Data' and 'Simulation' annotations on side of plot."
function annotate_data_simulation_side!(
    text1 = "Data",
    text2 = "Simulation";
    y1 = 0.77,
    y2 = 0.27,
)
    annotate(
        text1,
        (0.01, y1);
        fontsize = 16,
        fontweight = "bold",
        rotation = 90,
        xycoords = "figure fraction",
        verticalalignment = "center",
    )
    annotate(
        text2,
        (0.01, y2);
        fontsize = 16,
        fontweight = "bold",
        rotation = 90,
        xycoords = "figure fraction",
        verticalalignment = "center",
    )
end

"Add 'Data' and 'Simulation' annotations on top of plot."
function annotate_data_simulation_top!(
    text1 = "Data",
    text2 = "Simulation";
    x1 = 0.26,
    x2 = 0.76,
)
    annotate(
        text1,
        (x1, 0.99);
        fontsize = 16,
        fontweight = "bold",
        xycoords = "figure fraction",
        horizontalalignment = "center",
    )
    annotate(
        text2,
        (x2, 0.99);
        fontsize = 16,
        fontweight = "bold",
        xycoords = "figure fraction",
        horizontalalignment = "center",
    )
end


;