import sys
import string
import json
import itertools as it
from operator import itemgetter
import collections
import sys
figpath = "dist.png"

traces = collections.defaultdict(list)
chroms = collections.OrderedDict()

for f in sys.argv[1:]:
    sample = f.replace(".mosdepth.dist.txt", "")
    gen = (x.rstrip().split("\t") for x in open(f))
    for chrom, data in it.groupby(gen, itemgetter(0)):
        chroms[chrom] = True
        xs, ys = [], []
        v50 = 0
        found = False
        for _, x, y in data:
            y = float(y)
            if y < 0.01: continue
            if not found and y > 0.5:
                v50 = x
                found = True
                print "%s\t%s\t%s\t%.3f" % (sample, chrom, x, y)

            xs.append(float(x))
            ys.append(y)

        traces[chrom].append({
               'x': [round(x, 3) for x in xs],
               'y': [round(y, 3) for y in ys],
               'mode': 'lines',
               'name': sample + (" (%.1f)" % float(v50))
        })

tmpl = """<html>
<head>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>$plot_divs</div>
<script>
var layout = {
    hovermode: 'closest',
    width: 400,
    height: 400,
    xaxis: {title: 'Coverage'},
    yaxis: {title: 'Proportion of bases at coverage', domain: [0, 1], dtick: 0.25},
    showlegend: $showlegend,
    legend: {
        x: 0.1,
        y: 0.1
    },
}
"""
footer = """
</script>
</body>
</html>"""

chr_tmpl = """
Plotly.newPlot('plot-div-$chrom', $data, layout);
"""


tmpl = string.Template(tmpl)
with open("dist.html", "w") as html:
    divs = "\n".join("<h4>{chrom}</h4><div id='plot-div-{chrom}'></div>".format(chrom=c) for c in chroms)
    html.write(tmpl.substitute(showlegend="true" if len(sys.argv[1:]) < 20 else "false", plot_divs=divs))
    for chrom in chroms:
        html.write(string.Template(chr_tmpl).substitute(chrom=chrom, data=json.dumps(traces[chrom])))
    html.write(footer)
