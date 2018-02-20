import sys
import string
import json
import itertools as it
from operator import itemgetter
import collections
import sys

traces = collections.defaultdict(list)
chroms = collections.OrderedDict()
chroms["total"] = True

for f in sys.argv[1:]:
    sample = f.replace(".mosdepth.dist.txt", "")
    gen = (x.rstrip().split("\t") for x in open(f))
    for chrom, data in it.groupby(gen, itemgetter(0)):
        if chrom.startswith("GL"): continue
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
                print("{}\t{}\t{}\t{:.3f}".format(sample, chrom, x, y))

            xs.append(float(x))
            ys.append(y)

        if len(xs) > 100:
            xs = [x for i, x in enumerate(xs) if ys[i] > 0.02]
            ys = [y for y in ys if y > 0.02]
            if len(xs) > 100:
                xs = xs[::2]
                ys = ys[::2]

        traces[chrom].append({
               'x': [round(x, 3) for x in xs],
               'y': [round(y, 3) for y in ys],
               'mode': 'lines',
               'name': sample + (" (%.1f)" % float(v50))
        })

tmpl = """<html>
<head>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <style>
div {
  width: 800px;
  height: 420px;
}
</style
</head>
<body>$plot_divs</div>
<script>
var layout = {
    hovermode: 'closest',
    xaxis: {title: 'Coverage'},
    yaxis: {title: 'Proportion of bases at coverage', domain: [0, 1], dtick: 0.25},
    showlegend: $showlegend,
    autosize: true,
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
Plotly.newPlot('plot-div-$chrom', $data, layout, {displayModeBar: false, displaylogo: false, fillFrame: false, autosizeable: true});
"""

tmpl = string.Template(tmpl)
with open("dist.html", "w") as html:
    divs = "\n".join("<{div}>{chrom}</{div}><div id='plot-div-{chrom}'></div><hr/>".format(
        chrom=c, div="h2" if c == "total" else "b") for c in chroms)
    html.write(tmpl.substitute(showlegend="true" if len(sys.argv[1:]) < 20 else "false", plot_divs=divs))
    for chrom in chroms:
        html.write(string.Template(chr_tmpl).substitute(chrom=chrom, data=json.dumps(traces[chrom])))
    html.write(footer)
