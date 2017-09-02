import sys
import string
import json
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
import sys
figpath = "dist.png"

sns.set_style('whitegrid')
sns.set_palette('Set1', 13)

traces = []

for f in sys.argv[1]:
    xs, ys = [], []
    v50 = 0
    sample = f.split(".")[0]
    found = False
    for x, y in (x.rstrip().split("\t") for x in open(f)):
        y = float(y)
        if y < 0.01: continue
        if not found and y > 0.5:
            v50 = x
            found = True
            print f, x, y

        xs.append(float(x))
        ys.append(y)

    traces.append({
           'x': list(np.round(xs, 3)),
           'y': list(np.round(ys, 3)),
           'mode': 'lines',
           'name': sample + (" (%.1f)" % float(v50))
    })
    plt.plot(xs, ys, label=sample)


plt.xlabel("Coverage")
plt.ylabel("Proportion of bases at coverage")
plt.xlim(xmin=0)
plt.ylim(ymin=0, ymax=1)
plt.savefig(figpath)

tmpl = """<html>
<head>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body><div id="plot-div"></div>
<script>
var layout = {
    hovermode: 'closest',
    width: 900,
    height: 900,
    xaxis: {title: 'Coverage', domain: [0, 1]},
    yaxis: {title: 'Proportion of bases at coverage', domain: [0, 1]},
    legend: {
        x: 0.1,
        y: 0.1
    },
}
Plotly.newPlot('plot-div', $data, layout);
</script>
</body>
</html>"""

with open("dist.html", "w") as html:
    tmpl = string.Template(tmpl)
    html.write(tmpl.substitute(data=json.dumps(traces)))
