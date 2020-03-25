import sys
import string
import json
import itertools as it
from operator import itemgetter
import collections
from argparse import ArgumentParser


def main():
    args = get_args()
    traces = collections.defaultdict(list)
    chroms = collections.OrderedDict()
    chroms["total"] = True

    for f in args.input:
        sample = f.replace(".mosdepth.global.dist.txt", "")
        gen = (x.rstrip().split("\t") for x in open(f))
        for chrom, data in it.groupby(gen, itemgetter(0)):
            if chrom.startswith("GL"):
                continue
            if "Un" in chrom: continue
            if "random" in chrom or "HLA" in chrom: continue
            if chrom.endswith("alt"): continue
            chroms[chrom] = True
            xs, ys = [], []
            v50 = 0
            found = False
            for _, x, y in data:
                y = float(y)
                if y < 0.01:
                    continue
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
    try:
        with open(args.output, "w") as html:
            divs = "\n".join("<{div}>{chrom}</{div}><div id='plot-div-{chrom}'></div><hr/>".format(
                chrom=c, div="h2" if c == "total" else "b") for c in chroms)
            html.write(tmpl.substitute(showlegend="true" if len(
                sys.argv[1:]) < 20 else "false", plot_divs=divs))
            for chrom in chroms:
                html.write(string.Template(chr_tmpl).substitute(
                    chrom=chrom, data=json.dumps(traces[chrom])))
            html.write(footer)
    except FileNotFoundError:
        sys.exit("ERROR: failed creating output file, does the path exist?")


def get_args():
    parser = ArgumentParser(description="Creates html plots from mosdepth results.")
    parser.add_argument("-o", "--output",
                        default="dist.html",
                        help="path and name of output file. Directories must exist.")
    parser.add_argument("input",
                        nargs='+',
                        help="the dist file(s) to use for plotting")
    return parser.parse_args()


if __name__ == '__main__':
    main()
