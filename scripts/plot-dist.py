import sys
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
import sys
figpath = "dist.png"

sns.set_style('whitegrid')
sns.set_palette('Set1', 13)

plt.axhline(y=0.5, zorder=0, color='#aaaaaa')
for f in sys.argv[1:]:
    xs, ys = [], []
    sample = f.split(".")[0]
    found = False
    for x, y in (x.rstrip().split("\t") for x in open(f)):
        y = float(y)
        if y < 0.01: continue
        if not found and y > 0.5:
            found = True
            print f, x, y

        xs.append(float(x))
        ys.append(y)

    plt.plot(xs, ys, label=sample)
plt.xlabel("Coverage")
plt.ylabel("Proportion of bases at coverage")
plt.xlim(xmin=0)
plt.ylim(ymin=0, ymax=1)
plt.savefig(figpath)
