import seaborn as sns
from matplotlib import pyplot as plt
sns.set_style('whitegrid')


stimes = """27:09.25
23:05.32
15:00.63
13:05.17
11:26.52
12:06.25
11:37.36
11:35.57
11:38.39"""

def to_secs(t):
    m, s = map(float, t.split(":"))
    return s + 60 * m

times = [to_secs(t) / 60. for t in stimes.strip().split()]

ideal = [(i, times[0] / (i+1)) for i, x in enumerate(times)]

plt.figure(figsize=(6,4))
plt.plot(times, label="actual")
xs, ys = zip(*ideal)
plt.plot(xs, ys, label="ideal", ls="--", color="0.5")
plt.xlabel("additional decompression threads")
plt.ylabel("time (minutes)")
plt.legend()
plt.ylim(ymin=0)
plt.xlim(xmin=0, xmax=8)
plt.title("time to calculate coverage")
plt.savefig("mosdepth-scaling.eps")
plt.savefig("mosdepth-scaling.png")
plt.show()

