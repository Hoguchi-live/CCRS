import json
import matplotlib.pyplot as plt

DIR = "./files/"
FILENAME = "timings.json"

with open(DIR+FILENAME) as f:
    timings = json.load(f)

x = sorted([int(l) for l in timings.keys()])[::]
print(x)
y = [timings[str(l)]["forward"] for l in x]
print(y)
plt.scatter(x, y)
plt.yscale("log")
plt.xscale("log")
plt.show()
