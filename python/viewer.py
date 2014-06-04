#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np

def to_object(h, d):
    class Snapshot(object):
        pass

    out = Snapshot
    setattr(out, "time", float(h[0][2]))
    setattr(out, "niter", int(h[0][4]))
    setattr(out, "header", h)
    setattr(out, "data", d)
    for column, title in zip(d,h[2][1:]):
        setattr(out, title, column)

    return out

def chunkify(path = "output.log"):
    f = open(path, "r")
    data = []
    infos = []
    cols = []
    header = []
    # Iterate over the lines until we reach two empty lines
    for l in f.readlines():
        line = l.split()
        
        if line == []:
            if len(data) > 0:
                yield to_object(np.array(data[:3]),
                                np.array(data[3:]).transpose().astype(float))
                data = []
            else:
                chunkify()
        else:
            data.append(line)

if __name__ == "__main__":
    for snap in chunkify():
        plt.plot(snap.k, snap.T_0, label="$T_0$ @ t= %f" % snap.time)

    plt.legend()
    plt.show()
