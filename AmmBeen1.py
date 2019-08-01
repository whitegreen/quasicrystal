import matplotlib.pyplot as plt
import numpy as np

# projection method: Chapter 3, Grimm & Schreiber, 2002, Fig3.9
k = 0.5 * np.sqrt(2)
mat = [[k, 0.5, 0, -0.5], [0, 0.5, k, 0.5], [k, -0.5, 0, 0.5], [0, 0.5, -k, 0.5]]

def polytope():
    ps = []
    for i in range(2):
        for j in range(2):
            for m in range(2):
                for n in range(2):
                    ps.append([i - 0.5, j - 0.5, m - 0.5, n - 0.5])
    m = np.array(mat[2:])
    arr = []
    for i in (9, 13, 5, 4, 6, 2, 10, 11):
        v = m.dot(np.array(ps[i]))
        arr.append(v)
    return np.array(arr)


def inside(p, vs): # check if point p is inside polygon vs
    i = len(vs) - 1
    j = i
    oddNodes = False
    for i in range(len(vs)):
        if (vs[i][1] < p[1] and vs[j][1] >= p[1] or vs[j][1] < p[1] and vs[i][1] >= p[1]) and (
                vs[i][0] <= p[0] or vs[j][0] <= p[0]):
            if vs[i][0] + (p[1] - vs[i][1]) / (vs[j][1] - vs[i][1]) * (vs[j][0] - vs[i][0]) < p[0]:
                oddNodes = not oddNodes
        j = i
    return oddNodes


poly = polytope() # thispolytope will cut the entire 4D lattice
en = 5   # using 7 is very slow, poor python
s = en // 2
lps = []
for i in range(en):
    for j in range(en):
        for m in range(en):
            for n in range(en):
                lps.append([i - s, j - s, m - s, n - s])
ps = np.array(lps)  # lattice points (4D)
nodesize = ps.shape[0]
edges = []
for i in range(nodesize):
    for j in range(i + 1, nodesize):
        if np.linalg.norm(ps[i] - ps[j]) < 1.1:
            edges.append([ps[i], ps[j]])
m = np.array(mat[2:])
pro = np.array(mat[:2])
for pa, pb in edges:  # project all inside points
    if inside(m.dot(pa), poly) and inside(m.dot(pb), poly):
        a = pro.dot(pa)
        b = pro.dot(pb)
        plt.plot([a[0], b[0]], [a[1], b[1]], color='C0')
plt.axis('scaled')
plt.show()
