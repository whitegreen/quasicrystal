import matplotlib.pyplot as plt
import numpy as np

# projection method: Chapter 3, Grimm & Schreiber, 2002, Fig3.9
k = 0.5 * np.sqrt(2)
mat = np.array([[k, 0.5, 0, -0.5], [0, 0.5, k, 0.5], [k, -0.5, 0, 0.5], [0, 0.5, -k, 0.5]])

def lattice(en=2):  # 4D lattic
    s = (en - 1) / 2
    lps = []
    for i in range(en):
        for j in range(en):
            for m in range(en):
                for n in range(en):
                    lps.append([i - s, j - s, m - s, n - s])
    return np.array(lps)


def polytope():
    ps = lattice()
    m = mat[2:]
    h = ps[[9, 13, 5, 4, 6, 2, 10, 11]]  # fansy indexing
    return h.dot(m.transpose())


def inside(p, vs):  # check if point p is inside polygon vs
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


poly = polytope()  # this polytope will cut the entire 4D lattice
ps = lattice(5)  # 5 is fast, 7 is slow,  lattice points (4D)
nodesize = ps.shape[0]
edges = []
for i in range(nodesize):
    for j in range(i + 1, nodesize):
        if np.linalg.norm(ps[i] - ps[j]) < 1.1:
            edges.append([ps[i], ps[j]])
m = mat[2:]
pro = mat[:2]
for pa, pb in edges:  # project all inside points & windowing the 4D lattice
    if inside(m.dot(pa), poly) and inside(m.dot(pb), poly):
        a = pro.dot(pa)
        b = pro.dot(pb)
        plt.plot([a[0], b[0]], [a[1], b[1]], color='C0')
plt.axis('scaled')
plt.show()
