import matplotlib.pyplot as plt
import numpy as np

# projection method: Chapter 3, Grimm & Schreiber, 2002
basis = []
for i in range(5):
    a = np.pi * 0.4 * i
    basis.append([np.cos(a), np.sin(a), np.cos(2 * a), np.sin(2 * a), np.sqrt(0.5)])
basis = np.transpose(basis)


def lattice(en=2):  # 5D lattice
    s = (en - 1) / 2
    lps = []
    for i in range(en):
        for j in range(en):
            for m in range(en):
                for n in range(en):
                    for l in range(en):
                       lps.append([i - s, j - s, m - s, n - s, l-s])
    return np.array(lps)

def polytope():  # of 20 faces
    c0 = [27, 30, 23, 29, 15]
    c1 = [1, 4, 16, 2, 8]
    c2 = [26, 18, 22, 20, 21, 5, 13, 9, 11, 10]
    fs = []
    for i in range(5):
        fs.append([31, c0[i], c0[(i + 1) % 5]])
        fs.append([0, c1[i], c1[(i + 1) % 5]])
    for i in range(10):
        fs.append([c2[i], c2[(i + 1) % 10], c2[(i + 2) % 10]])
    return fs


def normalize(v):
    norm = np.linalg.norm(v)
    return v / norm


def window(allpoints):
    m_orth = np.array(basis[2:])
    pps= lattice().dot(m_orth.transpose())
    fs = polytope()
    face_ps = []
    face_ns = []    # each (face_ps[i], face_ns[i]) pair describes a face of the 3D polytope
    for f in fs:
        p0 = pps[f[0]]
        p1 = pps[f[1]]
        p2 = pps[f[2]]
        face_ps.append(p1)
        nor = normalize(np.cross(p0 - p1, p2 - p1))
        if np.dot(p1, nor) > 0:
            nor = nor* -1
        face_ns.append(nor)
    inside = []
    for p in allpoints:
        pp = np.dot(m_orth, p)
        flag = True
        for j in range(len(fs)):
            if np.dot(pp - face_ps[j], face_ns[j]) < 0:
                flag = False
                break
        inside.append(flag)
    return inside


poly = polytope()  # this polytope will cut the entire 5D lattice
ps = lattice(3)  #3 is very fast, 5 takes 10 seconds, lattice points (5D)
nodesize = ps.shape[0]
edges = []
for i in range(nodesize):
    for j in range(i + 1, nodesize):
        dis = np.linalg.norm(ps[i] - ps[j])
        if 0.99 < dis < 1.01:
            edges.append([i, j])

inside = window(ps)
pro = basis[:2]
for ida, idb in edges:  # project all inside points & windowing the 5D lattice
    if inside[ida] and inside[idb]:
        a = pro.dot(ps[ida])
        b = pro.dot(ps[idb])
        plt.plot([a[0], b[0]], [a[1], b[1]], color='C0')
plt.axis('scaled')
plt.show()
