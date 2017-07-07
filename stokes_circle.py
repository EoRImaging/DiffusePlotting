# creates a circle of points.
import numpy as np
import matplotlib.pyplot as plt
T = [100]
R = [1]
def rtpairs(r, n):

    for i in range(len(r)):
       for j in range(n[i]):
        yield r[i], j*(np.pi / n[i])

def coords(x, y):
    ax = []
    ay = []
    for r, t in rtpairs(R, T):
        ax.append(r*np.cos(t))
        ay.append(r*np.sin(t))
        plt.plot(ax, ay, 'bo')
    return [ax, ay]
plt.show()
ax, ay = coords('ax', 'ay')
ax = np.array(ax)
ay = np.array(ay)
# print np.power(ax, 2), np.power(ay, 2)
Q = np.power(ax, 2) - np.power(ay,2)
U = 2 * ax * ay
ay2 = -.5 + (np.sqrt(np.power(Q, 2) + np.power(U, 2)) / 2)
arc = plt.plot(Q, U, 'bo')
plt.show(arc)
