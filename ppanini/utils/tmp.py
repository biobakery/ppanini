import matplotlib
from matplotlib import pyplot

x = [1,1,2,2,3,3,4,4]
y = [1,1,2,2,3,3,4,4]
color = ['red','b','red','b','red','b','red','b']
pyplot.figure()
pyplot.scatter([x[0],x[2],x[5]],y[:3], c=color[:3], zorder=1)
pyplot.scatter(x[3:],y[3:], c=color[3:], zorder=2)
pyplot.savefig('tmp.pdf')