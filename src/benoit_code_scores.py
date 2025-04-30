###################################
###################################
import numpy as np
import matplotlib.pyplot as plt
clear('essentialdomain')
for ii in np.arange(1,len(genes.coordinates)+1).reshape(-1):
    tnindex = find(tncoordinates_copy(:,2) > np.logical_and(genes.coordinates(ii,1),tncoordinates_copy(:,2)) < genes.coordinates(ii,2))
    tns(ii).data = np.array([[genes.coordinates(ii,1)],[tncoordinates_copy(tnindex,2)],[genes.coordinates(ii,2)]])
    tnnumber[ii] = len(tns(ii).data)
    genelength[ii] = genes.coordinates(ii,2) - genes.coordinates(ii,1)

mingap = 200
maxlength = 0.9
minlength = 0.1
mintnnumber = 20
skiptn = 5
for ii in np.arange(1,len(genes.coordinates)+1).reshape(-1):
    if tnnumber(ii) < skiptn + 1:
        tndif[ii] = 0
    else:
        tndif[ii] = np.amax(tns(ii).data(np.arange(skiptn + 1,len(tns(ii).data)+1)) - tns(ii).data(np.arange(1,len(tns(ii).data) - skiptn+1)))

tndif = double(tndif)
essentialdomain = np.multiply(double(tndif),tnnumber) / genelength ** 1.5
essentialdomain[tndif < np.logical_or[mingap,tndif / genelength] > np.logical_or[maxlength,tndif / genelength] < np.logical_or[minlength,tnnumber] < mintnnumber] = 0
dump,index = __builtint__.sorted(essentialdomain,'descend')
clear('dump')
plt.plot(essentialdomain(index),'.')
scipy.io.loadmat('names.mat')
scipy.io.loadmat('essentialgenesindex.mat')
names(index(np.arange(1,200+1)))
heatplot = np.zeros((len(essentialdomain),1))
heatplot[essentialgenes] = 1
plt.figure(2)
for ii in np.arange(1,6+1).reshape(-1):
    jj = ii * 1000
    subplot(21,1,ii * 3 - 2)
    image(rot90(256 * heatplot(index(np.arange(jj - 999,jj+1)))))
    colormap('gray')

plt.figure(3)
for ii in np.arange(1,6+1).reshape(-1):
    jj = ii * 1000
    subplot(21,1,ii * 3 - 1)
    image(50 * (essentialdomain(index(np.arange(jj - 999,jj+1)))))
    colormap('jet')
