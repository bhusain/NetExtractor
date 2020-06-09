import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import itertools

import numpy as np
from scipy import linalg
'''
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
'''
import pandas as pd
from random import randint
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import normalized_mutual_info_score

from matplotlib import pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import matplotlib as mpl
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.stats import gaussian_kde
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

from scipy.stats import kde
from mpl_toolkits.mplot3d import Axes3D


from sklearn import mixture

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from sklearn import metrics
import pandas as pd
import multiprocessing as mp


color_iter = ['navy', 'darkorange']

#Read input GEM

data = pd.read_csv("./GTEx_v7_brain_subGEM-log-no.txt", sep='\t')
data = data.iloc[:, 1:]
data = data.transpose()
data = data.fillna(0)
data = data.fillna(0)
data[data < 0] = 0
global_max = data.max().max()
global_min = data.min().min()

print(data.shape, global_max, global_min)
data = data.fillna(0)
#bins = 2 * int(math.ceil(global_max))
#print(bins)

count = 0

threshold = 1400

f = open('NetExtractor_output.txt', 'w')

index_list = []

#Create list of genes that satisfy threshold criteria
for i in range(0, 56200):
    count_i = sum(m > 0 for m in data.ix[:, i])
    if count_i > threshold:
        index_list.append(i)


def my_func(x):
    for j in index_list:
        if j > x:
N = np.c_[data.ix[:, x], data.ix[:, j]]
            #Calculate Mututual information of all samples
            value = normalized_mutual_info_score(N[:, 0], N[:, 1])
            #Calculate mixture models. Specifiy the number of clusters: n_components
            dpgmm = mixture.BayesianGaussianMixture(n_components=5, covariance_type='full', max_iter=100, ).fit(N)
            Y_ = dpgmm.predict(N)
            unique, counts = np.unique(Y_, return_counts=True)
            S = metrics.silhouette_score(N, Y_)
            x_1 = []
            y_1 = []
            x_2 = []
            y_2 = []
            #Ensure enough samples fall into each cluster and have a minimum intercluster score
            if len(counts) == 2 and counts[0] > 30 and counts[1] > 30 and S > 0.5:
                count_cluster_1 = 0
                count_cluster_2 = 0
                covar = dpgmm.covariances_
                mean = dpgmm.means_
                v, w = linalg.eigh(covar[0])
                v = 2. * np.sqrt(2.) * np.sqrt(v)
                u = w[0] / linalg.norm(w[0])
                angle = np.arctan(u[1] / u[0])
                angle = 180. * angle / np.pi  # convert to degrees
                cos_angle = np.cos(np.radians(180. + angle))
                sin_angle = np.sin(np.radians(180. + angle))
                xc = mean[1][0] - mean[0][0]
                yc = mean[1][1] - mean[0][1]
                xct = xc * cos_angle - yc * sin_angle
                yct = xc * sin_angle + yc * cos_angle
                rad_cc = (xct ** 2 / (v[0] / 2.) ** 2) + (yct ** 2 / (v[1] / 2.) ** 2)
                #Ensure clusters are mutually separate using their mean and variance
                if rad_cc > 1:
                    #Write the edge, MI value, and intercluster score
                    string_val = '%s,%s,%0.3f,%0.3f\n' % (str(data.columns.values[x]), str(data.columns.values[j]), value, S)
                    f.write(string_val)
                    f.flush()


#Calling multiprocessing function. Specify number of threads to launch the job.
def main():
    pool = mp.Pool(20)
    result = pool.map(my_func, index_list)





if __name__ == "__main__":
  main()
