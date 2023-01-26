# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3.8.10 ('satay-dev')
#     language: python
#     name: python3
# ---

# +
# Plot a histogram of the fitness differences between both methods 

## Normalize the SGA values to the value of HO in the SGA dataset (HO is not in the SGA dataset)
# -

keys=['dbem3_b',
 'WT_1-Benoit',
 'dnrp1_2',
 'bem1-aid_a',
 'dbem1dbem3_b',
 'wt_merged',
 'dbem3_a-trimmed',
 'dbem1dbem3_a',
 'bem1-aid-dbem3_a',
 'bem1-aid_merged',
 'bem1-aid-dbem3_b',
 'dnrp1_1',
 'wt_b',
 'wt_a',
 'dnrp1_merged',
 'WT_2-Benoit',
 'bem1-aid_b',
 'dbem3_merged',
 'dbem3_a']

# +
import pandas as pd
sga_fitness=pd.read_excel('../data/fitness-SGD.xlsx',index_col=0)
import pickle
with open("../postprocessed-data/fitness_models_all_backgrounds", "rb") as fp:   # Unpickling
    b = pickle.load(fp)

satay_fitness=pd.concat(b,axis=0,keys=keys)
# -

satay_wt=satay_fitness.loc["wt_merged"]


# +
# satay_wt=satay_fitness[satay_fitness.loc[:,"background"]=="wt_merged"]
# satay_wt.index=satay_wt.loc[:,"Gene name"]

# #satay_wt_ho=satay_wt[satay_wt.loc[:,"Gene name"]=="HO"]
# satay_wt_ho=satay_wt.loc["HO"]

#values2ho=satay_wt["fitness"]/satay_wt_ho["fitness"].values[0]

# +
#satay_wt["fitness2HO"]=satay_wt["fitness"]/satay_wt_ho["fitness"]

# +
from collections import defaultdict
import numpy as np

diff_sga_satay=defaultdict(dict)
for i in sga_fitness.index:
    if i.upper() in satay_wt.index and satay_wt.loc[i.upper(),"fitness_gene"]!="Not enough flanking regions":
        if type(sga_fitness.loc[i,"fitness"])!=np.float64:
            diff_sga_satay[i]["difference"]=(sga_fitness.loc[i,"fitness"].unique()[0]-satay_wt.loc[i.upper(),"fitness_gene"])
            diff_sga_satay[i]["sga"]=sga_fitness.loc[i,"fitness"].unique()[0]
            diff_sga_satay[i]["satay"]=satay_wt.loc[i.upper(),"fitness_gene"]
        else:
            diff_sga_satay[i]["difference"]=(sga_fitness.loc[i,"fitness"]-satay_wt.loc[i.upper(),"fitness_gene"])
            diff_sga_satay[i]["sga"]=sga_fitness.loc[i,"fitness"]
            diff_sga_satay[i]["satay"]=satay_wt.loc[i.upper(),"fitness_gene"]
    else:
        diff_sga_satay[i]["difference"]=np.nan

# +
from collections import defaultdict
import numpy as np

diff_sga_satay=defaultdict(dict)
for i in sga_fitness.index:
    if i.upper() in satay_wt.index and satay_wt.loc[i.upper(),"fitness_gene"]!="Not enough flanking regions":
        if type(sga_fitness.loc[i,"fitness"])!=np.float64:
            diff_sga_satay[i]["difference"]=(sga_fitness.loc[i,"fitness"].unique()[0]-satay_wt.loc[i.upper(),"fitness_domains_corrected"])
            diff_sga_satay[i]["sga"]=sga_fitness.loc[i,"fitness"].unique()[0]
            diff_sga_satay[i]["satay"]=satay_wt.loc[i.upper(),"fitness_gene"]
        else:
            diff_sga_satay[i]["difference"]=(sga_fitness.loc[i,"fitness"]-satay_wt.loc[i.upper(),"fitness_domains_corrected"])
            diff_sga_satay[i]["sga"]=sga_fitness.loc[i,"fitness"]
            diff_sga_satay[i]["satay"]=satay_wt.loc[i.upper(),"fitness_domains_corrected"]
    else:
        diff_sga_satay[i]["difference"]=np.nan
# -

diff_sga_satay_pd=pd.DataFrame(diff_sga_satay).T

diff_sga_satay_pd

# +
from cmath import isfinite


data=diff_sga_satay_pd.dropna(how='all')
data=data[~data.isin([np.nan, np.inf, -np.inf]).any(1)]


# +
data_satay=data["satay"].astype(float)
data_sga=data["sga"].astype(float)
data_difference=data["difference"].astype(float)

# remove nan and inf .values()
data_satay=data_satay[data_satay.apply(lambda x: isfinite(x))]
data_sga=data_sga[data_sga.apply(lambda x: isfinite(x))]
data_difference=data_difference[data_difference.apply(lambda x: isfinite(x))]

# remove genes that are not in both datasets
for i in data_satay.index:
    if i not in data_sga.index:
        data_satay.drop(i,inplace=True)
        

# +
from scipy.stats import pearsonr

corr, p = pearsonr(data_sga, data_satay)
corr,p

# +
import seaborn as sns
from matplotlib import pyplot as plt

plt.figure(figsize=(10,10))
sns.heatmap(data_sga[200:300].to_frame().join(data_satay[200:300].to_frame()),
cmap="RdBu_r",vmin=0,vmax=1.4)


# +

# make a 2d histogram for the scatter plot !!!

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from scipy.stats import norm

##data2fit = data.loc[:,"difference"]
data2fit=data_difference
mu, std = norm.fit(data2fit)

fig,ax=plt.subplots(1,3,figsize=(12,3))
# plt.subplots_adjust(wspace=0.5)


# Plot the PDF.
#ax[2].hist(data.loc[:,"difference"],bins=100,alpha=0.4,density=True,color="gray");
ax[2].hist(data_difference,bins=100,alpha=0.4,density=True,color="gray");
xmin, xmax = ax[2].get_xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)

ax[2].plot(x, p, 'k', linewidth=2)


ax[2].text(-0.8, 1, r'$\mu=%.2f,\ \sigma=%.2f$' % (mu, std),fontsize=12)

ax[2].set_xlabel("Fitness SGA - Fitness SATAY in WT")
ax[2].set_ylabel("Density")
ax[2].set_xlim(-1,1.2)
#ax[0].set_title("Distribution of fitness differences")

ax[1].hist2d(data_sga,data_satay,bins=(100,100),vmin=0,
vmax=1.4,cmap=plt.cm.Blues,density=True);

fig.colorbar(ax[1].collections[0],ax=ax[1])

ax[1].plot([0,1.4],[0,1.4],color="red",alpha=0.2)

#ax[1].scatter(data.loc[:,"sga"],data.loc[:,"satay"],s=10,alpha=0.2,color="gray");
#sns.jointplot(data.loc[:,"sga"],data.loc[:,"satay"],kind="scatter",s=10,alpha=0.2,ax=ax[1]);
ax[1].set_xlabel("SGA fitness")
ax[1].set_ylabel("SATAY fitness")

ax[1].set_xlim(0,1.4)
ax[1].set_ylim(0,1.4)
ax[1].text(0,1.2,"Pearson corr={:.2f}".format(corr))
#ax[1].set_title("How much are the fitness values different?")



ax[0].hist(data_sga,bins=100,alpha=0.6,density=True,label="SGA");
ax[0].hist(data_satay,bins=100,alpha=0.6,density=True,label="SATAY");

axins0 = ax[0].inset_axes([0.6, 0.2, 0.3, 0.3])
# axins0.errorbar(["SGA","SATAY"],[data.loc[:,"sga"].mean(),
# data.loc[:,"satay"].mean()],yerr=[data.loc[:,"sga"].std(),
# data.loc[:,"satay"].std()],fmt='o',alpha=0.6,
#               capsize=5,capthick=2,elinewidth=2);
axins0.errorbar(["SGA"],data_sga.mean(),yerr=data_sga.std(),fmt='o',alpha=0.6,
              capsize=5,capthick=2,elinewidth=2);
axins0.errorbar(["SATAY"],data_satay.mean(),yerr=data_satay.std(),fmt='o',alpha=0.6,
                capsize=5,capthick=2,elinewidth=2);
#axins0.set_ylabel("Normalized fitness to WT")
axins0.set_xlim(-0.5,1.5)
axins0.set_xticks([0,1])
axins0.set_xticklabels(["",""])
axins0.set_ylabel("$\mu \pm \sigma$")
ax[0].set_xlabel("Normalized fitness to WT")
ax[0].set_ylabel("Density")
ax[0].legend()
#plt.tight_layout()

# for axes in ax:
#     axes.set_aspect('equal', adjustable='box')

#fig.savefig("../figures/fig_Fitness_SGA_vs_SATAY.png",dpi=300,transparent=True)

# +
import plotly.express as px

df = px.data.tips()

fig = px.density_heatmap(data, x="sga", y="satay", color_continuous_scale="Viridis",
width=400, height=400)
fig.show()
# -


