
import matplotlib.pyplot as plt
import numpy as np


from data_import_tools import per_gene, wiggle,profile,correct


## Importing pergene files 

pergene_files=["WT_merged_pergene_insertions.txt"]

wig_files=["WT_merged.wig"]

# Import data

gene_inserts=[]
wig_inserts=[]
for i in range(len(pergene_files)):
    gene_inserts.append(per_gene(pergene_files[i]))
    wig_inserts.append(wiggle(wig_files[i]))




# +
## Plots for the figure of the average transposon density based on the missing sites 

# +
maxdist = 200000
chr_file="../data/chromosome_sizes.txt"
cen_file="../data/centromere_locations.txt"
# Compute the profile and correct it 

c_dis, tn_csum, cen = profile(wig_inserts,chr_file=chr_file,cen_file=cen_file)
poly_fit = correct(c_dis, tn_csum,method="poly",maxdist=1000000)
poly_der = np.polyder(poly_fit, 1)

# +
 ## Linear fit of the centromere bias 

for i in np.arange(0,tn_csum.shape[0]):
    plt.plot(c_dis,tn_csum[i],color="grey",alpha=0.3)

mean_tn_csum=np.nanmean(tn_csum,axis=0)
plt.plot(c_dis,mean_tn_csum,color="blue",linewidth=1,label="Average")

x=np.arange(0,400000,10000)
x2fit=np.arange(200000,400000,10000)
model = np.polyfit(x2fit, mean_tn_csum[20:40], 1)
plt.plot(x,model[0]*x+model[1],color="red",label="Linear fit",alpha=0.8)
#plt.plot(model[0]*x+model[1],color="blue",label="Linear fit")
plt.text(200000, 5000, 'y=%.3fx+%.1f' % (model[0], model[1]), fontsize=10)
xticks_values = [100000, 200000, 300000, 400000]  # Specify the x-axis tick values
xticks_labels = np.round(np.divide(xticks_values,1000),0)  # Specify the labels for the tick values

# Set the x-axis tick locations and labels
plt.xticks(xticks_values, xticks_labels)
plt.xlim(0,400000)
plt.ylim(0,25000)
plt.legend()

# plt.savefig("../figures/centromere_bias_linear_fit.png",dpi=300)

# +
## average cumulative insertions vs distance to centromere 

mean_tn_csum=np.nanmean(tn_csum,axis=0)
plt.plot(c_dis,mean_tn_csum,color="blue",linewidth=1,label="Average")
plt.plot(c_dis,poly_fit(c_dis),color="red",linewidth=1,label="Polynomial fit")
plt.xlim(0,400000)
plt.ylim(0,25000)

plt.xlabel("Distance to centromere (bp)")
plt.ylabel("Cumulative insertions")
plt.legend()

#plt.savefig("../figures/centromere_bias_polynomial_fit.png",dpi=300)

# +
## Plot insertion rate vs distance to centromere and polynomial fit

plt.plot(c_dis[0:-1],np.diff(mean_tn_csum)/10000,color="blue",linewidth=1,label="Average insertion rate")
plt.plot(c_dis,poly_der(c_dis),color="red",linewidth=1,label="Polynomial fit")
plt.xlim(0,400000)
plt.ylim(0,0.1)
plt.ylabel("Insertion rate ($bp^{-1}$)")
plt.xlabel("Distance to centromere (bp)")

plt.legend()



# +
## Expected insertion density minus the observed insertion density over essential and not essential genes

standard_essentials=np.loadtxt("standard_essentials.txt",dtype=str) 
diff_essentials=[]
diff=[]

for g in gene_inserts.index:
    
        # Get the info of the current gene
        gene_info = gene_inserts.loc[g]
        gene_size = gene_info['End location'] - gene_info['Start location']
        
        reads_b =  ( wig_inserts["Chrom"] == "chr"+gene_info["Chromosome"] ) & \
                ( wig_inserts["Position"] > gene_info['Start location']+gene_size*0.1 ) & \
                ( wig_inserts["Position"] < gene_info['Start location']+gene_size*0.9 )
        
        reads = wig_inserts["Reads"][reads_b]
        
        # Calculate smallest distance from the gene to the centromere
        curr_cen = cen[ "chr"+gene_info["Chromosome"] == cen["chrom"] ].reset_index()
        
        if len(curr_cen):
            
            gene_center = gene_info['Start location'] + gene_size/2
            d_cen = min( abs( gene_center - curr_cen["start"][0] ), 
                        abs( gene_center - curr_cen["stop"][0] ) )
            
            # Determine the expected insertion rate
            if d_cen <= maxdist:           
                E_tn = np.floor( poly_der(d_cen)*gene_size*0.8 ) 
                
            else:
                E_tn = np.floor( poly_der(maxdist)*gene_size*0.8  )
            if g in standard_essentials:
            # compute the differences from the expected insertion density and the observed insertion density
                diff_essentials.append(E_tn-len(reads))
            else:
                diff.append(E_tn-len(reads))

# +
plt.hist(diff,bins=600,histtype="stepfilled",color="gray",label="Non essential");
plt.hist(diff_essentials,bins=600,histtype="stepfilled",color="pink",label="Essential");

plt.xlabel("E(x)-O(x)")
plt.ylabel("Count")

plt.legend()
plt.tight_layout()

