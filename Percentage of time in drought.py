# First import all necessary packages
from netCDF4 import Dataset
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import time
%matplotlib inline

# Define moving average function
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

## Define Changing Parameters ##
VAR=['pr','s0','qtot']
var_names=['Precipitation','Soil Moisture','Runoff']

GCM=['CSIRO-BOM-ACCESS1-0','CNRM-CERFACS-CNRM-CM5','NOAA-GFDL-GFDL-ESM2M','MIROC-MIROC5']
gcm_names=['ACCESS1-0','CNRM-CM5','GFDL-ESM2M','MIROC-MIROC5']

BC_compliant=["ISIMIP2b","QME","CCAM",'MRNBC']

MONTH=['1-month',"24-month"]

## Define constant parameters
rcp='rcp85'

## Define file paths
data_path= "/g/data/eg3/st9855/drought_metrics"
mask_path= "/g/data/wj02/MISC/NRM_CLUSTERS/NRM_clusters.nc"

## Define cluster information
cluster_regions=[1,2,4,5,6,7,8,9]
cluster_index=[0,1,2,3,4,5,6,7]
cluster_titles=["Central Slopes",'East Coast','Murray Basin','Monsoonal North','Rangelands','Southern Slopes','Southern and South-Western Flatlands','Wet Tropics']
cluster_names=['CS','EC','MB','MN','R','SS','SSWF','WT']

# Define masks
mask      = xr.open_dataset(mask_path)["NRM_cluster"]
mask2      = np.tile(mask,(1,1)) ## Applied to 30 year chunks
mask1      = np.tile(mask,(1,1,1)) ## Applied to 30 year chunks

# For Burdekin case
cluster_index=[0]
# mask  = xr.open_dataset('/g/data/er4/st9855/mth3000/Burdekin_extent.nc')['Band1']
# mask2  = np.tile(mask,(1,1))
# mask1 = np.tile(mask,(1,1,1))

## Define year ranges
sel_years=[["1976-01-01","2005-12-31"],["2016-01-01","2045-12-31"],["2036-01-01","2065-12-31"],["2056-01-01","2085-12-31"],["2070-01-01","2099-12-31"]]
center_years=[1990,2030,2050,2070,2085]

# Define plot features
colors=['deepskyblue','sandybrown','yellowgreen','royalblue','peru','green']
symbols=['x','o','^','x','o','^']

# Loops through all clusters
for cluster in cluster_index:
    # Intialise figure and grid
    fig =plt.figure(figsize=(10,6),constrained_layout=True)
    grid = plt.GridSpec(6, 11)
    
    main_plot=fig.add_subplot(grid[:,:7])
    plot1=fig.add_subplot(grid[:3,7:])
    plot24=fig.add_subplot(grid[3:,7:])
    
    color_index=0 #keeps track of colors
    
    for month in MONTH:
        print('starting new month')
        for var in VAR:
            print('starting new var')
            
            #List to hold data for each 30 year period
            data_to_plot=[[],[],[],[],[]]
            
            for gcm in GCM:
                
                for bc in BC_compliant:
                    # Load dataset (loads yearly drought frequnecy: no. of months in drought)
                    # 1 and 24 month have slightly different paths/file names
                    if month!='1-month':
                        file_path=data_path+"/"+month+'/'+bc+"/"+gcm+"/frequency_"+month+"_"+bc+"_"+gcm+"_"+var+"_"+rcp+".nc"
                    else:
                        file_path=data_path+"/"+month+'/'+bc+"/"+gcm+"/frequency_"+bc+"_"+gcm+"_"+var+"_"+rcp+".nc"
                
                    dataset=xr.open_dataset(file_path)
        
                    # Loop through all 30 year periods
            
                    for y in range(len(sel_years)):
                        # Select desired time slice
                        data=dataset.sel(time=slice(sel_years[y][0],sel_years[y][1]))
                        data=data['timing']
                        
                        # Apply mask to select desired cluster
                        data=xr.where(mask1==cluster_regions[cluster],data,np.nan)
                        
                        # Sum over all 30 years (this removes mask)
                        data=data.sum(axis=0,skipna=True)
                        
                        #Rapply mask
                        
                        data=xr.where(mask2==cluster_regions[cluster],data,np.nan)
                        data_to_plot[y].append(data)
                        
                        # To do this for all of Australia make so both mask just look for values >0
                        #data=xr.where(mask1>0,data,np.nan)
                        #data=xr.where(mask2>0,data,np.nan)
                    
            # Initialise lists to hold stats
            historical=data_to_plot[0]
            medians=[]
            lowers=[]
            uppers=[]
            
            # Loops through all periods to calculate relevant stats
            # Divide by 360 (number of months in 30 years) and * by 100 to get percentage then - 50/3 (historical average)
            # Adding median, 10th quantile and 90th quantile
            for i in data_to_plot:
                medians.append(np.nanmedian(i)/360*100-(50/3))
                lowers.append(np.nanquantile(i,0.1)/360*100-(50/3))
                uppers.append(np.nanquantile(i,0.9)/360*100-(50/3))
                
            # Define label names    
            if month=="1-month":
                label_name='Short-term '
            if month=="24-month":
                label_name='Long-term '
            if var=="pr":
                label_name=label_name+"rainfall extreme dry"
            if var=="qtot":
                label_name=label_name+"runoff extreme dry"
            if var=="s0":
                label_name=label_name+"soil moisture extreme dry"
            
            # Plot main plot
            # Only plotting 
            main_plot.scatter(center_years[1:],medians[1:],label=label_name,color=colors[color_index],marker=symbols[color_index])
            main_plot.plot(center_years[1:],medians[1:],linestyle='-',alpha=0.7,color=colors[color_index])
            
            # Plot 1-month 
            if month=="1-month":
                plot1.plot(center_years[1:],medians[1:],linestyle='-',alpha=0.7,color=colors[color_index],marker=symbols[color_index])
                plot1.fill_between(center_years[1:],uppers[1:],lowers[1:],color=colors[color_index],alpha=0.1)
                plot1.scatter(center_years[1:],medians[1:],color=colors[color_index],marker=symbols[color_index])
                
            # Plot 24-month
            if month=="24-month":
                plot24.plot(center_years[1:],medians[1:],linestyle='-',alpha=0.7,color=colors[color_index],marker=symbols[color_index])
                plot24.fill_between(center_years[1:],uppers[1:],lowers[1:],color=colors[color_index],alpha=0.1)
                plot24.scatter(center_years[1:],medians[1:],color=colors[color_index],marker=symbols[color_index])
                
            # Update colour index    
            color_index+=1

    plot1.set_title('Short-term')
    plot1.set_xticks([2030,2050,2070,2085])
    
    plot24.set_title('Long-term')
    plot24.set_xticks([2030,2050,2070,2085])
    
    main_plot.set_xticks([2030,2050,2070,2085])
    main_plot.set_ylabel("Change in % of time in extreme dry conditions")
    
    fig.suptitle(f'Change in % of time in extreme dry conditions from 1976-2005',fontsize=16)
    fig.tight_layout()
    fig.legend(bbox_to_anchor=(0.85,-0.05),ncol=2)
    
    plt.savefig(f'/g/data/er4/st9855/mth3000/plots/example_bomplots/{cluster_names[cluster]}_updated_timeline2_{rcp}.png',bbox_inches='tight',transparent=False)
    #plt.savefig(f'/g/data/er4/st9855/mth3000/plots/cluster_plots/region_figures/{cluster_names[cluster]}/{cluster_names[cluster]}_updated_timeline2_{rcp}.png',bbox_inches='tight')
    #plt.savefig(f'/g/data/er4/st9855/mth3000/plots/burdekin_timeline_{rcp}.png',bbox_inches='tight',transparent=False)
