import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Files to be read...
infile1='calculated_d18O_2015'          #model output
infile2='calculated_d18O_2006.nc'       #from site
infile3='file1'

# Parameters
im = 360
jm = 180
lm = 33
sizeTitle = 80
offSet = 2
sizeData = im*jm
d1 = np.float32
test=2

d18o_model=np.zeros((360,180,33))
d18o_site_2=np.zeros((360,180,33))
d18o_diff=np.zeros((360,180))

sizeData3=259208

# import files...
f1=open(infile1,"rb")
f2=Dataset(infile2, mode='r')
f3=open(infile3,"rb")


# Time for some data analysis! let's start with the binary model output...
f1.read(offSet)                    #skip first 2 bytes of nonsense
header = f1.read(sizeTitle+2) #read the first 80+2 bytes for title 
#print header
record = np.fromfile(f1,dtype=d1,count=sizeData) #read the first record from the Fortran file
recordT = np.reshape(record,(180,360)).transpose()  #turn the 1D input into two dimesnions and transpose for lon,lat
d18o_model[:,:,0] = recordT
#print d18o_model[:,:,0]

# The next records are more predictable and can be read with a forloop
for k in range(1,3):
    f1.read(8)
    header = f1.read(sizeTitle)
    record = np.fromfile(f1,dtype=d1,count=sizeData)
    recordT=np.reshape(record,(180,360)).transpose()
    d18o_model[:,:,k]=recordT
#    print d18o_model[:,:,k]

if(test==2): #test for estimated D18
    infile4='estimated_d18O_2015'
    f4=open(infile4,"rb")
    d18o_est=np.zeros((360,180,33))
    d18o_est_calc=np.zeros((360,180))
    d18o_est_site=np.zeros((360,180))
    
    print("opening estimated d180..")    
    
    f4.read(offSet)                    #skip first 2 bytes of nonsense
    header = f4.read(sizeTitle+2) #read the first 80+2 bytes for title 
    print header
    record = np.fromfile(f4,dtype=d1,count=sizeData) #read the first record from the Fortran file
    recordT = np.reshape(record,(180,360)).transpose()  #turn the 1D input into two dimesnions and transpose for lon,lat
    d18o_est[:,:,0] = recordT
    print d18o_model[:,:,0]
    
    # The next records are more predictable and can be read with a forloop
    for k in range(1,3):
        f4.read(8)
        header = f4.read(sizeTitle)
        record = np.fromfile(f4,dtype=d1,count=sizeData)
        recordT=np.reshape(record,(180,360)).transpose()
        d18o_est[:,:,k]=recordT
    #    print d18o_model[:,:,k]



# Next the siteloaded netCDF... easy!
lons = f2.variables['lon'][:]
lats = f2.variables['lat'][:]
d18o_site = f2.variables['d18o'][:]

# site output...
f3.read(offSet+2)                    #skip first 2 bytes of nonsense
header = f3.read(sizeTitle) #read the first 80+2 bytes for title 
#print header
record = np.fromfile(f3,dtype=d1,count=sizeData) #read the first record from the Fortran file
recordT = np.reshape(record,(180,360)).transpose()  #turn the 1D input into two dimesnions and transpose for lon,lat
d18o_site_2[:,:,0] = recordT
#print d18o_site_2[:,:,0]

for k in range(1,lm):
    f3.read(8)
#    print header
    header = f3.read(sizeTitle)
    record = np.fromfile(f3,dtype=d1,count=sizeData)
    recordT=np.reshape(record,(180,360)).transpose()
    d18o_site_2[:,:,k]=recordT
#    print d18o_site_2[:,:,k]



    
# Scrub out -1e30 for NaN
for j in range(0,jm):
    for i in range(0,im):
        for k in range(0,lm):
            # Scrub out -1e30 for NaN for model output
            if(d18o_model[i,j,k]<=-1.000000e+30):
                d18o_model[i,j,k]=np.nan
            if(d18o_site_2[i,j,k]<=-1.000000e+30):
                d18o_site_2[i,j,k]=np.nan

if (test==2):
    for j in range(0,jm):
        for i in range(0,im):
            for k in range(0,lm):
                # Scrub out -1e30 for NaN for model output
                if(d18o_est[i,j,k]<=-1.000000e+30):
                    d18o_est[i,j,k]=np.nan
    d18o_est=np.ma.masked_invalid(d18o_est)
    print("testing difference between est and calc...")
    for j in range(0,jm):
        for i in range(0,im):
            # Between calc and est   
            if(np.isnan(d18o_model[i,j,0]) or np.isnan(d18o_est[i,j,0])):
                d18o_est_calc[i,j]=np.nan  # ignore NaN
            else:
                d18o_est_calc[i,j]=d18o_model[i,j,0]-d18o_est[i,j,0]
                if(d18o_est_calc[i,j]>-.01) and (d18o_est_calc[i,j]<.01):
                    d18o_est_calc[i,j]=np.nan
            # Between site and est        
            if(np.isnan(d18o_est[i,j,0]) or np.isnan(d18o_site[0,j,i])):
                d18o_est_site[i,j]=np.nan  # ignore NaN
            else:
                d18o_diff[i,j]=d18o_est[i,j,0]-d18o_site[0,j,i]
                if(d18o_est_site[i,j]>-.01) and (d18o_diff[i,j]<.01):
                    d18o_est_site[i,j]=np.nan        
            
      # Mask the NaN values
    d18o_est_calc=np.ma.masked_invalid(d18o_est_calc)
    

# Mask the NaN values
d18o_model=np.ma.masked_invalid(d18o_model)
d18o_site_2=np.ma.masked_invalid(d18o_site_2)


# Calculate difference between outputs
for j in range(0,jm):
    for i in range(0,im):
        if(np.isnan(d18o_model[i,j,0]) or np.isnan(d18o_site[0,j,i])):
            d18o_diff[i,j]=np.nan  # ignore NaN
        else:
            d18o_diff[i,j]=d18o_model[i,j,0]-d18o_site[0,j,i]
            if(d18o_diff[i,j]>-.01) and (d18o_diff[i,j]<.01):
                d18o_diff[i,j]=np.nan
  # Mask the NaN values
d18o_diff=np.ma.masked_invalid(d18o_diff) 
             
print "Artcic Ocean... ",d18o_model[78,63,0]


#plotting nonsese
lon = np.arange(-180.0,180.0,1)
lat = np.arange(-90.0,90.0,1)
m1 = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,
             llcrnrlon=-180,urcrnrlon=180,resolution='c')
x,y = m1(*np.meshgrid(lon, lat))
colors=[-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-.8,-.6,-.4,
        -.2,.2,.4,.6,.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.7]
fig1 = plt.figure()

# Make map for model output
m1.contourf(x, y, np.transpose(d18o_model[:,:,0]),colors)
m1.drawcoastlines()
m1.fillcontinents()
m1.drawmapboundary()
plt.title('Calculated Surface d18O Seawater')
cbar = plt.colorbar(orientation='horizontal', extend='both')
cbar.ax.set_xlabel('O18/O16 permilli')
plt.savefig('model_out.png', dpi = 300)
plt.show()

# Make map for site output
m1.contourf(x, y, d18o_site[1,:,:],colors)
m1.drawcoastlines()
m1.fillcontinents()
m1.drawmapboundary()
plt.title('Surface d18O Seawater from Site')
cbar = plt.colorbar(orientation='horizontal', extend='both')
cbar.ax.set_xlabel('O18/O16 permilli')
plt.savefig('site_out.png', dpi = 300)
plt.show()

if(test==2):
    # Make map for estimate-calc dif
    m1.contourf(x, y, np.transpose(d18o_est[:,:,0]),colors)
    m1.drawcoastlines()
    m1.fillcontinents()
    m1.drawmapboundary()
    plt.title('Estimated Surface d18O Seawater')
    cbar = plt.colorbar(orientation='horizontal', extend='both')
    cbar.ax.set_xlabel('O18/O16 permilli')
    plt.savefig('Estimate', dpi = 300)
    plt.show()
    # Make map for estimate
    m1.contourf(x, y, np.transpose(d18o_est_calc[:,:]),colors)
    m1.drawcoastlines()
    m1.fillcontinents()
    m1.drawmapboundary()
    plt.title('Surface d18O Seawater Est-Cal discrepencies ')
    cbar = plt.colorbar(orientation='horizontal', extend='both')
    cbar.ax.set_xlabel('O18/O16 permilli')
    plt.savefig('EstimateCalcAnomoly', dpi = 300)
    plt.show()    
    # Make map for estimate-calc dif
    m1.contourf(x, y, np.transpose(d18o_est_site[:,:]),colors)
    m1.drawcoastlines()
    m1.fillcontinents()
    m1.drawmapboundary()
    plt.title('Surface d18O Seawater est-site difference')
    cbar = plt.colorbar(orientation='horizontal', extend='both')
    cbar.ax.set_xlabel('O18/O16 permilli')
    plt.savefig('EstimateSiteAnomoly', dpi = 300)
    plt.show()
    
    


# Make map for difference
m1.contourf(x, y,np.transpose(d18o_diff[:,:]),colors)
m1.drawcoastlines()
m1.fillcontinents()
m1.drawmapboundary()
plt.title('Surface d18O Seawater Anomoly')
cbar = plt.colorbar(orientation='horizontal', extend='both')
cbar.ax.set_xlabel('O18/O16 permilli')
plt.savefig('difference.png', dpi = 300)
plt.show()
