import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os
from mpl_toolkits.basemap import Basemap
import colorsys
from scipy.interpolate import spline

"""
Pt I
This is a program that takes in a list of properly formatted CSV files(toAdd), 
to append to the Global Seawater d18O database (infile) using Pandas.

This thanks to Pandas, there are dynamic database (DB) querrying for any of 
the 11 fields-- including seasonality, reference, and location.

The database is converted to Pandas dataframe(DF) data structure, 
The DF can be accessed by typing 'df' into the command line after running.

Pandas has thorough documentation on their website below:
--> http://pandas.pydata.org/pandas-docs/stable/tutorials.html

The querry results can be plotted globally using Basemaps. 
There is an example available by turning the 'seasons' switch to 1.

the integer seasons is a logic switch(much like a boolean), used to access
specific functions in this program for writing out files, searching the DF

At this time, this program is intended to update these files...

Ref                 File Name               # of Recs       Compiler    Dates
Strauss 2015        'Strauss2015.csv'       271  records    L           7.17.15
Cases 2012          'Cases2012.csv'         485  records    C   
COX 2010,2001       'COX2010_2001.csv'      88   records    C
Conroy 2012         'Conroy2012.csv'        191  records    S
Tiwari 2013         'Tiwari2013.csv'        76   records    S           7.20.15
Brown 2014          'BROWN2014.csv'         1207 records    C           7.25.15
Geotraces 2014...
Rijkenberg 2011     'GeoRijkenberg2011.csv' 432 records     C
Jenkins 2010        'GeoJenkins2010.csv'    192 records     C  

total:                                      2942

Pt II
The Global Gridded Dataset records the coordinates for each water mass according
to WAOMASK prn files. These coordinates are from using LeonArcMask.prn in grid_d18o_7_9_15
It makes a dataframe for each water mass, plots d18O-S relationship, and jackknife slope analysis.

Written by Leon Yin for NASA GISS NYCRI 2015
Please consult NASA GISS or the author (leon.yin@nyu.edu) for questions

v - 0.1.2
"""
# Files to be read/written in the parent directory
filein      = 'dbO18.txt'       # this file is downloaded from the d18o website.
fileout     = 'dbO18_2015.txt'  # output of the entire DB as a fixed-width txt/ASCII format
csvout      = 'dbO18_2015.csv'  # replace 2015 with specifics!

# an example of a list of properly formatted csv files to add to the DB.
toAdd       = ('Strauss2015.csv','Cases2012.csv','COX2010_2001.csv',
              'Conroy2012.csv','Tiwari2013.csv','BROWN2014.csv',
              'GeoRijkenberg2011.csv','GeoJenkins2010.csv')    
# Lon and Lat coordinates (as integers) derived from LeonArcMask.prn in grid_d18o_7_9_15
beauLonIN   =   'beauLon.txt'
beauLatIN   =   'beauLat.txt'
canLonIN    =   'canLon.txt'
canLatIN    =   'canLat.txt'
sibLonIN    =   'sibLon.txt'
sibLatIN    =   'sibLat.txt'
eurLonIN    =   'eurLon.txt'
eurLatIN    =   'eurLat.txt'
karLonIN    =   'karLon.txt'
karLatIN    =   'karLat.txt'
lapLonIN    =   'lapLon.txt'
lapLatIN    =   'lapLat.txt'
lomLonIN    =   'lomLon.txt'
lomLatIN    =   'lomLat.txt'
norLonIN    =   'norLon.txt'
norLatIN    =   'norLat.txt'
barLonIN    =   'barLon.txt'
barLatIN    =   'barLat.txt'
chukLonIN   =   'chukLon.txt'
chukLatIN   =   'chukLat.txt'

# Logic switches -- use these guys to turn functions on-and-off.
outlierELIM     =  -1
writeOutIce     =  -1
save            =  -1
updates         =  -1

"""****************************************************************************
These three functions are from Stackoverflow posts
****************************************************************************"""
def _get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness   = (50 + np.random.rand() * 10) / 100.
        saturation  = (90 + np.random.rand() * 10) / 100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
    return colors

def mad_based_outlier(points, thresh=3.5):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

def percentile_based_outlier(data, threshold=95):
    diff = (100 - threshold) / 2.0
    minval, maxval = np.percentile(data, [diff, 100 - diff])
    return (data < minval) | (data > maxval)

"""****************************************************************************
Use this function to convert DF to to CVS, and write it out as a csv.
****************************************************************************"""
def convert2CSV(df,csvout):
    df2csv = 'dbo18_2015.csv'
    toCSV=pd.DataFrame.to_csv(df2csv,sep=',',line_terminator='\n')
    csvout = open(csvout, 'w')      # opens the file to write into
    csvout.write(toCSV)             # writes df to csv... 
    print "Database transferred to csv..."
    csvout.close()

"""****************************************************************************
Use this function to write the database into a fixed-width text file.
****************************************************************************"""
def writeOut(df):
    df=df.fillna('')                #replace NaN values with spaces
    temp = 'temp.txt'
    print("******************************************************************")        
    dbOut = open(temp, 'w')
    v = df.values.T.tolist()        # The DF must be converted to a list first.
    for i in range(0,dfRows):       # Right each row into the file.
        dbOut.write(( \
        '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.1f}{:>15}{:>60}'\
        .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
        v[9][i],v[10][i]) ))
        dbOut.write("\n")
    dbOut.close
    # swap -99.90 with -99.9
    with open(fileout, "wt") as fout:
        with open(temp, "rt") as fin:
            for line in fin:
                fout.write(line.replace('-99.90', ' -99.9'))
    print "Database transferred to txt..."
    fout.close
    fin.close
    os.remove(temp)                 # Remove the temporary file.
"""****************************************************************************
Plots all the measurements according to the d18O value
****************************************************************************"""
def pltVals(df):
    plt.figure(figsize=(7,7))  
    map = Basemap(projection='npstere',boundinglat=68,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,4.)
    meridians = np.arange(-180.,181.,10.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    map.drawmeridians(meridians,labels=[True,False,False,False])
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    x1 = df.lon.values.T.tolist()
    y1 = df.lat.values.T.tolist()
    z1 = df.d18o.values.T.tolist()
    x, y = map(x1, y1)
    map.scatter(x,y,c=z1,marker='o',s=16,linewidth=.2)
    plt.title('Arctic Observed $\delta^{18}$O',fontsize=14)
    plt.clim(-5.6,2.8)
    cbar = plt.colorbar(orientation='vertical')
    cbar.ax.set_xlabel('$\delta^{18}$O')
    if(save == 1):
        plt.savefig('valDistro', dpi = 300)
    plt.show()

"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following section contains search functions and several examples.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
"""****************************************************************************
Use this function to search within one field (column) several values.
Make sure args are: search(dataframe, string of field, list of search results)
ex:
    fromCol     = 'month'           # What column are you interested in?(ln:72)
    values      = (12,1,2)          # Put what you're looking for here!
    search(df,fromCol,values)
****************************************************************************"""
def search(df,fromCol,values):
     # How many automations/searches are required?
    uniqueVal   = len(values)      
    # Collect values from specified col into new DF.    
    querry = df.loc[df[fromCol] == values[0]]
    for i in range(1,uniqueVal):
        querry=querry.append(df.loc[df[fromCol] == values[i]])
    print"found",len(querry),"instances of",fromCol,"with values of:",values
    
    # Book keeping...
    refList = querry.ref.unique()   # List of search-relevant references.
    refSize = refList.size
    print"Found",refSize,"references containing search results listed below...\n",refList
    
    # plot the results! Below is optional
    plt.figure(figsize=(9,9))
    map = Basemap(projection='mill', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    x1 = querry.lon.values.T.tolist()
    y1 = querry.lat.values.T.tolist() 
    x, y = map(x1, y1)
    
    map.scatter(x, y, marker='o',color='k',s=9)
    
    plt.title(r'Winter Distribution of Observed $\delta^{18}$O',size=16)
    plt.show
"""****************************************************************************
Use this function to section off the observed datapoints seasonally.
This is an extension of the above function.
****************************************************************************"""
def pltSeason(df):
    # Seasonal split!
    fromCol     = 'month'           # What column/field are we interested in?
    winter      = (12,1,2)          # Winter Months
    spring      = (3,4,5)           # Spring Months
    summer      = (6,7,8)           # Summer Months
    fall        = (9,10,11)         # Fall Months (my favorite...)
    uniqueVal   = len(winter)       # This is for automation'ssake
    
    #   Section off winter months  
    querryW = df.loc[df[fromCol] == winter[0]]
    for i in range(1,uniqueVal):
        querryW=querryW.append(df.loc[df[fromCol] == winter[i]])
    refListW = querryW.ref.unique()
    refSize = refListW.size
    print"found",len(querryW),"instances of",fromCol,"with values of:",winter,\
    "from",refSize,"unique references"
    
    # Section off spring months  
    querryS = df.loc[df[fromCol] == spring[0]]
    for i in range(1,uniqueVal):
        querryS=querryS.append(df.loc[df[fromCol] == spring[i]])
    refListS = querryS.ref.unique()
    refSize = refListS.size
    print"found",len(querryS),"instances of",fromCol,"with values of:",spring,\
    "from",refSize,"unique references"
    
    # Section off summer months
    querrySM = df.loc[df[fromCol] == summer[0]]
    for i in range(1,uniqueVal):
        querrySM=querrySM.append(df.loc[df[fromCol] == summer[i]])
    refListSM = querrySM.ref.unique()
    refSize = refListSM.size    
    print"found",len(querrySM),"instances of",fromCol,"with values of:",summer,\
    "from",refSize,"unique references"
    
    # section off fall months
    querryF = df.loc[df[fromCol] == fall[0]]
    for i in range(1,uniqueVal):
        querryF=querryF.append(df.loc[df[fromCol] == fall[i]])
    refListF = querryF.ref.unique()
    refSize = refListF.size
    print"found",len(querryF),"instances of",fromCol,"with values of:",fall,\
    "from",refSize,"unique references"   
    
    # plot the dictributions...
    plt.figure(figsize=(7,7))
    # choose one of the projections--
    # global projections
#    map = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90,
#    llcrnrlon=-180,urcrnrlon=180,resolution='c')
    # Arctic projection    
    map = Basemap(projection='npstere',boundinglat=68,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,4.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    meridians = np.arange(-180.,181.,10.)
    map.drawmeridians(meridians,labels=[True,False,False,False])
    
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    #establish lon and Lat per season...    
    x1 = querryS.lon.values.T.tolist()
    y1 = querryS.lat.values.T.tolist() 
    xS, yS = map(x1, y1)
    map.scatter(xS, yS, marker='o',color='g',label='Mar Apr May',s=16)

    x1 = querrySM.lon.values.T.tolist()
    y1 = querrySM.lat.values.T.tolist() 
    xSM, ySM = map(x1, y1)
    map.scatter(xSM, ySM, marker='o',color='r',label='Jun Jul Aug',s=16)

    x1 = querryF.lon.values.T.tolist()
    y1 = querryF.lat.values.T.tolist() 
    xF, yF = map(x1, y1)
    map.scatter(xF, yF, marker='o',color='orange',label='Sept Oct Nov',s=16)

    x1 = querryW.lon.values.T.tolist()
    y1 = querryW.lat.values.T.tolist() 
    xW, yW = map(x1, y1)
    map.scatter(xW, yW, marker='o',color='k',label='Dec Jan Feb',s=16)
    # stuff so the for legend is saved w/ the fig!
    art = []
    lgd = plt.legend(loc='upper center',scatterpoints=1,ncol=4,
         bbox_to_anchor=(0.5, -0.05),prop={'size':12})
    art.append(lgd)      
    
    plt.title(r'Tempoal Distribution of Observed $\delta^{18}$O',size=16)
    plt.show
    if (save ==1):
        plt.savefig('seasonality.png', bbox_inches="tight",
                    additional_artists=art,format='png',dpi = 300)
    # Bookkeeping...    
    seasonDF = (querryF,querryW,querryS,querrySM)
"""****************************************************************************
Inter-annual distribution...
Something to be explored -LY
****************************************************************************"""
def pltAnnual(df):
    plt.figure(figsize=(9,7))  
    map = Basemap(projection='npstere',boundinglat=68,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,4.)
    meridians = np.arange(-180.,181.,10.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    map.drawmeridians(meridians,labels=[True,False,False,False])
#    map = Basemap(projection='mill', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    df=df[df.year!=-999]
    annie = df.year.unique()
    annie.sort()
    annieHall = len(annie)
    
    colorSwatch=_get_colors(annieHall)
    mark = ['o','v','s','d']        # pool of shapes for scatter    
    
    for i in range(annieHall):
#        print annie[i]
        annualDF = df.loc[df.year==annie[i]]
        annualDF = annualDF[annualDF.depth<=10]
        

        x1 = annualDF.lon.values.T.tolist()
        y1 = annualDF.lat.values.T.tolist()
        x, y = map(x1, y1)
        
        map.scatter(x, y, marker=mark[i%4],color=colorSwatch[i],
                    label=annie[i],s=8)
    art = []
    lgd = plt.legend(loc='upper center',scatterpoints=1,ncol=6,
         bbox_to_anchor=(0.5, -0.05),prop={'size':10})
    plt.title('Interannual Distribution of $\delta^{18}$O',fontsize=14)
    if(save == 1):
        plt.savefig('interannual', dpi = 300)
    plt.show()    
"""****************************************************************************
Use this to generate a scatter plot for all points by unique reference.
****************************************************************************"""
def pltRef(df):
    refList = df.ref.unique()       # creates list of unique references from DB
    refSize = refList.size          # how large is the list?
    
    colorSwatch=_get_colors(refSize)# get enough colors for each point. 
    mark = ['o','v','s','d']        # pool of shapes for scatter
    
    plt.figure(figsize=(10,10))     # initiate plot w/ dimensions.
    map = Basemap(projection='robin', lat_0=0, lon_0=0)
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    for i in range(0,refSize):
        #breaks dataframe into subdataframes of unique references    
        uniqueRef=refList[i]
        dfSub=df.loc[df['ref']==uniqueRef]
        x1 = dfSub.lon.values.T.tolist()    # convert longitude values to list
        y1 = dfSub.lat.values.T.tolist()    # convert latitude valus to list
        x, y = map(x1, y1)                  # convert the lon and lat into map
        map.scatter(x, y, marker=mark[i%4],color=colorSwatch[i],
                    label=str(refList[i]),s=3)
    art = []
    lgd = plt.legend(scatterpoints=1,bbox_to_anchor=(1.05, 1), 
                     loc=2, borderaxespad=0.0,prop={'size':8})
#    art.append(lgd)     
   
    plt.title(r'Global Observed $\delta18O$ data')
    plt.show
    if (save == 1):
        plt.savefig('observedD18O.png', bbox_inches="tight",
                    additional_artists=art,format='png',dpi = 300)

"""****************************************************************************
Use this to generate a scatter plot for all arctic points by unique reference.
****************************************************************************""" 
def pltArcRef(df):
    # Define the Arctic...
    newArc = df.loc[df.lat>=65.8]
    refList = newArc.ref.unique()
    refList.sort()
    refSize = refList.size

    # Various Schemes of getting colors
    colorSwatch=_get_colors(refSize)
    
    mark = ['o','v','s','d']
    
    plt.figure(figsize=(10,10))
    map = Basemap(projection='npstere',boundinglat=66,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,2.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    meridians = np.arange(-180.,181.,5.)
    map.drawmeridians(meridians,labels=[True,False,False,False])
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    # Make a plot according to reference
    for i in range(0,refSize):
        #breaks dataframe into subdataframes of unique references    
        uniqueRef=refList[i]
        dfSub=newArc.loc[newArc['ref']==uniqueRef]
        x1 = dfSub.lon.values.T.tolist()    # convert longitude values to list
        y1 = dfSub.lat.values.T.tolist()    # convert latitude valus to list
        x, y = map(x1, y1)
        map.scatter(x, y, marker=mark[i%4],color=colorSwatch[i],label=str(uniqueRef),s=7)
    art = []
    lgd = plt.legend(scatterpoints=1,bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0,prop={'size':9})
    art.append(lgd)    
    plt.title(r'Spatial Distribution of Observed Arctic $\delta^{18}$O')
    if (save == 1):
        plt.savefig('arcDist.png', bbox_inches="tight" ,additional_artists=art,format='png', dpi=200)
#       plt.savefig('arcDist.eps', format='eps', dpi=1000) #for vectors
    plt.show
"""****************************************************************************
Separates DB into two DF's one for melt and one for frozen NOTE: fix this...
****************************************************************************"""
def pltMelt(df):
    print("******************************************************************")        
    fromCol = 'month'           # What column/field are we interested in?
    ice     = (1,2,3,4,10,11,12)    # Winter Months
    water   = (5,6,7,8,9)           # Spring Months may-sept
    iLen    = len(ice)
    wLen    = len(water)            # This is for automation'ssake
    
    #   Section off frozen months  
    iceDF = df.loc[df[fromCol] == ice[0]]
    for i in range(1,iLen):
        iceDF=iceDF.append(df.loc[df[fromCol] == ice[i]])
    refListIce = iceDF.ref.unique()
    refSize = refListIce.size
    print"found",len(iceDF),"instances of",fromCol,"with values of:",ice,\
    "from",refSize,"unique references"
    
    # Section off melt months  
    waterDF = df.loc[df[fromCol] == water[0]]
    for i in range(1,wLen):
        waterDF=waterDF.append(df.loc[df[fromCol] == water[i]])
    refListWater = waterDF.ref.unique()
    refSize = refListWater.size
    print"found",len(waterDF),"instances of",fromCol,"with values of:",water,\
    "from",refSize,"unique references"
   
   # write the DB for ice and melt into .dat files for global gridded dataset.
    if (writeOutIce == 1):
        # Ice season
        iceDf=iceDF.fillna('')
        temp = 'temp.txt'
        fileout = 'ice.dat'
        print("******************************************************************")        
        dbOut = open(temp, 'w')
        v = iceDf.values.T.tolist()        # The DF must be converted to a list first.
        for i in range(0,len(iceDF)):       # Right each row into the file.
            dbOut.write(( \
            '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.1f}{:>15}{:>60}'\
            .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
            v[9][i],v[10][i]) ))
            dbOut.write("\n")
        dbOut.close
        # swap -99.90 with -99.9
        with open(fileout, "wt") as fout:
            with open(temp, "rt") as fin:
                for line in fin:
                    fout.write(line.replace('-99.90', ' -99.9'))
        print "Database transferred to txt..."
        fout.close
        fin.close
        os.remove(temp) 
        
        # Time for meltseason
        waterDf=waterDF.fillna('')
    
        fileout = 'melt.dat'
        dbOut = open(temp, 'w')
        v = waterDf.values.T.tolist()        # The DF must be converted to a list first.
        for i in range(0,len(waterDF)):       # Right each row into the file.
            dbOut.write(( \
            '{:7.2f}{:>6.2f}{:>2.0f}{:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:6.1f}{:>15}{:>60}'\
            .format(v[0][i],v[1][i],v[2][i],v[3][i],v[4][i],v[5][i],v[6][i],v[7][i],v[8][i],\
            v[9][i],v[10][i]) ))
            dbOut.write("\n")
        dbOut.close
        # swap -99.90 with -99.9
        with open(fileout, "wt") as fout:
            with open(temp, "rt") as fin:
                for line in fin:
                    fout.write(line.replace('-99.90', ' -99.9'))
        print "Database transferred to txt..."
        fout.close
        fin.close
        os.remove(temp) 

    #surface only
    surfDep     =  10
    iceDF       =  iceDF.loc[iceDF.depth<=surfDep]
    waterDF     =  waterDF.loc[waterDF.depth<=surfDep]
    
    # Plot 1
    plt.figure(figsize=(7,7))  
    map = Basemap(projection='npstere',boundinglat=68,lon_0=330,resolution='l')
    parallels = np.arange(-90.,90.,4.)
    meridians = np.arange(-180.,181.,10.)
    map.drawparallels(parallels,labels=[False,False,False,True])
    map.drawmeridians(meridians,labels=[True,False,False,False])
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    x1 = iceDF.lon.values.T.tolist()
    y1 = iceDF.lat.values.T.tolist()
    z1 = iceDF.d18o.values.T.tolist()
    x2, y2 = map(x1, y1)
    map.scatter(x2,y2,c=z1,marker='o',s=16,linewidth=.2)
    plt.title('Arctic Frozen Season Observed $\delta^{18}$O',fontsize=14)
    plt.clim(-5.6,2.8)
    cbar = plt.colorbar(orientation='vertical')
    cbar.ax.set_xlabel('$\delta^{18}$O')
    if(save == 1):
        plt.savefig('frozenObs', dpi = 300)
    plt.show()
    
    # Plot 2
    plt.figure(figsize=(7,7))  
    map = Basemap(projection='npstere',boundinglat=68,lon_0=330,resolution='l')
    map.drawparallels(parallels,labels=[False,False,False,True])
    map.drawmeridians(meridians,labels=[True,False,False,False])
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()
    
    x1 = waterDF.lon.values.T.tolist()
    y1 = waterDF.lat.values.T.tolist()
    z1 = waterDF.d18o.values.T.tolist()
    x2, y2 = map(x1, y1)
    map.scatter(x2,y2,c=z1,marker='o',s=16,linewidth=.2)
    plt.title('Arctic Melt Season Observed $\delta^{18}$O',fontsize=14)
    plt.clim(-5.6,2.8)
    cbar = plt.colorbar(orientation='vertical')
    cbar.ax.set_xlabel('$\delta^{18}$O')
    if(save == 1):
        plt.savefig('meltObs', dpi = 300)
    plt.show()
"""****************************************************************************
Monthly averages of salinty calclulated from WOA and comapred to sal and d18O observations
****************************************************************************"""
def pltSal():
    if(os.path.isfile(WOAF)):
        WOA    =   Dataset(WOAF, mode='r')
    else:
        print chukLonIN,"is either unreadable or DNE."
        raise SystemExit
    # Looking at the netCDF surface lvl 0 where: sal[month][depth][Lat][Lon]
    sal =WOA.variables['salinity'][:,0,:,:]
    
    beauSalAv   =   (np.zeros(12))
    canSalAv    =   (np.zeros(12))
    sibSalAv    =   (np.zeros(12))
    eurSalAv    =   (np.zeros(12))
    lapSalAv    =   (np.zeros(12))
    lomSalAv    =   (np.zeros(12))
    norSalAv    =   (np.zeros(12))
    karSalAv    =   (np.zeros(12))
    barSalAv    =   (np.zeros(12))
    chukSalAv   =   (np.zeros(12))
    
    print"\n*******************************************************************"
    print"Calculating ave sainity of Arctic reigons from Wolrd Ocean Atlas..."
    print"*******************************************************************"

    for month in range(0,12):
        beauSal=np.zeros(beauPts)
        for i in range(0,beauPts):
            beauSal[i]=sal[month][beauLat[i]][beauLon[i]]
        #mask the nan values and calculated the average of sal
        beauSal=np.ma.masked_invalid(beauSal)
        beauSalAv[month]=round(np.ma.mean(beauSal),2)
        print"the average salinity in the beaufort sea during month",month+1,\
        "is",beauSalAv[month]
        
        canSal=np.zeros(canPts)
        for i in range(0,canPts):
            canSal[i]=sal[month][canLat[i]][canLon[i]]
        #mask the nan values and calculated the average of sal
        canSal=np.ma.masked_invalid(canSal)
        canSalAv[month]=round(np.ma.mean(canSal),2)
#        print"the average salinity in the Canadian basin during month",month+1,\
#        "is",canSalAv[month]
        
        sibSal=np.zeros(sibPts)
        for i in range(0,sibPts):
            sibSal[i]=sal[month][sibLat[i]][sibLon[i]]
        #mask the nan values and calculated the average of sal
        sibSal=np.ma.masked_invalid(sibSal)
        sibSalAv[month]=round(np.ma.mean(sibSal),2)
#        print"the average salinity in the sibadian basin during month",month+1,\
#        "is",sibSalAv[month]
        
        eurSal=np.zeros(eurPts)
        for i in range(0,eurPts):
            eurSal[i]=sal[month][eurLat[i]][eurLon[i]]
        #mask the nan values and calculated the average of sal
        eurSal=np.ma.masked_invalid(eurSal)
        eurSalAv[month]=round(np.ma.mean(eurSal),2)
#        print"the average salinity in the euraadian basin during month",month+1,\
#        "is",euraSalAv[month]

        lapSal=np.zeros(lapPts)
        for i in range(0,lapPts):
            lapSal[i]=sal[month][lapLat[i]][lapLon[i]]
        #mask the nan values and calculated the average of sal
        lapSal=np.ma.masked_invalid(lapSal)
        lapSalAv[month]=round(np.ma.mean(lapSal),2)
        
        lomSal=np.zeros(lomPts)
        for i in range(0,lomPts):
            lomSal[i]=sal[month][lomLat[i]][lomLon[i]]
        #mask the nan values and calculated the average of sal
        lomSal=np.ma.masked_invalid(lomSal)
        lomSalAv[month]=round(np.ma.mean(lomSal),2)
#        print"the average salinity in the lomadian basin during month",month+1,\
#        "is",lomSalAv[month]
        
        norSal=np.zeros(norPts)
        for i in range(0,norPts):
            norSal[i]=sal[month][norLat[i]][norLon[i]]
        #mask the nan values and calculated the average of sal
        norSal=np.ma.masked_invalid(norSal)
        norSalAv[month]=round(np.ma.mean(norSal),2)
#        print"the average salinity in the noradian basin during month",month+1,\
#        "is",norSalAv[month]
        
        karSal=np.zeros(karPts)
        for i in range(0,karPts):
            karSal[i]=sal[month][karLat[i]][karLon[i]]
        #mask the nan values and calculated the average of sal
        karSal=np.ma.masked_invalid(karSal)
        karSalAv[month]=round(np.ma.mean(karSal),2)
#        print"the average salinity in the karaadian basin during month",month+1,\
#        "is",karaSalAv[month]
        
        barSal=np.zeros(barPts)
        for i in range(0,barPts):
            barSal[i]=sal[month][barLat[i]][barLon[i]]
        #mask the nan values and calculated the average of sal
        barSal=np.ma.masked_invalid(barSal)
        barSalAv[month]=round(np.ma.mean(barSal),2)
        
        chukSal=np.zeros(chukPts)
        for i in range(0,chukPts):
            chukSal[i]=sal[month][chukLat[i]][chukLon[i]]
        #mask the nan values and calculated the average of sal
        chukSal=np.ma.masked_invalid(chukSal)
        chukSalAv[month]=round(np.ma.mean(chukSal),2)
        
        
    y1List      =   (eurSalAv,beauSalAv,canSalAv,sibSalAv,lapSalAv,lomSalAv,
                     norSalAv,karSalAv,chukSalAv,barSalAv)

    fig = plt.figure(figsize=(10,6))     
    ax1 = fig.add_subplot(111)                              
    
    ax2=ax1.twinx()
    ax1.set_ylabel('salinity (psu)',size=14)
    ax1.set_xlabel('month',size=14)
    ax2.set_ylabel('$\delta^{18O}$',size=14)
    
    ax1.set_xlim((0,13))
    ax1.set_ylim(10,36)
    ax2.set_xlim((0,13))
    ax2.set_ylim((-20,3))
    
    for i in range(0,reigons):
        # points on the surface of the first (surfDep) meters
        arcSurf=DFList[i].loc[DFList[i].depth<surfDep]        
        # remove the -99.90 (skip values) from the Db salinity.
        arcSurf=arcSurf[arcSurf.sal!=skip]
        arcSurf=arcSurf[arcSurf.d18o!=skip]
        arcSurf=arcSurf[arcSurf.month!=13]
        
        # plotting parameters
        ySal    =   arcSurf.sal
        xArc    =   arcSurf.month # month coordiantes per reigon
        y1      =   y1List[i]
        x1      =   np.array([1,2,3,4,5,6,7,8,9,10,11,12])
        xMonSp  =   np.linspace(x1.min(),x1.max(),300)
        ySalSp  =   spline(x1,y1, xMonSp)
        yD18O   =   arcSurf.d18o
        # plot the monthly average salinity
        ax1.plot(xMonSp,ySalSp ,linewidth=linWid,color=colorSwatch[i],zorder=1,alpha=.8)
        # plot the observed salinity
        ax1.scatter(xArc ,ySal ,marker=markSal[i%4],s=sizeSal[i%4],
                    color=colorSwatch[i],label=lblSal[i],alpha=.7,zorder=10)
        # plot the observed d18O
        ax2.scatter(xArc ,yD18O,marker=markD18[i%4],s=sizeD18[i%4],
                    color=colorSwatch[i],label=lblD18[i],alpha=.7,zorder=5)
    
    art = []
    lg2 = ax2.legend(scatterpoints=1,loc='center left',markerscale=2,
                     bbox_to_anchor=(1.1,.5),ncol=1,prop={'size':12})
    lgd = ax1.legend(scatterpoints=1,loc='upper center',markerscale=2,
                     bbox_to_anchor=(0.5, -0.08),ncol=2,prop={'size':12})
    art.append(lgd) 
    art.append(lg2)
#    plt.title('Temporal Distribution of Salinity in Arctic Surface')

    plt.title('Monthly Variability of Salinity in Arcitic Ocean',size=16)
    if (save == 1):
        plt.savefig('salDistro.png',additional_artists=art,bbox_inches="tight", dpi=300)
    plt.show
"""****************************************************************************
Scatterplot and linear regression of reigonal d18o-sal relationships
****************************************************************************"""
def pltD18O_sal():
    fig = plt.figure(figsize=(10,6))     
    ax = fig.add_subplot(111)    
    ax.set_ylabel('$\delta^{18}$O (permilli)',size=14)
    ax.set_xlabel('Salinity (psu)',size=14)
    ax.set_xlim((1,36))
    ax.set_ylim(-20,1)   
    colorSwatch=_get_colors(reigons)
    print"\n*******************************************************************"
    print"The d18O-S slope and y-intercept of Arctic reigons..."
    print"*******************************************************************"
    order = (1,2,0,4,7,5,8,9,3,6)

    for i in range(len(order)):
        # points on the surface of the first (surfDep) meters
        arcSurf=DFList[order[i]].loc[DFList[order[i]].depth<surfDep]
    
        # remove the -99.90 (skip values) from the Db salinity.
        arcSurf=arcSurf[arcSurf.sal!=skip]
        arcSurf=arcSurf[arcSurf.d18o!=skip]
        arcSurf=arcSurf[arcSurf.month!=13]
    
        # Calculate slope and y-int
        x = arcSurf.sal  #sal
        y = arcSurf.d18o #d18O        
    
        par = np.polyfit(x, y, 1, full=True)
        m=par[0][0]     # SLOPE
        b=par[0][1]     # Y-INTERCEPT
        # Put the average WOA salinity in smooth plot-able format.
        yp = np.polyval([m,b],x)
        #plotting
        ax.plot(x,yp,color=colorSwatch[order[i]],linewidth=4.5,zorder=1)   
        print reigonName[order[i]],len(arcSurf)," points\nY-intercept:",\
                    round(b,2),'d18O-S Slope:',round(m,2)
        ax.scatter(x ,y ,marker=markSal[order[i]%4],s=sizeSal[order[i]%4],\
                   color=colorSwatch[order[i]],label=reigonName[order[i]],
                    alpha=.7,zorder=10)
        equation=str(round(m,2))+"x "+str(round(b,1))  
        
#        ax.text(-0.21, 0.75-(i*.055),equation,horizontalalignment='center',
#                color=colorSwatch[order[i]], verticalalignment='center',
#                 transform=ax.transAxes,size=12)
        ax.text(.92, 0.45-(i*.045),equation,horizontalalignment='center',
                color=colorSwatch[order[i]], verticalalignment='center',
                 transform=ax.transAxes,size=10)
    art = []
#    lgd = ax.legend(scatterpoints=1,markerscale=2,loc='center left', 
#                    bbox_to_anchor=(-0.55,.5),ncol=1,prop={'size':12},frameon=False)
    lgd = ax.legend(scatterpoints=1,markerscale=2, loc='center',
                    bbox_to_anchor=(.75,.25),ncol=1,prop={'size':10},frameon=False)
    art.append(lgd) 
    plt.title('$\delta^{18}$O-Salinity Relationships of Arctic Water Masses',size=16)
    if (save == 1):
        plt.savefig('slopes.png',additional_artists=art,bbox_inches="tight", dpi=300)
    plt.show
    print 'Note that the y-intercept is d18 of freshwater inputs'
"""****************************************************************************
Jackknife slope analysis with outlier flagging
****************************************************************************"""
def jack():
    import seaborn as sns
    # Initialize some empty lists...        
    outDF = pd.DataFrame(columns=col)
    indieList = [[] for i in range(reigons)]    
    jack = [[] for i in range(reigons)]
    outliers = [[] for i in range(reigons)]
    for i in range(reigons):
        arcSurf=DFList[i].loc[DFList[i].depth<surfDep]
        
        # remove the -99.90 (skip values) from the Db salinity.
        arcSurf=arcSurf[arcSurf.sal!=skip]
        arcSurf=arcSurf[arcSurf.d18o!=skip]
        arcSurf=arcSurf[arcSurf.month!=13]
        
        reigLen = len(arcSurf)
        jack[i] = np.zeros(reigLen)
        indieList[i] = arcSurf.index.values.T.tolist() 
        for j in range(reigLen): #jacknife
            temp = arcSurf[arcSurf.index!=indieList[i][j]] #remove on ID at a time
            x3 = temp.sal
            y3 = temp.d18o
            par2 = np.polyfit(x3, y3, 1, full=True)
            m1=par2[0][0]           # SLOPE
            b1=par2[0][1]           # Y-INTERCEPT
            jack[i][j]=round(m1,3)  # Record jackknife with 3-decimal percision
#            print jack[i][j]
            
        # Density distribution of jack knife and outliers
        fig, axes = plt.subplots(nrows=2)
        for ax, func in zip(axes, [percentile_based_outlier, mad_based_outlier]):
            sns.distplot(jack[i], ax=ax, rug=True, hist=False)
            outliers[i] = jack[i][func(jack[i])]
            ax.plot(outliers[i], np.zeros_like(outliers[i]), 'ro', clip_on=False)
            fig.suptitle('Jack Knife Analysis of '+reigonName[i]+' using {} points'.format(len(jack[i])), size=14)
        kwargs = dict(y=0.95, x=0.05, ha='left', va='top')
        axes[0].set_title('Percentile-based Outliers', **kwargs)
        axes[1].set_title('MAD-based Outliers', **kwargs)
        if (save == 1):
            plt.savefig(reigonName[i]+'outliers.png') 
        plt.show
    # find outliers
    toDestroy = [0,1,2,3,4,5,6,7,8,9]
#    toDestroy = [0,4,6,8]
    for i in range(len(toDestroy)):
        print reigonName[toDestroy[i]]+"~"
        # points on the surface of the first (surfDep) meters
        arcSurf=DFList[toDestroy[i]].loc[DFList[toDestroy[i]].depth<surfDep]
        
        # remove the -99.90 (skip values) from the Db salinity.
        arcSurf=arcSurf[arcSurf.sal!=skip]
        arcSurf=arcSurf[arcSurf.d18o!=skip]
        arcSurf=arcSurf[arcSurf.month!=13]
        print ( '{:>6}{:>5}{:>7}{:>6} {:>2} {:>4}{:>5}{:>6}{:>6}{:>6}{:>12}{:>40}'\
            .format('id','sig','lat','lon','m','y','dep','temp','sal','d18o','notes','ref'))
        for k in range(len(outliers[toDestroy[i]])):
            qq=np.where(jack[toDestroy[i]]==outliers[toDestroy[i]][k])
            outRow = arcSurf.loc[arcSurf.index==indieList[toDestroy[i]][qq[0][0]]]
            v = outRow.values.T.tolist()
#            print indieList[toDestroy[i]][qq[0][0]],outliers[toDestroy[i]][k]
            print ( '{:>6.0f}{:>5.2f}{:>7.2f}{:>6.2f} {:>2.0f} {:>4.0f}{:>5.0f}{:6.2f}{:6.2f}{:6.2f}{:>12}{:>40}'\
            .format(indieList[toDestroy[i]][qq[0][0]],outliers[toDestroy[i]][k],v[0][0],v[1][0],v[2][0],v[3][0],v[4][0],v[5][0],v[6][0],v[7][0],\
            'OUTLIER!',v[10][0]) )
            outDF=outDF.append(outRow)
        killList=outDF.index.unique() #print this out if a list of indicies of outliers is desired...
    # PLot outliers
    plt.figure(figsize=(7,7)) 
    parallels = np.arange(-90.,91.,4)
    meridians = np.arange(-180.,181.,10.)
    map = Basemap(projection='npstere',boundinglat=68,lon_0=330,resolution='l')
    map.drawparallels(parallels,labels=[False,False,False,True])
    map.drawmeridians(meridians,labels=[True,False,False,False])
    map.drawcoastlines()
    map.fillcontinents()
    map.drawmapboundary()

    x1 = df.lon.values.T.tolist()
    y1 = df.lat.values.T.tolist()
    z1 = df.d18o.values.T.tolist()
    x, y = map(x1, y1)
    map.scatter(x,y,c=z1,marker='o',s=16,linewidth=.2)
    plt.clim(-5.6,2.8)
    x2 = outDF.lon.values.T.tolist()
    y2 = outDF.lat.values.T.tolist()
    z2 = outDF.d18o.values.T.tolist()
    x3, y3 = map(x2, y2)
    # Outliers location shown as larger square
    map.scatter(x3,y3,c=z2,marker='s',s=40,linewidth=.6)
    plt.title('Arctic Observed $\delta^{18}$O Outliers',fontsize=14)
    plt.clim(-5.6,2.8)
    cbar = plt.colorbar(orientation='vertical')
    cbar.ax.set_xlabel('$\delta^{18}$O')
    if(save == 1):
        plt.savefig('arcticOutliers.png',dpi = 300)
    plt.show()


"""~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PT I This is the main function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"""
# Global Seawater d18O database format.
col=\
  ['lon','lat','month','year','depth','temp', 'sal', 'd18o','dD','notes','ref'] 
fWid=( 7,    6,      2,     4,      5,     6,     6,      6,   6,    15,   60)

# Read in Database into a Pandas Dataframe data structure
if(os.path.isfile(filein)):
    df=pd.read_fwf(filein,widths=fWid,header=None,names=col,skip_blank_lines=True)
    ogLen     = len(df.index-1) # number of records in DF from filein
    print ogLen,"records from Global Seawater Oxygen-18 Database imported to dataframe df."
else:
    print filein,"is either unreadable or DNE."
    raise SystemExit(0)

# some variables...
dfRows    = len(df.index-1) # number of records in current DF
ogLen     = len(df.index-1) # number of records in DF from filein
files2Add = len(toAdd)      # number of sources to add

"""****************************************************************************
Automates reading each CSV, and appending to existing DB.
To troubleshoot adjust the range and double-check 'csv' in the shell.
****************************************************************************"""
files2Add = len(toAdd)      # number of sources to add
for i in range(files2Add):
    if(os.path.isfile(toAdd[i])==False):
        print toAdd[i],"is either unreadable or DNE."
        raise SystemExit(0)
print("******************************************************************")        
print"Integrating CSV into Global D18O Seawater Database..."    
for i in range(0,files2Add):
#    print "opening",toAdd[i],""
    csv2DF=pd.read_csv(toAdd[i],sep=',',skiprows=1,names=col,na_filter=True)
    #              File Name,Delimiter, skip header
    dfRows= len(df.index-1)
    df = df.append(csv2DF)  # add contents of new csv into DB.
    dfRows= len(df.index-1) # update the number of records in DB.
print dfRows-ogLen,"records added to the database."
print("******************************************************************")        

"""****************************************************************************
Seek vengence on outliers from variable high latitutes.
ID's determined using Jack knife statstical analysis in osSlopes.py...
Known troublsome data...

id     lon  lat  month  year  depth   temp    sal  d18o    dD                ref 
18067    3   80      6  1984      0  -1.20  20.00  0.14 -99.9  Ostlund and Grall (1993)  

replace hard coded idOdeath with killList from osSlopes.py for easy(but mindless)
removal of outliers.
****************************************************************************"""
if (outlierELIM == 1):
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
    print "Destroying outliers..."
    print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")

    idOdeath    = [18067,519,825,613,18383]
    lenOdeath   = len(idOdeath)
    for toKill in range(lenOdeath):     # filters DF by natural index.
        df = df[df.index != idOdeath[toKill]] 
    print lenOdeath,"outliers ELIMINATED!"
    dfRows    = len(df.index-1)
    print("\m/ -___- \m/")
"""****************************************************************************
PT II
****************************************************************************"""
if(os.path.isfile(beauLonIN)):
    beauLonF    =   open(beauLonIN,"r")
else:
    print beauLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(beauLatIN)):
    beauLatF    =   open(beauLatIN,"r")
else:
    print beauLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(canLonIN)):
    canLonF    =   open(canLonIN,"r")
else:
    print canLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(canLatIN)):
    canLatF    =   open(canLatIN,"r")
else:
    print canLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(sibLonIN)):
    sibLonF    =   open(sibLonIN,"r")
else:
    print sibLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(sibLatIN)):
    sibLatF    =   open(sibLatIN,"r")
else:
    print sibLonIN,"is either unreadable or DNE."
    raise SystemExit
if(os.path.isfile(eurLonIN)):
    eurLonF    =   open(eurLonIN,"r")
else:
    print eurLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(eurLatIN)):
    eurLatF    =   open(eurLatIN,"r")
else:
    print eurLonIN,"is either unreadable or DNE."
    raise SystemExit
if(os.path.isfile(karLonIN)):
    karLonF    =   open(karLonIN,"r")
else:
    print karLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(karLatIN)):
    karLatF    =   open(karLatIN,"r")
else:
    print karLonIN,"is either unreadable or DNE."
    raise SystemExit
if(os.path.isfile(lapLonIN)):
    lapLonF    =   open(lapLonIN,"r")
else:
    print lapLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(lapLatIN)):
    lapLatF    =   open(lapLatIN,"r")
else:
    print lapLonIN,"is either unreadable or DNE."
    raise SystemExit
if(os.path.isfile(lomLonIN)):
    lomLonF    =   open(lomLonIN,"r")
else:
    print lomLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(lomLatIN)):
    lomLatF    =   open(lomLatIN,"r")
else:
    print lomLonIN,"is either unreadable or DNE."
    raise SystemExit
if(os.path.isfile(norLonIN)):
    norLonF    =   open(norLonIN,"r")
else:
    print norLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(norLatIN)):
    norLatF    =   open(norLatIN,"r")
else:
    print norLonIN,"is either unreadable or DNE."
    raise SystemExit
if(os.path.isfile(barLonIN)):
    barLonF    =   open(barLonIN,"r")
else:
    print barLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(barLatIN)):
    barLatF    =   open(barLatIN,"r")
else:
    print barLonIN,"is either unreadable or DNE."
    raise SystemExit
if(os.path.isfile(chukLonIN)):
    chukLonF    =   open(chukLonIN,"r")
else:
    print chukLonIN,"is either unreadable or DNE."
    raise SystemExit(0)
if(os.path.isfile(chukLatIN)):
    chukLatF    =   open(chukLatIN,"r")
else:
    print chukLonIN,"is either unreadable or DNE."
    raise SystemExit

# Parse file into a list
beauLat     =   beauLatF.read().split('\n')
beauLon     =   beauLonF.read().split('\n')
canLat      =   canLatF.read().split('\n')
canLon      =   canLonF.read().split('\n')
sibLat      =   sibLatF.read().split('\n')
sibLon      =   sibLonF.read().split('\n')
eurLat      =   eurLatF.read().split('\n')
eurLon      =   eurLonF.read().split('\n')
lapLat      =   lapLatF.read().split('\n')
lapLon      =   lapLonF.read().split('\n')
lomLat      =   lomLatF.read().split('\n')
lomLon      =   lomLonF.read().split('\n')
norLat      =   norLatF.read().split('\n')
norLon      =   norLonF.read().split('\n')
karLat      =   karLatF.read().split('\n')
karLon      =   karLonF.read().split('\n')
barLat      =   barLatF.read().split('\n')
barLon      =   barLonF.read().split('\n')
chukLat     =   chukLatF.read().split('\n')
chukLon     =   chukLonF.read().split('\n')

# Trim each array for the inevitable empty last index
beauLat     =   beauLat[:len(beauLat)-1]
beauLon     =   beauLon[:len(beauLon)-1]
canLat      =   canLat[:len(canLat)-1]
canLon      =   canLon[:len(canLon)-1]
sibLat      =   sibLat[:len(sibLat)-1]
sibLon      =   sibLon[:len(sibLon)-1]
eurLat      =   eurLat[:len(eurLat)-1]
eurLon      =   eurLon[:len(eurLon)-1]
lapLat      =   lapLat[:len(lapLat)-1]
lapLon      =   lapLon[:len(lapLon)-1]
lomLat      =   lomLat[:len(lomLat)-1]
lomLon      =   lomLon[:len(lomLon)-1]
norLat      =   norLat[:len(norLat)-1]
norLon      =   norLon[:len(norLon)-1]
karLat      =   karLat[:len(karLat)-1]
karLon      =   karLon[:len(karLon)-1]
barLat      =   barLat[:len(barLat)-1]
barLon      =   barLon[:len(barLon)-1]
chukLat     =   chukLat[:len(chukLat)-1]
chukLon     =   chukLon[:len(chukLon)-1]

# get the size of the list
beauPts     =   len(beauLon)
canPts      =   len(canLon)
sibPts      =   len(sibLon)
eurPts      =   len(eurLon)
lapPts      =   len(lapLon)
lomPts      =   len(lomLon)
norPts      =   len(norLon)
karPts      =   len(karLon)
barPts      =   len(barLon)
chukPts     =   len(chukLon)
arcticSum   =   beauPts+canPts+sibPts+eurPts+lapPts+lomPts+norPts+karPts+barPts+chukPts

"""****************************************************************************
Convert the read coordinates to integers that have the 
off one index for Python.
****************************************************************************"""
for i in range(0,beauPts):
    beauLat[i]=int(float(beauLat[i].strip()))-1
    beauLon[i]=int(float(beauLon[i].strip()))-1
for i in range(0,canPts):
    canLat[i]=int(float(canLat[i].strip()))-1
    canLon[i]=int(float(canLon[i].strip()))-1
for i in range(0,sibPts):
    sibLat[i]=int(float(sibLat[i].strip()))-1
    sibLon[i]=int(float(sibLon[i].strip()))-1
for i in range(0,eurPts):
    eurLat[i]=int(float(eurLat[i].strip()))-1
    eurLon[i]=int(float(eurLon[i].strip()))-1
for i in range(0,lapPts):
    lapLat[i]=int(float(lapLat[i].strip()))-1
    lapLon[i]=int(float(lapLon[i].strip()))-1
for i in range(0,lomPts):
    lomLat[i]=int(float(lomLat[i].strip()))-1
    lomLon[i]=int(float(lomLon[i].strip()))-1
for i in range(0,norPts):
    norLat[i]=int(float(norLat[i].strip()))-1
    norLon[i]=int(float(norLon[i].strip()))-1
for i in range(0,karPts):
    karLat[i]=int(float(karLat[i].strip()))-1
    karLon[i]=int(float(karLon[i].strip()))-1
for i in range(0,barPts):
    barLat[i]=int(float(barLat[i].strip()))-1
    barLon[i]=int(float(barLon[i].strip()))-1
for i in range(0,chukPts):
    chukLat[i]=int(float(chukLat[i].strip()))-1
    chukLon[i]=int(float(chukLon[i].strip()))-1
"""****************************************************************************
Search the database for coordinates in each reigon and parase those records 
into subDFs. Lists of references are avaiable to view and download.
****************************************************************************"""
# Convert the DB lat and lon to integers...
df.lon = df.lon.astype(int)
df.lat = df.lat.astype(int)

# Search the DB for Beaufort sea coordiantes, starting with lat then lon.
beauDF = df.loc[df.lat==beauLat[0]-89]           # create a subDF for each reigon starting with Lat.
beauDF = beauDF.loc[beauDF.lon==beauLon[0]-179]  # hone in with Lon!
for i in range(1,beauPts):
    # Make each coordinate of the arctic a unique dataframe.
    beauT  = df.loc[df.lat==beauLat[i]-89]         # create a temp DF for specific reigonal Lat.
    beauT  = beauT.loc[beauT.lon==beauLon[i]-179]  # within temp DF hone in on specific reigonal Lon.
    beauDF = beauDF.append(beauT)                  # Add the reigonal points to the Master.      
# Unique count and list of arctic refernces
beauRefList = beauDF.ref.unique()
beauRefSize = beauRefList.size

# We're going to do the exact same routine for each arctic reigon...
# Search the DB for Canadian basin coordiantes, starting with lat then lon.
canDF = df.loc[df.lat==canLat[0]-89]
canDF = canDF.loc[canDF.lon==canLon[0]-179]     
for i in range(1,canPts):
    canT  = df.loc[df.lat==canLat[i]-89]        
    canT  = canT.loc[canT.lon==canLon[i]-179]   
    canDF = canDF.append(canT)                        
# Unique count and list of arctic refernces
canRefList = canDF.ref.unique()
canRefSize = canRefList.size

# Search the DB for E siberian shelf coordiantes, starting with lat then lon.
sibDF = df.loc[df.lat==sibLat[0]-89]           
sibDF = sibDF.loc[sibDF.lon==sibLon[0]-179]      
for i in range(1,sibPts):
    sibT  = df.loc[df.lat==sibLat[i]-89]          
    sibT  = sibT.loc[sibT.lon==sibLon[i]-179]      
    sibDF = sibDF.append(sibT)                        
# Unique count and list of arctic refernces
sibRefList = sibDF.ref.unique()
sibRefSize = sibRefList.size

# Search the DB for eursian basin coordiantes, starting with lat then lon.
eurDF = df.loc[df.lat==eurLat[0]-89]         
eurDF = eurDF.loc[eurDF.lon==eurLon[0]-179]  
for i in range(1,eurPts):
    eurT  = df.loc[df.lat==eurLat[i]-89]         
    eurT  = eurT.loc[eurT.lon==eurLon[i]-179]  
    eurDF = eurDF.append(eurT)                     
# Unique count and list of arctic refernces
eurRefList = eurDF.ref.unique()
eurRefSize = eurRefList.size

lapDF = df.loc[df.lat==lapLat[0]-89]         # create a subDF for each reigon starting with Lat.
lapDF = lapDF.loc[lapDF.lon==lapLon[0]-179]  # hone in with Lon!
for i in range(1,lapPts):
    # Make each coordinate of the arctic a unique dataframe.
    lapT  = df.loc[df.lat==lapLat[i]-89]         # create a temp DF for specific reigonal Lat.
    lapT  = lapT.loc[lapT.lon==lapLon[i]-179]  # within temp DF hone in on specific reigonal Lon.
    lapDF = lapDF.append(lapT)                  # Add the reigonal points to the Master.      
# Unique count and list of arctic refernces
lapRefList = lapDF.ref.unique()
lapRefSize = lapRefList.size

# Search the DB for lomadian basin coordiantes, starting with lat then lon.
lomDF = df.loc[df.lat==lomLat[0]-89]
lomDF = lomDF.loc[lomDF.lon==lomLon[0]-179]     
for i in range(1,lomPts):
    lomT  = df.loc[df.lat==lomLat[i]-89]        
    lomT  = lomT.loc[lomT.lon==lomLon[i]-179]   
    lomDF = lomDF.append(lomT)                        
# Unique count and list of arctic refernces
lomRefList = lomDF.ref.unique()
lomRefSize = lomRefList.size

# Search the DB for E norerian shelf coordiantes, starting with lat then lon.
norDF = df.loc[df.lat==norLat[0]-89]           
norDF = norDF.loc[norDF.lon==norLon[0]-179]      
for i in range(1,norPts):
    norT  = df.loc[df.lat==norLat[i]-89]          
    norT  = norT.loc[norT.lon==norLon[i]-179]      
    norDF = norDF.append(norT)                        
# Unique count and list of arctic refernces
norRefList = norDF.ref.unique()
norRefSize = norRefList.size

# Search the DB for karsian basin coordiantes, starting with lat then lon.
karDF = df.loc[df.lat==karLat[0]-89]         
karDF = karDF.loc[karDF.lon==karLon[0]-179]  
for i in range(1,karPts):
    karT  = df.loc[df.lat==karLat[i]-89]         
    karT  = karT.loc[karT.lon==karLon[i]-179]  
    karDF = karDF.append(karT)                     
# Unique count and list of arctic refernces
karRefList = karDF.ref.unique()
karRefSize = karRefList.size

# Search the DB for E chukerian shelf coordiantes, starting with lat then lon.
chukDF = df.loc[df.lat==chukLat[0]-89]           
chukDF = chukDF.loc[chukDF.lon==chukLon[0]-179]      
for i in range(1,chukPts):
    chukT  = df.loc[df.lat==chukLat[i]-89]          
    chukT  = chukT.loc[chukT.lon==chukLon[i]-179]      
    chukDF = chukDF.append(chukT)                        
# Unique count and list of arctic refernces
chukRefList = chukDF.ref.unique()
chukRefSize = chukRefList.size

# Search the DB for barsian basin coordiantes, starting with lat then lon.
barDF = df.loc[df.lat==barLat[0]-89]         
barDF = barDF.loc[barDF.lon==barLon[0]-179]  
for i in range(1,barPts):
    barT  = df.loc[df.lat==barLat[i]-89]         
    barT  = barT.loc[barT.lon==barLon[i]-179]  
    barDF = barDF.append(barT)                     
# Unique count and list of arctic refernces
barRefList = barDF.ref.unique()
barRefSize = barRefList.size

arcMaster   =   eurDF.append(sibDF).append(canDF).append(beauDF).append(lapDF)\
                .append(lomDF).append(norDF).append(karDF).append(chukDF).append(barDF)
arcRef      =   arcMaster.ref.unique()
arcRefSize  =   arcRef.size

print"\nThere are",len(arcMaster.index),"unique coordinate points in the arctic."
print"There are",arcRefSize,"refs that can be viewed by typing arcRef into the shell."

# Revive DB! back to .00 percision from 0. percision   
df=pd.read_fwf(filein,widths=fWid,header=None,names=col,skip_blank_lines=True) 
for i in range(0,files2Add):
    print "opening",toAdd[i],""
    csv2DF=pd.read_csv(toAdd[i],sep=',',skiprows=1,names=col,na_filter=True)
    #              File Name,Delimiter, skip header
    dfRows= len(df.index-1)
    
    df = df.append(csv2DF)     # add the contents of new csv into DB.
    dfRows= len(df.index-1) # update the number of records in the DB.
print dfRows-ogLen,"records added to the database."
"""****************************************************************************
Statistics for each reigon
****************************************************************************"""
if (updates == 1):
    print"\n*******************************************************************"
    print"Separating Global Seawater O18 Database Arctic Points by Reigon...\n"
    print"*******************************************************************"
    print "There are",beauRefSize,"references to the Beaufort Sea from",\
    len(beauDF.index),"unique points.\
    \nWith a temporal distribution of",len(beauDF.month.unique()),"months and",len(beauDF.year.unique())\
    ,"years...\nIn the shell type beauRefSize for a list of all reigonal refs.\n"
    print "There are",canRefSize,"references to the Canadian Basin from",\
    len(canDF.index),"unique points.\
    \nWith a temporal distribution of",len(canDF.month.unique()),"months and",len(canDF.year.unique())\
    ,"years...\nIn the shell type canRefSize for a list of all reigonal refs.\n"
    print "There are",sibRefSize,"references to the E Siberian Shelf from",\
    len(sibDF.index),"unique points.\
    \nWith a temporal distribution of",len(sibDF.month.unique()),"months and",len(sibDF.year.unique())\
    ,"years...\nIn the shell type sibRefSize for a list of all reigonal refs.\n"
    print "There are",eurRefSize,"references to the Eursian Basin from",\
    len(eurDF.index),"unique points.\
    \nWith a temporal distribution of",len(eurDF.month.unique()),"months and",len(eurDF.year.unique())\
    ,"years...\nIn the shell type eurRefSize for a list of all reigonal refs.\n"
    print "There are",lapRefSize,"references to the Laplev Sea from",\
    len(lapDF.index),"unique points.\
    \nWith a temporal distribution of",len(lapDF.month.unique()),"months and",len(lapDF.year.unique())\
    ,"years...\nIn the shell type lapRefSize for a list of all reigonal refs.\n"
    print "There are",lomRefSize,"references to the Lomonosov Ridge from",\
    len(lomDF.index),"unique points.\
    \nWith a temporal distribution of",len(lomDF.month.unique()),"months and",len(lomDF.year.unique())\
    ,"years...\nIn the shell type lomRefSize for a list of all reigonal refs.\n"
    print "There are",norRefSize,"references to the Norweigian Sea from",\
    len(norDF.index),"unique points.\
    \nWith a temporal distribution of",len(norDF.month.unique()),"months and",len(norDF.year.unique())\
    ,"years...\nIn the shell type norRefSize for a list of all reigonal refs.\n"
    print "There are",karRefSize,"references to the Kara Sea from",\
    len(karDF.index),"unique points.\
    \nWith a temporal distribution of",len(karDF.month.unique()),"months and",len(karDF.year.unique())\
    ,"years...\nIn the shell type karRefSize for a list of all reigonal refs.\n"
    print "There are",chukRefSize,"references to the Chukchi Sea from",\
    len(chukDF.index),"unique points.\
    \nWith a temporal distribution of",len(chukDF.month.unique()),"months and",len(chukDF.year.unique())\
    ,"years...\nIn the shell type chukRefSize for a list of all reigonal refs.\n"
    print "There are",barRefSize,"references to the Barents Sea from",\
    len(barDF.index),"unique points.\
    \nWith a temporal distribution of",len(barDF.month.unique()),"months and",len(barDF.year.unique())\
    ,"years...\nIn the shell type barRefSize for a list of all reigonal refs."

"""****************************************************************************
Search each subDF for records of surfance depth, and plot the d18O and salinity
with the WOA-reigonal salinity averages.
****************************************************************************"""
# List of Params to be iterated through...
# Files
WOAF          =   'WOAsalMonthly.nc' # Average monthly salinity from WOA2001 in NetCDF
DFList      =   (eurDF,beauDF,canDF,sibDF,lapDF,lomDF,norDF,karDF,chukDF,barDF)

reigonName  =   ('Eurasian Basin','Beaufort Sea','Canadian Basin','E Siberian Shelf','Laplev Sea'
                ,'Lomonosov Ridge','Norweigian Sea','Kara Sea','Chukchi Sea','Barents Sea')

reigons     =   len(DFList)

lblSal      =     ('Eursian Shelf Observed Salinity ','Beaufort Sea Observed Salinity ',
                   'Canadian Basin Observed Salinity','E Siberian Shelf Observed Salinity '
                 ,'Laplev Sea Observed Salinity ','Lomonosov Ridge Observed Salinity ',
                 'Norweigian Sea Observed Salinity','Kara Sea Observed Salinity '
                 ,'Chuckchi Sea Observed Salinity ','Barents Sea Observed Salinity ')
                          
lblD18      =   ('Eursian Shelf $\delta$18O','Beaufort Sea $\delta$18O','Canadian Basin $\delta$18O',
                 'E Siberian Shelf $\delta$18O','Laptev Sea $\delta$18O','Lomonosov Ridge $\delta$18O',
                 'Norweigian Sea $\delta$18O','Kara Sea $\delta$18O','Chuckchi Sea $\delta$18O','Barents Sea $\delta$18O')

colorSwatch =   _get_colors(reigons)

markSal     =   ('1','2','x','+')
markD18     =   ('o','p','d','v')
sizeSal     =   (20,20,16,20)
sizeD18     =   (10,14,10,8)
linWid      =   3.5
skip        =   -99.9
surfDep     =   10      # depth of the surface (meters).
