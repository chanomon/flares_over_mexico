#!/usr/bin/env python
# coding: utf-8

# In[15]:


#!/usr/bin/env python
# By vdelaluz@igeofisica.unam.mx
#; https://www.ngdc.noaa.gov/stp/space-weather/solar-data/solar-features/solar-flares/x-rays/goes/xrs/

#;Output File Specification
#;   Column  Format  Description   
#;
#; 0    1- 2     I2    Data code: always 31 for x-ray events
#; 1    3- 5     I3    Station Code, 777 for GOES
#; 2    6- 7     I2    Year
#; 3    8- 9     I2    Month
#; 4   10-11     I2    Day
#; 5   12-13     A2    Astrisks mark record with unconfirmed change (What does this mean?)
#; 6   14-17     I4    Start time of x-ray event - SEE NOTE 1
#;     18        1X    <space>
#; 7   19-22     I4    End time
#;     23        1X    <space>
#; 8   24-27     I4    Max time
#;     28        1X    <space>
#; 9   29        A1    N or S for north or south latitude of xray flare if known
#;10   30-31     I2    Latitude of xray flare, if known
#;11   32        A1    E or W for east or west of longitude of xray flare, in known
#;12   33-34     I2    Central meridian distance of x-ray flare, if known
#;13   35-37     A3    SXI if data are from SXI imagery, blank otherwise
#;     38-59    22X    <space>
#;14   60        A1    X-ray class: C,M,X code - SEE NOTE 2
#;     61        1X    <space>
#;15   62-63     I2    X-ray intensity 10-99 for 1.0-9.9 x xray class
#;     64-67     4X    <space>
#;16   68-71     A4    Station ame abbreviation - "Gxx " for GOES
#;     72        1X    <space>
#;17   73-80   E7.1    Integrated flux (units = J/m**2)
#;18   81-85     I5    NOAA/USAF sunspot region number
#;     86        1X    <space>
#;19   87-88     I2    Year - central meridian passage (CMP)
#;20   89-90     I2    Month - central meridian passage (CMP)
#;21   91-94   F4.1    Day - central meridian passage (CMP)
#;     95        1X    <space>
#;22   96-102  F7.1    Total region area in squared arc seconds
#;    103        1X    <space>
#;23  104-110  F7.2    Total intensity (units - TBD) from SXI, if available
#;  ---------------------------------------------------------------------
#;  Note 1: Prior to 1997 if x-ray event could be corrolated to an optical
#;  event, then the time of the optical event was used.
#;  Note 2: X-ray class are classified according to the order of magnitude
#;  of the peak burst intensity (I) within the 0.1 - 0.8 nm band. The
#;  following apply:
#;  Class  +------Watt/m**2-----+
#;    B               I  <  10E-6         
#;    C    10E-6  <=  I  <  10E-5
#;    M    10E-5  <=  I  <  10E-4
#;    X               I  >  10E-4           
#;
#;------------------------------------------------------------------------
#C    2.1    G15    0.0019 = 1.9 * 10^-3
#C    2.1    G15    0.0015   = 1.5 * 10^-3
#C    1.8    G15    0.00039    = 3.9 * 10^-4
#B    7      G15    1e-04   
#B    9.3    G15    0.00015    = 15 * 10^-5
#B    9.7    G15    2e-04     

#;
#

import fortranformat as ff
from datetime import tzinfo, timedelta, datetime
from pysolar.solar import *
import glob
import os
from tqdm import tqdm


# In[16]:


ZERO = timedelta(0)
PATH= '/home/elizandro/MEGAsync/MAESTRIA/flares-sobre-mexico/'

def flux_by_class(flare_class,flare_subclass):###get the class, return the flux
    switcher = {
        'A': 1e-8,
        'B': 1e-7,
        'C': 1e-6,
        'M': 1e-5,
        'X': 1e-4,
    }
    scale = float(switcher.get(flare_class, 0))
    subclass= float(flare_subclass)
    #print subclass
    flux = subclass*scale  #scale+math.log10(subclass)*(scale*9.0)
    #print flare_class+str(flare_subclass)+'\t'+ str(flux)
    return flux

def get_date(data,time):
    if len(time.strip()) == 0:
        return True, datetime.datetime(1900,1,1,0,0,0,0, tzinfo=utc)
    #time = data[6]        
    year=int(data[2])
    if 0 <= year and year <= 70:  #so this function can endure far from 2017
        year = 2000+year
        hh= int(time[0:2])
        mm= int(time[2:4])
    else:
        year = 1900+year        
        hh= int(time[0:2])
        mm= int(time[2:4])

    if hh >= 24:
        hh = hh - 24
        time_event =  datetime.datetime(year,int(data[3]),int(data[4]),hh,mm,0,0, tzinfo=utc)
        time_event = time_event + datetime.timedelta(days=1)
    else:
        time_event =  datetime.datetime(year,int(data[3]),int(data[4]),hh,mm,0,0, tzinfo=utc)
    return False, time_event



class UTCtzinfo(tzinfo):
    def utcoffset(self, dt):
        return ZERO
    def tzname(self, dt):
        return "UTC"
    def dst(self, dt):
        return ZERO


# In[ ]:





# In[43]:



utc = UTCtzinfo()
#Coordenadas de coeneo
#LAT = 19.81 #deg
#LON = -100.3 #deg
           
#Coordenadas de isla guadalupe
LATG = 28.9
LONG = -118.3
#Coordenadas de isla mujeres
LATM = 21.2
LONM = -86.7


files = glob.glob(PATH+"goes*.txt")
#files = glob.glob(PATH+"goes-xrs-report_2017-lance-07-10.txt")
numofevents = 0
n=0
n_errors=0
eventTags = []
for filename in tqdm(files):
    print(filename)
    with open(filename, "r") as ins:
        for line in ins:
            try:
                #ffline = ff.FortranRecordReader('(I2,I3,I2,I2,I2,A2,A4,1X,A4,1X,A4,1X,A1,I2,A1,I2,A3,22X,A1,1X,I2,4X,A4,1X,E7.1,1X,I5,1X,I2,I2,F4.1,1X,F7.1,1X,F7.2)')
                ffline = ff.FortranRecordReader('(I2,I3,I2,I2,I2,A2,A4,1X,A4,1X,A4,1X,A1,I2,A1,I2,A3,22X,A1,I3,4X,A4,1X,E7.1,1X,I5,1X,I2,I2,F4.1,1X,F7.1,1X,F7.2)')
                data=ffline.read(line)
                T_value=False##################porque se hace esto?????=>creo que para llevar una bandera de que se leyo bien el tiempo
            except:
                try:
                    #Fixing GOES T
                    #ffline = ff.FortranRecordReader('(I2,I3,I2,I2,I2,A2,A4,1X,A4,1X,A4,1X,A1,I2,A1,I2,A3,22X,A1,1X,I2,4X,A4,1X,A1,E7.1,1X,I5,1X,I2,I2,F4.1,1X,F7.1,1X,F7.2)')
                    ffline = ff.FortranRecordReader('(I2,I3,I2,I2,I2,A2,A4,1X,A4,1X,A4,1X,A1,I2,A1,I2,A3,22X,A1,I3,4X,A4,1X,A1,E7.1,1X,I5,1X,I2,I2,F4.1,1X,F7.1,1X,F7.2)')
                    data=ffline.read(line)
                    T_value=True
                except:
                    break
            try:
                time = data[6]
                day = data[4]
                month = data[3]
                year = int(data[2])
            
                if 0 <= year and year <= 17:
                    year = 2000+year
                    directory=PATH+'database-mx/'+str(year)
                    #print directory
                    #if not os.path.exists(directory):
                    #    os.makedirs(directory)             #segun yo esto no es necesario
                    
                    error, t_start = get_date(data,time)
                
                
                else:
                    year = 1900+year        

                    directory=PATH+'database-mx/'+str(year)
                    #print directory
                    #if not os.path.exists(directory):
                    #    os.makedirs(directory)             #segun yo esto no es necesario
                    
                    error, t_start = get_date(data,time)
                
                #year=int(data[2])
                #if 0 <= year and year <= 17:
                #    year = 2000+year
                #else:
                #    year = 1900+year        
                #hh= int(time[0:2])
                #mm= int(time[2:4])
                #t_start = datetime.datetime(year,int(data[3]),int(data[4]),hh,mm,0,0, tzinfo=utc)
                time = data[8]
                error, t_max = get_date(data,time)

                #year=int(data[2])                
                    #if 0 <= year and year <= 17:
                #    year = 2000+year
                #else:
                #    year = 1900+year        
                #hh= int(time[0:2])
                #mm= int(time[2:4])
                    #t_max = datetime.datetime(year,int(data[3]),int(data[4]),hh,mm,0,0, tzinfo=utc)          	           	
                time = data[7]        
                error, t_end = get_date(data,time)

                if error:
                    t_end = t_max
                    
                    #year=int(data[2])
                #if 0 <= year and year <= 17:
                #    year = 2000+year
                #else:
                #    year = 1900+year
                    #
                    #month=int(data[3])
                #hh= int(time[0:2])
                    #
            	##plus1=False
            	##if hh>24          	
                #mm= int(time[2:4])
                #t_end = datetime.datetime(year,month,int(data[4]),hh,mm,0,0, tzinfo=utc)
            	    	    	
              
            
            
            
	    	
                #Si t_max es un objeto datetime.datetime las siguientes dos lineas funcionan
                altitudeG = get_altitude(LATG, LONG, t_max) #esta linea podria ser get_altitude(LAT,LON, t_max) 
                azimuthG = get_azimuth(LATG,LONG,t_max) #esta linea podria ser get_altitude(LAT,LON, t_max) si tienes la ultima version de pysolar
                # la mia es la 0.6
                altitudeM = get_altitude(LATM, LONM, t_max) #esta linea podria ser get_altitude(LAT,LON, t_max) 
                azimuthM = get_azimuth(LATM,LONM,t_max)
                
    
                if T_value:
                    integrated_flux = float(data[18])
                    sunspot =  data[19]                    
                else:
                    integrated_flux = float(data[17])
                    sunspot =  data[18]
    
                if sunspot is None:
                    sunspot = '0'
                #print t_start,t_max,t_end,  data[14], data[15], data[16], data[17],sunspot,altitude,azimuth
                flare_class = str(data[14])
                flare_subclass = str(float(data[15])/10.0)
                    
                try:
                    peak_flux=flux_by_class(flare_class,flare_subclass)
                except Exception as e:                
                    print(line)
                    n_errors=n_errors+1
                    print("error ni geting flux, ",str(e))
                    peak_flux=0.0
                if peak_flux > 1e-5 and (altitudeM >= 0 or altitudeG >= 0): #evento en alguno de los extremos longitudinales
                    numofevents += 1
                    text= str(t_start)+'\t'+str(t_max)+'\t'+str(t_end)+'\t'+flare_class+'\t'+flare_subclass+'\t'+str(data[16])+'\t'+'%1.5e'%peak_flux+'\t'+'%1.2e'%integrated_flux+'\t'+str(sunspot)+'\t'+ '%.2f'%altitude+'\t'+'%.2f'%azimuth+'\r\n'    
                    eventTag = str(flare_class)+str(flare_subclass)+'\t'+str(year)+'\t'+str(month)+'\t'+str(day)+'\t'+str(time)
                    eventTags.append(eventTag)
                    #print(eventTag)                    
                    #with open(PATH+"database-mx/"+str(year)+"/"+str(year)+'-'+str(month)+'-goes-xray.dat', "a") as myfile:
                    #    myfile.write(text)
                    with open(PATH+"database-mx/goes-xray-"+str(year)+".dat", "a") as myfile:
                        myfile.write(text)
                    n=n+1
            except Exception as e:                
                print(line)
                n_errors=n_errors+1
                print("error in asigning data, ",str(e))
    pass


# In[44]:


print('Processed:                    '+str(n))
print('Errors:                       '+str(n_errors))
print('Total:                        '+str(n+n_errors))
print('Number of events since 1975:  ' +str(numofevents))


# In[45]:


for event in eventTags:
    print(event)

