import numpy as np
import sys
from haversine import haversine
import os


"""

The routine computes the gas mass (no H20) from the outcomes of the simulation done by HYSPLIT

"""

filename = "partial_mass.gas"

file=open(filename,'w')

AIR=[]

#print ngas, n_levels, H_LEVELS
print ' '
print '*** MASS IN THE AIR ***'
print ' '
# Check mass deposited on the ground

fname = 'CON2ASC.AIR'

with open(fname) as f:

    for line in f:
         total_mass = 0
         filename = line.strip() 
         time = line.strip()[-4:]
         day = line.strip()[-8:-5]
           
         #print ' ---> day and time ',day,' ',time,' '
         file.writelines(' ---> day and time '+str(day)+' '+str(time)+' \n')

         f = open(filename)
         data = f.read()
         first_line = data.split('\n', 1)[0]
        
         line_split = first_line.split()

         m = []        

         for j in range(99):

             h_new = []

             to_find = 'CL'+str(int(j)).zfill(2)

             occurrence = 0
             
             for i in range(len(line_split)):

                 if to_find in  line_split[i]:

                     occurrence = occurrence + 1

                     h_new.append(int(line_split[i][4:]))


             if occurrence > 0 :

                 h_old=np.asarray(h_new)

                 h_old=h_old.reshape((-1,1))

                 m.append([j,occurrence])

                 n_height = h_old.shape[0]

             if occurrence == 0 :
 
                 h = h_old
                                 
                 break                 

         m = np.asarray(m)
         m=m.reshape((-1,2))
        
         ngas = m.shape[0]
         n_levels = h.shape[0]
         H_LEVELS = h

         #print 'Number of particle classes :'ngas
         #print 'Heights :',H_LEVELS

         a = np.loadtxt(filename, skiprows = 1)
         
         if a.shape[0] == 0 :
         
             #print 'No in the atmosphere at ',time
             file.writelines('No mass deposited at '+str(time)+'\n')   
        
         else:

             a = a.reshape((-1,(ngas * n_levels + 4)))

             lat = a[:,2]

             lat = lat.reshape((-1,1))

             lon = a[:,3]

             lon = lon.reshape((-1,1))

             lon_unique = np.unique(lon)

             lat_unique = np.unique(lat)

             lon_unique = lon_unique.reshape((-1,1))
             lat_unique = lat_unique.reshape((-1,1))

             spacing_lat = np.absolute(lat_unique[0,0] - lat_unique[1,0])

             spacing_lon = np.absolute(lon_unique[0,0] - lon_unique[1,0])

             dist = []

             for i in range(a.shape[0]):

                 point1 = ((lat[i,0]-spacing_lat/2),lon[i,0])
                 point2 = ((lat[i,0]+spacing_lat/2),lon[i,0])
                  
                 d_lat = haversine(point1, point2) 


                 point3 = (lat[i,0],(lon[i,0]-spacing_lon/2))
                 point4 = (lat[i,0],(lon[i,0]+spacing_lon/2))
                  
                 d_lon = haversine(point3, point4)       
                 
                 dist.append([d_lat*1000,d_lon*1000])

                 #print 'd_lat, d_lon ',d_lat,d_lon

             dist = np.asarray(dist)
             dist = dist.reshape((-1,2))
 
             
             a = a[:,4:]

             a = a.reshape((-1,(ngas * n_levels)))

             column = 0

             for i in range(ngas):

                     for j in range(n_levels-1):

                         #print 'j ',j

                         #print  column + j + 1
  
                         conc = a[:, column + j + 1]
            

                         conc = conc.reshape((-1,1))

                         mass_in_the_air = np.sum(conc[:,0] * dist[:,0] * dist[:,1] * (int(H_LEVELS[j+1,0])-int(H_LEVELS[j,0])))

                        

                         total_mass = total_mass + mass_in_the_air

                         #total_mass = total_mass +  mass_on_the_ground 

                         #print 'class CL',str(i+1).zfill(2),' level ',H_LEVELS[j+1,0],'  mass ',mass_in_the_air,' kg'   

                     #print 'class CL',str(i+1).zfill(2),'  mass ','%.1e'%mass_in_the_air,' kg'   
                     file.writelines("component G %d mass %.1e kg\n"%(i+1,mass_in_the_air))
                         

                     column = column + n_levels

         file.writelines('Total mass in the air %.1e kg \n'%(total_mass))              
         #print 'Total mass in the air ', '%.1e'%total_mass ,' kg'
         #print ' '
          
         AIR.append([day,time,total_mass])
         os.remove(filename)

file.close()

AIR = np.asarray(AIR)
AIR = AIR.reshape((-1,3))

file_mass=open('total_mass.gas','w')

file_mass.writelines("Day    Time    Mass in the Air[kg]\n")

print "--> Gas mass from HYSPLIT"

for i in range(AIR.shape[0]):

    #print '*** day', int(GROUND[i,1]),' time ',str(int(GROUND[i,0])).zfill(4),' ***'
    #print 'Mass deposited ', '%.1e'%GROUND[i,2]
    #print 'Mass in the air ',AIR[i,0]
    #print 'Mass in the domain ',GROUND[i,2] + AIR[i,0]
    print ' - day %d time %d : %.1e'%(int(AIR[i,0]),int(AIR[i,1]),float(AIR[i,2]))    
    
    file_mass.writelines("%d    %04d    %.1e\n"%(int(AIR[i,0]),int(AIR[i,1]),float(AIR[i,2])))
    
file_mass.close()
print "***------***"
          


