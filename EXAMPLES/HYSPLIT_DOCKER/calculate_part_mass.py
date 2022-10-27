import numpy as np
import sys
from haversine import haversine
import os
import datetime

"""

The routine computes the solid mass from the outcomes of the simulation done by HYSPLIT

"""


filename = "partial_mass.part"

file=open(filename,'w')

GROUND=[]
AIR=[]

#print npart, n_levels, H_LEVELS
file.writelines(' \n')
file.writelines('*** MASS ON THE GROUND *** \n')
file.writelines(' \n')
#print ' '
#print '*** MASS ON THE GROUND ***'
#print ' '
# Check mass deposited on the ground

fname = 'CON2ASC.GROUND'

with open(fname) as f:

    for line in f:
         total_mass = 0
         filename = line.strip() 
         time = line.strip()[-4:]
         day = line.strip()[-8:-5]
           
         #print ' ---> day and time ',day,' ',time,' '
         file.writelines(' ---> day and time '+str(day)+' '+str(time)+' \n')

         count = 0
         header = []
         with open(filename) as f:
             lines = f.readlines()
             for line in lines:
                 to_check= line[0:3].strip()
                 if to_check.isdigit() == False:
                     header = header + line.split()
                     #print header 
                     count = count + 1

                 else:

                     break

         f.close()    
         
         m = []        

         for j in range(99):

             h_new = []

             to_find = 'CL'+str(int(j)).zfill(2)

             occurrence = 0
             
             for i in range(len(header)):

                 if to_find in header[i]:

                     occurrence = occurrence + 1

                     h_new.append(int(header[i][4:]))


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
        
         npart = m.shape[0]
         n_levels = h.shape[0]
         H_LEVELS = h

         #print 'Number of particle classes :',npart
         #print 'Heights :',H_LEVELS

         a = np.loadtxt(filename, skiprows = int(count))
         
         if a.shape[0] == 0 :
         
             #print 'No mass deposited at ',time
             file.writelines('No mass deposited at '+str(time)+'\n')        


         else:

             a = np.asarray(a) 

             a = a.reshape((-1,(npart * n_levels + 4)))

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

                 point1 = ((lat[i,0]-spacing_lat/float(2)),lon[i,0])
                 point2 = ((lat[i,0]+spacing_lat/float(2)),lon[i,0])
                  
                 d_lat = haversine(point1, point2) 


                 point3 = (lat[i,0],(lon[i,0]-spacing_lon/float(2)))
                 point4 = (lat[i,0],(lon[i,0]+spacing_lon/float(2)))
                  
                 d_lon = haversine(point3, point4)       
                 
                 dist.append([d_lat*1000,d_lon*1000])

             dist = np.asarray(dist)
             dist = dist.reshape((-1,2))

             a = a[:,4:]

             a = a.reshape((-1,(npart * n_levels)))

             column = 0

             for i in range(npart):

                     conc = a[:, column]

                     conc = conc.reshape((-1,1))


                     mass_on_the_ground = np.sum(conc[:,0] * dist[:,0] * dist[:,1])
                     #mass_on_the_ground = np.sum(conc[:,0])


                     total_mass = total_mass +  mass_on_the_ground 

                     #print 'class CL',str(i+1).zfill(2),' mass ','%.1e'%mass_on_the_ground,' kg'                  
                     #file.writelines('class CL'+str(i+1).zfill(2)+' mass ','%.1e'%mass_on_the_ground,' kg \n')

                     file.writelines("CL %d mass %.1e kg\n"%(i+1,mass_on_the_ground))
                     column = column + n_levels

         #print 'Total mass deposited ','%.1e'%total_mass,' kg'
         file.writelines('Total mass deposited %.1e kg \n'%(total_mass))

         GROUND.append([int(time),int(day),total_mass])
         os.remove(filename)
         
         #print ' '
         file.writelines(' \n')

#print '*** MASS IN THE AIR ***'

file.writelines('*** MASS IN THE AIR *** \n')
#print ' '
file.writelines(' \n')
# Check mass still in the atmosphere

fname = 'CON2ASC.AIR'

with open(fname) as f:

    for line in f:


         total_mass = 0
         filename = line.strip() 
         time = line.strip()[-4:]
         day = line.strip()[-8:-5]


         #print ' ---> day and time ',day,' ',time,' '
         file.writelines(' ---> day and time '+str(day)+' '+str(time)+' \n')

         a = np.loadtxt(filename, skiprows = int(count))

         file.writelines('filename '+str(filename)+'\n')

         
         if a.shape[0] == 0 :
         
             #print 'No in the atmosphere at ',time
             file.writelines('No in the atmosphere at ',time,'\n')  
        
         else:

             a = a.reshape((-1,(npart * n_levels + 4)))

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
             
             #spacing_lat = 0.1

             #spacing_lon = 0.1

             dist = []

             for i in range(a.shape[0]):

                 point1 = ((lat[i,0]-spacing_lat/float(2)),lon[i,0])
                 point2 = ((lat[i,0]+spacing_lat/float(2)),lon[i,0])
                  
                 d_lat = haversine(point1, point2) 


                 point3 = (lat[i,0],(lon[i,0]-spacing_lon/float(2)))
                 point4 = (lat[i,0],(lon[i,0]+spacing_lon/float(2)))
                  
                 d_lon = haversine(point3, point4)  

                 #d_lat = 10 
                 #d_lon = 10     
                 
                 dist.append([d_lat*1000,d_lon*1000])

                 

             dist = np.asarray(dist)
             dist = dist.reshape((-1,2))
 
             
             a = a[:,4:]

             a = a.reshape((-1,(npart * n_levels)))

             column = 0


             for i in range(npart):

                     mass_part=0

                     for j in range(n_levels-1):
  
                         conc = a[:, column + j + 1]
            

                         conc = conc.reshape((-1,1))

                         mass_in_the_air = np.sum(conc[:,0] * dist[:,0] * dist[:,1] * (int(H_LEVELS[j+1,0])-int(H_LEVELS[j,0])))

                         #mass_in_the_air = np.sum(conc[:,0])

                         mass_part = mass_part + mass_in_the_air


                         total_mass = mass_in_the_air +  total_mass


                         #print total_mass
                         #print 'class CL',str(i+1).zfill(2),' level ',H_LEVELS[j+1,0],'  mass ',mass_in_the_air,' kg'   

                         #print 'class CL',str(i+1).zfill(2),'  mass ','%.1e'%mass_in_the_air,' kg'   
                         #file.writelines('class CL',str(i+1).zfill(2),'  mass ','%.1e'%mass_in_the_air,' kg \n')
                     file.writelines('CL %d mass %.1e kg \n'%(i+1,mass_part))    

                     column = column + n_levels

              
         #print 'Total mass in the air ', '%.1e'%total_mass ,' kg'
         file.writelines('Total mass in the air %.1e kg \n'%(total_mass))
         #print ' '
          
         AIR.append([total_mass])
         os.remove(filename)

file.close()

GROUND = np.asarray(GROUND)
AIR = np.asarray(AIR)

file_mass=open('total_mass.part','w')

file_mass.writelines("Day    Time    Mass Deposited[kg]    Mass in the Air[kg]        Tot Mass[kg]\n")

print ( "***------***" )

print ( "--> Solid particle mass from HYSPLIT" )

for i in range(GROUND.shape[0]):

    tot = GROUND[i,2] + AIR[i,0]
    print ( ' - day %d time %d : %.1e'%(int(GROUND[i,1]),int(GROUND[i,0]),tot) )
    #print 'Mass deposited ', '%.1e'%GROUND[i,2]
    #print 'Mass in the air ',AIR[i,0]
    #TOT = GROUND[i,2] + AIR[i,0]
    #print 'Mass in the domain %.1e'%(TOT)

    
    file_mass.writelines("%d    %04d    %.1e               %.1e                    %.1e\n"%(GROUND[i,1],GROUND[i,0],GROUND[i,2],AIR[i,0],tot))

    
file_mass.close()



print ( "***------***" )


