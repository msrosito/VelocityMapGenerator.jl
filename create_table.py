import eagleSqlTools as sql
import numpy as np
import matplotlib.pyplot as plt

# Creates a table from information in the EAGLE public database
# https://arxiv.org/pdf/1510.01320.pdf

# Array of chosen simulations. Entries refer to the simulation name and comoving box length.
Simu = np.array([('RefL0100N1504', 100.)])

# This uses the eagleSqlTools module to connect to the database with your username and password.
# If the password is not given, the module will prompt for it.

con = sql.connect("<username>", password="<password>")

# Construct and execute query for the simulation. This query returns the GalaxyID
# the group and subgroup numbers, the star hmr, the stellar mass and the virial radius
# for galaxies whose stellar masses are greater than 1e10.
            
myQuery = "SELECT \
    SH.GalaxyID as galID, \
	SH.GroupNumber as gr, \
	SH.SubGroupNumber as sub, \
	SH.HalfMassRad_Star as hmr, \
	SH.MassType_Star as Mstar, \
	FOF.Group_R_Crit200 as Rvir \
           FROM \
	RefL0100N1504_SubHalo as SH, \
	RefL0100N1504_FOF as FOF \
           WHERE \
	SH.SnapNum = 28 and \
	SH.MassType_Star > 1.e10 and \
	SH.SnapNum = FOF.SnapNum and \
	SH.GroupID = FOF.GroupID \
          GROUP BY \
   SH.GalaxyID, \
   SH.GroupNumber, \
   SH.SubGroupNumber, \
   SH.HalfMassRad_Star, \
   SH.MassType_Star, \
   FOF.Group_R_Crit200 \
          ORDER BY \
   galID"
    
# Execute query.
myData 	= sql.execute_query(con, myQuery)

n = len(myData)

# Construct the table
table = np.zeros((n, 6))

for i in range(n):
    table[i][0] = myData['galID'][i]
    table[i][1] = myData['gr'][i]
    table[i][2] = myData['sub'][i]
    table[i][3] = myData['hmr'][i]
    table[i][4] = myData['Mstar'][i]
    table[i][5] = myData['Rvir'][i]

np.savetxt('table_galaxies.dat', table)
