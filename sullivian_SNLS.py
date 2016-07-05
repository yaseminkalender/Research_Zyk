# SULLIVIAN 2010 PAPER
# Table 1 is used to work on the Star Formation Rate and Stellar Mass relation
#    ***************************************************
#    23 JUNE 2016
#    AUTHOR : Z. Yasemin Kalender

#Importing Packages
from astropy.io import ascii
from astropy.io import fits
import numpy as  np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from astropy.units import Unit, Quantity
import math

# Show the table that we use
t1=Table.read('sullivian2.fit')
print t1


def pecz_to_muerrsq(z, zerr, pecz=300):
    """
        Convert redshift uncertainty to equivalent uncertainty in mu
        Expects peculiar cz in km/s
        This is valid for almost any smooth cosmology
        """
    # Simple would be:
    #   c = 2.99792458e5 # km/s
    #   pez = pecz / c  # unitless fractional "redshift"
    # Formally using astropy constants and units
    from astropy.constants import c  # MKS
    # Define that pecz is in km/s
    # 'pez' is in unitless fractional "redshift"
    pez = ((Quantity(value=pecz, unit=Unit("km")/Unit("s")) / c).decompose()).value
    
    zfacsq = (5.0/np.log(10.0))**2
    dzerrsq = pez*pez + zerr*zerr
    dzerrsq *= zfacsq/(z*z)
    
    return dzerrsq


#def pecz_to_muerr(z, zerr, pecz=300):
def pecz_to_muerr(*args, **kwargs):
    return np.sqrt(pecz_to_muerrsq(*args, **kwargs))


i_mag_err=np.empty(231)
i_mag_err=i_mag_err.fill(0.01)
cos=FlatLambdaCDM(72,0.28)
z_values=t1["zcmb"] #importing redshifts from file
z_model=np.linspace(0,0.1,1000)+0.001
mu=cos.distmod(z_values).value # distance modulus for data
mu_model=cos.distmod(z_model).value # distance modulus for model
mu_model_err=pecz_to_muerr(z_model,0,pecz=300)

print "pec velocity:" , pecz_to_muerrsq(z_values,0,pecz=300)

intrinsic_dispersion = 0.08
#errsq=i_mag_err*i_mag_err+pecz_to_muerrsq(z_values,0,pecz=300)+intrinsic_dispersion*intrinsic_dispersion
absmagH=np.average(t1["imag"]-mu),#1/errsq)
print "Absolute M: " , absmagH

residuals=t1["imag"]-(mu+absmagH)


plt.errorbar(t1["logM_"],residuals,t1["e_logM_1"]*0,t1["e_logM_1"],'o')
plt.xlabel('Stellar Mass')
plt.ylabel('Residuals')
plt.title('Host Galaxy Sullivian')
plt.show()



'''

slope,intercept=np.polyfit(lum_dist_modulus,t1["imag"],1)
print "y-intercept: " ,intercept
print "slope of the function: " ,slope
plt.figure(2)
plt.plot(lum_dist_modulus,t1["imag"],'o')
plt.xlabel("Distance Modulus")
plt.ylabel("Apparent Magnitude")
plt.title('Sullivian 2010 Papaer')
#plt.show()



Apparent_model=intercept+lum_dist_modulus
Residuals=Apparent_model-t1["imag"]
Total_Residuals=(Residuals*Residuals)/(150*150)
#print ("Hubble Residuals: ") ,Residuals
#print ("Total Residuals: ") , Total_Residuals


Total_Apparent_Mag= np.sum(t1["imag"])
#print "Sum of apparent mag times their errors: ", Total_Apparent_Mag
Average_Apparent_Mag=Total_Apparent_Mag/231
#print "Average apparent Magnitude: ",Average_Apparent_Mag

Total_luminosity=np.sum(lum_dist_modulus)
print "Sum of all distance modulues: " , Total_luminosity
Average_lum_dist=Total_luminosity/231
print "Average value for distance modulus: " ,Average_lum_dist

Ave_Abs_Mag=Average_Apparent_Mag-Average_lum_dist
print "LCDM model Absolute Magnitude: " ,Ave_Abs_Mag

Model_Apparent=lum_dist_modulus+Ave_Abs_Mag
print "Model Apparent Magnitude" , Model_Apparent
Residuals=t1["imag"]-Model_Apparent
print Residuals

'''





#Plot for Stellar Mass versus Star Formation

plt.figure(1)
plt.errorbar(t1["logM_"],t1["logSFR"],[t1["e_logSFR1"],t1["E_logSFR"]],[t1["e_logM_1"],t1["E_logM_"]],fmt='o',ecolor='g',label="SNLS")
plt.xlabel("Log host galaxy M_stellar")
plt.ylabel("LOG host galaxy SFR")
plt.title("SULLIVIAN et al 2010, SNLS data")
plt.axis([7, 12, -3.5, 2.5])
plt.legend(loc=1)
plt.savefig('SFR_vs_StellarMass')
plt.show()





# the histogram of the data
n, bins, patches = plt.hist(t1["logM_"], 100 , normed=1, facecolor='green')

plt.xlabel('Log M_Stellar')
plt.ylabel('Relative Number')
plt.title(r'$\mathrm{Histogram\ of\ Log\ host\ Galaxy}\ $')
plt.axis([6, 13, 0, 1])
plt.grid(True)
plt.savefig('Stellar_Mass_Histogram')
plt.show()


plt.figure(3)
n, bins, patches = plt.hist(t1["logSFR"], 50 , normed=1, facecolor='red')
plt.xlabel('LogSFR')
plt.ylabel('Relative Number')
plt.title(r'$\mathrm{Histogram\ of\ Log\ host\ Galaxy}\ $')
plt.axis([-4,4, 0, 1])
plt.grid(True)
plt.savefig('SFR_Histogram')
plt.show()

 