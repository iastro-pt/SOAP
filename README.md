# SOAP

SOAP version 2.0, downloaded from 
[this page](http://astro.up.pt/resources/soap2/).
The original readme is reproduced below.




      ######          ###########            #####         ###########        #######
    #        #      #             #        #       #       #          #     #         #
   #               #               #     #           #     #           #               #
  #                #        ##     #    #             #    #           #               #
   #               #        ##     #    #             #    #           #              #
     #             #               #    #             #    #          #             #
        #          #               #    #             #    ###########            #
           #       #               #    #             #    #                    #
             #     #               #    # ########### #    #                  #
              #    #               #    #             #    #                #
             #     #               #    #             #    #               #
   #        #       #             #     #             #    #               #
     ######           ###########       #             #    #               ###########
 

   |                               #
   |                           #        #
   |                         #            #
   |                       #                #
   |                     #                    #
   | ###################                        #                       #################
   |                                              #                   #
   |                                                #               #
   | WWW.ASTRO.UP.PT/RESOURCES/SOAP2                  #           #
   | Version 1.0 - Xavier Dumusque                      #       #
   | xdumusque@cfa.harvard.edu                              #
   |_____________________________________________________________________________________
   

#########################################
INTRODUCTION
#########################################

SOAP 2.0 is a code that estimates the effects of active regions, spots or plages, on radial-velocity and photometric measurements. SOAP 2.0 gives access to the radial-velocity, the bisector span, the full width at half maximum variations as defined with the cross-correlation technique and optimized to reproduce the HARPS (HARPS-N) observations. The photometric variation induced by active regions (flux at 5293 Angstrom) is also returned by SOAP 2.0.

The code is published in Dumusque et al. (2014, ApJ, 796, 132) and you are advised to look at this paper and reference therein for more information about the code. To summarize, the code is based on SOAP (Boise et al. 2012, A&A, 545, 109) with some modifications to include the effect of the inhibition of convection in magnetic regions due to strong local magnetic fields. This effect has been shown to be significant for stars similar to the Sun for which plages are dominating (Meunier et al 2010, A&A, 512, 39). To include this effect, we used observed spectra of the Sun taken with the Kitt Peak Fourrier Transform Spectrograph (FTS). Including this inhibition of the convective blueshift effect using solar spectra allows reproducing better the observations of slow rotators (see Dumusque 2014, ApJ, 796, 133).


#########################################
INSTALLATION
#########################################

The SOAP 2.0 code is written in C, with a python interface for plotting the results.
The code has been tested on python 2.6 and 2.7 and though it is expected to work properly on more recent versions, the user is advised to use this version in case of troubles.

The following standard python libraries are required:
- numpy 1.6.2 or higher
- scipy 0.11.0 or higher
- matplotlib 1.1.1 or higher
- pyfits 3.3 or higher

The following C libraries are required:
- gsl

You can download the latest version of SOAP 2.0 here:
www.astro.up.pt/ressources/soap2

Once downloaded, unzip the SOAP_2.zip file
$ unzip SOAP_2.zip

To run SOAP 2.0, you first have to compile the C code to be usable by your version of python.
To do so go on to the StarSpot folder and compile the C code using the setup.py file. Then copy the created file in the StarSpot folder:
$ cd SOAP_2/StarSpot                     (where SOAP_PATH is where you have installed SOAP 2.0)
$ python setup.py build                     (this create the wrapped file in build/lib.***) 
$ cp build/lib.****/starspot.so .           (where **** corresponds to your operating system)

You might have a warning "Using deprecated NumPy API...” when compiling the C file using the command “python setup.py build”, however this will not affect SOAP 2.


#########################################
RUNNING SOAP 2.0
#########################################

CONFIGURATION:
The configuration of SOAP 2.0 is done using the config.cfg file. Inside this ascii file, you can configure the resolution of the grid used for the simulation (grid) and the resolution of the active region circumference (nrho). The instrumental resolution can also be changed.  It is of course possible to configure the stellar properties (inclination, rotational period, radius..) and the active regions properties (active region type (spot, plage), longitude, latitude, size). Finally, it is possible to configure at which phase the output is generated. Everything is commented in the config.cfg file and it should be easy to configure it yourself.

Up to four active regions, spots or plages, can be created with SOAP 2.0, however it is easy to configure the program to accept more active regions.

Once the config.cfg file is configured to your needs, run the soap2.py file with the command  (you must be in the folder where you installed SOAP 2.0):

$ python soap2.py
or
$ ipython soap2.py (for the one that have installed ipython (much more flexible than python))


DISPLAY:
To work with python in interactive mode, which display interactive figures, you should configure the correct backend for you operating system. This can be done by opening the matplotlibrc file in the main SOAP_2 folder and change the value of the backend variable. For MacOSX, choose MacOSX, for Linux, choose belong GTKAgg, TkAGG, WX, WXAgg, QtAgg, Qt4Agg and test if it works (see here for more information http://matplotlib.org/faq/usage_faq.html#what-is-a-backend). You can check in your original matplotlibrc file what is the value for the backend, as this one should be correctly configured (this original file should be in the matplotlib directory that can be find in the python site-package folder: /Library/Python/2.7/site-packages/matplotlib-1.3.1-py2.7-macosx-10.9-intel.egg/matplotlib/mpl-data for my MacOSX installation)

If you are in interactive mode two plots appear. The first plot shows the two CCF derived from the quiet photosphere spectrum and from the spot spectrum of the Kitt Peak FTS solar observations. The second plot shows as a function of phase the effect in flux, radial-velocity (RV), bisector span (Bis Span) and full width at half maximum (Fwhm) of the disc-integrated CCF. Note that for each observable, the flux effect (only due to the contrast difference of the active regions), the convective blueshift effect (only due to the inhibition of the convective blueshift inside active regions) and the combined effect are returned.

If you are not in interactive mode (Agg backend), the figures are not displayed.

Note that if the wrong interactive backend is selected, python might crash (not always the case) because it cannot find the backend, and therefore no output will be generated. If you cannot find any interactive backend that works, you should select the Agg backend. You will not be able to have interactive figures, but they will be saved in .pdf format anyway.


OUTPUTS:
The following files will be created in a specific folder that can be found in the ‘outputs’ directory. The name of this folder depends on the active region configuration that you chose and will contain:
- A .pdf file for each generated figure
- a .rdb file that gives the value of all the observables for all calculated phase
- a ‘fits’ folder where .fits files for each phase is created. The headers of these files contain all the information about the configuration used and the calculated observables. The data of the fits files contain the disc-integrated CCF from which the RV, Bis Span and Fwhm are calculated.


CODE RAPIDITY:
The rapidity of the code can be increased by reducing the number of phase at which the outputs are calculated. In addition the efficiency of the code can be improved by reducing  the stellar grid resolution (GRID = 300) or the spot circumference resolution (NRHO = 20). I would not select a grid resolution smaller than 100 and a circumference resolution smaller than 10. For big active regions, the decrease in resolution will not influence strongly the observable variations, but be careful for small active regions. No serious test have been done to check the effect of small resolution and I would advise to always compare the results with the standard resolution (GRID = 300, NRHO = 20) to be sure that the observed variations are induced by active regions and not by numerical imprecision.


#########################################
TROUBLESHOOTINGS
#########################################

Did you make sure you are using python 2.6 or 2.7, and have installed the libraries specified above?

If you do not have one of the python libraries, you can install them that way:
$ sudo easy_install library_name (where library_name should be replaced by numpy, scipy...).
Another way is to use pip (https://pypi.python.org/pypi/pip)

If you do not have the C gsl library, you can install it easily on MacOSX with MacPort (https://www.macports.org/install.php). Once MacPort is installed, just type in a terminal the following command:
$ sudo port install gsl
Under Linux you can for example use yum to install gsl. Under Linux you will also need the gsl development package:
$ sudo yum install gsl
$ sudo yum install gsl-devel


If compiling the C code does not work because it can not find the numpy/arrayobject.h file
fatal error: ‘numpy/arrayobject.h’ file not found,
you have to manually update the path to where the file is. Open the setup.py file an update the  python_include_dir variable. For MacOSX, this should be something similar to:
‘/System/Library/Frameworks/Python.framework/Versions/2.7/Extras/lib/python/numpy/core/include’
or
‘/Library/Python/2.7/site-packages/numpy/core/include/’

Similar problem can happen for the gsl/gsl_errno.h file
‘gsl/gsl_errno.h' file not found.
In this case, update the  gsl_include_dir variable

If it cannot find the gsl libraries, it can also return an error message
‘library not found for -lgslcblas’.
In this case, update the gsl_lib_dir variable.

If still you cannot install and/or give the correct path to your libraries, please go and see you system administrator that will help you. This is not a problem of SOAP 2 and I will not reply to library problem requests.


#########################################
CONTACT
#########################################

Xavier Dumusque

xdumusque@cfa.harvard.edu

Harvard Smithsonian Center for Astrophysics
60 Garden Street
02138 Cambridge, MA
USA



