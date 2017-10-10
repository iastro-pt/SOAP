#!/usr/bin/env python

import numpy, os, sys
import pyfits
from pylab import *
import ConfigParser
import shutil
import string
sys.path.append('StarSpot/')
import PyStarSpot_soap2 as StSp

######################################################################
###           Functions            ###################################
######################################################################

# To read rdb file
def read_rdb(filename):
    
    f = open(filename, 'r')
    data = f.readlines()
    f.close()
    
    z=0
    while data[z][:2] == '# ' or data[z][:2] == ' #':
        z += 1

    key = string.split(data[z+0][:-1],'\t')
    output = {}
    for i in range(len(key)): output[key[i]] = []
    
    for line in data[z+2:]:
        qq = string.split(line[:-1],'\t')
        for i in range(len(key)):
            try: value = float(qq[i])
            except ValueError: value = qq[i]
            output[key[i]].append(value)

    return output

def write_rdb(filename,data,keys,format):
    
    f = open(filename, 'w')
        
    head1 = string.join(keys,'\t')
    head2 = ''
    for i in head1:
        if i=='\t': head2 = head2+'\t'
        else: head2 = head2+'-'

    f.write(head1+'\n')
    f.write(head2+'\n')
    
    if len(data.values()) > 0:
        for i in range(len(data.values()[0])):
            line = []
            for j in keys: line.append(data[j][i])
            f.write(format % tuple(line))
                            
    f.close()

####################################################################################
#                    INITITIALIZATION                     ##########################
####################################################################################

config = ConfigParser.ConfigParser()
config.read("config.cfg")

GRID              = int(config.get('main','grid' ))
NRHO              = int(config.get('main','nrho' ))
INST_RESO         = int(config.get('main','instrument_reso' ))
RAD_Sun           = int(config.get('star','radius_sun' ))
RAD               = float(config.get('star','radius' )) * RAD_Sun
PROT              = float(config.get('star','prot' ))
STAR              = StSp.Star(prot           = PROT,\
                              vrot           = (2.*pi*RAD)/(PROT*86400.),\
                              incl           = float(config.get('star','I')),\
                              limba1         = float(config.get('star','limb1')),\
                              limba2         = float(config.get('star','limb2')),\
                              psi            = float(config.get('star','psi')),\
                              rad_sun        = RAD_Sun,\
                              rad            = RAD,\
                              Temp           = int(config.get('star','Tstar')),\
                              Temp_diff_spot = int(config.get('star','Tdiff_spot')))

######################################################################
###           Read solar CCF       ###################################
######################################################################

if STAR.vrot < 10.:
    file1 = 'solar_CCFs/CCF_solar_spectrum_G2_FTS_reso_not_evenly_sampled_in_freq.rdb'
else:
    file1 = 'solar_CCFs/CCF_solar_spectrum_G2_FTS_reso_not_evenly_sampled_in_freq_extra_large_low_resolution.rdb'

data_ccf = read_rdb(file1)
rv_ccf = array(data_ccf['vrad'])
rv_ccf_magn_region = array(data_ccf['vrad'])
rv_ccf_sampling = rv_ccf[1]-rv_ccf[0]
intensity_ccf = array(data_ccf['CCF'])
rv_ccf_magn_region = array(data_ccf['vrad'])
intensity_ccf_magn_region = array(data_ccf['CCF_spot'])

####################################################################################
#                    INITITIALIZATION                     ##########################
####################################################################################

CCF               = StSp.Ccf(rv        = rv_ccf,\
                             intensity = intensity_ccf,\
                             width     = max(rv_ccf),\
                             step      = rv_ccf_sampling,\
                             star      = STAR)
               
CCF_active_region = StSp.Ccf(rv        = rv_ccf_magn_region,\
                             intensity = intensity_ccf_magn_region,\
                             width     = max(rv_ccf_magn_region),\
                             step      = rv_ccf_sampling,\
                             star      = STAR)

active_region1    = StSp.Active_region(check              = int(config.get('active_regions','check1')),\
                                       long               = float(config.get('active_regions','long1')),\
                                       lat                = float(config.get('active_regions','lat1')),\
                                       size               = float(config.get('active_regions','size1')),\
                                       active_region_type = int(config.get('active_regions','act_reg_type1')))

active_region2    = StSp.Active_region(check              = int(config.get('active_regions','check2')),\
                                       long               = float(config.get('active_regions','long2')),\
                                       lat                = float(config.get('active_regions','lat2')),\
                                       size               = float(config.get('active_regions','size2')),\
                                       active_region_type = int(config.get('active_regions','act_reg_type2')))

active_region3    = StSp.Active_region(check              = int(config.get('active_regions','check3')),\
                                       long               = float(config.get('active_regions','long3')),\
                                       lat                = float(config.get('active_regions','lat3')),\
                                       size               = float(config.get('active_regions','size3')),\
                                       active_region_type = int(config.get('active_regions','act_reg_type3')))

active_region4    = StSp.Active_region(check              = int(config.get('active_regions','check4')),\
                                       long               = float(config.get('active_regions','long4')),\
                                       lat                = float(config.get('active_regions','lat4')),\
                                       size               = float(config.get('active_regions','size4')),\
                                       active_region_type = int(config.get('active_regions','act_reg_type4')))

ACTIVE_REGIONS    = [active_region1,active_region2,active_region3,active_region4]

#calculate the visibility of the active regions
active_region1.calc_maps(STAR, GRID, NRHO)
active_region2.calc_maps(STAR, GRID, NRHO)
active_region3.calc_maps(STAR, GRID, NRHO)
active_region4.calc_maps(STAR, GRID, NRHO)

# Phase at which the signal we be generated
phase_step  = float(config.get('output','ph_step'))
phase_in    = config.get('output','ph_in')

# Read the phase in the file given in the "phase_in" argument in "config.cfg"
if phase_in!='None':
    f   = open(phase_in)
    PSI = array([float(l.split()[0]) for l in f.readlines()])
    f.close()
# Otherwise create a evenly sampled phase between 0 and 1 with the step given
# by the "ph_step" argument in "config.cfg"
elif phase_step>0.:
    PSI = arange(0.,1.,phase_step)
else:
    print
    print "WARNING: put a positive phase step ('ph_step' argument in 'config.cfg') or give a file name that contain the phases at which you want to compute the model ('phase_in' argument in 'config.cfg')"
    print


##########################################
#Calculate photometry and spectroscopy
##########################################

# Test if at least one active region has been activated in the config.cfg file
try:
    print
    print "****************************************************"
    print "Estimating the photometric and spectroscopic effects"
    print "****************************************************"
    print
    FLUXstar_quiet, CCFstar_quiet, flux, CCFstar_flux, CCFstar_bconv, CCFstar_tot, rv_flux, rv_bconv, rv_tot, span_flux, span_bconv, span_tot, fwhm_flux, fwhm_bconv, fwhm_tot, depth_flux, depth_bconv, depth_tot =\
                            StSp.Calculate_activity_signal(NRHO, GRID, STAR, PSI, CCF, CCF_active_region, ACTIVE_REGIONS,INST_RESO)

    # Calculates where rv_flux = rv_bconv, which corresponds to the phases where the active region is not visible
    index_equal_rv = where((rv_flux-rv_bconv)==0)[0]
    if len(index_equal_rv) != 0:
        zero_velocity = rv_flux[index_equal_rv][0] # velocity when the spot is not visible
        rv_flux  -= zero_velocity # Put the velocity when the spot is not visible to 0
        rv_bconv -= zero_velocity # Put the velocity when the spot is not visible to 0
        rv_tot   -= zero_velocity # Put the velocity when the spot is not visible to 0

    ########################################################################
    ################ Write to file #########################################
    ########################################################################

    print
    print "****************************************************"
    print "Writing to file"
    print "****************************************************"
    print
    CCF_folder_outputs = 'outputs/CCF_PROT=%.2f_i=%.2f_lon=(%.1f,%.1f,%.1f,%.1f)_lat=(%.1f,%.1f,%.1f,%.1f)_size=(%.4f,%.4f,%.4f,%.4f)/' % (STAR.prot,STAR.incl,active_region1.long,active_region2.long,active_region3.long,active_region4.long,active_region1.lat,active_region2.lat,active_region3.lat,active_region4.lat,active_region1.s,active_region2.s,active_region3.s,active_region4.s)
    CCF_folder_fits = CCF_folder_outputs + 'fits/'

    if os.path.exists(CCF_folder_outputs):
        shutil.rmtree(CCF_folder_outputs)
        os.mkdir(CCF_folder_outputs)
    else:
        os.mkdir(CCF_folder_outputs)

    if os.path.exists(CCF_folder_fits):
        shutil.rmtree(CCF_folder_fits)
        os.mkdir(CCF_folder_fits)
    else:
        os.mkdir(CCF_folder_fits)

    filename = "CCF_PROT=%.2f_i=%.2f_lon=(%.1f,%.1f,%.1f,%.1f)_lat=(%.1f,%.1f,%.1f,%.1f)_size=(%.4f,%.4f,%.4f,%.4f)" % (STAR.prot,STAR.incl,active_region1.long,active_region2.long,active_region3.long,active_region4.long,active_region1.lat,active_region2.lat,active_region3.lat,active_region4.lat,active_region1.s,active_region2.s,active_region3.s,active_region4.s)

    # Writing to fits file
    for i in range(len(PSI)):
        hdu = pyfits.PrimaryHDU(CCFstar_tot[i])
        hdulist = pyfits.HDUList([hdu])
        hdulist[0].header.set('STEP',CCF.step,'CCF Step [km/s]')
        hdulist[0].header.set('GRID',GRID,'Grid resolution')
        hdulist[0].header.set('NRHO',NRHO,'Spot circonference resolution')
        hdulist[0].header.set('RESO',NRHO,'Instrument resolution')
        hdulist[0].header.set('RADSUN',RAD_Sun,'Radius of the sun [km]')
        hdulist[0].header.set('RAD',STAR.rad,'Stellar radius [Rsun]')
        hdulist[0].header.set('PROT',STAR.prot,'Stellar Rotation Period [day]')
        hdulist[0].header.set('INCL',STAR.incl,'Stellar Inclination [degree]')
        hdulist[0].header.set('TSTAR',STAR.Temp,'Stellar effective temperature [K]')
        hdulist[0].header.set('TSPOTDIF',STAR.Temp_diff_spot,'Temp diff between photosphere and spot')
        hdulist[0].header.set('LIMBA1',STAR.limba1,'Stellar linear limb darkening coefficient')
        hdulist[0].header.set('LIMBA2',STAR.limba2,'Stellar linear limb darkening coefficient')
        hdulist[0].header.set('VROT',STAR.vrot,'Stellar Rotation [km/s] (calculated)')
        hdulist[0].header.set('TYPE1',active_region1.active_region_type_str,'Active region 1 type')
        hdulist[0].header.set('LONG1',active_region1.long,'Active region 1 longitude [degree]')
        hdulist[0].header.set('LAT1',active_region1.lat,'Active region 1 latitude [degree]')
        hdulist[0].header.set('SIZE1',active_region1.s,'Active region size [fraction of stellar radius]')
        hdulist[0].header.set('TYPE2',active_region2.active_region_type_str,'Active region 2 type')
        hdulist[0].header.set('LONG2',active_region2.long,'Active region 2 longitude [degree]')
        hdulist[0].header.set('LAT2',active_region2.lat,'Active region 2 latitude [degree]')
        hdulist[0].header.set('SIZE2',active_region2.s,'Active region size [fraction of stellar radius]')
        hdulist[0].header.set('TYPE3',active_region3.active_region_type_str,'Active region 3 type')
        hdulist[0].header.set('LONG3',active_region3.long,'Active region 3 longitude [degree]')
        hdulist[0].header.set('LAT3',active_region3.lat,'Active region 3 latitude [degree]')
        hdulist[0].header.set('SIZE3',active_region3.s,'Active region size [fraction of stellar radius]')
        hdulist[0].header.set('TYPE4',active_region4.active_region_type_str,'Active region 4 type')
        hdulist[0].header.set('LONG4',active_region4.long,'Active region 4 longitude [degree]')
        hdulist[0].header.set('LAT4',active_region4.lat,'Active region 4 latitude [degree]')
        hdulist[0].header.set('SIZE4',active_region4.s,'Active region size [fraction of stellar radius]')
        hdulist[0].header.set('RV_F',rv_flux[i],'RV for flux effect [km/s] (measured)')
        hdulist[0].header.set('SPAN_F',span_flux[i],'BIS for flux effect [km/s] (measured)')
        hdulist[0].header.set('FWHM_F',fwhm_flux[i],'FWHM for flux effect [km/s] (measured)')
        hdulist[0].header.set('RV_BC',rv_bconv[i],'RV for conv blue [km/s] (measured)')
        hdulist[0].header.set('SPAN_BC',span_bconv[i],'BIS for conv blue [km/s] (measured)')
        hdulist[0].header.set('FWHM_BC',fwhm_bconv[i],'FWHM for conv blue [km/s] (measured)')
        hdulist[0].header.set('RV',rv_tot[i],'RV [km/s] (measured)')
        hdulist[0].header.set('SPAN',span_tot[i],'BIS [km/s] (measured)')
        hdulist[0].header.set('FWHM',fwhm_tot[i],'FWHM [km/s] (measured)')
        hdulist[0].header.set('FLUX',flux[i],'Norm. Flux (measured)')
        hdulist[0].header.set('PSI',PSI[i],'Phase [0-1]')
        hdulist.writeto(CCF_folder_fits+filename+"_PSI=%.3f.fits" % PSI[i])

    # Writing to rdb file
    keys = ["Phase","Flux","RV_tot","Bis_Span_tot","Fwhm_tot","RV_flux","Bis_Span_flux","Fwhm_flux","RV_conv_blue","Bis_Span_conv_blue","Fwhm_conv_blue"]
    value = [PSI,flux,rv_tot,span_tot,fwhm_tot,rv_flux,span_flux,fwhm_flux,rv_bconv,span_bconv,fwhm_bconv]
    data = dict([(keys[i],value[i]) for i in arange(len(value))])
    format = "%f\t"*(len(value)-1) + "%f\n"
    write_rdb(CCF_folder_outputs+filename+'.rdb',data,keys,format)

    ########################################################################
    ################## Figures #############################################
    ########################################################################

    # Plot the solar CCFs used in the simulation
    figure(0)
    title('solar CCFs')
    plot(rv_ccf,intensity_ccf,color='b',lw=3,label='Quiet photosphere')
    plot(rv_ccf_magn_region,intensity_ccf_magn_region,color='r',ls='--',lw=3,label='Spot')
    ylim(0.4,1.05)
    xlim(-20,20)
    ylabel('Normalized flux')
    xlabel('RV [km/s]')
    legend(loc=4)
    subplots_adjust(top=0.93,left=0.1,right=0.96,bottom=0.13)
    savefig(CCF_folder_outputs+'FTS_solar_CCFs.pdf')

    majorFormatter = ScalarFormatter()
    majorFormatter.set_powerlimits((-8,8))
    majorFormatter.set_useOffset(0)

    # Plot the Flux, RV, BIS SPAN and FWHM variations induced by the active regions defined in "config.cfg"
    figure(1,[16,9])
    title('')
    subplot(511)
    plot(PSI,flux,color='r')
    ylabel('Norm. Flux',size=16)
    ax=gca()
    setp(ax.get_xticklabels(), visible=False)
    ax.yaxis.set_major_formatter(majorFormatter)
    subplot(512,sharex=ax)
    plot(PSI,rv_flux*1000,color='b',label='flux')
    plot(PSI,rv_bconv*1000,color='g',label='conv. blue.')
    plot(PSI,rv_tot*1000,color='r',label='tot')
    legend()
    ylabel('RV [m/s]',size=16)
    ax2=gca()
    setp(ax2.get_xticklabels(), visible=False)
    subplot(513,sharex=ax)
    plot(PSI,span_flux*1000,color='b')
    plot(PSI,span_bconv*1000,color='g')
    plot(PSI,span_tot*1000,color='r')
    ylabel('Bis Span [m/s]',size=16)
    ax3=gca()
    setp(ax3.get_xticklabels(), visible=False)
    subplot(514,sharex=ax)
    plot(PSI,fwhm_flux*1000,color='b',label='flux')
    plot(PSI,fwhm_bconv*1000,color='g',label='conv. blue.')
    plot(PSI,fwhm_tot*1000,color='r',label='tot')
    ylabel('Fwhm [m/s]',size=16)
    ax = gca()
    ax.yaxis.set_major_formatter(majorFormatter)
    xlabel('Phase',size=16)
    subplots_adjust(top=0.98,left=0.07,right=0.98,bottom=0.07,hspace=0.15)

    text1 =   '{0:19} = {1:<8d}, {2:23} = {3:<8d}\n'.format('Grid reso', GRID, 'Circumference reso', NRHO) + \
              '{0:19} = {1:<8d}, {2:23} = {3:<8d}\n'.format('Instr reso', INST_RESO, 'Radius of the sun [km]', RAD_Sun) + \
              '{0:19} = {1:<8.2f}, {2:23} = {3:<8.2f}\n'.format('Star Rot Period [d]', STAR.prot, 'Stellar radius [Rsun]', STAR.rad/RAD_Sun) + \
              '{0:19} = {1:<8.1f}, {2:23} = {3:<8.2f}\n'.format('Star incli', STAR.incl, 'Stellar vsini [m/s]', STAR.vrot) + \
              '{0:19} = {1:<8d}, {2:23} = {3:<8d}\n'.format('Stellar Teff [K]', STAR.Temp, 'Tdiff spot-photo [K]', STAR.Temp_diff_spot) + \
              '{0:19} = {1:<8.3f}, {2:23} = {3:<8.3f}'.format('Limb-dark lin', STAR.limba1, 'Limb-dark quad', STAR.limba2)

    text2 = 'Act Reg {0}: {1:<6}, lon = {2:<5.1f}, lat = {3:<4.1f}, size = {4:<5.3f} [Rsun] (= {5:.2f}%)\n'.format('1',active_region1.active_region_type_str,active_region1.long,active_region1.lat,active_region1.s,active_region1.s**2/2.*100) +\
            'Act Reg {0}: {1:<6}, lon = {2:<5.1f}, lat = {3:<4.1f}, size = {4:<5.3f} [Rsun] (= {5:.2f}%)\n'.format('2',active_region2.active_region_type_str,active_region2.long,active_region2.lat,active_region2.s,active_region2.s**2/2.*100) +\
            'Act Reg {0}: {1:<6}, lon = {2:<5.1f}, lat = {3:<4.1f}, size = {4:<5.3f} [Rsun] (= {5:.2f}%)\n'.format('3',active_region3.active_region_type_str,active_region3.long,active_region3.lat,active_region3.s,active_region3.s**2/2.*100) +\
            'Act Reg {0}: {1:<6}, lon = {2:<5.1f}, lat = {3:<4.1f}, size = {4:<5.3f} [Rsun] (= {5:.2f}%)'.format('4',active_region4.active_region_type_str,active_region4.long,active_region4.lat,active_region4.s,active_region4.s**2/2.*100)

    figtext(0.015,0.03,text1, name='Bitstream Vera Sans Mono',size=13) #name='Courier' very important because monotype font (all the characters have the same space)
    figtext(0.485,0.076,text2, name='Bitstream Vera Sans Mono',size=13) #name='Courier' very important because monotype font (all the characters have the same space)

    savefig(CCF_folder_outputs+filename+'.pdf')
#
## If no active region has been activated in the config.cfg file
except IndexError:
    print
    print 'ERROR: There is no active region selected for the star, please select at least one in the config.cfg file'
    print

show()