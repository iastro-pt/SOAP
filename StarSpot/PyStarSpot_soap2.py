from numpy import array, zeros, where, arange, searchsorted, r_, concatenate, mean, pi, sqrt, exp, log, ones, append
import starspot as stsp
from scipy import signal

#################################################
######               Classes              #######
#################################################

class Star:
    def __init__(self, prot, vrot, incl, limba1, limba2, psi, rad_sun, rad, Temp, Temp_diff_spot):
        self.prot           = prot           # Rotational period [days]
        self.vrot           = vrot           # Rotational velocity
        self.incl           = incl           # Inclination
        self.limba1         = limba1         # linear limb-darkening coeff.
        self.limba2         = limba2         # quadratic limb-darkening coeff.
        self.psi            = psi            # Phase
        self.rad_sun        = rad_sun        # Sun radius [km]
        self.rad            = rad            # Radius [Msol]
        self.Temp           = Temp           # Effective temperature of the star
        self.Temp_diff_spot = Temp_diff_spot # difference in temperature between the spot temperature and the effective temperature (<0 for spots, > 1 for plages)

class Active_region:
    def __init__(self, long, lat, size, active_region_type, check=False):
        self.long = long    # Longitude
        self.lat  = lat     # Latitude
        self.s    = size    # Size
        self.active_region_type = active_region_type # 0 if spot, 1 if plage
        if self.active_region_type == 0:
            self.active_region_type_str = 'spot' #active region type in string
        else:
            self.active_region_type_str = 'plage'  #active region type in string
        self.check = check                           # 1 if the active region is considered, 0 if no active region

    def calc_maps(self, star, grid, nrho):
        xyz  = stsp.spot_init(self.s, self.long, self.lat, star.incl, nrho) # Calculates the initial coordinates of the active region circumference
        xyz2 = stsp.spot_phase(xyz, star.incl, nrho, star.psi) # Calculates the coordinates of the active region circumference at phase psi
        self.visible, self.iminy, self.iminz, self.imaxy, self.imaxz = stsp.spot_area(xyz2, nrho, grid) # Check is the active region is visible at phase psi
                                                                                                        # and return the min and max y and z values for the circumference

class Ccf:
    def __init__(self, rv, intensity, width, step, star):
        self.rv         = rv             # Radial velocity of each point of the CCF
        self.intensity  = intensity      # Intensity of each point of the CCF
        self.width      = width          # Width of the CCF
        self.step       = step           # Sampling of the CCF
        self.n          = len(self.rv)   # Number of points of the CCF
        # The CCF has a width of self.width. This CCF has a rotational velocity of 0 because taken in the stellar disk center.
        # To consider rotation when simulating the emerging spectrum on the limb of the star (velocity different from 0), We have to extrapolate
        # outside of the range [-self.width;self.width]. This is done by calculating self.n_v which is the number of point of the CCF required
        # to consider a stellar rotation of star.vrot.
        # self.n_v must by an odd integer so that we have as many values on the positive side of the CCF than on the negative one (because 0 is present).
        if round(((1.1*star.vrot)/self.step)*2)%2 == 0:
            self.n_v        = int(self.n + round(((1.1*star.vrot)/self.step)*2))
        else:
            self.n_v        = int(self.n + round(((1.1*star.vrot)/self.step)*2)+1)
        # self.v_interval gives the difference in radial velocity between the minimum and maximum values of the CCF, once the extrapolation is done
        # self.n_v gives the number of points of the extrapolated CCF. (self.n_v-1) gives the number of intervals, which is then multiplied by the
        # sampling of the CCF self.step to give self.v_interval
        self.v_interval = self.step*(self.n_v-1) / 2.

#################################################
######            Functions               #######
#################################################

def compute_bis(RV,CCF):
    err = ones(len(CCF),'d')
    [mod,c,k,v0,fwhm,sig_c,sig_k,sig_v0,sig_fwhm,span,bis,depth] = stsp.gauss_bis(RV,CCF,err)
    depth=-k;
    return depth,bis,span,v0,fwhm


def Calculate_activity_signal(nrho, grid, star, psi, ccf, ccf_active_region, active_regions, inst_reso):

    # Calculates the flux and CCF in each cell of the grid and integrate
    # over the entire stellar disc to have the integrated flux (FLUXstar) and CCF (CCFstar).
    # Calculate the CCF (CCFstar) and the total flux (FLUXstar) for the quiet star
    CCFstar_quiet, FLUXstar_quiet = stsp.itot(star.vrot, star.incl, star.limba1, star.limba2, 0., 0., 0., grid, ccf.rv, ccf.intensity, ccf.v_interval, ccf.n_v, ccf.n)
    
    # initialization of the variables
    FLUXstar      = FLUXstar_quiet
    CCFstar_flux  = CCFstar_quiet
    CCFstar_bconv = CCFstar_quiet
    CCFstar_tot   = CCFstar_quiet

    #Loop over the different active regions
    for active_region in active_regions:
        if not active_region.check: continue # active_region.check = 1 if the active region is activated, 0 otherwise
        
        # Calculates the position of the spot initialized at the disc center
        xyz = stsp.spot_init(active_region.s, active_region.long, active_region.lat, star.incl, nrho)
        
        # Scans the yz-area where the spot is for different phases (psi) and
        # returns the spot's "non-contribution" to the total flux and its
        # "non-contribution" to the ccf, for each phase.
        # Thus the result is to be subtracted to the output of the itot() function.
        CCFactive_region_flux,CCFactive_region_bconv,CCFactive_region_tot,FLUXactive_region,\
                    = stsp.spot_scan_npsi(xyz, nrho, psi, len(psi), star.vrot, star.incl, star.limba1, star.limba2, 0., 0., 0., grid, ccf.rv, ccf.intensity,\
                                            ccf_active_region.intensity, ccf.v_interval, ccf.n_v, ccf.n, active_region.s, active_region.long,\
                                            active_region.lat, active_region.active_region_type, star.Temp, star.Temp_diff_spot)

        FLUXstar      = FLUXstar-FLUXactive_region           # Calculate the total flux of the star affected by active regions
        CCFstar_flux  = CCFstar_flux-CCFactive_region_flux   # Calculate the CCF of the star affected by the flux effect of active regions
        CCFstar_bconv = CCFstar_bconv-CCFactive_region_bconv # Calculate the CCF of the star affected by the convective blueshift effect of active regions
        CCFstar_tot   = CCFstar_tot-CCFactive_region_tot     # Calculate the CCF of the star affected by the total effect of active regions

    #Calculate where the extrapolated CCF corresponds to the boundaries of the non-extrapolated CCF.
    istart = (ccf.n_v-len(ccf.rv))/2 # ccf.n_v is odd given our definition, like len(ccf.rv), therefore the difference can be divided by 2. The CCF has
                                     # been extrapolated the same way on each boundary, so dividing the difference between the extrapolated and non
                                     # extrapolated CCF by 2 gives by how many points the CCF was extrapolated on each side
    iend   = istart+len(ccf.rv)

    # truncate the extrapolated rotating CCF to the same interval as the non-rotating CCF
    CCFstar_flux  = CCFstar_flux[:,istart:iend]
    CCFstar_bconv = CCFstar_bconv[:,istart:iend]
    CCFstar_tot   = CCFstar_tot[:,istart:iend]
    CCFstar_quiet = CCFstar_quiet[istart:iend]

    # Normalization
    CCFstar_flux  = array([f/max(f) for f in CCFstar_flux], dtype='d')
    CCFstar_bconv = array([f/max(f) for f in CCFstar_bconv], dtype='d')
    CCFstar_tot   = array([f/max(f) for f in CCFstar_tot], dtype='d')
    CCFstar_quiet = CCFstar_quiet/max(CCFstar_quiet)
    
    #Convolution with a Gaussian instrumental profile of a given resolution, given by the "instrument_reso" variable in the "config.cfg" file
    if inst_reso != 0:
        
        c = 299792458. # speed of light in m/s
        HARPS_resolution         = inst_reso # resolution R = lambda / Delta(lambda)
        HARPS_inst_profile_FWHM  = c/HARPS_resolution/1000. # Resolution = c/Delta_v -> Delta_v = c/R
        HARPS_inst_profile_sigma = HARPS_inst_profile_FWHM/(2*sqrt(2*log(2)))
        Gaussian_low_reso        = exp(-ccf.rv**2/(2*(HARPS_inst_profile_sigma)**2))
        
        CCFstar_quiet_tmp = signal.convolve(-CCFstar_quiet+1,Gaussian_low_reso,'same')
        CCFstar_quiet     = 1-CCFstar_quiet_tmp*(1-min(CCFstar_quiet))/max(CCFstar_quiet_tmp) # normalization
        
        #Convolution of the CCF with the Gaussian instrumental profile to reduce the resolution of the CCF
        for i in arange(len(CCFstar_flux)):
            CCFstar_flux_tmp  = signal.convolve(-CCFstar_flux[i]+1,Gaussian_low_reso,'same')
            CCFstar_flux[i]   = 1-CCFstar_flux_tmp*(1-min(CCFstar_flux[i]))/max(CCFstar_flux_tmp)
            CCFstar_bconv_tmp = signal.convolve(-CCFstar_bconv[i]+1,Gaussian_low_reso,'same')
            CCFstar_bconv[i]  = 1-CCFstar_bconv_tmp*(1-min(CCFstar_bconv[i]))/max(CCFstar_bconv_tmp)
            CCFstar_tot_tmp   = signal.convolve(-CCFstar_tot[i]+1,Gaussian_low_reso,'same')
            CCFstar_tot[i]    = 1-CCFstar_tot_tmp*(1-min(CCFstar_tot[i]))/max(CCFstar_tot_tmp)

    # Calculate the CCF parameters depth, BIS SPAN, RV and FWHM
    depth_flux,span_flux,rv_flux,fwhm_flux     = [],[],[],[]
    depth_bconv,span_bconv,rv_bconv,fwhm_bconv = [],[],[],[]
    depth_tot,span_tot,rv_tot,fwhm_tot         = [],[],[],[]
    for i in arange(len(CCFstar_tot)):
        #Flux effect
        CCF_prop_flux = compute_bis(ccf.rv,CCFstar_flux[i])
        depth_flux = append(depth_flux,CCF_prop_flux[0])
        span_flux = append(span_flux,CCF_prop_flux[2])
        rv_flux = append(rv_flux,CCF_prop_flux[3])
        fwhm_flux = append(fwhm_flux,abs(CCF_prop_flux[4]))
        
        #Convective Blueshift effect
        CCF_prop_bconv = compute_bis(ccf.rv,CCFstar_bconv[i])
        depth_bconv = append(depth_bconv,CCF_prop_bconv[0])
        span_bconv = append(span_bconv,CCF_prop_bconv[2])
        rv_bconv = append(rv_bconv,CCF_prop_bconv[3])
        fwhm_bconv = append(fwhm_bconv,abs(CCF_prop_bconv[4]))
        
        #Total effect
        CCF_prop_tot = compute_bis(ccf.rv,CCFstar_tot[i])
        depth_tot = append(depth_tot,CCF_prop_tot[0])
        span_tot = append(span_tot,CCF_prop_tot[2])
        rv_tot = append(rv_tot,CCF_prop_tot[3])
        fwhm_tot = append(fwhm_tot,abs(CCF_prop_tot[4]))

    #Normalize the flux of the star
    FLUXstar /= FLUXstar_quiet
    
    #return the variables
    return FLUXstar_quiet, CCFstar_quiet, FLUXstar, CCFstar_flux, CCFstar_bconv, CCFstar_tot, rv_flux, rv_bconv, rv_tot, span_flux, span_bconv, span_tot, fwhm_flux, fwhm_bconv, fwhm_tot, depth_flux, depth_bconv, depth_tot

    
