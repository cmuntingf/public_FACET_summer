# E320 Laser and Beam Parameters
import numpy as np

# Laser parameters
laserw0 = 20 # laser waist in mu meter
pulselengthFWHM = 48e-3 #[ps] FWHM measured in power
LASER_PARAMETERS = {
        'angle_deg':28.1,                # in degree
        'pulseE':0.1,                # laser pulse energy [J]
        'sigLrx':laserw0/2,                 #/2; % given in [mu m] micro meter like 2 waist w0=28;
        'sigLry':laserw0/2,                 #/2; % given in [mu m] micro meter like 2 waist w0=28;
        'laserwl':800,                # laser wavelenth [nm] nano meters
        'sigt': pulselengthFWHM / np.sqrt(8*np.log(2)),      #RMS pulse length [ps]
        'shifting_laser_x' : 0,  #
        'shifting_laser_y' : 0,  #
        'shifting_laser_s' : 0,  #
        'shifting_laser_t' : 0,  #shifting_laser_t;  %

        'NPH':0, #is Maximum number of laser photons to be absorbed in one process NPH=0 linear, NPH>= 1, use nonlinear formula.

        'N_t_steps':250, #Number of time steps for linear can be 250-300 for non-linear bigger

        #STOKES:  linear [0 0 1]; circular [0 1 0]
        'STOKES_1':0,
        'STOKES_2':0,
        'STOKES_3':1,
    }
    
# Electron Beam Parameters
n_macro=5e3;               # Number of macroparticles
chargebunch = 1.6e-9;     # Charge per electrons bunch [C] pico->10^-12 femto->10^-15
n_electrons = chargebunch / (1.6e-19)
beam_energy_MeV=10e3;      # initial beam energy MeV
energy_spread=0.01;      # initial relative energy spread (not in [%])
norm_emit_x=164e-6;        # Normalised emittance x [m rad]
norm_emit_y=8e-6;        # Normalised emittance н [m rad]
sigma_e_x=46e-6;         # horizontal electron beam size [m]
sigma_e_y=23e-6;         # vertical electron beam size [m]
bunch_length=20e-6;    # initial bunch length [m]


Emass = 0.511e6;      #Electrons energy of rest
gamma=beam_energy_MeV*1e6/Emass
momentum_GeV = np.sqrt((beam_energy_MeV)**2 - (Emass/1e6)**2) / 1000 # Relativistic momentum in GeV/c

emit_x=norm_emit_x/np.sqrt(gamma**2-1)
emit_y=norm_emit_y/np.sqrt(gamma**2-1)

Betax=sigma_e_x**2/emit_x
Betay=sigma_e_y**2/emit_y


sigma_xp=(sigma_e_x/Betax)
sigma_yp=(sigma_e_y/Betay)
