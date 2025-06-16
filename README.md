*** nHCal_VM_Analysis.C reads one or more MC files and produces a root file out.strang-nHCal_VM.root

- usage:
  
root -l -q 'nHCal_VM_Analysis.C(tracking,  mode, "strang")'
  
tracking=
	
1: tracking available (ReconstructedChargedParticles is filled) = default
			
0: no tracking available
			
mode=
			
1: streaming runlist (SDCC via Jlab) = default
			
2: local runlist
			
strang is the data production string, for example=
			
nhcal_only_tile5cm_absorber4cm_scintillator0.4cm_10layers_mu-_p1gev_phi45_theta170_10events


- EXAMPLE of an output file: out.nhcal_only_tile5cm_absorber4cm_scintillator0.4cm_10layers_mu-_p1gev_phi45_theta170_10events-nHCal_VM.root

*** nHCal_VM_Plotting.C takes the passed "strang", picks up the according out.strang-nHCal_VM.root file, and produces plots in the subdirectory "strang". 

- usage:

root -l -q 'nHCal_VM_Plotting.C("strang")'

EXAMPLE: root -l -q 'nHCal_VM_Plotting.C("nhcal_only_tile5cm_absorber4cm_scintillator0.4cm_10layers_mu-_p1gev_phi45_theta170_10events")'
 

LOGS:


CKR 2024-08-20: clean up strang definitions 

CKR 2024-10-03: continue cleaning up; remove option of streaming 1 file - use runlist always (it's easy to generate a runlist with 1 file)

CKR 2024-10-28: continue cleaning up

CKR 2024-11-12: set up this code on github and set it up both on my laptop and on SDCC

CKR 2024-12-03: calculate decay length 

CKR June 2025 - read cluster and hit information
