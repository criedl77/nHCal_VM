nHCal_VM_Analysis.C reads one or more MC files and produces a root file output*.root

nHCal_VM_Plotting.C then picks up this output*.root file and produces plots. 

TODO: add description how to stream from a runlist, etc.

TODO: introduce flags to easily switch between local and streaming runlist; between different MC productions; etc.

LOGS:

CKR 2024-08-20: clean up strang definitions (not complete yet)

CKR 2024-10-03: continue cleaning up; remove option of streaming 1 file - use runlist always (it's easy to generate a runlist with 1 file)

CKR 2024-10-28: continue cleaning up

CKR 2024-11-12: set up this code on github and set it up both on my laptop and on SDCC

CKR 2024-12-03: calculate decay length 
