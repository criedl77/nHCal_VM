#ifndef MY_CONSTANTS_H
#define MY_CONSTANTS_H
#include <TMath.h>

// define HCal acceptances (2024-10-28) 
const double eta_min_nhcal = -4.05; // old 2024-07-16: -4.14
const double eta_max_nhcal = -1.2; // old 2024-07-16: -1.18
const double eta_min_bhcal = -1.2;
const double eta_max_bhcal = 1.18;
const double eta_min_lfhcal = 1.18;			    
const double eta_max_lfhcal = 4.2;

// nHCal design from which z_thickness_nhcal is built below
const int nlayers_nhcal=10; // 2025-11-14
const double steel_thickness_nhcal=4; // [cm] 2025-11-14
const double scintillator_thickness_nhcal=2.4; // [cm] 2025-11-14
const double extra_z_nhcal = 1.; // [cm]
const double tilesize_nhcal=10; // [cm] 2025-11-19
const double tilesize_lfhcal=5; // [cm] 2025-11-19

// from https://github.com/eic/epic/blob/main/compact/definitions.xml :
const double BackwardServiceGap_zmin = 316.; // [cm]
const double BackwardServiceGap_length = 15.; // [cm]
const double GapToBackwardServiceGap = 44.; // [cm], (called "gap between service gap" in definitions.xml,) includes 28.5cm oculus

const double z_min_nhcal = -1*(BackwardServiceGap_zmin + BackwardServiceGap_length + GapToBackwardServiceGap); // start of nHCal in z-direction [cm]
const double z_thickness_nhcal = nlayers_nhcal * (steel_thickness_nhcal + scintillator_thickness_nhcal) + extra_z_nhcal; // nHCal thickness in z [cm]
const double z_max_nhcal = z_min_nhcal - z_thickness_nhcal; // what should be done is that z_max_nhcal is fixed (the pivot) and z_min_nhcal adapted accordingly. Also the gap to the gap, GapToBackwardServiceGap, shouldn't be hardcoded like this.

// define max and min for histograms (need mm):
const double hx_min_nhcal = -2800.;
const double hx_max_nhcal = 2800.;
const double hy_min_nhcal = -2800.;
const double hy_max_nhcal = 2800.;
const double hz_min_nhcal = (z_min_nhcal*10)*0.95; //extends to more positive by 5% (0.95 because negative number) than the actual calo to display the bin edge in a well-defined way
const double hz_max_nhcal = (z_max_nhcal*10)*1.05; // extends to more negative by 5% than the actual calo to display the bin edge in a well-defined way

// define physical constants:
const double speedoflight = 299792458;        // speed of light in m/s
const double kpmlifetime = 0.00000001238;     // charged kaon life time in s
const double kpmmass = 0.493677;              // charged kaon mass in GeV
const double kPi = TMath::Pi();               // Value of Pi
const double kElectronMass = 9.10938356e-31;  // Electron mass in kg
const double kProtonMass = 1.6726219e-27;     // Proton mass in kg

const char* nhcal_name = "nhcal";
const char* bhcal_name = "bhcal";
const char* lfhcal_name = "lfhcal";

//From Dhruv - accepts a string cal_name containing the desired calorimeter name and returns true if an eta particle_eta is within acceptance for the given calorimeter
bool calo_eta_acceptance(const char* cal_name, float particle_eta) {
    if (cal_name == nhcal_name) {
        return (particle_eta >=  eta_min_nhcal && particle_eta <  eta_max_nhcal);
    } else if (cal_name == bhcal_name) {
        return (particle_eta >= eta_min_bhcal && particle_eta < eta_max_bhcal);
    } else if (cal_name == lfhcal_name) {
        return (particle_eta >= eta_min_lfhcal && particle_eta < eta_max_lfhcal);
    }
    return false;
}

// Example: Define some useful functions
inline double ConvertToGeV(double massInKg) {
  return massInKg * (TMath::Power(10, 9) / speedoflight); // Convert kg to GeV
}


#endif // MY_CONSTANTS_H
