#ifndef MY_CONSTANTS_H
#define MY_CONSTANTS_H

// Include necessary ROOT headers
#include <TMath.h>

// define HCal acceptances (2024-10-28) 
const double eta_min_nhcal = -4.05; // old 2024-07-16: -4.14
const double eta_max_nhcal = -1.2; // old 2024-07-16: -1.18
const double eta_min_bhcal = -1.2;
const double eta_max_bhcal = 1.18;
const double eta_min_lfhcal = 1.18;			    
const double eta_max_lfhcal = 4.2;

// from $DETECTOR_PATH/compact/definitions.xml:
const double z_nhcal_min = -395; // (2024-12-03) start of nHCal in z-direction [cm]
const double z_nhcal_thickness = 45; // nHCal thickness in z [cm]
const double z_nhcal_max = z_nhcal_min - z_nhcal_thickness;

// define physical constants:
const double speedoflight = 299792458; // speed of light in m/s
const double kpmlifetime = 0.00000001238; // charged kaon life time in s
const double kpmmass = 0.493677; // charged kaon mass in GeV
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
