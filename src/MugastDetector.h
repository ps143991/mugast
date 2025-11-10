#ifndef MugastDetector_h
#define MugastDetector_h

#include <map>
#include <stdlib.h>
#include <vector>

#include "Math/Vector3D.h"
#include "MugastData.h"
#include "MugastPhysics.h"
#include "NPApplication.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"

// ROOT
#include "TGraphErrors.h"
#include "TH1.h"
#include "TObject.h"
#include "TVector2.h"
#include "TVector3.h"

#include "TRandom3.h"
// include xyz vector from root
using ROOT::Math::XYZVector;

// MFM
#ifdef MFM_FOUND
#include "MFMAllFrames.h"
#endif

#ifdef Geant4_FOUND
#include "NPG4VDetector.h"
#endif

namespace mugast {

class MugastSpectra;
class MugastDetector : public nptool::VDetector {
public: // Constructor and Destructor
  MugastDetector();
  ~MugastDetector();

public: // Data member
  mugast::MugastData*                    m_RawData;
  mugast::MugastData*                    m_CalData;
  mugast::MugastPhysics*                 m_PhysicsData;
  nptool::CalibrationManager             m_Cal;
  std::shared_ptr<mugast::MugastSpectra> m_Spectra;

public: // inherrited from nptool::VPlugin
  std::vector<std::string> GetDependencies() { return {"root"}; };

public: // inherrited from nptool::VDetector
        //  Read stream at ConfigFile to pick-up parameters of detector
        //  (Position,...) using Token
  void ReadConfiguration(nptool::InputParser);

  //   Initialize the standard parameter for analysis
  //   ie: all channel enable, maximum multiplicity for strip = number of
  //   telescope
  void InitializeStandardParameter();

  //   Read the user configuration file; if no file found, load standard one
  void ReadAnalysisConfig();

  //  Add Parameter to the CalibrationManger
  void AddParameterToCalibrationManager() {};

  void InitializeDataInputConversion(std::shared_ptr<nptool::VDataInput>){};

  //  Activated associated Branches and link it to the private member
  //  DetectorData address In this method mother Branches (Detector) AND
  //  daughter leaf (fDetector_parameter) have to be activated
  void InitializeDataInputRaw(std::shared_ptr<nptool::VDataInput>);

  //  Activated associated Branches and link it to the private member
  //  DetectorPhysics address In this method mother Branches (Detector) AND
  //  daughter leaf (parameter) have to be activated
  void InitializeDataInputPhysics(std::shared_ptr<nptool::VDataInput>);

  //  Create associated branches and associated private member DetectorPhysics
  //  address
  void InitializeDataOutputRaw(std::shared_ptr<nptool::VDataOutput>);

  //  Create associated branches and associated private member DetectorPhysics
  //  address
  void InitializeDataOutputPhysics(std::shared_ptr<nptool::VDataOutput>);

  //  This method is called at each event of a conversion to fill in the data
  //  class
  void BuildRawEvent(const std::string& daq_name, const std::string& label,
                     void* commonframe);
  void TreatFrame(void* commonframe);
  void ReadActionFile(std::string action_file);
  bool faction_file_initialized;
  std::map<int, std::string> fLabelMap;

  //  This method is called at each event read from the Input Tree. Aim is to
  //  build treat Raw dat in order to extract physical parameter.
  void BuildPhysicalEvent();
  //   Remove bad channel, calibrate the data and apply threshold
  void PreTreat();
  void ClearPreTreatedData() { m_CalData->Clear(); }


  //  Those two method all to clear the Event Physics or Data
  void ClearEventPhysics() { m_PhysicsData->Clear(); };
  void ClearEventData() { m_RawData->Clear(); };
  // Method related to the TSpectra classes, aimed at providing a framework for
  // online applications Instantiate the Spectra class and the histogramm
  // throught it
  void InitSpectra();
  // Fill the spectra hold by the spectra class
  void FillSpectra();
  // Write the spectra to a file
  void WriteSpectra();
  // Used for Online mainly, perform check on the histo and for example change
  // their color if issues are found
  void CheckSpectra();
  // Used for Online only, clear all the spectra hold by the Spectra class
  void ClearSpectra();
  // Used for interoperability with other framework
  void SetRawDataPointer(void*) {};


/******************************************Functions related to MUGAST now*****************************************/

 private:
  map<int, MG_DetectorType> DetectorType; //!


  //   Return false if the channel is disabled by user
  //   Frist argument is either 0 for X,1 Y,2 SecondLayer 3
  bool IsValidChannel(const int& Type, const int& telescope, const int& channel);

  //   Initialize the standard parameter for analysis
  //   ie: all channel enable, maximum multiplicity for strip = number of
  //   telescope

  //   Read the user configuration file; if no file found, load standard one

  //   Add a Telescope using Corner Coordinate information
  void AddTelescope(MG_DetectorType type, XYZVector C_X1_Y1, XYZVector C_X128_Y1, XYZVector C_X1_Y128,
                    XYZVector C_X128_Y128);

  //   Add a Telescope using R Theta Phi of Si center information
  void AddTelescope(MG_DetectorType type, double theta, double phi, double distance, double beta_u, double beta_v,
                    double beta_w);

  //   Special Method for Annular S1
  void AddTelescope(MG_DetectorType type, XYZVector C_Center);

  // Give and external TMustData object to TMugastPhysics. Needed for online
  // analysis for example.



  // Use to access the strip position
  double GetStripPositionX(const int N, const int X, const int Y) {
    // if (N==9)
    // cout << N << " " << X << " " << Y << " " << m_DetectorNumberIndex[N] << " " << m_StripPositionX[
    // m_DetectorNumberIndex[N] - 1][X - 1][Y - 1] << endl;
    return m_StripPositionX[m_DetectorNumberIndex[N] - 1][X - 1][Y - 1];
  };
  double GetStripPositionY(const int N, const int X, const int Y) {
    return m_StripPositionY[m_DetectorNumberIndex[N] - 1][X - 1][Y - 1];
  };
  double GetStripPositionZ(const int N, const int X, const int Y) {
    return m_StripPositionZ[m_DetectorNumberIndex[N] - 1][X - 1][Y - 1];
  };

  double GetNumberOfTelescope() const { return m_NumberOfTelescope; };

  // To be called after a build Physical Even
  int GetEventMultiplicity() const { return m_PhysicsData->EventMultiplicity; };

  double GetEnergyDeposit(const int i) const { return m_PhysicsData->TotalEnergy[i]; };

  TVector3 GetPositionOfInteraction(const int i, bool random = false);
  TVector3 GetTelescopeNormal(const int i);


 private: //   Parameter used in the analysis
  // Shape of the detector Trapezoid or Square
  map<int, int> m_DetectorNumberIndex;

  // By default take EX and TY.
  bool m_Take_E_Y; //!
  bool m_Take_T_Y; //!

  //   Event over this value after pre-treatment are not treated / avoid long
  //   treatment time on spurious event
  unsigned int m_MaximumStripMultiplicityAllowed; //!
  //   Give the allowance in percent of the difference in energy between X and Y
  double m_StripEnergyMatching; //!

  // Raw Threshold
  int m_DSSD_X_E_RAW_Threshold;      //!
  int m_DSSD_Y_E_RAW_Threshold;      //!
  int m_SecondLayer_E_RAW_Threshold; //!

  // Calibrated Threshold
  double m_DSSD_X_E_Threshold;      //!
  double m_DSSD_Y_E_Threshold;      //!
  double m_SecondLayer_E_Threshold; //!


 private:                                            //   Map of activated channel
  map<int, vector<bool>> m_XChannelStatus;           //!
  map<int, vector<bool>> m_YChannelStatus;           //!
  map<int, vector<bool>> m_SecondLayerChannelStatus; //!

 private:                  // Spatial Position of Strip Calculated on bases of detector position
  int m_NumberOfTelescope; //!

  vector<vector<vector<double>>> m_StripPositionX; //!
  vector<vector<vector<double>>> m_StripPositionY; //!
  vector<vector<vector<double>>> m_StripPositionZ; //!

 private:
  map<int, int> m_HitDSSDX; //!
  map<int, int> m_HitDSSDY; //!


  //   DSSD
  //   X
  double fDSSD_X_E(/* const TMugastData* Data, */ const int& i);
  double fDSSD_X_T(/* const TMugastData* Data, */ const int& i);

  //   Y
  double fDSSD_Y_E(/* const TMugastData* Data, */ const int& i);
  double fDSSD_Y_T(/* const TMugastData* Data, */ const int& i);

  //  Second Layer
  double fSecondLayer_E(/* const TMugastData* Data, */ const int& i);
  double fSecondLayer_T(/* const TMugastData* Data, */ const int& i);


public:
  TRandom3* m_random; //!

  vector<TVector2> Match_X_Y();
  bool Match_SecondLayer(int X, int Y, int StripNbr);
  bool ResolvePseudoEvent();
  int CheckEvent();
  bool m_reader = true;
  double GetNumberOfTelescopes() const { return m_NumberOfTelescope; };









  // simulation
public:
  void        InitSimulation(std::string simtype);
  void        ConstructGeometry();
  std::string m_simtype;

#ifdef Geant4_FOUND
private:
  std::shared_ptr<nptool::geant4::VDetector> m_Geant4;
#endif
};

} // namespace mugast
#endif
