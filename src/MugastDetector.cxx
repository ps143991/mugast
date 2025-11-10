#include "MugastDetector.h"
#include "MugastSpectra.h"
#include "TRandom3.h"
#include "NPFunction.h"
#include "NPRootPlugin.h"
#include <iostream>

// For Geant4 Simulation
// #ifdef Geant4_FOUND
// #include "MugastGeant4.h"
// #endif

// std
#include <dlfcn.h>

using namespace mugast;
using namespace std;
using namespace ROOT::Math;

////////////////////////////////////////////////////////////////////////////////
MugastDetector::MugastDetector() {
  m_RawData     = new mugast::MugastData();
  m_CalData     = new mugast::MugastData();
  m_PhysicsData = new mugast::MugastPhysics();
  m_Cal.InitCalibration();
  m_PhysicsData->EventMultiplicity = 0;

  m_random = new TRandom3();
  m_Spectra = NULL;
  m_NumberOfTelescope = 0;
  m_MaximumStripMultiplicityAllowed = 10;
  m_StripEnergyMatching = 0.050;
  // Raw Threshold
  m_DSSD_X_E_RAW_Threshold = 8200;
  m_DSSD_Y_E_RAW_Threshold = 8200;
  m_SecondLayer_E_RAW_Threshold = 8200;
  // Calibrated Threshold
  m_DSSD_X_E_Threshold = 0;
  m_DSSD_Y_E_Threshold = 0;
  m_SecondLayer_E_Threshold = 0;

  m_Take_E_Y = false;
  m_Take_T_Y = true;
}

MugastDetector::~MugastDetector() {}



//////////////////////////////////////////////////////////////////////////
void MugastDetector::PreTreat() {
  ClearPreTreatedData();
  static unsigned int DSSDX_EMult, DSSDY_EMult, SecondLayer_EMult;
  static unsigned int DSSDX_TMult, DSSDY_TMult, SecondLayer_TMult;
  DSSDX_EMult = m_RawData->GetDSSDXEMult();
  DSSDY_EMult = m_RawData->GetDSSDYEMult();
  DSSDX_TMult = m_RawData->GetDSSDXTMult();
  DSSDY_TMult = m_RawData->GetDSSDYTMult();
  SecondLayer_EMult = m_RawData->GetSecondLayerEMult();
  SecondLayer_TMult = m_RawData->GetSecondLayerTMult();
  MG_DetectorType type = MG_NOCHANGE;
  //   X
  //   E
  for (unsigned int i = 0; i < DSSDX_EMult; ++i) {
    type = DetectorType[m_RawData->GetDSSDXEDetectorNbr(i)];
    if (m_RawData->GetDSSDXEEnergy(i) > m_DSSD_X_E_RAW_Threshold &&
        IsValidChannel(0, m_RawData->GetDSSDXEDetectorNbr(i), m_RawData->GetDSSDXEStripNbr(i))) {
      double EX = fDSSD_X_E(/* m_RawData, */ i);
      if (EX > m_DSSD_X_E_Threshold)
        m_CalData->SetDSSDXE(type, m_RawData->GetDSSDXEDetectorNbr(i), m_RawData->GetDSSDXEStripNbr(i), EX);
    }
  }

  //   T
  for (unsigned int i = 0; i < DSSDX_TMult; ++i) {
    type = DetectorType[m_RawData->GetDSSDXTDetectorNbr(i)];
    if (IsValidChannel(0, m_RawData->GetDSSDXTDetectorNbr(i), m_RawData->GetDSSDXTStripNbr(i)))
      m_CalData->SetDSSDXT(type, m_RawData->GetDSSDXTDetectorNbr(i), m_RawData->GetDSSDXTStripNbr(i),
                                  fDSSD_X_T(/* m_RawData, */ i));
  }

  //   Y
  //   E
  for (unsigned int i = 0; i < DSSDY_EMult; ++i) {
    type = DetectorType[m_RawData->GetDSSDYEDetectorNbr(i)];
    if (m_RawData->GetDSSDYEEnergy(i) < m_DSSD_Y_E_RAW_Threshold &&
        IsValidChannel(1, m_RawData->GetDSSDYEDetectorNbr(i), m_RawData->GetDSSDYEStripNbr(i))) {
      double EY = fDSSD_Y_E(/* m_RawData, */ i);
      if (EY > m_DSSD_Y_E_Threshold)
        m_CalData->SetDSSDYE(type, m_RawData->GetDSSDYEDetectorNbr(i), m_RawData->GetDSSDYEStripNbr(i), EY);
    }
  }

  //   T
  for (unsigned int i = 0; i < DSSDY_TMult; ++i) {
    type = DetectorType[m_RawData->GetDSSDYTDetectorNbr(i)];
    if (IsValidChannel(1, m_RawData->GetDSSDYTDetectorNbr(i), m_RawData->GetDSSDYTStripNbr(i)))
      m_CalData->SetDSSDYT(type, m_RawData->GetDSSDYTDetectorNbr(i), m_RawData->GetDSSDYTStripNbr(i),
                                  fDSSD_Y_T(/* m_RawData, */ i));
  }

  //   SecondLayer
  //   E
  for (unsigned int i = 0; i < SecondLayer_EMult; ++i) {
    if (m_RawData->GetSecondLayerEEnergy(i) > m_SecondLayer_E_RAW_Threshold &&
        IsValidChannel(2, m_RawData->GetSecondLayerEDetectorNbr(i), m_RawData->GetSecondLayerEStripNbr(i))) {
      double ESecondLayer = fSecondLayer_E(/* m_RawData, */ i);
      if (ESecondLayer > m_SecondLayer_E_Threshold)
        m_CalData->SetSecondLayerE(m_RawData->GetSecondLayerEDetectorNbr(i),
                                          m_RawData->GetSecondLayerEStripNbr(i), ESecondLayer);
    }
  }

  //   T
  for (unsigned int i = 0; i < SecondLayer_TMult; ++i) {
    if (IsValidChannel(2, m_RawData->GetSecondLayerTDetectorNbr(i), m_RawData->GetSecondLayerTStripNbr(i)))
      m_CalData->SetSecondLayerT(m_RawData->GetSecondLayerTDetectorNbr(i),
                                        m_RawData->GetSecondLayerTStripNbr(i), fSecondLayer_T(/* m_RawData, */ i));
  }
  return;
}

///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////

void MugastDetector::BuildPhysicalEvent() {


  PreTreat();

  bool check_SecondLayer = false;
  static unsigned int DSSDXEMult, DSSDYEMult, DSSDXTMult, DSSDYTMult, SecondLayerEMult, SecondLayerTMult;
  DSSDXEMult = m_CalData->GetDSSDXEMult();
  DSSDYEMult = m_CalData->GetDSSDYEMult();
  DSSDXTMult = m_CalData->GetDSSDXTMult();
  DSSDYTMult = m_CalData->GetDSSDYTMult();
  SecondLayerEMult = m_CalData->GetSecondLayerEMult();
  SecondLayerTMult = m_CalData->GetSecondLayerTMult();

  // random->SetSeed(42);

  // srand(time(NULL));

  if (1 /*CheckEvent() == 1*/) {
    vector<TVector2> couple = Match_X_Y();

    m_PhysicsData->EventMultiplicity = couple.size();
    for (unsigned int i = 0; i < couple.size(); ++i) {
      check_SecondLayer = false;

      int N = m_CalData->GetDSSDXEDetectorNbr(couple[i].X());

      int X = m_CalData->GetDSSDXEStripNbr(couple[i].X());
      int Y = m_CalData->GetDSSDYEStripNbr(couple[i].Y());

      double DSSD_X_E = m_CalData->GetDSSDXEEnergy(couple[i].X());
      double DSSD_Y_E = m_CalData->GetDSSDYEEnergy(couple[i].Y());

      //  Search for associate Time
      double DSSD_X_T = -1000;
      for (unsigned int t = 0; t < DSSDXTMult; ++t) {
        if (m_CalData->GetDSSDXTStripNbr(couple[i].X()) == m_CalData->GetDSSDXTStripNbr(t) &&
            m_CalData->GetDSSDXTDetectorNbr(couple[i].X()) == m_CalData->GetDSSDXTDetectorNbr(t)) {
          DSSD_X_T = m_CalData->GetDSSDXTTime(t);
          break;
        }
      }

      double DSSD_Y_T = -1000;
      for (unsigned int t = 0; t < DSSDYTMult; ++t) {
        if (m_CalData->GetDSSDYTStripNbr(couple[i].Y()) == m_CalData->GetDSSDYTStripNbr(t) &&
            m_CalData->GetDSSDYTDetectorNbr(couple[i].Y()) == m_CalData->GetDSSDYTDetectorNbr(t)) {
          DSSD_Y_T = m_CalData->GetDSSDYTTime(t);
          break;
        }
      }

      m_PhysicsData->DSSD_X.push_back(X);
      m_PhysicsData->DSSD_Y.push_back(Y);
      m_PhysicsData->TelescopeNumber.push_back(N);

      // Randomize annular detector in Phi
      if (DetectorType[N] == MG_ANNULAR) {
        TVector3 Inter = GetPositionOfInteraction(i);
        double Phi = Inter.Phi();
        double Perp = Inter.Perp();
        double rPhi = m_random->Uniform(-11.25 / 180. * M_PI, 11.25 / 180. * M_PI);
        double rPerp = m_random->Uniform(-0.75, 0.75);

        Inter.SetPhi(Phi + rPhi);
        Inter.SetPerp(Perp + rPerp);
        m_PhysicsData->PosX.push_back(Inter.X());
        m_PhysicsData->PosY.push_back(Inter.Y());
        m_PhysicsData->PosZ.push_back(Inter.Z());
      }

      // No random for Trapezoid and Square
      else {
        m_PhysicsData->PosX.push_back(GetPositionOfInteraction(i).x());
        m_PhysicsData->PosY.push_back(GetPositionOfInteraction(i).y());
        m_PhysicsData->PosZ.push_back(GetPositionOfInteraction(i).z());
      }

      if (m_Take_E_Y) {
        m_PhysicsData->DSSD_E.push_back(DSSD_Y_E);
        m_PhysicsData->TotalEnergy.push_back(DSSD_Y_E);
      }
      else {
        m_PhysicsData->DSSD_E.push_back(DSSD_X_E);
        m_PhysicsData->TotalEnergy.push_back(DSSD_X_E);
      }

      if (m_Take_T_Y)
        m_PhysicsData->DSSD_T.push_back(DSSD_Y_T);
      else
        m_PhysicsData->DSSD_T.push_back(DSSD_X_T);

      for (unsigned int j = 0; j < SecondLayerEMult; ++j) {
        if (m_CalData->GetSecondLayerEDetectorNbr(j) == N) {
          if (Match_SecondLayer(X, Y, m_CalData->GetSecondLayerEStripNbr(j))) {
            m_PhysicsData->SecondLayer_N.push_back(m_CalData->GetSecondLayerEStripNbr(j));
            m_PhysicsData->SecondLayer_E.push_back(m_CalData->GetSecondLayerEEnergy(j));
            m_PhysicsData->SecondLayer_T.push_back(-1000);
            //   Look for associate Time
            for (unsigned int k = 0; k < SecondLayerTMult; ++k) {
              // Same DSSD, Same Detector
              if (m_CalData->GetSecondLayerEStripNbr(j) == m_CalData->GetSecondLayerTStripNbr(k) &&
                  m_CalData->GetSecondLayerEDetectorNbr(j) == m_CalData->GetSecondLayerTDetectorNbr(k)) {
                m_PhysicsData->SecondLayer_T[m_PhysicsData->SecondLayer_T.size() - 1] = m_CalData->GetSecondLayerTTime(j);
                break;
              }
            }

            check_SecondLayer = true;
          }
        }
      }

      if (!check_SecondLayer) {
        m_PhysicsData->SecondLayer_N.push_back(0);
        m_PhysicsData->SecondLayer_E.push_back(-1000);
        m_PhysicsData->SecondLayer_T.push_back(-1000);
      }
    } // loop on couples
  }   // if (CheckEvent)
  return;
}


///////////////////////////////////////////////////////////////////////////
int MugastDetector::CheckEvent() {
  // Check the size of the different elements
  if (m_CalData->GetDSSDXEMult() == m_CalData->GetDSSDYEMult())
    return 1; // Regular Event

  // INterstrip management is not coded, so waste of time to make this test
  /*  else if(   m_CalData->GetMMStripXEMult() ==
      m_CalData->GetMMStripYEMult()+1
      || m_CalData->GetMMStripXEMult() ==
      m_CalData->GetMMStripYEMult()-1  )
      return 2 ; // Pseudo Event, potentially interstrip*/

  else
    return -1; // Rejected Event
}

///////////////////////////////////////////////////////////////////////////
bool MugastDetector::ResolvePseudoEvent() { return false; }

///////////////////////////////////////////////////////////////////////////
vector<TVector2> MugastDetector::Match_X_Y() {
  vector<TVector2> ArrayOfGoodCouple;
  static unsigned int m_DSSDXEMult, m_DSSDYEMult;
  m_DSSDXEMult = m_CalData->GetDSSDXEMult();
  m_DSSDYEMult = m_CalData->GetDSSDYEMult();

  // Prevent code from treating very high multiplicity Event
  // Those event are not physical anyway and that improve speed.
  if (m_DSSDXEMult > m_MaximumStripMultiplicityAllowed || m_DSSDYEMult > m_MaximumStripMultiplicityAllowed) {
    return ArrayOfGoodCouple;
  }

  for (unsigned int i = 0; i < m_DSSDXEMult; i++) {
    for (unsigned int j = 0; j < m_DSSDYEMult; j++) {

      // Declaration of variable for clarity
      double DSSDXDetNbr = m_CalData->GetDSSDXEDetectorNbr(i);
      double DSSDYDetNbr = m_CalData->GetDSSDYEDetectorNbr(j);

      //   if same detector check energy
      if (DSSDXDetNbr == DSSDYDetNbr) {

        // Declaration of variable for clarity
        double DSSDXEnergy = m_CalData->GetDSSDXEEnergy(i);
        double DSSDXNbr = m_CalData->GetDSSDXEStripNbr(i);
        double DSSDYEnergy = m_CalData->GetDSSDYEEnergy(j);
        double DSSDYNbr = m_CalData->GetDSSDYEStripNbr(j);

        //   Look if energy match
        if (abs((DSSDXEnergy - DSSDYEnergy) / 2.) < m_StripEnergyMatching) {
          // Gives a unique ID for every telescope and strip combination
          int IDX = m_NumberOfTelescope * DSSDXNbr + DSSDXDetNbr;
          int IDY = m_NumberOfTelescope * DSSDYNbr + DSSDYDetNbr;

          m_HitDSSDX[IDX]++;
          m_HitDSSDY[IDY]++;

          ArrayOfGoodCouple.push_back(TVector2(i, j));
        }
      }
    }
  }

  // Prevent to treat event with ambiguous matching beetween X and Y
  map<int, int>::iterator itX = m_HitDSSDX.begin();
  for (; itX != m_HitDSSDX.end(); itX++) {
    if (itX->second > 1) {
      ArrayOfGoodCouple.clear();
    }
  }

  map<int, int>::iterator itY = m_HitDSSDY.begin();
  for (; itY != m_HitDSSDY.end(); itY++) {
    if (itY->second > 1) {
      ArrayOfGoodCouple.clear();
    }
  }

  m_HitDSSDX.clear();
  m_HitDSSDY.clear();

  return ArrayOfGoodCouple;
}

////////////////////////////////////////////////////////////////////////////
bool MugastDetector::IsValidChannel(const int& Type,
                                    // Uses raw channel number
                                    const int& telescope, const int& channel) {
  if (Type == 0) {
    return *(m_XChannelStatus[telescope].begin() + channel - 1);
  }
  else if (Type == 1) {
    return *(m_YChannelStatus[telescope].begin() + channel - 1);
  }

  else if (Type == 2)
    return *(m_SecondLayerChannelStatus[telescope].begin() + channel - 1);

  else
    return false;
}


///////////////////////////////////////////////////////////////////////////
void MugastDetector::ReadAnalysisConfig() {

  nptool::InputParser parser("./configs/ConfigMugast.dat", false);
  // vector<nptool::InputBlock*> blocks = parser.GetAllBlocksWithToken("ConfigMugast");
  auto blocks = parser.GetAllBlocksWithToken("MugastAnalysis");             // PS: changed to "auto"

  cout << endl << "//// Read MUGAST analysis configuration" << endl;

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasToken("MAX_STRIP_MULTIPLICITY"))
      m_MaximumStripMultiplicityAllowed = blocks[i]->GetInt("MAX_STRIP_MULTIPLICITY");

    if (blocks[i]->HasToken("STRIP_ENERGY_MATCHING"))
      m_StripEnergyMatching = blocks[i]->GetDouble("STRIP_ENERGY_MATCHING", "MeV");

    if (blocks[i]->HasToken("DISABLE_CHANNEL_X")) {
      vector<int> v = blocks[i]->GetVectorInt("DISABLE_CHANNEL_X");
      *(m_XChannelStatus[v[0]].begin() + v[1] - 1) = false;
    }

    if (blocks[i]->HasToken("DISABLE_CHANNEL_Y")) {
      vector<int> v = blocks[i]->GetVectorInt("DISABLE_CHANNEL_Y");
      *(m_YChannelStatus[v[0]].begin() + v[1] - 1) = false;
    }

    if (blocks[i]->HasToken("DISABLE_ALL")) {
      int telescope = blocks[i]->GetInt("DISABLE_ALL");
      vector<bool> ChannelStatus;
      ChannelStatus.resize(128, false);
      m_XChannelStatus[telescope] = ChannelStatus;
      m_YChannelStatus[telescope] = ChannelStatus;
      ChannelStatus.resize(16, false);
      m_SecondLayerChannelStatus[telescope] = ChannelStatus;
    }

    if (blocks[i]->HasToken("TAKE_E_Y"))
      m_Take_E_Y = blocks[i]->GetInt("TAKE_E_Y");

    if (blocks[i]->HasToken("TAKE_T_Y"))
      m_Take_T_Y = blocks[i]->GetInt("TAKE_T_Y");

    if (blocks[i]->HasToken("TAKE_E_X"))
      m_Take_E_Y = !(blocks[i]->GetInt("TAKE_E_X"));

    if (blocks[i]->HasToken("TAKE_T_X"))
      m_Take_T_Y = !(blocks[i]->GetInt("TAKE_T_X"));

    if (blocks[i]->HasToken("DSSD_X_E_RAW_THRESHOLD"))
      m_DSSD_X_E_RAW_Threshold = blocks[i]->GetInt("DSSD_X_E_RAW_THRESHOLD");

    if (blocks[i]->HasToken("DSSD_Y_E_RAW_THRESHOLD"))
      m_DSSD_Y_E_RAW_Threshold = blocks[i]->GetInt("DSSD_Y_E_RAW_THRESHOLD");

    if (blocks[i]->HasToken("SECONDLAYER_E_RAW_THRESHOLD"))
      m_SecondLayer_E_RAW_Threshold = blocks[i]->GetInt("SECONDLAYER_E_RAW_THRESHOLD");

    if (blocks[i]->HasToken("DSSD_X_E_THRESHOLD"))
      m_DSSD_X_E_Threshold = blocks[i]->GetDouble("DSSD_X_E_THRESHOLD", "MeV");

    if (blocks[i]->HasToken("DSSD_Y_E_THRESHOLD"))
      m_DSSD_Y_E_Threshold = blocks[i]->GetDouble("DSSD_Y_E_THRESHOLD", "MeV");

    if (blocks[i]->HasToken("SECONDLAYER_E_THRESHOLD"))
      m_SecondLayer_E_Threshold = blocks[i]->GetDouble("SECONDLAYER_E_THRESHOLD", "MeV");
  }
}

///////////////////////////////////////////////////////////////////////////
bool MugastDetector::Match_SecondLayer(int X, int Y, int StripNbr) {
  /*
     if (abs(m_CsI_MatchingX[CristalNbr - 1] - X) < (double)m_CsI_Size / 2.
     && abs(m_CsI_MatchingY[CristalNbr - 1] - Y) < (double)m_CsI_Size / 2.)
     return true;

     else
     return false;*/
  return true;
}



///////////////////////////////////////////////////////////////////////////
void MugastDetector::ReadConfiguration(nptool::InputParser parser) {
  // vector<nptool::InputBlock*> blocks = parser.GetAllBlocksWithToken("Mugast");
  auto blocks = parser.GetAllBlocksWithToken("Mugast");  //PS: becasue definition of GetAllBlocksWithToken changed in NP4

  // Cartesian Case
  vector<string> cart = {"DetectorNumber", "X001_Y001", "X001_Y128", "X128_Y001", "X128_Y128"};
  // Spherical Case
  vector<string> sphe = {"DetectorNumber", "R", "THETA", "PHI", "BETA"};
  // Annular Case
  vector<string> annular = {"DetectorNumber", "Center"};
  string Type;

  int number = 0;
  for (const auto& block : blocks) {
    if (block->HasTokenList(cart)) {

        Type = block->GetMainValue();
      nptool::message("green","","","Muagst Telescope: "+Type);   //PS: correct syntax??

      int detectorNbr = block->GetInt("DetectorNumber");
      if (Type == "Square")
        DetectorType[detectorNbr] = MG_SQUARE;
      else if (Type == "Trapezoid")
        DetectorType[detectorNbr] = MG_TRAPEZE;
      else {
        cout << "ERROR bad Annular token" << endl;
        exit(1);
      }

      number++;
      m_DetectorNumberIndex[detectorNbr] = number;
      auto A = nptool::ConvertXYZVector(block->GetVector3("X001_Y001", "mm"));    //PS: changed the variable names
      auto B = nptool::ConvertXYZVector(block->GetVector3("X128_Y001", "mm"));         // to match MUST2
      auto C = nptool::ConvertXYZVector(block->GetVector3("X001_Y128", "mm"));
      auto D = nptool::ConvertXYZVector(block->GetVector3("X128_Y128", "mm"));

      AddTelescope(DetectorType[detectorNbr], A, B, C, D);
    }

    else if (block->HasTokenList(annular)) {
        Type = block->GetMainValue();
      nptool::message("green","","","Muagst Telescope: "+Type);   //PS: correct syntax??

      int detectorNbr = block->GetInt("DetectorNumber");
      if (Type == "Annular")
        DetectorType[detectorNbr] = MG_ANNULAR;
      else {
        cout << "ERROR: Using Mugast Annular Token for Square or Trapezoid detector " << endl;
        exit(1);
      }

      number++;
      m_DetectorNumberIndex[detectorNbr] = number;
      auto Center = nptool::ConvertXYZVector(block->GetVector3("Center", "mm"));
      AddTelescope(DetectorType[detectorNbr], Center);
    }

    else if (block->HasTokenList(sphe)) {
        Type = block->GetMainValue();
      nptool::message("green","","","Muagst Telescope: "+Type);   //PS: correct syntax??

      int detectorNbr = block->GetInt("DetectorNumber");
      if (Type == "Square")
        DetectorType[detectorNbr] = MG_SQUARE;
      else if (Type == "Trapezoid")
        DetectorType[detectorNbr] = MG_TRAPEZE;
      else {
        cout << "ERROR bad Annular token" << endl;
        exit(1);
      }

      number++;
      m_DetectorNumberIndex[detectorNbr] = number;
      double Theta = block->GetDouble("THETA", "deg");
      double Phi = block->GetDouble("PHI", "deg");
      double R = block->GetDouble("R", "mm");
      vector<double> beta = block->GetVectorDouble("BETA", "deg");
      AddTelescope(DetectorType[detectorNbr], Theta, Phi, R, beta[0], beta[1], beta[2]);
    }

    else {
      cout << "ERROR: Missing token for Mugast, check your detector definition (yaml)"
              "file"
           << endl;
      exit(1);
    }
  }

  InitializeStandardParameter();
  // Create a file to be read by Ganil2Root telling which detector
  // is which shape
  std::ofstream shapeFile(".MugastShape");
  for (auto& it : DetectorType) {
    shapeFile << it.first << " " << it.second << endl;
  }
  shapeFile.close();
  ReadAnalysisConfig();
}
//////////////////////////////////////////////////////////////////////////

/* 
///////////////////////////////////////////////////////////////////////////
void MugastDetector::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  // Good for simulation, close to typical values
  vector<double> standardX = {-63, 63. / 8192.};
  vector<double> standardY = {63, -63. / 8192.};
  vector<double> standardSecondLayer = {-63, 63. / 8192.};
  vector<double> standardT = {-1000, 1000. / 8192.};

  map<int, int>::iterator it = m_DetectorNumberIndex.begin();

  for (; it != m_DetectorNumberIndex.end(); it++) {

    for (int j = 0; j < 128; j++) {
      Cal->AddParameter("Mugast", "T" + nptool::itoa(it->first) + "_DSSD_X" + nptool::itoa(j + 1) + "_E",
                        "Mugast_T" + nptool::itoa(it->first) + "_DSSD_X" + nptool::itoa(j + 1) + "_E", standardX);
      Cal->AddParameter("Mugast", "T" + nptool::itoa(it->first) + "_DSSD_Y" + nptool::itoa(j + 1) + "_E",
                        "Mugast_T" + nptool::itoa(it->first) + "_DSSD_Y" + nptool::itoa(j + 1) + "_E", standardY);
      Cal->AddParameter("Mugast", "T" + nptool::itoa(it->first) + "_DSSD_X" + nptool::itoa(j + 1) + "_T",
                        "Mugast_T" + nptool::itoa(it->first) + "_DSSD_X" + nptool::itoa(j + 1) + "_T", standardT);
      Cal->AddParameter("Mugast", "T" + nptool::itoa(it->first) + "_DSSD_Y" + nptool::itoa(j + 1) + "_T",
                        "Mugast_T" + nptool::itoa(it->first) + "_DSSD_Y" + nptool::itoa(j + 1) + "_T", standardT);
    }

    for (int j = 0; j < 16; ++j) {
      Cal->AddParameter("Mugast", "T" + nptool::itoa(it->first) + "_SecondLayer" + nptool::itoa(j + 1) + "_E",
                        "Mugast_T" + nptool::itoa(it->first) + "_SecondLayer" + nptool::itoa(j + 1) + "_E",
                        standardSecondLayer);
      Cal->AddParameter("Mugast", "T" + nptool::itoa(it->first) + "_SecondLayer" + nptool::itoa(j + 1) + "_T",
                        "Mugast_T" + nptool::itoa(it->first) + "_SecondLayer" + nptool::itoa(j + 1) + "_T", standardT);
    }
  }

  return;
}


 */
////////////////////////////////////////////////////////////////////////////////
void MugastDetector::ReadActionFile(std::string ActionFile) {
  faction_file_initialized = true;
  std::string LabelName;
  ifstream    file(ActionFile, std::ifstream::in);
  int         Label, c, d;
  if (!file)
    exit(1);
  while (file >> LabelName >> Label >> c >> d) {
    fLabelMap[Label] = LabelName;
  }
}

////////////////////////////////////////////////////////////////////////////////

/////   Specific to MugastArray   ////
void MugastDetector::AddTelescope(MG_DetectorType type, XYZVector C_X1_Y1, XYZVector C_X128_Y1, XYZVector C_X1_Y128,
                                  XYZVector C_X128_Y128) {
  // To avoid warning
  C_X128_Y128 *= 1;

  m_NumberOfTelescope++;

  // Vector U parallel to BaseLarge
  XYZVector U = C_X128_Y1 - C_X1_Y1;
  U = U.Unit();

  // Vector V parallel to height
  XYZVector V = 0.5 * (C_X1_Y128 + C_X128_Y128 - C_X1_Y1 - C_X128_Y1);
  V = V.Unit();

  //   Position Vector of Strip Center
  XYZVector StripCenter = XYZVector(0, 0, 0);
  //   Position Vector of X=1 Y=1 Strip
  XYZVector Strip_1_1;

  //   Geometry Parameter
  double Base, Height;
  if (type == MG_TRAPEZE) {
    Base = 91.48;     // mm
    Height = 104.688; // mm
  }

  if (type == MG_SQUARE) {
    Base = 91.716;   // mm
    Height = 94.916; // mm
                     //
                     //   Height        = 194.916; // mm
  }
  // double Face          = 98; // mm
  double NumberOfStrip = 128;
  double StripPitchBase = Base / NumberOfStrip;     // mm
  double StripPitchHeight = Height / NumberOfStrip; // mm
  //   Buffer object to fill Position Array
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneTelescopeStripPositionX;
  vector<vector<double>> OneTelescopeStripPositionY;
  vector<vector<double>> OneTelescopeStripPositionZ;

  //   Moving StripCenter to 1.1 corner:
  // Strip_1_1 = C_X1_Y1 + U  * (StripPitchBase / 2.) + V * (StripPitchHeight / 2.);
  // This calculation recenter the strip around the detector center.
  // This account for cases where the provided corner coordinates
  // does not match the detector size
  XYZVector Center = 0.25 * (C_X1_Y128 + C_X128_Y128 + C_X1_Y1 + C_X128_Y1);
  Strip_1_1 = Center - (0.5 * Base * U + 0.5 * Height * V) + U * (StripPitchBase / 2.) + V * (StripPitchHeight / 2.);

  for (int i = 0; i < 128; ++i) {
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for (int j = 0; j < 128; ++j) {
      // StripCenter = Strip_1_1 + StripPitch * (i * U + j * V);
      StripCenter = Strip_1_1 + i * U * StripPitchBase + j * V * StripPitchHeight;
      lineX.push_back(StripCenter.X());
      lineY.push_back(StripCenter.Y());
      lineZ.push_back(StripCenter.Z());
    }

    OneTelescopeStripPositionX.push_back(lineX);
    OneTelescopeStripPositionY.push_back(lineY);
    OneTelescopeStripPositionZ.push_back(lineZ);
  }
  m_StripPositionX.push_back(OneTelescopeStripPositionX);
  m_StripPositionY.push_back(OneTelescopeStripPositionY);
  m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
}


///////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void MugastDetector::AddTelescope(MG_DetectorType type, XYZVector C_Center) {
  // To avoid warning
  m_NumberOfTelescope++;
  double Z = C_Center.Z();

  double R_Min = 24;
  double R_Max = 48;

  double Phi_Min = 0;
  double Phi_Max = 360;

  int NumberOfQuadrant = 4;
  int NumberofRing = 16;   // Per Quadrant
  int NumberofSector = 16; // Per detector, ( 4 in each Quad)

  double StripPitchSector = (Phi_Max - Phi_Min) / (NumberofSector); // radial strip spacing in rad
  double StripPitchRing = (R_Max - R_Min) / NumberofRing;           // ring strip spacing in mm

  // double Phi_0 = 8*StripPitchSector; // Phi Offset: 1st sector starts at 180 degrees and ends at 180-22.5 degrees in
  // the lab frame, numbering goes clockwise
  double Phi_0 = 90;
  TVector3 Strip_1_1 = TVector3(0, 0, Z);
  TVector3 StripCenter = Strip_1_1;

  //   Buffer object to fill Position Array
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;
  vector<vector<double>> OneStripPositionX;
  vector<vector<double>> OneStripPositionY;
  vector<vector<double>> OneStripPositionZ;

  for (int iQuad = 0; iQuad < NumberOfQuadrant; iQuad++) {
    for (int iRing = 0; iRing < NumberofRing; iRing++) {

      lineX.clear();
      lineY.clear();
      lineZ.clear();

      for (int iSector = 0; iSector < NumberofSector; iSector++) {

        // Build vector
        StripCenter = TVector3(C_Center.X() + R_Min + (iRing + 0.5) * StripPitchRing, C_Center.Y(), Z);
        StripCenter.RotateZ((Phi_0 + (iSector + 0.5) * StripPitchSector) * M_PI / 180.);

        // these vectors will contain 16x4 = 64 elements
        lineX.push_back(StripCenter.X());
        lineY.push_back(StripCenter.Y());
        lineZ.push_back(StripCenter.Z());
      }
      OneStripPositionX.push_back(lineX);
      OneStripPositionY.push_back(lineY);
      OneStripPositionZ.push_back(lineZ);
    }
  }

  // Increase the size of the Position array to 128 to avoid seg fault
  // in case of connecting a trapezoid to an annular
  vector<double> defaultLine;
  defaultLine.resize(128, -1000);
  OneStripPositionX.resize(128, defaultLine);
  OneStripPositionY.resize(128, defaultLine);
  OneStripPositionZ.resize(128, defaultLine);

  m_StripPositionX.push_back(OneStripPositionX);
  m_StripPositionY.push_back(OneStripPositionY);
  m_StripPositionZ.push_back(OneStripPositionZ);

  return;
}

////////////////////////////////////////////////////////////////////////////////
void MugastDetector::AddTelescope(MG_DetectorType type, double theta, double phi, double distance, double beta_u,
                                  double beta_v, double beta_w) {

  m_NumberOfTelescope++;

  double Pi = 3.141592654;

  // convert from degree to radian:
  theta = theta * Pi / 180.;
  phi = phi * Pi / 180.;

  // Vector U on Telescope Face (paralelle to Y Strip) (NB: remember that Y
  // strip are allong X axis)
  TVector3 U;
  // Vector V on Telescope Face (parallele to X Strip)
  TVector3 V;
  // Vector W normal to Telescope Face (pointing CsI)
  TVector3 W;
  // Vector position of Telescope Face center
  TVector3 C;

  C = TVector3(distance * sin(theta) * cos(phi), distance * sin(theta) * sin(phi), distance * cos(theta));

  TVector3 P = TVector3(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));

  W = C.Unit();
  U = W.Cross(P);
  V = W.Cross(U);

  U = U.Unit();
  V = V.Unit();

  U.Rotate(beta_u * Pi / 180., U);
  V.Rotate(beta_u * Pi / 180., U);

  U.Rotate(beta_v * Pi / 180., V);
  V.Rotate(beta_v * Pi / 180., V);

  U.Rotate(beta_w * Pi / 180., W);
  V.Rotate(beta_w * Pi / 180., W);

  double Face = 98; // mm
  double NumberOfStrip = 128;
  double StripPitch = Face / NumberOfStrip; // mm

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneTelescopeStripPositionX;
  vector<vector<double>> OneTelescopeStripPositionY;
  vector<vector<double>> OneTelescopeStripPositionZ;

  double X, Y, Z;

  // Moving C to the 1.1 corner:
  C.SetX(C.X() - (Face / 2 - StripPitch / 2) * (V.X() + U.X()));
  C.SetY(C.Y() - (Face / 2 - StripPitch / 2) * (V.Y() + U.Y()));
  C.SetZ(C.Z() - (Face / 2 - StripPitch / 2) * (V.Z() + U.Z()));

  for (int i = 0; i < 128; ++i) {

    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for (int j = 0; j < 128; ++j) {
      X = C.X() + StripPitch * (U.X() * i + V.X() * j);
      Y = C.Y() + StripPitch * (U.Y() * i + V.Y() * j);
      Z = C.Z() + StripPitch * (U.Z() * i + V.Z() * j);

      lineX.push_back(X);
      lineY.push_back(Y);
      lineZ.push_back(Z);
    }

    OneTelescopeStripPositionX.push_back(lineX);
    OneTelescopeStripPositionY.push_back(lineY);
    OneTelescopeStripPositionZ.push_back(lineZ);
  }
  m_StripPositionX.push_back(OneTelescopeStripPositionX);
  m_StripPositionY.push_back(OneTelescopeStripPositionY);
  m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
}

///////////////////////////////////////////////////////////////////////////////

void MugastDetector::InitializeStandardParameter() {
  //   Enable all channel
  vector<bool> ChannelStatus;
  m_XChannelStatus.clear();
  m_YChannelStatus.clear();
  m_SecondLayerChannelStatus.clear();

  ChannelStatus.resize(128, true);
  for (int i = 0; i < m_NumberOfTelescope; ++i) {
    auto it = m_DetectorNumberIndex.begin();
    for (unsigned int j = 0; j < i; j++) {
      it++;
    }
    int det = it->first;
    m_XChannelStatus[det] = ChannelStatus;
    m_YChannelStatus[det] = ChannelStatus;
  }

  ChannelStatus.resize(16, true);
  for (int i = 0; i < m_NumberOfTelescope; ++i) {
    auto it = m_DetectorNumberIndex.begin();
    for (unsigned int j = 0; j < i; j++) {
      it++;
    }

    int det = it->first;
    m_SecondLayerChannelStatus[det] = ChannelStatus;
  }

  m_MaximumStripMultiplicityAllowed = m_NumberOfTelescope;

  return;
}

///////////////////////////////////////////////////////////////////////////////
TVector3 MugastDetector::GetPositionOfInteraction(const int i, bool random) {
  TVector3 Position = TVector3(GetStripPositionX(m_PhysicsData->TelescopeNumber[i], m_PhysicsData->DSSD_X[i], m_PhysicsData->DSSD_Y[i]),
                               GetStripPositionY(m_PhysicsData->TelescopeNumber[i], m_PhysicsData->DSSD_X[i], m_PhysicsData->DSSD_Y[i]),
                               GetStripPositionZ(m_PhysicsData->TelescopeNumber[i], m_PhysicsData->DSSD_X[i], m_PhysicsData->DSSD_Y[i]));

  return Position;
}

///////////////////////////////////////////////////////////////////////////////
TVector3 MugastDetector::GetTelescopeNormal(const int i) {
  TVector3 U = TVector3(GetStripPositionX(m_PhysicsData->TelescopeNumber[i], 128, 1), GetStripPositionY(m_PhysicsData->TelescopeNumber[i], 128, 1),
                        GetStripPositionZ(m_PhysicsData->TelescopeNumber[i], 128, 1))

               - TVector3(GetStripPositionX(m_PhysicsData->TelescopeNumber[i], 1, 1), GetStripPositionY(m_PhysicsData->TelescopeNumber[i], 1, 1),
                          GetStripPositionZ(m_PhysicsData->TelescopeNumber[i], 1, 1));

  TVector3 V =
      TVector3(GetStripPositionX(m_PhysicsData->TelescopeNumber[i], 128, 128), GetStripPositionY(m_PhysicsData->TelescopeNumber[i], 128, 128),
               GetStripPositionZ(m_PhysicsData->TelescopeNumber[i], 128, 128))

      - TVector3(GetStripPositionX(m_PhysicsData->TelescopeNumber[i], 128, 1), GetStripPositionY(m_PhysicsData->TelescopeNumber[i], 128, 1),
                 GetStripPositionZ(m_PhysicsData->TelescopeNumber[i], 128, 1));

  TVector3 Normal = U.Cross(V);

  return (Normal.Unit());
}



  //   DSSD
  //   X
  double MugastDetector::fDSSD_X_E(const int& i) {
    static string name;
    name = "Mugast/T";
    name += nptool::itoa(m_RawData->GetDSSDXEDetectorNbr(i));
    name += "_DSSD_X";
    name += nptool::itoa(m_RawData->GetDSSDXEStripNbr(i));
    name += "_E";
    return m_Cal.ApplyCalibration(name, m_RawData->GetDSSDXEEnergy(i), 1);
  }

  double MugastDetector::fDSSD_X_T(const int& i) {
    static string name;
    name = "Mugast/T";
    name += nptool::itoa(m_RawData->GetDSSDXTDetectorNbr(i));
    name += "_DSSD_X";
    name += nptool::itoa(m_RawData->GetDSSDXTStripNbr(i));
    name += "_T";
    return m_Cal.ApplyCalibration(name, m_RawData->GetDSSDXTTime(i), 1);
  }

  //   Y
  double MugastDetector::fDSSD_Y_E(const int& i) {
    static string name;
    name = "Mugast/T";
    name += nptool::itoa(m_RawData->GetDSSDYEDetectorNbr(i));
    name += "_DSSD_Y";
    name += nptool::itoa(m_RawData->GetDSSDYEStripNbr(i));
    name += "_E";
    return m_Cal.ApplyCalibration(name, m_RawData->GetDSSDYEEnergy(i), 1);
  }

  double MugastDetector::fDSSD_Y_T(const int& i) {
    static string name;
    name = "Mugast/T";
    name += nptool::itoa(m_RawData->GetDSSDYTDetectorNbr(i));
    name += "_DSSD_Y";
    name += nptool::itoa(m_RawData->GetDSSDYTStripNbr(i));
    name += "_T";
    return m_Cal.ApplyCalibration(name, m_RawData->GetDSSDYTTime(i), 1);
  }

  //   SecondLayer
  double MugastDetector::fSecondLayer_E(const int& i) {
    static string name;
    name = "Mugast/T";
    name += nptool::itoa(m_RawData->GetSecondLayerEDetectorNbr(i));
    name += "_SecondLayer";
    name += nptool::itoa(m_RawData->GetSecondLayerEStripNbr(i));
    name += "_E";
    return m_Cal.ApplyCalibration(name, m_RawData->GetSecondLayerEEnergy(i), 1);
  }

  double MugastDetector::fSecondLayer_T(const int& i) {
    static string name;
    name = "Mugast/T";
    name += nptool::itoa(m_RawData->GetSecondLayerTDetectorNbr(i));
    name += "_SecondLayer";
    name += nptool::itoa(m_RawData->GetSecondLayerTStripNbr(i));
    name += "_T";
    return m_Cal.ApplyCalibration(name, m_RawData->GetSecondLayerTTime(i), 1);
  }



///////////////////////////////////******************************STOPPPPPPPPPPPPPPPPP PPPPPPPPPPPSSSSSSSSSSSSSSSSS */


void MugastDetector::BuildRawEvent(const std::string& daq_name,
                                  const std::string& st_type_key,
                                  void*              commonframe) {
#ifdef MFM_FOUND
  int type_key = ((MFMCommonFrame*)commonframe)->GetFrameType();
  if (type_key == MFM_EBY_EN_TS_FRAME_TYPE) {
    TreatFrame((MFMCommonFrame*)commonframe);
  }
#endif
}

void MugastDetector::TreatFrame(void* commonframe) {
#ifdef MFM_FOUND
  unsigned short     NItems = 0;
  unsigned short     value, label_id;
  std::string        label;
  unsigned short     LblData[2];
  unsigned long long TS;
  MG_DetectorType type = MG_NOCHANGE ;

  std::shared_ptr<MFMEbyedatFrame> EbyEframe
      = std::make_shared<MFMEbyedatFrame>();
  EbyEframe->SetAttributs(((MFMCommonFrame*)commonframe)->GetPointHeader());

  NItems = EbyEframe->GetNbItems();
  TS     = EbyEframe->GetTimeStamp();
  if (TS == 0) {
    std::cout << "/////////////////////////////" << std::endl;
    std::cout << "   MUGAST WARNING TS IS 0    " << std::endl;
    std::cout << "/////////////////////////////" << std::endl;
  }

  for (unsigned short i = 0; i < NItems; i++) {
    EbyEframe->EbyedatGetParameters(i, &label_id, &value);
    LblData[0] = label_id;
    LblData[1] = value;

    label = fLabelMap[label_id];

    if (label.compare(0, 2, "MG") == 0) {
      double det = atoi(label.substr(2, 1).c_str());
      double channel;
      if (label.compare(5, 6, "STRX_E") == 0) {
        channel = atoi(label.substr(11).c_str());
        m_RawData->SetDSSDXE(type, det, channel, value);
      }

      else if (label.compare(5, 6, "STRX_T") == 0) {
        channel = atoi(label.substr(11).c_str());    
        m_RawData->SetDSSDXT(type, det, channel, value);
      }

      else if (label.compare(5, 6, "STRY_E") == 0) {
        channel = atoi(label.substr(11).c_str());
        m_RawData->SetDSSDYE(type, det, channel, value);
      }

      else if (label.compare(5, 6, "STRY_T") == 0) {
        channel = atoi(label.substr(11).c_str());
        m_RawData->SetDSSDYT(type, det, channel, value);
      }
      // m_RawData->SetMMTS(TS);
    }
  }
#endif
}







////////////////////////////////////////////////////////////////////////////////
void MugastDetector::InitializeDataInputRaw(
    std::shared_ptr<nptool::VDataInput> input) {
  input->Attach("mugast", "mugast::MugastData", &m_RawData);
}
////////////////////////////////////////////////////////////////////////////////
void MugastDetector::InitializeDataInputPhysics(
    std::shared_ptr<nptool::VDataInput> input) {
  input->Attach("mugast", "mugast::MugastPhysics", &m_PhysicsData);
}
////////////////////////////////////////////////////////////////////////////////
void MugastDetector::InitializeDataOutputRaw(
    std::shared_ptr<nptool::VDataOutput> output) {
  output->Attach("mugast", "mugast::MugastData", &m_RawData);
}
////////////////////////////////////////////////////////////////////////////////
void MugastDetector::InitializeDataOutputPhysics(
    std::shared_ptr<nptool::VDataOutput> output) {
  output->Attach("mugast", "mugast::MugastPhysics", &m_PhysicsData);
}





////////////////////////////////////////////////////////////////////////////////
void MugastDetector::InitSimulation(std::string simtype) {
  // store the loaded simulation
  m_simtype = simtype;

  // load the plugin
  auto handle = nptool::Application::GetApplication()->LoadPlugin(
      "mugast-" + m_simtype, true);
  // build the class
#ifdef Geant4_FOUND
  auto func = (std::shared_ptr<nptool::geant4::VDetector>(*)())dlsym(
      handle, "ConstructDetectorSimulation");
  if (func)
    m_Geant4 = (*func)();
  else
    throw(nptool::Error("MugastDetector", "Fail to load Geant4 module"));
#endif
}
////////////////////////////////////////////////////////////////////////////////

void MugastDetector::ConstructGeometry() {
#ifdef Geant4_FOUND
  if (m_Geant4) {
    // should load the library here and find the external constructor
    std::cout << "I CONSTRUCT GEOMETRY" << std::endl;
    m_Geant4->ConstructDetector();
  }
#endif
}




////////////////////////////////////////////////////////////////////////////////
void MugastDetector::InitSpectra() {
  m_Spectra = std::make_shared<mugast::MugastSpectra>();
};
////////////////////////////////////////////////////////////////////////////////
void MugastDetector::FillSpectra() {
  // m_Spectra->FillRaw();
  // m_Spectra->FillPhy();
};
////////////////////////////////////////////////////////////////////////////////
void MugastDetector::WriteSpectra() {};
////////////////////////////////////////////////////////////////////////////////
void MugastDetector::CheckSpectra() {};
////////////////////////////////////////////////////////////////////////////////
void MugastDetector::ClearSpectra() { m_Spectra->Clear(); };
////////////////////////////////////////////////////////////////////////////////
extern "C" {
shared_ptr<nptool::VDetector> ConstructDetector() {
  return make_shared<mugast::MugastDetector>();
};
}
