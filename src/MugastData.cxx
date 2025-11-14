#include <iostream>

/*****************************************************************************
 * Original Author: V. Girard-Alcindor                                       *
 * contact address: girard-alcindor@ijclab.in2p3.fr                          *
 *                                                                           *
 * Creation Date  : 08/03/24                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * This class hold Mugast raw data                                            *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This class is heavily based on the nptool v3 Mugast detector               *
 *                                                                           *
 *****************************************************************************/

using namespace std;

#include "MugastData.h"
#include "MugastMap.h"

mugast::MugastData::MugastData() {
      // Init the correspondace table
  for (unsigned int i = 0; i < 128; i++) {
    fMG_MapTrapezeX[i + 1] = MUGAST_MAP::TrapezeX[i];
    fMG_MapTrapezeY[i + 1] = MUGAST_MAP::TrapezeY[i];
    fMG_MapSquareX[i + 1] = MUGAST_MAP::SquareX[i];
    fMG_MapSquareY[i + 1] = MUGAST_MAP::SquareY[i];
    fMG_MapAnnularX[i + 1] = MUGAST_MAP::AnnularX[i];
    fMG_MapAnnularY[i + 1] = MUGAST_MAP::AnnularY[i];
  }
}

mugast::MugastData::~MugastData() {}

void mugast::MugastData::Clear() {

  fMG_DSSDXE_DetectorNbr.clear();
  fMG_DSSDXE_StripNbr.clear();
  fMG_DSSDXE_Energy.clear();
  fMG_DSSDXT_DetectorNbr.clear();
  fMG_DSSDXT_StripNbr.clear();
  fMG_DSSDXT_Time.clear();
  fMG_DSSDYE_DetectorNbr.clear();
  fMG_DSSDYE_StripNbr.clear();
  fMG_DSSDYE_Energy.clear();
  fMG_DSSDYT_DetectorNbr.clear();
  fMG_DSSDYT_StripNbr.clear();
  fMG_DSSDYT_Time.clear();
  fMG_SecondLayerE_DetectorNbr.clear();
  fMG_SecondLayerE_StripNbr.clear();
  fMG_SecondLayerE_Energy.clear();
  fMG_SecondLayerT_DetectorNbr.clear();
  fMG_SecondLayerT_StripNbr.clear();
  fMG_SecondLayerT_Time.clear();
}

void mugast::MugastData::Dump() const {
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX Mugast Event XXXXXXXXXXXXXXXXX" << endl;

  cout << "// First Layer " << endl;
  // (X,E)
  cout << " DSSDXE_Mult = " << fMG_DSSDXE_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fMG_DSSDXE_DetectorNbr.size(); i++)
    cout << "  DetNbr: " << fMG_DSSDXE_DetectorNbr[i] << " DSSD: " << fMG_DSSDXE_StripNbr[i]
         << " Energy: " << fMG_DSSDXE_Energy[i] << endl;
  // (X,T)
  cout << " DSSDXT_Mult = " << fMG_DSSDXT_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fMG_DSSDXT_DetectorNbr.size(); i++)
    cout << "  DetNbr: " << fMG_DSSDXT_DetectorNbr[i] << " DSSD: " << fMG_DSSDXT_StripNbr[i]
         << " Time: " << fMG_DSSDXT_Time[i] << endl;
  // (Y,E)
  cout << " DSSDYE_Mult = " << fMG_DSSDYE_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fMG_DSSDYE_DetectorNbr.size(); i++)
    cout << "  DetNbr: " << fMG_DSSDYE_DetectorNbr[i] << " DSSD: " << fMG_DSSDYE_StripNbr[i]
         << " Energy: " << fMG_DSSDYE_Energy[i] << endl;
  // (Y,T)
  cout << " DSSDYT_Mult = " << fMG_DSSDYT_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fMG_DSSDYT_DetectorNbr.size(); i++)
    cout << "  DetNbr: " << fMG_DSSDYT_DetectorNbr[i] << " DSSD: " << fMG_DSSDYT_StripNbr[i]
         << " Time: " << fMG_DSSDYT_Time[i] << endl;

  // SecondLayer
  // Energy
  cout << "// Second Layer " << endl;
  cout << " SecondLayerE_Mult = " << fMG_SecondLayerE_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fMG_SecondLayerE_DetectorNbr.size(); i++)
    cout << "  Det: " << fMG_SecondLayerE_DetectorNbr[i] << " DSSD: " << fMG_SecondLayerE_StripNbr[i]
         << " Energy: " << fMG_SecondLayerE_Energy[i] << endl;
  // Time
  cout << " SecondLayerT_Mult = " << fMG_SecondLayerT_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fMG_SecondLayerT_DetectorNbr.size(); i++)
    cout << "  Det: " << fMG_SecondLayerT_DetectorNbr[i] << " DSSD: " << fMG_SecondLayerT_StripNbr[i]
         << " Time: " << fMG_SecondLayerT_Time[i] << endl;
}
