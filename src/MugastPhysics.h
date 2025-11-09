#ifndef NPMugastPhysics_H
#define NPMugastPhysics_H

/*****************************************************************************
 * Original Author: V. Girard-Alcindor                                       *
 * contact address: girard-alcindor@ijclab.in2p3.fr                          *
 *                                                                           *
 * Creation Date  : 08/03/24                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * This class hold MUGAST Physics data                                        *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This class is heavily based on the nptool v3 MUGAST detector               *
 *                                                                           *
 *****************************************************************************/
using namespace std;
#include <vector>
#include <string>
#include "TObject.h"
namespace mugast {
class MugastPhysics {
public:
  MugastPhysics(){};
  ~MugastPhysics(){};

public:
  void Clear();
  void Dump() const;


public:
  //   Provide Physical Multiplicity
  Int_t EventMultiplicity;

  //   Provide a Classification of Event
  vector<int> EventType;

  // Telescope
  vector<int> TelescopeNumber;
  //   DSSD
  vector<double> DSSD_E;
  vector<double> DSSD_T;
  vector<int> DSSD_X;
  vector<int> DSSD_Y;

  vector<double> PosX;
  vector<double> PosY;
  vector<double> PosZ;
  // vector<double> Theta;        //PS: not used in np3

  //   Second Layer
  vector<double> SecondLayer_E;
  vector<double> SecondLayer_T;
  vector<int> SecondLayer_N;

  // Physical Value
  vector<double> TotalEnergy;

};
} // namespace must2
#endif
