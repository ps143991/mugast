#include "MugastPhysics.h"

/*****************************************************************************
 * Original Author: V. Girard-Alcindor                                       *
 * contact address: girard-alcindor@ijclab.in2p3.fr                          *
 *                                                                           *
 * Creation Date  : 08/03/24                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * This class hold MUST2 Physics data                                        *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This class is heavily based on the nptool v3 MUST2 detector               *
 *                                                                           *
 *****************************************************************************/

void mugast::MugastPhysics::Clear() {
  EventMultiplicity = 0;

  TelescopeNumber.clear();
  EventType.clear();
  TotalEnergy.clear();

  PosX.clear();
  PosY.clear();
  PosZ.clear();

  // DSSD
  DSSD_E.clear();
  DSSD_T.clear();
  DSSD_X.clear();
  DSSD_Y.clear();

  // SecondLayer
  SecondLayer_E.clear();
  SecondLayer_T.clear();
  SecondLayer_N.clear();
}
