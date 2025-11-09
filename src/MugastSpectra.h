#ifndef MugastSpectra_h
#define MugastSpectra_h

/*****************************************************************************
 * Original Author: V. Girard-Alcindor                                       *
 * contact address: girard-alcindor@ijclab.in2p3.fr                          *
 *                                                                           *
 * Creation Date  : 08/03/24                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * This class hold MUST2 Spectra definitions                                 * 
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This class is heavily based on the nptool v3 MUST2 detector               *
 *                                                                           *
 *****************************************************************************/

#include "MugastDetector.h"
// root
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
// std
#include <map>
#include <memory>
namespace mugast {

  // forward declaration is necessary
  // class MugastDetector;
  class MugastSpectra {
   public:
    MugastSpectra();
    ~MugastSpectra(){};

   private:
    std::shared_ptr<mugast::MugastDetector> m_detector;
    mugast::MugastData* m_RawData;
    mugast::MugastPhysics* m_PhysicsData;
    std::map<std::string, TH1*> m_raw_hist;
    std::map<std::string, TCanvas*> m_raw_canvas;
    std::map<std::string, TH1*> m_phy_hist;
    std::map<std::string, TCanvas*> m_phy_canvas;
    std::shared_ptr<nptool::Application> m_app;
    // Boolean used to display more spectra
    bool is_expert;

   public:
    void FillRaw();
    void FillPhy();
    void Clear();
  };

} // namespace mugast
#endif
