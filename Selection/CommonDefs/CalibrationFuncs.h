#ifndef CALIBRATIONFUNCS_H
#define CALIBRATIONFUNCS_H

#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibProvider.h"
#include "ubevt/Database/TPCEnergyCalib/TPCEnergyCalibService.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "SCECorrections.h"

namespace searchingfornues
{

    void getCali(
        std::vector<art::Ptr<recob::SpacePoint>> spcpnts,
        art::FindManyP<recob::Hit> hits_per_spcpnts,
        std::vector<float> &cali_corr)
    {
        const lariov::TPCEnergyCalibProvider &_energy_calib_provider = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();

        cali_corr.resize(3);
        std::vector<float> total_charge(3, 0);

        for (auto _sps : spcpnts)
        {
            std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(_sps.key());
            const double *xyz = _sps->XYZ();

            for (auto &hit : hits)
            {
                auto plane_nr = hit->View();

                if (plane_nr > 2 || plane_nr < 0)
                    continue;

                total_charge[plane_nr] += hit->Integral();
                float yzcorrection = _energy_calib_provider.YZdqdxCorrection(plane_nr, xyz[1], xyz[2]);
                float xcorrection = _energy_calib_provider.XdqdxCorrection(plane_nr, xyz[0]);

                if (!yzcorrection)
                    yzcorrection = 1.0;
                if (!xcorrection)
                    xcorrection = 1.0;

                cali_corr[plane_nr] += yzcorrection * xcorrection * hit->Integral();
            }
        }

        for (unsigned short i = 0; i < 3; ++i)
        {
            if (total_charge[i] > 0)
            {
                cali_corr[i] /= total_charge[i];
            }
            else
            {
                cali_corr[i] = 1;
            }
        }
    }

    void getDQdxCali(art::Ptr <recob::Shower> shower_obj,
                     std::vector<float> &dqdx_cali)
    {
        TVector3 pfp_dir;
        const lariov::TPCEnergyCalibProvider &_energy_calib_provider = art::ServiceHandle<lariov::TPCEnergyCalibService>()->GetProvider();

        // Field needed for calibration factor
        float x_start, y_start, z_start;
        float x_middle, y_middle, z_middle;
        float x_end, y_end, z_end;
        float start_corr, middle_corr, end_corr;

        pfp_dir.SetX(shower_obj->Direction().X());
        pfp_dir.SetY(shower_obj->Direction().Y());
        pfp_dir.SetZ(shower_obj->Direction().Z());

        x_start = shower_obj->ShowerStart().X();
        y_start = shower_obj->ShowerStart().Y();
        z_start = shower_obj->ShowerStart().Z();

        float _dQdx_rectangle_length = 4;
        pfp_dir.SetMag(_dQdx_rectangle_length / 2.); //Go 2cm along the direction of the object.
        x_middle = x_start + pfp_dir.X();
        y_middle = y_start + pfp_dir.Y();
        z_middle = z_start + pfp_dir.Z();
        x_end = x_middle + pfp_dir.X();
        y_end = y_middle + pfp_dir.Y();
        z_end = z_middle + pfp_dir.Z();
        pfp_dir.SetMag(1.); //Normalise again for safety (not needed).

        for (int plane_nr = 0; plane_nr < 3; ++plane_nr)
        {
            float yzcorrection_start = _energy_calib_provider.YZdqdxCorrection(plane_nr, y_start, z_start);
            float xcorrection_start = _energy_calib_provider.XdqdxCorrection(plane_nr, x_start);
            if (!yzcorrection_start)
                yzcorrection_start = 1.0;
            if (!xcorrection_start)
                xcorrection_start = 1.0;
            start_corr = yzcorrection_start * xcorrection_start;

            float yzcorrection_middle = _energy_calib_provider.YZdqdxCorrection(plane_nr, y_middle, z_middle);
            float xcorrection_middle = _energy_calib_provider.XdqdxCorrection(plane_nr, x_middle);
            if (!yzcorrection_middle)
                yzcorrection_middle = 1.0;
            if (!xcorrection_middle)
                xcorrection_middle = 1.0;
            middle_corr = yzcorrection_middle * xcorrection_middle;

            float yzcorrection_end = _energy_calib_provider.YZdqdxCorrection(plane_nr, y_end, z_end);
            float xcorrection_end = _energy_calib_provider.XdqdxCorrection(plane_nr, x_end);
            if (!yzcorrection_end)
                yzcorrection_end = 1.0;
            if (!xcorrection_end)
                xcorrection_end = 1.0;
            end_corr = yzcorrection_end * xcorrection_end;
            //std::cout << "[EnergyHelper] dqdx_cali " << start_corr << middle_corr << end_corr << std::endl;
            dqdx_cali[plane_nr] = (start_corr + middle_corr + end_corr) / 3;
        }
    }
    
    //----------------------------------------------------------------------------------                                                                           
    // Modified Box model correction                                                                                                                                     
    double ModBoxCorrection(const double dQdx, const float x, const float y, const float z) {
      
      double rho = 1.383;//detprop->Density();            // LAr density in g/cm^3
      double Wion = 23.6/1e6;//util::kGeVToElectrons;  // 23.6 eV = 1e, Wion in MeV/e
      
      auto E_field = GetLocalEFieldMag(x,y,z); // kV / cm
      
      double fModBoxA = 0.930;
      double fModBoxB = 0.212;
      
      double Beta = fModBoxB / (rho * E_field);
      double Alpha = fModBoxA;
      double dEdx = (exp(Beta * Wion * dQdx ) - Alpha) / Beta;
      
      return dEdx;
      
    }
    
} // namespace searchingfornues

#endif
