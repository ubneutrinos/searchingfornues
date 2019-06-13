#ifndef PIDFUNCS_H
#define PIDFUNCS_H
#include "ubana/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

namespace searchingfornues
{

    double PID(art::Ptr<anab::ParticleID> selected_pid,
                                std::string AlgName,
                                anab::kVariableType VariableType,
                                anab::kTrackDir TrackDirection,
                                int pdgCode,
                                int selectedPlane)
    {

    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = selected_pid->ParticleIDAlgScores();
    for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
    {
        anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
        int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneMask);

        if (selectedPlane != planeid) {
        continue;
        }

        if (AlgScore.fAlgName == AlgName)
        {
            if (anab::kVariableType(AlgScore.fVariableType) == VariableType && anab::kTrackDir(AlgScore.fTrackDir) == TrackDirection)
            {
                if (AlgScore.fAssumedPdg == pdgCode)
                {
                    double alg_value = AlgScore.fValue;
                    return alg_value;
                }
            }
        }
    }
    return std::numeric_limits<double>::lowest();
    }
} // namespace searchingfornues

#endif