#ifndef PIDFUNCS_H
#define PIDFUNCS_H


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
        int planeid = -1;
        // std::cout << "AlgScore index " << i_algscore << " plane " << AlgScore.fPlaneMask<< std::endl;
        if (AlgScore.fPlaneMask.none() || AlgScore.fPlaneMask.count() > 1 || (AlgScore.fPlaneMask.count() == 1 && !(AlgScore.fPlaneMask.test(0) || AlgScore.fPlaneMask.test(1) || AlgScore.fPlaneMask.test(2))))
        {
        // std::cout << "[uB_PlaneIDBitsetHelper] Cannot return a single MicroBooNE plane for bitset " << AlgScore.fPlaneMask << ". Returning -1 (invalid planeID)." << std::endl;
        continue;
        }
        else if (AlgScore.fPlaneMask.test(0))
        planeid = 2;
        else if (AlgScore.fPlaneMask.test(1))
        planeid = 1;
        else if (AlgScore.fPlaneMask.test(2))
        planeid = 0;
        else
        planeid = -1;

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