#include <math.h>
#include <stdio.h>  /* printf */

#include "header.hxx"
#include "Log.hxx"
#include "BeamManager.hxx"

BeamManager* BeamManager::fBeamManager = NULL;

BeamManager::BeamManager():
    beamType(kSPring8),
    beamSlz(0),
    beamSlx(0),
    beamInz(0),
    beamInx(0),
    beamSlzRange(0),
    beamSlxRange(0),
    beamInzRange(0),
    beamInxRange(0)
{
}

bool BeamManager::Initialize(BeamType theBeamType){
    beamType = theBeamType;
    if (beamType == kSPring8 || beamType == kSPring8Tilted){
        // FIXME: currently set a broader range. Need further investigation
        beamSlz    = 0.02;
        beamSlx    = 0;
        beamInz    = 150; // mm
        beamInx    = 0;
        beamSlzRange = 0.3;
        beamSlxRange = 0.1;
        beamInzRange = chamberHL;
        beamInxRange = chamberHH;
        if (beamType == kSPring8Tilted){
            beamSlx += tan(-18.4*M_PI/180);
        }
    }
    else if (beamType == kCosmic){
        MyError("Cosmic type of beam is not supported yet!");
        return false;
    }
    return true;
}

void BeamManager::Print(){
    printf("#########################BeamManager###########################\n");
    printf("  beam type %d:\n",beamType);
    printf("  beamSlz       = %.3e\n",beamSlz);
    printf("  beamSlx       = %.3e\n",beamSlx);
    printf("  beamInz       = %.3e\n",beamInz);
    printf("  beamInx       = %.3e\n",beamInx);
    printf("  beamSlzRange  = %.3e\n",beamSlzRange);
    printf("  beamSlxRange  = %.3e\n",beamSlxRange);
    printf("  beamInzRange  = %.3e\n",beamInzRange);
    printf("  beamInxRange  = %.3e\n",beamInxRange);
}
