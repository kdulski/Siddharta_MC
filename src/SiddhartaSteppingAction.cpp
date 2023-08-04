#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaSteppingAction.h"
#include "../include/SiddhartaHisto.h"

#include <G4SteppingManager.hh>
#include <G4RunManager.hh>
#include <G4Event.hh>

SiddhartaSteppingAction::SiddhartaSteppingAction() {}

void SiddhartaSteppingAction::UserSteppingAction(const G4Step* gst1)
{
  SiddhartaAnalysisManager* analysis = SiddhartaAnalysisManager::getInstance();
  G4String matname = "";
  G4String volname = "";
  G4String partname = gst1->GetTrack()->GetDefinition()->GetParticleName();

  if (partname == "kaon-") {
    G4TrackStatus kStatus = gst1->GetTrack()->GetTrackStatus();
    G4double Kinener = gst1->GetTrack()->GetKineticEnergy();
    G4ThreeVector pos_km_bp = gst1->GetTrack()->GetPosition();
    G4double sinethetaKm = sqrt(pos_km_bp[0]*pos_km_bp[0]/(pos_km_bp[0]*pos_km_bp[0] + pos_km_bp[1]*pos_km_bp[1]));
    volname = gst1->GetPreStepPoint()->GetPhysicalVolume()->GetName();

    if (volname == "BeamPipeIpVacuum") {
      if(analysis->YesHistos)
        analysis->histo->fillHisto("4", Kinener, 1.);
    }

    G4int evid2 = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    if (evid2 != prev_evid2) {
      analysis->histo->fillHistogram("kaonParticles", 1, doubleCheck(1), doubleCheck(1));
      prev_evid2 = evid2;

      bool initialMomentumConditionOK = analysis->CheckKaonIniMom(pos_km_bp, gst1->GetTrack()->GetMomentum());

      if (!initialMomentumConditionOK)
        G4RunManager::GetRunManager()->AbortEvent();
    }

    if (Kinener <= 1001.*eV) {
      G4int evid = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
      if (evid != prev_evid) {
        G4ThreeVector position = gst1->GetTrack()->GetPosition();
        if(analysis->YesHistos)
          analysis->histo->fillHisto("5", position[1], 1.);

        matname = gst1->GetPreStepPoint()->GetMaterial()->GetName();
        volname = gst1->GetPreStepPoint()->GetPhysicalVolume()->GetName();
        if(analysis->YesHistos) {
          if (volname == "Target") {
            analysis->histo->fillHisto("14", position[1], 1.);
            analysis->histo->fillHisto3("1", position[0], position[2], position[1], 1.);
            analysis->histo->ntuData.XYZstopK[0] = position[0];
            analysis->histo->ntuData.XYZstopK[1] = position[1];
            analysis->histo->ntuData.XYZstopK[2] = position[2];
            analysis->histo->fillHisto3("5", position[0], position[2], position[1], 1.);

            if (analysis->histo->ntuData.XYZKMTop[0] > 0 && analysis->histo->ntuData.XYZKMTop[0] < 20) {
              analysis->histo->fillHisto("1401", position[1], 1.);
            } else if (analysis->histo->ntuData.XYZKMTop[0]>20 && analysis->histo->ntuData.XYZKMTop[0] < 40) {
              analysis->histo->fillHisto("1402", position[1], 1.);
            } else if (analysis->histo->ntuData.XYZKMTop[0]>40 && analysis->histo->ntuData.XYZKMTop[0] < 60) {
              analysis->histo->fillHisto("1403", position[1], 1.);
            } else if (analysis->histo->ntuData.XYZKMTop[0]>-20 && analysis->histo->ntuData.XYZKMTop[0] < 0) {
              analysis->histo->fillHisto("1404", position[1], 1.);
            } else if (analysis->histo->ntuData.XYZKMTop[0]>-40 && analysis->histo->ntuData.XYZKMTop[0] < -20) {
              analysis->histo->fillHisto("1405", position[1], 1.);
            } else if (analysis->histo->ntuData.XYZKMTop[0]>-60 && analysis->histo->ntuData.XYZKMTop[0] < -40) {
              analysis->histo->fillHisto("1406", position[1], 1.);
            }

            if (analysis->histo->ntuData.XYZKMTop[2] > -30 && analysis->histo->ntuData.XYZKMTop[2] < 30) {
              if (analysis->histo->ntuData.XYZKMTop[0] > 0 && analysis->histo->ntuData.XYZKMTop[0] < 20) {
                analysis->histo->fillHisto("140001", position[1], 1.);
              } else if (analysis->histo->ntuData.XYZKMTop[0] > 20 && analysis->histo->ntuData.XYZKMTop[0] < 40) {
                analysis->histo->fillHisto("140002", position[1], 1.);
              } else if (analysis->histo->ntuData.XYZKMTop[0] > 40 && analysis->histo->ntuData.XYZKMTop[0] < 60) {
                analysis->histo->fillHisto("140003", position[1], 1.);
              } else if (analysis->histo->ntuData.XYZKMTop[0] > -20 && analysis->histo->ntuData.XYZKMTop[0] < 0) {
                analysis->histo->fillHisto("140004", position[1], 1.);
              } else if (analysis->histo->ntuData.XYZKMTop[0] > -40 && analysis->histo->ntuData.XYZKMTop[0] < -20) {
                analysis->histo->fillHisto("140005", position[1], 1.);
              } else if (analysis->histo->ntuData.XYZKMTop[0] > -60 && analysis->histo->ntuData.XYZKMTop[0] < -40) {
                analysis->histo->fillHisto("140006", position[1], 1.);
              }
            }
          } else {
            analysis->histo->ntuData.XYZstopK[0] = -900;
            analysis->histo->ntuData.XYZstopK[1] = -900;
            analysis->histo->ntuData.XYZstopK[2] = -900;
          }

          if (matname == "deuterium") {
            analysis->histo->fillHisto("15", position[1], 1.);
            analysis->histo->fillHisto3("2", position[0], position[2], position[1], 1.);
          } else if (matname == "mylar") {
            analysis->histo->fillHisto("16", position[1], 1.);
            analysis->histo->fillHisto3("3", position[0], position[2], position[1], 1.);
          }
        }
      }
      prev_evid = evid;
    }
  } else if (partname == "kaon+") {
    G4double KeneKP = gst1->GetTrack()->GetKineticEnergy();
    volname = gst1->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    if (KeneKP <= 1001.*eV) {
      G4int evid_KP = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
      if (evid_KP != prev_evid_KP) {
        G4ThreeVector posKPstop = gst1->GetTrack()->GetPosition();
        if (analysis->YesHistos) {
          if (volname == "Target") {
            analysis->histo->ntuData.XYZstopKP[0] = posKPstop[0];
            analysis->histo->ntuData.XYZstopKP[1] = posKPstop[1];
            analysis->histo->ntuData.XYZstopKP[2] = posKPstop[2];
          } else {
            analysis->histo->ntuData.XYZstopKP[0] = -900;
            analysis->histo->ntuData.XYZstopKP[1] = -900;
            analysis->histo->ntuData.XYZstopKP[2] = -900;
          }
        }
      }
      prev_evid_KP = evid_KP;
    }
  }
}
