#include "../include/SiddhartaAnalysisManager.h"
#include "../include/SiddhartaHisto.h"

#include <G4TrajectoryContainer.hh>
#include <G4UnitsTable.hh>
#include <G4Trajectory.hh>

#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>

SiddhartaAnalysisManager* SiddhartaAnalysisManager::fManager = 0;

SiddhartaAnalysisManager* SiddhartaAnalysisManager::getInstance()
{
  if(!fManager) {
    fManager = new SiddhartaAnalysisManager();
  }
  return fManager;
}

void SiddhartaAnalysisManager::dispose()
{
  delete fManager;
  fManager = 0;
}

SiddhartaAnalysisManager::SiddhartaAnalysisManager()
{
  YesHistos = 1;
  verbose = 0;
  nEvt1 = -1;
  nEvt2 = -1;
  targetThresE = 10*keV;
  detectorThresE = 10*keV;
  pulseWidth = 1.*microsecond;
  histo = new SiddhartaHisto();
  mycard = SiddhartaCard::getInstance();
  iniKaonMomCond = mycard->variables["initialKaonMomentumCheck"];
  bookHisto();
}

SiddhartaAnalysisManager::~SiddhartaAnalysisManager()
{
//??????????????????? -> Histo always are created in constructor, therefore does not need to check those

#ifdef G4ANALYSIS_USE
#define HISTDELETE
#endif

#ifdef G4ANALYSIS_USE_ROOT
#define HISTDELETE
#endif

#ifdef HISTDELETE
  delete histo;
#endif
}

/*
Histograms Part to separate file -> HistoCreator managed by HistoManager. AnalysisManager removed or renamed
 */
void SiddhartaAnalysisManager::bookHisto()
{
  histEMax = 15.0*MeV;
  histEMin = .0*MeV;
  histNBin = 100;

  if (YesHistos) {
//  histo->add1D("0","Energy deposit (MeV) in the traget",histNBin,histEMin,histEMax,MeV);
    histo->add1D("0","Vertex X (mm)",10000,-100.*mm,100.*mm,mm);
    histo->add1D("1","Vertex Y (um)",1000,-1500.*um,1500.*um,um);
    histo->add1D("2","Vertex Z (cm)",1000,-500.*cm,500.*cm,cm);
    histo->add3D("4","X vs Y vs Z vertex ",100,-500.*cm,500.*cm,1,1000,-100.*mm,100.*mm,1,100,-50.*mm,50.*mm,1);
    histo->add1D("3","Phi Kinetic Energy",2000,0.,20.*MeV,MeV);
    histo->add1D("4","Al BP hit Kaon- Kinetic Energy",300,0.*MeV,30.*MeV,MeV);
    histo->add1D("5","Kaon- stopped Y",1000,0.*cm,100.*cm,cm);
    histo->add1D("6","SDD Deposit Energy",200000,0.,2000000.*eV,eV);
    histo->add1D("7","KaonMonitorTop Deposit Energy",2000,0.,20000000.*eV,eV);
    histo->add1D("71","KaonMonitorTop Deposit Energy",2000,0.,20000000.*eV,eV);
    histo->add1D("72","KaonMonitorTop Deposit Energy > 3MeV/c",2000,0.,20000000.*eV,eV);
    histo->add1D("73","KaonMonitorTop Deposit Energy > 3MeV/c && EnergyBottom > 3 MeV/c",2000,0.,20000000.*eV,eV);
    histo->add1D("74","KaonMonitorTop Deposit Energy > 3MeV/c && EnergyBottom > 3 MeV/c && sqrt(XKM**2+YKM**2)<45",2000,0.,20000000.*eV,eV);
    histo->add1D("75","KaonMonitorTop Deposit Energy > 3MeV/c && EnergyBottom > 3 MeV/c && sqrt(XKM**2+YKM**2)<40",2000,0.,20000000.*eV,eV);
    histo->add1D("76","KaonMonitorTop Deposit Energy > 3MeV/c && EnergyBottom > 3 MeV/c && sqrt(XKM**2+YKM**2)<35",2000,0.,20000000.*eV,eV);
    histo->add1D("77","KaonMonitorTop Deposit Energy > 3MeV/c && EnergyBottom > 3 MeV/c && sqrt(XKM**2+YKM**2)<30",2000,0.,20000000.*eV,eV);
    histo->add1D("78","KaonMonitorTop Deposit Energy > 3MeV/c && EnergyBottom > 3 MeV/c && sqrt(XKM**2+YKM**2)<25",2000,0.,20000000.*eV,eV);
    histo->add1D("8","KaonMonitorBottom Deposit Energy",2000,0.,20000000.*eV,eV);
    histo->add1D("81","KaonMonitorBottom Deposit Energy",2000,0.,20000000.*eV,eV);
    histo->add1D("82","KaonMonitorBottom Deposit Energy > 3MeV/c",2000,0.,20000000.*eV,eV);
    histo->add1D("9", "KaonMonitorTop Time",100,0.,4.,ns);
    histo->add1D("10", "KaonMonitorTop EnergyDep",100,0.,10000000.,1.);
    histo->add1D("11", "KaonMonitorTop Time",100,0.,4.,ns);
    histo->add1D("12","KaonMonitorTop Y Position",160,-80.,80.,mm);
    histo->add1D("13", "KaonMonitorTop Time",100,0.,4.,ns);
    histo->add1D("14","Kaon- stopped Y in Target",10000,0.*mm,500.*mm,mm);
    histo->add1D("1401","Kaon- stopped Y in Hydrogen for X in [0,20]mm",10000,0.*mm,500.*mm,mm);
    histo->add1D("1402","Kaon- stopped Y in Hydrogen for X in [20,40]mm",10000,0.*mm,500.*mm,mm);
    histo->add1D("1403","Kaon- stopped Y in Hydrogen for X in [40,60]mm",10000,0.*mm,500.*mm,mm);
    histo->add1D("1404","Kaon- stopped Y in Hydrogen for X in [-20,0]mm",10000,0.*mm,500.*mm,mm);
    histo->add1D("1405","Kaon- stopped Y in Hydrogen for X in [-40,-20]mm",10000,0.*mm,500.*mm,mm);
    histo->add1D("1406","Kaon- stopped Y in Hydrogen for X in [-60,-40]mm",10000,0.*mm,500.*mm,mm);
    histo->add1D("140001","Kaon- stopped Y in Hydrogen for X in [0,20]mm and Z in [-30,30]mm",1000,0.*cm,100.*cm,cm);
    histo->add1D("140002","Kaon- stopped Y in Hydrogen for X in [20,40]mm and Z in [-30,30]mm",1000,0.*cm,100.*cm,cm);
    histo->add1D("140003","Kaon- stopped Y in Hydrogen for X in [40,60]mm and Z in [-30,30]mm",1000,0.*cm,100.*cm,cm);
    histo->add1D("140004","Kaon- stopped Y in Hydrogen for X in [-20,0]mm and Z in [-30,30]mm",1000,0.*cm,100.*cm,cm);
    histo->add1D("140005","Kaon- stopped Y in Hydrogen for X in [-40,-20]mm and Z in [-30,30]mm",1000,0.*cm,100.*cm,cm);
    histo->add1D("140006","Kaon- stopped Y in Hydrogen for X in [-60,-40]mm and Z in [-30,30]mm",1000,0.*cm,100.*cm,cm);
    histo->add1D("15","Kaon- stopped Y in deuterium",1000,0.*cm,100.*cm,cm);
    histo->add1D("16","Kaon- stopped Y in mylar",1000,0.*cm,100.*cm,cm);
    histo->add1D("20","Kminus sinetheta at BP", 100,0.,1.,1.);
    histo->add3D("1","Kaon- gas target stop",100,-100.,100.,1,100,-100.,100.,1,300,100.,400.,1);
    histo->add3D("8","Kaon- gas target stop w kcoin",100,-100.,100.,1,100,-100.,100.,1,300,100.,400.,1);
    histo->add3D("9","Kaon+ gas target stop w kcoin",100,-100.,100.,1,100,-100.,100.,1,300,100.,400.,1);
    histo->add3D("5","Kaon- solid target stop",    48, -60., 60.,1, 48, -60., 60.,1,2600,170.,300.,1); // for mylar target check
    histo->add3D("2","Kaon- stopped deuterium",100,-100.,100.,1,100,-100.,100.,1,300,100.,400.,1);
    histo->add3D("3","Kaon- stopped mylar    ",100,-100.,100.,1,100,-100.,100.,1,300,100.,400.,1);
    histo->addTuple("0", "siddharta", "XYZvertex[3]/D:nHitSDD/I:EnergySDD[nHitSDD]:D");

    histo->createHistogramWithAxes(new TH2D("kaonAbsorptionEnergies_vs_Nucleus", "Energy of the generated X rays vs the nucleus atomic number",
                                            2600, 4000., 30000.*eV, 70, 0.5, 70.5), "Energy [eV]", "Atomic number");
    histo->createHistogramWithAxes(new TH1D("kaonParticles", "Number of events with kaon produced", 2, 0.5, 2.5), "Counts", "Control");
    histo->createHistogramWithAxes(new TH1D("materialIDmap", "Map of the material ID for translation", 300, 0.5, 300.5), "Material name", "Counts");

    histo->createHistogramWithAxes(new TH1D("KaonPlusDetector_time", "Time of the hits in the kaon plus detector", 250, -5, 20), "Time [ns]", "Counts");
    histo->createHistogramWithAxes(new TH1D("Veto2_time", "Time of the hits in the Veto2", 100, 0.1, 100.1), "Time [ns]", "Counts");
  }
}

void SiddhartaAnalysisManager::BeginOfRun()
{
  histo->book();
  G4cout << "SiddhartaAnalysisManager: Histograms are booked and the run has been started" << G4endl;
}

void SiddhartaAnalysisManager::EndOfRun()
{
  histo->save();
}

void SiddhartaAnalysisManager::BeginOfEvent()
{
  histo->ntuData.nHitSDD = 0;
  histo->ntuData.NbAnti = 0;
  for (G4int i=0;i<MaxHits;i++) {
    for (G4int j=0;j<20;j++) {
      histo->ntuData.parentName[i][j] = '\0';
    }
  }
  for (G4int i=0;i<MaxAnti;i++) {
    for (G4int j=0;j<20;j++) {
      histo->ntuData.parentNameAnti[i][j] = '\0';
    }
  }
}



void SiddhartaAnalysisManager::EndOfEvent(const G4Event* evt)
{
  if ( histo->ntuData.EnergyDepKMTop > 0 ) histo->fillHisto("71",histo->ntuData.EnergyDepKMTop*eV,1.);
  if ( histo->ntuData.EnergyDepKMBottom > 0 ) histo->fillHisto("81",histo->ntuData.EnergyDepKMBottom*eV,1.);
  if ( histo->ntuData.EnergyDepKMTop > 3000000 ) {
    histo->fillHisto("72",histo->ntuData.EnergyDepKMTop*eV,1.);
    if ( histo->ntuData.EnergyDepKMBottom > 3000000 ) {
      histo->fillHisto("73",histo->ntuData.EnergyDepKMTop*eV,1.);
      if ( sqrt(histo->ntuData.XYZKMTop[0]*histo->ntuData.XYZKMTop[0]+histo->ntuData.XYZKMTop[2]*histo->ntuData.XYZKMTop[2]) < 45 ) {
        histo->fillHisto("74",histo->ntuData.EnergyDepKMTop*eV,1.);
      }
      if ( sqrt(histo->ntuData.XYZKMTop[0]*histo->ntuData.XYZKMTop[0]+histo->ntuData.XYZKMTop[2]*histo->ntuData.XYZKMTop[2]) < 40 ) {
        histo->fillHisto("75",histo->ntuData.EnergyDepKMTop*eV,1.);
      }
      if ( sqrt(histo->ntuData.XYZKMTop[0]*histo->ntuData.XYZKMTop[0]+histo->ntuData.XYZKMTop[2]*histo->ntuData.XYZKMTop[2]) < 35 ) {
        histo->fillHisto("76",histo->ntuData.EnergyDepKMTop*eV,1.);
      }
      if ( sqrt(histo->ntuData.XYZKMTop[0]*histo->ntuData.XYZKMTop[0]+histo->ntuData.XYZKMTop[2]*histo->ntuData.XYZKMTop[2]) < 30 ) {
        histo->fillHisto("77",histo->ntuData.EnergyDepKMTop*eV,1.);
      }
      if ( sqrt(histo->ntuData.XYZKMTop[0]*histo->ntuData.XYZKMTop[0]+histo->ntuData.XYZKMTop[2]*histo->ntuData.XYZKMTop[2]) < 25 ) {
        histo->fillHisto("78",histo->ntuData.EnergyDepKMTop*eV,1.);
      }
    }
  }
  if (histo->ntuData.EnergyDepKMTop > 3000000) {
    histo->fillHisto("82", histo->ntuData.EnergyDepKMTop*eV, 1.);
  }
  if (histo->ntuData.EnergyDepKMTop > 3e6 &&
      histo->ntuData.EnergyDepKMTop < 7e6 &&
      histo->ntuData.EnergyDepKMBottom > 3e6 &&
      histo->ntuData.EnergyDepKMBottom < 7e6 &&
      histo->ntuData.TimeKMTop > 1.2         &&
      histo->ntuData.TimeKMTop < 2.          &&
      histo->ntuData.TimeKMBottom > 1.2         &&
      histo->ntuData.TimeKMBottom < 2.          &&
      histo->ntuData.XYZKMBottom[0] > -50.    &&
      histo->ntuData.XYZKMBottom[0] <  50.)
  {
    histo->fillHisto3("8",histo->ntuData.XYZstopK[0],histo->ntuData.XYZstopK[2],histo->ntuData.XYZstopK[1],1.); // with kaon trigger cut
    histo->fillHisto3("9",histo->ntuData.XYZstopKP[0],histo->ntuData.XYZstopKP[2],histo->ntuData.XYZstopKP[1],1.); // with kaon trigger cut
  }

 //SiddhartaCard* mycard = SiddhartaCard::getInstance();
  int SiddhartaSetup = mycard->variables["SiddhartaSetupVersion"];

  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;

  if (trajectoryContainer)
    n_trajectories = trajectoryContainer->entries();

  
  bool sddhit = false;
  bool czthit = false;
  bool hpgehit = false;

  int ACTIVE_DET_ID = mycard->variables["activeDetMode"];

  if (histo->ntuData.nHitSDD > 0 || histo->ntuData.NbAnti > 0 || histo->ntuData.EnergyDepKMTop > 0 ||
      histo->ntuData.EnergyDepKMBottom > 0 || histo->ntuData.EnergyDepLMAntiBoost > 0 || histo->ntuData.EnergyDepLMBoost > 0) 
  {
    int Anti_size = 9;

    if (SiddhartaSetup == 2 || SiddhartaSetup == 6)
      Anti_size = 9;

    if (SiddhartaSetup == 7 || SiddhartaSetup == 8 || SiddhartaSetup == 2020 || SiddhartaSetup == 2023)
      Anti_size = 17;

    G4String pname ="";
    for (G4int j=0; j<histo->ntuData.nHitSDD; j++) {
      for (G4int i=0; i<n_trajectories; i++) {
        G4Trajectory* trj = (G4Trajectory*)(*(trajectoryContainer))[i];
        pname = trj->GetParticleName();
        int ilen = strlen(pname);

        if (trj->GetTrackID() == histo->ntuData.parentID[j]) {
          for(G4int k=0; k<ilen; k++) {
            histo->ntuData.parentName[j][k] = (pname.data())[k];
          }
        histo->ntuData.parentName[j][ilen] = '\0';
        }
      }
    }

    for (G4int i=0;i<n_trajectories;i++) {
      G4Trajectory* trj = (G4Trajectory*)(*(trajectoryContainer))[i];
      pname = trj->GetParticleName();
      int ilen = strlen(pname);

      for (G4int l=0; l<Anti_size; l++) {
        if (trj->GetTrackID() == histo->ntuData.parentIDAnti[l]) {
          for(G4int k=0;k<ilen;k++) {
            histo->ntuData.parentNameAnti[l][k] = (pname.data())[k];
          }
        histo->ntuData.parentNameAnti[l][ilen] = '\0';
        }
      }
    }

    sddhit = true;
  }

  if ((histo->ntuData.EnergyDepLMAntiBoost>-1) || (histo->ntuData.EnergyDepKLDAntiBoost>0) || (histo->ntuData.EnergyDepKLT1AntiBoost>0) || 
     (histo->ntuData.EnergyDepKLT2AntiBoost>0) || (histo->ntuData.EnergyDepKLBWAntiBoost>0) || (histo->ntuData.TotalEnergyDepKLCZTAntiBoost>0))
    czthit = true;

  if ((histo->ntuData.EnergyDepLMBoost>0) || (histo->ntuData.EnergyDepKLDBoost>0) || (histo->ntuData.EnergyDepKLTBoost>0) || (histo->ntuData.TotalEnergyDepKLHPGeBoost>0))
    hpgehit = true;

  switch(ACTIVE_DET_ID) {
    case 1:
      if(sddhit)
        histo->fillTuple("0");
        break;
    case 2:
      if (czthit)
        histo->fillTuple("0");
      break;
    case 3:
      if (hpgehit)
        histo->fillTuple("0");
      break;
    case 4:
      if (sddhit || czthit || hpgehit)
        histo->fillTuple("0");
      break;
    case 5:
      if (sddhit || czthit)
        histo->fillTuple("0");
      break;
    case 6:
      if (sddhit || hpgehit)
        histo->fillTuple("0");
      break;
    case 7:
      if (czthit || hpgehit)
        histo->fillTuple("0");
      break;
  }

  G4double cut1 = 10000.;
  G4double cut2 = 0.8;
  G4double cut3 = 1.05;
  G4double cut4 = 3000000;
  if (SiddhartaSetup == 1) {
    cut1 = 10000;
    cut2 = 0.75;
    cut3 = 1.20;
    cut4 = 3000000;
  } else if (SiddhartaSetup == 2) {
    cut1 = 45.; // Kaon detector on the top of the shielding
    cut2 = 1.6;
    cut3 = 2.10;
    cut4 = 3000000;
  } else if (SiddhartaSetup == 3) {
    cut1 = 45.;
    cut2 = 1.3;
    cut3 = 1.90;
    cut4 = 3000000;
  } else if (SiddhartaSetup == 6) {
    cut1 = 45.;
    cut2 = 1.6;
    cut3 = 2.25;
    cut4 = 3000000;
  } else {
    cut1 = 45.;
    cut2 = 1.3;
    cut3 = 1.90;
    cut4 = 3000000;
  }

  if (sqrt(histo->ntuData.XYZKMTop[0]*histo->ntuData.XYZKMTop[0] + histo->ntuData.XYZKMTop[2]*histo->ntuData.XYZKMTop[2]) < cut1) {
    histo->fillHisto("9", histo->ntuData.TimeKMTop, 1.);
    histo->fillHisto("10", histo->ntuData.EnergyDepKMTop, 1.);
    if (histo->ntuData.EnergyDepKMTop > cut4) {
      histo->fillHisto("11", histo->ntuData.TimeKMTop, 1.);
      if (histo->ntuData.TimeKMTop > cut2 && histo->ntuData.TimeKMTop < cut3) {
        histo->fillHisto("12", histo->ntuData.XYZKMTop[2], 1.);
        if (SiddhartaSetup == 6) {
          bool triggerAnti = 0;
          G4double cutT1 = 5.;
          G4double cutT2 = 8.;
          G4double cutT3 = 3000000.;

//Create some function that will simplify these code -> Cuts should be defined or loaded before or elsewhere
// Function to check the value in range
          if ((histo->ntuData.TimeAnti[0] > cutT1 && histo->ntuData.TimeAnti[0] < cutT2 && histo->ntuData.EnergyAnti[0] > cutT3) ||
              (histo->ntuData.TimeAnti[1] > cutT1 && histo->ntuData.TimeAnti[1] < cutT2 && histo->ntuData.EnergyAnti[1] > cutT3) ||
              (histo->ntuData.TimeAnti[2] > cutT1 && histo->ntuData.TimeAnti[2] < cutT2 && histo->ntuData.EnergyAnti[2] > cutT3) ||
              (histo->ntuData.TimeAnti[3] > cutT1 && histo->ntuData.TimeAnti[3] < cutT2 && histo->ntuData.EnergyAnti[3] > cutT3) ||
              (histo->ntuData.TimeAnti[4] > cutT1 && histo->ntuData.TimeAnti[4] < cutT2 && histo->ntuData.EnergyAnti[4] > cutT3) ||
              (histo->ntuData.TimeAnti[5] > cutT1 && histo->ntuData.TimeAnti[5] < cutT2 && histo->ntuData.EnergyAnti[5] > cutT3) ||
              (histo->ntuData.TimeAnti[6] > cutT1 && histo->ntuData.TimeAnti[6] < cutT2 && histo->ntuData.EnergyAnti[6] > cutT3) ||
              (histo->ntuData.TimeAnti[7] > cutT1 && histo->ntuData.TimeAnti[7] < cutT2 && histo->ntuData.EnergyAnti[7] > cutT3))
            triggerAnti = 1;

          if (triggerAnti) {
            histo->fillHisto("13", histo->ntuData.TimeKMTop, 1.);
          }
        }
      }
    }
  }

}

void SiddhartaAnalysisManager::AddEnergy(G4double edep, G4double weight, G4double time)
{
//??
/*
  if(1 < verbose) {
    G4cout << "SiddhartaAnalysisManager::AddEnergy: e(keV)= " << edep/keV
	   << " weight = " << weight << " time (s) = " <<  time/second
           << G4endl;
  }
  histo->fillTuple(2, 0, edep/MeV);
  histo->fillTuple(2,1,time/second);
  histo->fillTuple(2,2,weight);
  histo->addRow(2);
  //
  SiddhartaEnergyDeposition A(edep,time,weight);
  Edepo.push_back(A);
  */
}

void SiddhartaAnalysisManager::AddParticle(G4double pid, G4double energy, G4double weight, G4double time )
{
  if(1 < verbose) {
    G4cout << "SiddhartaAnalysisManager::AddParticle: " << pid << G4endl;
  }
  histo->fillTuple(0,0, pid);
  histo->fillTuple(0,1,energy/MeV);
  histo->fillTuple(0,2,time/second);
  histo->fillTuple(0,3,weight);
  histo->addRow(0);
}

void SiddhartaAnalysisManager::AddIsotope(G4double pid,G4double weight, G4double time )
{
  if(1 < verbose) {
    G4cout << "SiddhartaAnalysisManager::AddIsotope: " << pid << G4endl;
  }
  histo->fillTuple(1,0,pid);
  histo->fillTuple(1,1,time/second);
  histo->fillTuple(1,2,weight);
  histo->addRow(1);
}

bool SiddhartaAnalysisManager::CheckKaonIniMom(G4ThreeVector iniPos, G4ThreeVector momentum)
{
  int condition = true;
  if (iniKaonMomCond) {
    std::vector<G4ThreeVector> lumiPlane;
    if (iniKaonMomCond < 0) {
      if (momentum.x() > 0)
        return false;

      lumiPlane = mycard->GetLumiAntiBoostPosition();
    } else {
      if (momentum.x() < 0)
        return false;

      lumiPlane = mycard->GetLumiBoostPosition();
    }

    if (lumiPlane.size() == 2) {
      G4double lumiPlaneXMax = lumiPlane.at(0).x();
      G4double distanceFromIniPos = lumiPlaneXMax - iniPos.x();
      G4double mult = distanceFromIniPos/momentum.x();
      G4ThreeVector potentialKaonPos = iniPos + mult*momentum;

      if (!(potentialKaonPos.y() > lumiPlane.at(0).y() && potentialKaonPos.y() < lumiPlane.at(1).y() &&
          potentialKaonPos.z() > lumiPlane.at(0).z() && potentialKaonPos.z() < lumiPlane.at(1).z()))
        condition = false;
    }
  }
  return condition;
}
