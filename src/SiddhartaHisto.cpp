#include "../include/SiddhartaHisto.h"
#include "../include/SiddhartaCard.h"

#include <G4ParticleTable.hh>
#include <G4Tokenizer.hh>

#include <TApplication.h>
#include <TGClient.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TBranch.h>
#include <TTree.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>

SiddhartaHisto::SiddhartaHisto()
{
  verbose = 1;
  histName = "siddharta";
  histType = "root";
  nHisto = 0;
  nHisto3 = 0;
  nTuple = 0;
  defaultAct = 1;

  ROOThisto.clear();
  ROOThisto3.clear();
  ROOTntup.clear();
  Rarray.clear();
  Rcol.clear();
  active.clear();
  bins.clear();
  xmin.clear();
  xmax.clear();
  unit.clear();
  ids.clear();
  titles.clear();
  tupleName.clear();
  tupleId.clear();
  tupleList.clear();
  tupleListROOT.clear();
  active3.clear();
  bins31.clear();
  xmin31.clear();
  xmax31.clear();
  unit31.clear();
  bins32.clear();
  xmin32.clear();
  xmax32.clear();
  unit32.clear();
  bins33.clear();
  xmin33.clear();
  xmax33.clear();
  unit33.clear();
  ids3.clear();
  titles3.clear();
}

SiddhartaHisto::~SiddhartaHisto()
{
  //FIXME : G.Barrand : the below is crashy.
  //        In principle the TH are deleted
  //        when doing the TFile::Close !
  //         In fact the hfileROOT should
  //        be deleted in save(). And I am pretty
  //        sure that the TApplication is not needed.
  //
  // removed by F.Lei
  //  for(G4int i=0; i<nHisto; i++) {
  //   if(ROOThisto[i]) delete ROOThisto[i];
  // }
}

void SiddhartaHisto::book()
{
  SiddhartaCard* mycard = SiddhartaCard::getInstance();
  int SiddhartaSetup = mycard->variables["SiddhartaSetupVersion"];

  histName = histName + "_" + mycard->GetCardFileName();

  std::stringstream ss;
  ss << mycard->GetG4RandGaussSeed();
  G4String fileNameROOT = histName + ss.str() + G4String(".root");

  hfileROOT = new TFile(fileNameROOT.c_str() ,"RECREATE","ROOT file for Siddharta");
  G4cout << "Root file: " << fileNameROOT << G4endl;

  for (G4int i=0; i<nHisto; i++) {
    if (active[i]) {
      G4String id = G4String("h") + ids[i];
      ROOThisto[i] = new TH1D(id, titles[i], bins[i], xmin[i], xmax[i]);
      G4cout << "ROOT Histo " << ids[i] << " " << titles[i] << " booked " << G4endl;
    }
  }

  for (G4int i=0; i<nHisto3; i++) {
    if (active3[i]) {
      G4String id = G4String("h3d")+ids3[i];
      ROOThisto3[i] = new TH3D(id, titles3[i], bins31[i], xmin31[i], xmax31[i], bins32[i], xmin32[i], xmax32[i], bins33[i], xmin33[i], xmax33[i]);
      G4cout << "ROOT Histo3 " << ids3[i] << " " << titles3[i] << " booked " << G4endl;
    }
  }

  for (G4int i=0; i<nTuple; i++) {
    if (tupleListROOT[i] != "") {
      G4String id = G4String("t")+tupleId[i];
      G4cout << "Creating Ntuple " << tupleId[i] << " in ROOT file: " << tupleName[i] << G4endl;
      ROOTntup[i] = new TTree(id, tupleName[i]);
      ROOTntup[i]->Branch("Vertex", ntuData.XYZvertex, "XYZvertex[3]/D:VertexMomentum[3]");
      ROOTntup[i]->Branch("nHitSDD", &ntuData.nHitSDD, "nHitSDD/I");
      ROOTntup[i]->Branch("NoSDD", ntuData.NoSDD, "NoSDD[nHitSDD]/I");
      ROOTntup[i]->Branch("NbHitSDD", ntuData.NbHitSDD, "NbHitSDD[nHitSDD]/I");
      ROOTntup[i]->Branch("EnergySDD", ntuData.EnergySDD, "EnergySDD[nHitSDD]/D");
      ROOTntup[i]->Branch("TimeSDD", ntuData.TimeSDD, "TimeSDD[nHitSDD]/D");
      ROOTntup[i]->Branch("XYZvertexSDD", ntuData.XYZvertexSDD, "XYZvertexSDD[nHitSDD][3]/D");
      ROOTntup[i]->Branch("XYZSDD", ntuData.XYZSDD, "XYZSDD[nHitSDD][3]/D");
      ROOTntup[i]->Branch("XYZInSDD", ntuData.XYZInSDD, "XYZInSDD[nHitSDD][3]/D");
      ROOTntup[i]->Branch("KinvertexSDD", ntuData.KinvertexSDD, "KinvertexSDD[nHitSDD]/D");
      ROOTntup[i]->Branch("particleVertexPDGEncoding", ntuData.particleVertexPDGEncoding, "particleVertexPDGEncoding[nHitSDD]/I");
      ROOTntup[i]->Branch("matvert",ntuData.matvert,"matvert[nHitSDD]/l");
      ROOTntup[i]->Branch("matVertIDFromMap",ntuData.matVertIDFromMap,"matVertIDFromMap[nHitSDD]/l");
      ROOTntup[i]->Branch("materialVertex", ntuData.materialVertex, "materialVertex[nHitSDD][20]/C");
      ROOTntup[i]->Branch("parentName", ntuData.parentName, "parentName[nHitSDD][20]/C");
      ROOTntup[i]->Branch("TimeKMTop", &ntuData.TimeKMTop, "TimeKMTop/D");
      ROOTntup[i]->Branch("EnergyDepKMTop", &ntuData.EnergyDepKMTop, "EnergyDepKMTop/D");
      ROOTntup[i]->Branch("XYZKMTop", ntuData.XYZKMTop, "XYZKMTop[3]/D");
      ROOTntup[i]->Branch("particleNameKMTop", ntuData.particleNameKMTop, "particleNameKMTop[20]/C");
      ROOTntup[i]->Branch("TimeKMBottom", &ntuData.TimeKMBottom, "TimeKMBottom/D");
      ROOTntup[i]->Branch("EnergyKMBottom", &ntuData.EnergyDepKMBottom, "EnergyDepKMBottom/D");
      ROOTntup[i]->Branch("XYZKMBottom", ntuData.XYZKMBottom, "XYZKMBottom[3]/D");
      ROOTntup[i]->Branch("particleNameKMBottom", ntuData.particleNameKMBottom, "particleNameKMBottom[20]/C");

      ROOTntup[i]->Branch("TimeKPlusDetector", &ntuData.TimeKPlusDetector, "TimeKPlusDetector/D");
      ROOTntup[i]->Branch("EnergyKPlusDetector", &ntuData.EnergyDepKPlusDetector, "EnergyDepKPlusDetector/D");
      ROOTntup[i]->Branch("XYZKPlusDetector", ntuData.XYZKPlusDetector, "XYZKPlusDetector[3]/D");
      ROOTntup[i]->Branch("particleNameKPlusDetector", ntuData.particleNameKPlusDetector, "particleNameKPlusDetector[20]/C");

      ROOTntup[i]->Branch("TimeGhost", &ntuData.TimeGhost, "TimeGhost/D");
      ROOTntup[i]->Branch("EnergyDepGhost", &ntuData.EnergyDepGhost, "EnergyDepGhost/D");
      ROOTntup[i]->Branch("XYZGhost", ntuData.XYZGhost, "XYZGhost[3]/D");
      ROOTntup[i]->Branch("KineticGhost", &ntuData.KineticGhost, "KineticGhost/D");

      if (SiddhartaSetup == 2 || SiddhartaSetup == 6) {
        ROOTntup[i]->Branch("NbAnti", &ntuData.NbAnti, "NbAnti/I");
        ROOTntup[i]->Branch("TimeAnti", ntuData.TimeAnti, "TimeAnti[9]/D");
        ROOTntup[i]->Branch("EnergyAnti", ntuData.EnergyAnti, "EnergyAnti[9]/D");
        ROOTntup[i]->Branch("XYZAnti", ntuData.XYZAnti, "XYZAnti[9][3]/D");
        ROOTntup[i]->Branch("XYZVertexAnti", ntuData.XYZVertexAnti, "XYZVertexAnti[9][3]/D");
        ROOTntup[i]->Branch("KinVertexAnti", ntuData.KinVertexAnti, "KinVertexAnti[9]/D");
        ROOTntup[i]->Branch("particleAntiPDGEncoding", ntuData.particleAntiPDGEncoding, "particleAntiPDGEncoding[9]/I");
        ROOTntup[i]->Branch("parentNameAnti", ntuData.parentNameAnti, "parentNameAnti[9][20]/C");
        ROOTntup[i]->Branch("matvertAnti",ntuData.matvertAnti,"matvertAnti[17]/l");
        ROOTntup[i]->Branch("materialVertexAnti", ntuData.materialVertexAnti, "materialVertexAnti[17][20]/C");
        ROOTntup[i]->Branch("nHitSciAnti", &ntuData.nHitSciAnti, "nHitSciAnti/I");
        ROOTntup[i]->Branch("TimeScintAnti", ntuData.TimeScintAnti, "TimeScintAnti[nHitSciAnti]/D");
        ROOTntup[i]->Branch("EnergyScintAnti", ntuData.EnergyScintAnti, "EnergyScintAnti[nHitSciAnti]/D");
        ROOTntup[i]->Branch("copyScintAnti", ntuData.copyScintAnti, "copyScintAnti[nHitSciAnti]/I");
      } else if (SiddhartaSetup == 7 || SiddhartaSetup == 8) {
        ROOTntup[i]->Branch("NbAnti", &ntuData.NbAnti, "NbAnti/I");
        ROOTntup[i]->Branch("TimeAnti", ntuData.TimeAnti, "TimeAnti[17]/D");
        ROOTntup[i]->Branch("EnergyAnti", ntuData.EnergyAnti, "EnergyAnti[17]/D");
        ROOTntup[i]->Branch("XYZAnti", ntuData.XYZAnti, "XYZAnti[17][3]/D");
        ROOTntup[i]->Branch("XYZVertexAnti", ntuData.XYZVertexAnti, "XYZVertexAnti[17][3]/D");
        ROOTntup[i]->Branch("KinVertexAnti", ntuData.KinVertexAnti, "KinVertexAnti[17]/D");
        ROOTntup[i]->Branch("particleAntiPDGEncoding", ntuData.particleAntiPDGEncoding, "particleAntiPDGEncoding[17]/I");
        ROOTntup[i]->Branch("parentNameAnti", ntuData.parentNameAnti, "parentNameAnti[17][20]/C");
        ROOTntup[i]->Branch("matvertAnti",ntuData.matvertAnti,"matvertAnti[17]/l");
        ROOTntup[i]->Branch("materialVertexAnti", ntuData.materialVertexAnti, "materialVertexAnti[17][20]/C");
        ROOTntup[i]->Branch("nHitSciAnti", &ntuData.nHitSciAnti, "nHitSciAnti/I");
        ROOTntup[i]->Branch("TimeScintAnti", ntuData.TimeScintAnti, "TimeScintAnti[nHitSciAnti]/D");
        ROOTntup[i]->Branch("EnergyScintAnti", ntuData.EnergyScintAnti, "EnergyScintAnti[nHitSciAnti]/D");
        ROOTntup[i]->Branch("copyScintAnti", ntuData.copyScintAnti, "copyScintAnti[nHitSciAnti]/I");
      } else if (SiddhartaSetup == 2020 || SiddhartaSetup == 2023) { // 17 -> MaxAnti ??
	      ROOTntup[i]->Branch("NbAnti", &ntuData.NbAnti, "NbAnti/I");
	      ROOTntup[i]->Branch("TimeAnti", ntuData.TimeAnti, "TimeAnti[17]/D");
	      ROOTntup[i]->Branch("EnergyAnti", ntuData.EnergyAnti, "EnergyAnti[17]/D");
	      ROOTntup[i]->Branch("XYZAnti", ntuData.XYZAnti, "XYZAnti[17][3]/D");
	      ROOTntup[i]->Branch("XYZVertexAnti", ntuData.XYZVertexAnti, "XYZVertexAnti[17][3]/D");
	      ROOTntup[i]->Branch("KinVertexAnti", ntuData.KinVertexAnti, "KinVertexAnti[17]/D");
	      ROOTntup[i]->Branch("particleAntiPDGEncoding", ntuData.particleAntiPDGEncoding, "particleAntiPDGEncoding[17]/I");
	      ROOTntup[i]->Branch("parentNameAnti", ntuData.parentNameAnti, "parentNameAnti[17][20]/C");
	      ROOTntup[i]->Branch("matvertAnti",ntuData.matvertAnti,"matvertAnti[17]/l");
	      ROOTntup[i]->Branch("materialVertexAnti", ntuData.materialVertexAnti, "materialVertexAnti[17][20]/C");
	      ROOTntup[i]->Branch("matVertIDFromMap_Veto1", ntuData.matVertIDFromMap_Veto1, "matVertIDFromMap_Veto1[20]/l");
	      ROOTntup[i]->Branch("nHitSciAnti", &ntuData.nHitSciAnti, "nHitSciAnti/I");
	      ROOTntup[i]->Branch("TimeScintAnti", ntuData.TimeScintAnti, "TimeScintAnti[nHitSciAnti]/D");
	      ROOTntup[i]->Branch("EnergyScintAnti", ntuData.EnergyScintAnti, "EnergyScintAnti[nHitSciAnti]/D");
	      ROOTntup[i]->Branch("copyScintAnti", ntuData.copyScintAnti, "copyScintAnti[nHitSciAnti]/I");
	      ROOTntup[i]->Branch("particleScintAntiPDGEncoding", ntuData.particleScintAntiPDGEncoding, "particleScintAntiPDGEncoding[nHitSciAnti]/I");
	      ROOTntup[i]->Branch("matVertIDFromMap_Veto2", ntuData.matVertIDFromMap_Veto2, "matVertIDFromMap_Veto2[nHitSciAnti]/l");

	      ROOTntup[i]->Branch("TimeLMBoost", &ntuData.TimeLMBoost, "TimeLMBoost/D");
	      ROOTntup[i]->Branch("EnergyDepLMBoost", &ntuData.EnergyDepLMBoost, "EnergyDepLMBoost/D");
	      ROOTntup[i]->Branch("XYZLMBoost", ntuData.XYZLMBoost, "XYZLMBoost[3]/D");
	      ROOTntup[i]->Branch("XYZLMBoostKaonstop", ntuData.XYZLMBoostKaonstop, "XYZLMBoostKaonstop[3]/D");
	      ROOTntup[i]->Branch("particleNameLMBoost", ntuData.particleNameLMBoost, "particleNameLMBoost[20]/C");
	      ROOTntup[i]->Branch("pdgcodeLMBoost", &ntuData.pdgcodeLMBoost, "pdgcodeLMBoost/I");
	      ROOTntup[i]->Branch("kaonKinELMBoost", &ntuData.kaonKinELMBoost, "kaonKinELMBoost/D");
	      ROOTntup[i]->Branch("lastkaonKinELMBoost", &ntuData.lastkaonKinELMBoost, "lastkaonKinELMBoost/D");
	      ROOTntup[i]->Branch("TimeLMAntiBoost", &ntuData.TimeLMAntiBoost, "TimeLMAntiBoost/D");
	      ROOTntup[i]->Branch("EnergyDepLMAntiBoost", &ntuData.EnergyDepLMAntiBoost, "EnergyDepLMAntiBoost/D");
	      ROOTntup[i]->Branch("XYZLMAntiBoost", ntuData.XYZLMAntiBoost, "XYZLMAntiBoost[3]/D");
	      ROOTntup[i]->Branch("XYZLMAntiBoostKaonstop", ntuData.XYZLMAntiBoostKaonstop, "XYZLMAntiBoostKaonstop[3]/D");
	      ROOTntup[i]->Branch("particleNameLMAntiBoost", ntuData.particleNameLMAntiBoost, "particleNameLMAntiBoost[20]/C");
	      ROOTntup[i]->Branch("pdgcodeLMAntiBoost", &ntuData.pdgcodeLMAntiBoost, "pdgcodeLMAntiBoost/I");
	      ROOTntup[i]->Branch("kaonKinELMAntiBoost", &ntuData.kaonKinELMAntiBoost, "kaonKinELMAntiBoost/D");
	      ROOTntup[i]->Branch("lastkaonKinELMAntiBoost", &ntuData.lastkaonKinELMAntiBoost, "lastkaonKinELMAntiBoost/D");

	      ROOTntup[i]->Branch("TimeKLDAntiBoost",&ntuData.TimeKLDAntiBoost,"TimeKLDAntiBoost/D");
	      ROOTntup[i]->Branch("EnergyDepKLDAntiBoost",&ntuData.EnergyDepKLDAntiBoost,"EnergyDepKLDAntiBoost/D");
	      ROOTntup[i]->Branch("XYZKLDAntiBoost",ntuData.XYZKLDAntiBoost,"XYZKLDAntiBoost[3]/D");
	      ROOTntup[i]->Branch("KLDAntiBooststop",ntuData.KLDAntiBooststop,"KLDAntiBooststop[3]/D");
	      ROOTntup[i]->Branch("pdgcodeKLDAntiBoost",&ntuData.pdgcodeKLDAntiBoost,"pdgcodeKLDAntiBoost/I");
	      ROOTntup[i]->Branch("kaonKinEKLDAntiBoost",&ntuData.kaonKinEKLDAntiBoost,"kaonKinEKLDAntiBoost/D");
	      ROOTntup[i]->Branch("lastkaonKinEKLDAntiBoost",&ntuData.lastkaonKinEKLDAntiBoost,"lastkaonKinEKLDAntiBoost/D");

	      ROOTntup[i]->Branch("TimeKLT1AntiBoost",&ntuData.TimeKLT1AntiBoost,"TimeKLT1AntiBoost/D");
	      ROOTntup[i]->Branch("EnergyDepKLT1AntiBoost",&ntuData.EnergyDepKLT1AntiBoost,"EnergyDepKLT1AntiBoost/D");
	      ROOTntup[i]->Branch("XYZKLT1AntiBoost",ntuData.XYZKLT1AntiBoost,"XYZKLT1AntiBoost[3]/D");
	      ROOTntup[i]->Branch("KLT1AntiBooststop",ntuData.KLT1AntiBooststop,"KLT1AntiBooststop[3]/D");
	      ROOTntup[i]->Branch("KLT1AntiBoostKaonstop",ntuData.KLT1AntiBoostKaonstop,"KLT1AntiBoostKaonstop[3]/D");
	      ROOTntup[i]->Branch("pdgcodeKLT1AntiBoost",&ntuData.pdgcodeKLT1AntiBoost,"pdgcodeKLT1AntiBoost/I");
	      ROOTntup[i]->Branch("kaonKinEKLT1AntiBoost",&ntuData.kaonKinEKLT1AntiBoost,"kaonKinEKLT1AntiBoost/D");
          ROOTntup[i]->Branch("gammaKinEKLT1AntiBoost",&ntuData.gammaKinEKLT1AntiBoost,"gammaKinEKLT1AntiBoost/D");
	      ROOTntup[i]->Branch("lastkaonKinEKLT1AntiBoost",&ntuData.lastkaonKinEKLT1AntiBoost,"lastkaonKinEKLT1AntiBoost/D");

	      ROOTntup[i]->Branch("TimeKLT2AntiBoost",&ntuData.TimeKLT2AntiBoost,"TimeKLT2AntiBoost/D");
	      ROOTntup[i]->Branch("EnergyDepKLT2AntiBoost",&ntuData.EnergyDepKLT2AntiBoost,"EnergyDepKLT2AntiBoost/D");
	      ROOTntup[i]->Branch("XYZKLT2AntiBoost",ntuData.XYZKLT2AntiBoost,"XYZKLT2AntiBoost[3]/D");
	      ROOTntup[i]->Branch("KLT2AntiBooststop",ntuData.KLT2AntiBooststop,"KLT2AntiBooststop[3]/D");
	      ROOTntup[i]->Branch("KLT2AntiBoostKaonstop",ntuData.KLT2AntiBoostKaonstop,"KLT2AntiBoostKaonstop[3]/D");
	      ROOTntup[i]->Branch("pdgcodeKLT2AntiBoost",&ntuData.pdgcodeKLT2AntiBoost,"pdgcodeKLT2AntiBoost/I");
	      ROOTntup[i]->Branch("kaonKinEKLT2AntiBoost",&ntuData.kaonKinEKLT2AntiBoost,"kaonKinEKLT2AntiBoost/D");
          ROOTntup[i]->Branch("gammaKinEKLT2AntiBoost",&ntuData.gammaKinEKLT2AntiBoost,"gammaKinEKLT2AntiBoost/D");
	      ROOTntup[i]->Branch("lastkaonKinEKLT2AntiBoost",&ntuData.lastkaonKinEKLT2AntiBoost,"lastkaonKinEKLT2AntiBoost/D");

	      ROOTntup[i]->Branch("TimeKLBWAntiBoost",&ntuData.TimeKLBWAntiBoost,"TimeKLBWAntiBoost/D");
	      ROOTntup[i]->Branch("EnergyDepKLBWAntiBoost",&ntuData.EnergyDepKLBWAntiBoost,"EnergyDepKLBWAntiBoost/D");
	      ROOTntup[i]->Branch("XYZKLBWAntiBoost",ntuData.XYZKLBWAntiBoost,"XYZKLBWAntiBoost[3]/D");
	      ROOTntup[i]->Branch("KLBWAntiBooststop",ntuData.KLBWAntiBooststop,"KLBWAntiBooststop[3]/D");
	      ROOTntup[i]->Branch("pdgcodeKLBWAntiBoost",&ntuData.pdgcodeKLBWAntiBoost,"pdgcodeKLBWAntiBoost/I");
	      ROOTntup[i]->Branch("kaonKinEKLBWAntiBoost",&ntuData.kaonKinEKLBWAntiBoost,"kaonKinEKLBWAntiBoost/D");
	      ROOTntup[i]->Branch("lastkaonKinEKLBWAntiBoost",&ntuData.lastkaonKinEKLBWAntiBoost,"lastkaonKinEKLBWAntiBoost/D");

	      ROOTntup[i]->Branch("TotalEnergyDepKLCZTAntiBoost",&ntuData.TotalEnergyDepKLCZTAntiBoost,"TotalEnergyDepKLCZTAntiBoost/D");
	      ROOTntup[i]->Branch("gammaKinEKLCZTAntiBoost",&ntuData.gammaKinEKLCZTAntiBoost,"gammaKinEKLCZTAntiBoost/D");
	      ROOTntup[i]->Branch("nHitCZTAntiBoost", &ntuData.nHitCZTAntiBoost, "nHitCZTAntiBoost/I");
	      ROOTntup[i]->Branch("TimeKLCZTAntiBoost",ntuData.TimeKLCZTAntiBoost,"TimeKLCZTAntiBoost[nHitCZTAntiBoost]/D");
	      ROOTntup[i]->Branch("EnergyDepKLCZTAntiBoost",ntuData.EnergyDepKLCZTAntiBoost,"EnergyDepKLCZTAntiBoost[nHitCZTAntiBoost]/D");
	      ROOTntup[i]->Branch("XYZKLCZTAntiBoost",ntuData.XYZKLCZTAntiBoost,"XYZKLCZTAntiBoost[3][nHitCZTAntiBoost]/D");
	      ROOTntup[i]->Branch("pdgcodeKLCZTAntiBoost",ntuData.pdgcodeKLCZTAntiBoost,"pdgcodeKLCZTAntiBoost[nHitCZTAntiBoost]/I");
	      ROOTntup[i]->Branch("kaonKinEKLCZTAntiBoost",ntuData.kaonKinEKLCZTAntiBoost,"kaonKinEKLCZTAntiBoost[nHitCZTAntiBoost]/D");

	      ROOTntup[i]->Branch("TimeKLDBoost",&ntuData.TimeKLDBoost,"TimeKLDBoost/D");
	      ROOTntup[i]->Branch("EnergyDepKLDBoost",&ntuData.EnergyDepKLDBoost,"EnergyDepKLDBoost/D");
	      ROOTntup[i]->Branch("XYZKLDBoost",ntuData.XYZKLDBoost,"XYZKLDBoost[3]/D");
	      ROOTntup[i]->Branch("KLDBooststop",ntuData.KLDBooststop,"KLDBooststop[3]/D");
	      ROOTntup[i]->Branch("pdgcodeKLDBoost",&ntuData.pdgcodeKLDBoost,"pdgcodeKLDBoost/I");
	      ROOTntup[i]->Branch("kaonKinEKLDBoost",&ntuData.kaonKinEKLDBoost,"kaonKinEKLDBoost/D");
	      ROOTntup[i]->Branch("lastkaonKinEKLDBoost",&ntuData.lastkaonKinEKLDBoost,"lastkaonKinEKLDBoost/D");

	      ROOTntup[i]->Branch("TimeKLTBoost",&ntuData.TimeKLTBoost,"TimeKLTBoost/D");
	      ROOTntup[i]->Branch("EnergyDepKLTBoost",&ntuData.EnergyDepKLTBoost,"EnergyDepKLTBoost/D");
	      ROOTntup[i]->Branch("XYZKLTBoost",ntuData.XYZKLTBoost,"XYZKLTBoost[3]/D");
	      ROOTntup[i]->Branch("KLTBooststop",ntuData.KLTBooststop,"KLTBooststop[3]/D");
	      ROOTntup[i]->Branch("KLTBoostKaonstop",ntuData.KLTBoostKaonstop,"KLTBoostKaonstop[3]/D");
	      ROOTntup[i]->Branch("pdgcodeKLTBoost",&ntuData.pdgcodeKLTBoost,"pdgcodeKLTBoost/I");
	      ROOTntup[i]->Branch("kaonKinEKLTBoost",&ntuData.kaonKinEKLTBoost,"kaonKinEKLTBoost/D");
          ROOTntup[i]->Branch("gammaKinEKLTBoost",&ntuData.gammaKinEKLTBoost,"gammaKinEKLTBoost/D");
	      ROOTntup[i]->Branch("lastkaonKinEKLTBoost",&ntuData.lastkaonKinEKLTBoost,"lastkaonKinEKLTBoost/D");

	      ROOTntup[i]->Branch("TotalEnergyDepKLHPGeBoost",&ntuData.TotalEnergyDepKLHPGeBoost,"TotalEnergyDepKLHPGeBoost/D");
	      ROOTntup[i]->Branch("gammaKinEKLHPGeBoost",&ntuData.gammaKinEKLHPGeBoost,"gammaKinEKLHPGeBoost/D");
	      ROOTntup[i]->Branch("nHitHPGeBoost", &ntuData.nHitHPGeBoost, "nHitHPGeBoost/I");
	      ROOTntup[i]->Branch("TimeKLHPGeBoost",ntuData.TimeKLHPGeBoost,"TimeKLHPGeBoost[nHitHPGeBoost]/D");
	      ROOTntup[i]->Branch("EnergyDepKLHPGeBoost",ntuData.EnergyDepKLHPGeBoost,"EnergyDepKLHPGeBoost[nHitHPGeBoost]/D");
	      ROOTntup[i]->Branch("XYZKLHPGeBoost",ntuData.XYZKLHPGeBoost,"XYZKLHPGeBoost[3][nHitHPGeBoost]/D");
	      ROOTntup[i]->Branch("pdgcodeKLHPGeBoost",ntuData.pdgcodeKLHPGeBoost,"pdgcodeKLHPGeBoost[nHitHPGeBoost]/I");
	      ROOTntup[i]->Branch("kaonKinEKLHPGeBoost",ntuData.kaonKinEKLHPGeBoost,"kaonKinEKLHPGeBoost[nHitHPGeBoost]/D");
      }
      G4cout << "ROOT Ntuple " << id << " " << tupleName[i] << " "<< tupleListROOT[i] << " booked " << G4endl;
    }
  }
}

void SiddhartaHisto::save()
{
  G4cout << "ROOT: files writing..." << G4endl;
  hfileROOT->Write();
  hfileROOT->cd();
  fStats.Write();
  G4cout << "ROOT: files closing..." << G4endl;
  hfileROOT->Close();
  delete hfileROOT;
}

void SiddhartaHisto::add1D(const G4String& id, const G4String& name, G4int nb, G4double x1, G4double x2, G4double u)
{
  if (verbose > 0) {
    G4cout << "New histogram will be booked: #" << id << "  <" << name << "  " << nb << "  " << x1 << "  " << x2 << "  " << u << G4endl;
  }
  nHisto++;
  x1 /= u;
  x2 /= u;
  active.push_back(defaultAct);
  bins.push_back(nb);
  xmin.push_back(x1);
  xmax.push_back(x2);
  unit.push_back(u);
  ids.push_back(id);
  titles.push_back(name);
  ROOThisto.push_back(0);
  hashhisto[id] = nHisto;
}

void SiddhartaHisto::add3D(const G4String& id, const G4String& name, G4int nb1, G4double x11, G4double x21, G4double u1,
                           G4int nb2, G4double x12, G4double x22, G4double u2, G4int nb3, G4double x13, G4double x23, G4double u3)
{
  if (verbose > 0) {
    G4cout << "New histogram will be booked: #" << id << "  <" << name << "  " << nb1 << "  " << x11 << "  " << x21 << "  " << u1 << G4endl;
  }
  nHisto3++;
  x11 /= u1;
  x21 /= u1;
  active3.push_back(defaultAct);
  bins31.push_back(nb1);
  xmin31.push_back(x11);
  xmax31.push_back(x21);
  unit31.push_back(u1);
  bins32.push_back(nb2);
  xmin32.push_back(x12);
  xmax32.push_back(x22);
  unit32.push_back(u2);
  bins33.push_back(nb3);
  xmin33.push_back(x13);
  xmax33.push_back(x23);
  unit33.push_back(u3);
  ids3.push_back(id);
  titles3.push_back(name);
  ROOThisto3.push_back(0);
  hashhisto3[id] = nHisto3;
}

void SiddhartaHisto::setHisto1D(G4int i, G4int nb, G4double x1, G4double x2, G4double u)
{
  if (i>=0 && i<nHisto) {
    if(verbose > 0) {
      G4cout << "Update histogram: #" << i << "  " << nb << "  " << x1 << "  " << x2 << "  " << u << G4endl;
    }
    bins[i] = nb;
    xmin[i] = x1;
    xmax[i] = x2;
    unit[i] = u;
  } else {
    G4cout << "SiddhartaHisto::setSiddhartaHisto1D: WARNING! wrong histogram index " << i << G4endl;
  }
}

void SiddhartaHisto::fillHisto(G4int i, G4double x, G4double w)
{
  if(verbose > 1) {
    G4cout << "fill histogram: #" << i << " at x= " << x << "  weight= " << w << G4endl;
  }
  if(i>=0 && i<nHisto) {
    ROOThisto[i]->Fill(x/unit[i], w);
  } else {
    G4cout << "SiddhartaHisto::fill: WARNING! wrong ROOT histogram index " << i << G4endl;
  }
}

void SiddhartaHisto::fillHisto(G4String s, G4double x, G4double w)
{
  if( hashhisto[s] == 0) {
    G4cout << "Trying to fill Non-existing histogram h"<<hashhisto[s]<<G4endl;
  }
  G4int i = hashhisto[s]-1;
  SiddhartaHisto::fillHisto(i, x, w);
}

void SiddhartaHisto::fillHisto3(G4int i, G4double x, G4double y, G4double z, G4double w)
{
  if(verbose > 1) {
    G4cout << "fill histogram: #" << i << " at x= " << x << "  weight= " << w << G4endl;
  }
  if(i>=0 && i<nHisto3) {
    ROOThisto3[i]->Fill(x/unit31[i], y/unit32[i], z/unit33[i], w);
  } else {
    G4cout << "SiddhartaHisto3::fill: WARNING! wrong ROOT histogram index " << i << G4endl;
  }
}

void SiddhartaHisto::fillHisto3(G4String s, G4double x, G4double y, G4double z, G4double w)
{
  if( hashhisto3[s] == 0) {
    G4cout << "Trying to fill Non-existing histogram h3d" << hashhisto3[s] << G4endl;
  }
  G4int i = hashhisto3[s] - 1;
  SiddhartaHisto::fillHisto3(i, x, y, z, w);
}

void SiddhartaHisto::scaleHisto(G4int i, G4double x)
{
  if(verbose > 0) {
    G4cout << "Scale histogram: #" << i << " by factor " << x << G4endl;
  }
  if(i>=0 && i<nHisto) {
    ROOThisto[i]->Scale(x);
  } else {
    G4cout << "SiddhartaHisto::scale: WARNING! wrong ROOT histogram index " << i << G4endl;
  }
}

void SiddhartaHisto::addTuple(const G4String& w1, const G4String& w2, const G4String& w3)
{
  G4cout << w1 << " " << w2 << " " << w3 << G4endl;
  nTuple++;
  tupleId.push_back(w1);
  tupleName.push_back(w2) ;
  tupleListROOT.push_back(w3);
  ROOTntup.push_back(0);
  hashNtuple[w1] = nTuple;
}

void SiddhartaHisto::fillTuple(const G4String& id)
{
  ROOTntup[hashNtuple[id]-1]->Fill();
}

void SiddhartaHisto::fillTuple(G4int i, const G4String& parname, G4double x)
{
  G4cout << "fill tuple # " << i <<" with  parameter <" << parname << "> = " << x << G4endl;
}

void SiddhartaHisto::fillTuple(G4int i, G4int col, G4double x)
{
  if(verbose > 1) {
    G4cout << "fill tuple # " << i <<" in column < " << col << "> = " << x << G4endl;
  }
  if(ROOTntup[i])
    (Rarray[i])[col] = float(x);
}

void SiddhartaHisto::fillTuple(G4int i, const G4String& parname, G4String& x)
{
  if(verbose > 1) {
    G4cout << "fill tuple # " << i <<" with  parameter <" << parname << "> = " << x << G4endl;
  }
}

void SiddhartaHisto::addRow(G4int i)
{
  if(verbose > 1)
    G4cout << "Added a raw #" << i << " to tuple" << G4endl;
  float ar[4];
  for (G4int j=0; j < Rcol[i]; j++) {
    ar[j] = Rarray[i][j];
  }
}

void SiddhartaHisto::setFileName(const G4String& nam)
{
  histName = nam;
}

void SiddhartaHisto::setFileType(const G4String& nam)
{
  histType = nam;
}

const G4String& SiddhartaHisto::FileType() const
{
  return histType;
}

void SiddhartaHisto::createHistogramWithAxes(TObject* object, TString xAxisName, TString yAxisName, TString zAxisName)
{
  TClass *cl = object->IsA();
  if (cl->InheritsFrom("TH1D") || cl->InheritsFrom("TH1F")) {
    TH1D* tempHisto = new TH1D();
    if (cl->InheritsFrom("TH1F")) {
      TH1F* floatTemp = dynamic_cast<TH1F*>(object);
      floatTemp->Copy(*tempHisto);
      tempHisto->SetDirectory(nullptr);
      delete floatTemp;
    } else
      tempHisto = dynamic_cast<TH1D*>(object);
    tempHisto->GetXaxis()->SetTitle(xAxisName);
    tempHisto->GetYaxis()->SetTitle(yAxisName);
    fStats.Add(tempHisto);
  } else if (cl->InheritsFrom("TH2D") || cl->InheritsFrom("TH2F")) {
    TH2D* tempHisto = new TH2D();
    if (cl->InheritsFrom("TH2F")) {
      TH2F* floatTemp = dynamic_cast<TH2F*>(object);
      floatTemp->Copy(*tempHisto);
      tempHisto->SetDirectory(nullptr);
      delete floatTemp;
    } else
      tempHisto = dynamic_cast<TH2D*>(object);
    tempHisto->GetXaxis()->SetTitle(xAxisName);
    tempHisto->GetYaxis()->SetTitle(yAxisName);
    fStats.Add(tempHisto);
  } else if (cl->InheritsFrom("TH3D") || cl->InheritsFrom("TH3F")) {
    TH3D* tempHisto = new TH3D();
    if (cl->InheritsFrom("TH3F")) {
      TH3F* floatTemp = dynamic_cast<TH3F*>(object);
      floatTemp->Copy(*tempHisto);
      tempHisto->SetDirectory(nullptr);
      delete floatTemp;
    } else
      tempHisto = dynamic_cast<TH3D*>(object);
    tempHisto->GetXaxis()->SetTitle(xAxisName);
    tempHisto->GetYaxis()->SetTitle(yAxisName);
    tempHisto->GetZaxis()->SetTitle(zAxisName);
    fStats.Add(tempHisto);
  }
}

void SiddhartaHisto::fillHistogram(const char* name, double xValue, doubleCheck yValue, doubleCheck zValue)
{
  fillHistogramWithWeight(name, 1., xValue, yValue, zValue);
}

void SiddhartaHisto::fillHistogramWithWeight(const char* name, double weight, double xValue, doubleCheck yValue, doubleCheck zValue)
{
  TObject *tempObject = getObject<TObject>(name);
  if (!tempObject) {
    return;
  }
  TClass *cl = tempObject->IsA();
  if (cl->InheritsFrom("TH1D")) {
    TH1D* tempHisto = dynamic_cast<TH1D*>(tempObject);
    tempHisto->Fill(xValue, weight);
  } else if (cl->InheritsFrom("TH2D")) {
    TH2D* tempHisto = dynamic_cast<TH2D*>(tempObject);
    if (yValue.isChanged)
      tempHisto->Fill(xValue, yValue.value, weight);
  } else if (cl->InheritsFrom("TH3D")) {
    TH3D* tempHisto = dynamic_cast<TH3D*>(tempObject);
    if (zValue.isChanged)
      tempHisto->Fill(xValue, yValue.value, zValue.value, weight);
  }
}

int SiddhartaHisto::getMatIDFromString(std::string matName)
{
  if (matID.find(matName) == matID.end()) {
    matID.insert({matName, matID.size()});
    TH1D* tempHisto = dynamic_cast<TH1D*>(getObject<TObject>("materialIDmap"));
    tempHisto->Fill(matID.size());
    tempHisto->GetXaxis()->SetBinLabel(matID.size(), (matName).c_str());
    return matID.size();
  } else {
    TH1D* tempHisto = dynamic_cast<TH1D*>(getObject<TObject>("materialIDmap"));
    tempHisto->Fill(matID.at(matName));
    return matID.at(matName);
  }
}

void SiddhartaHisto::addMatIDFromString(std::string matName)
{
  if (matID.find(matName) == matID.end()) {
    matID.insert({matName, matID.size()});
    TH1D* tempHisto = dynamic_cast<TH1D*>(getObject<TObject>("materialIDmap"));
    tempHisto->GetXaxis()->SetBinLabel(matID.size(), (matName).c_str());
  }
}

