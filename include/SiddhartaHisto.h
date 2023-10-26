#ifndef SiddhartaHisto_h
#define SiddhartaHisto_h 1

#include <THashTable.h>

#include <globals.hh>

#include <vector>
#include <string>
#include <map>

//???
class TFile;
class TH1D;
class TH2D;
class TH3D;
class TTree;

class doubleCheck
{
public:
  bool isChanged = false;
  double value = 0.;
  doubleCheck(){};
  doubleCheck(double newValue) {value=newValue; isChanged=true;};
};

const unsigned int MaxHits = 700;
const unsigned int MaxCZTHits = 50;
const unsigned int MaxAnti = 40;

struct siddharta_t {
  double XYZvertex[3];
  double XYZstopK[3];   //18.20.2021
  double XYZstopKP[3];   //22.20.2021
  double VertexMomentum[3];
  int nHitSDD;
  int nHitCZTAntiBoost;
  int nHitHPGeBoost;
  int NoSDD[MaxHits];
  int NbHitSDD[MaxHits];
  double EnergySDD[MaxHits];
  double TimeSDD[MaxHits];
  double XYZvertexSDD[MaxHits][3];
  double XYZSDD[MaxHits][3];
  double XYZInSDD[MaxHits][3];
  double KinvertexSDD[MaxHits];
  char particleNameVertex[MaxHits][20];
  int particleVertexPDGEncoding[MaxHits];
  unsigned long int matvert[MaxHits];
  unsigned matVertIDFromMap[MaxHits];
  char materialVertex[MaxHits][20];
  int parentID[MaxHits];
  char parentName[MaxHits][20];
  int parentPDGEncoding[MaxHits];
  double TimeKMTop;
  double EnergyDepKMTop;
  double XYZKMTop[3];
  char particleNameKMTop[20];
  double TimeKMBottom;
  double EnergyDepKMBottom;
  double XYZKMBottom[3];
  char particleNameKMBottom[20];
  int particleKMBottomPDGEncoding;

  double TimeKPlusDetector;
  double EnergyDepKPlusDetector;
  double XYZKPlusDetector[3];
  char particleNameKPlusDetector[20];
  int particleKPlusDetectorPDGEncoding;

  double TimeLMBoost;
  double EnergyDepLMBoost;
  double XYZLMBoost[3];
  double XYZLMBoostKaonstop[3];
  char particleNameLMBoost[20];
  int pdgcodeLMBoost;
  double kaonKinELMBoost;
  double lastkaonKinELMBoost;
  double TimeLMAntiBoost;
  double EnergyDepLMAntiBoost;
  double XYZLMAntiBoost[3];
  double XYZLMAntiBoostKaonstop[3];
  char particleNameLMAntiBoost[20];
  int pdgcodeLMAntiBoost;
  double kaonKinELMAntiBoost;
  double lastkaonKinELMAntiBoost;

  double TimeKLDAntiBoost;
  double EnergyDepKLDAntiBoost;
  double kaonKinEKLDAntiBoost;
  double lastkaonKinEKLDAntiBoost;
  double XYZKLDAntiBoost[3];
  double KLDAntiBooststop[3];
  int pdgcodeKLDAntiBoost;

  double TimeKLT1AntiBoost;
  double EnergyDepKLT1AntiBoost;
  double kaonKinEKLT1AntiBoost;
  double gammaKinEKLT1AntiBoost;
  double lastkaonKinEKLT1AntiBoost;
  double XYZKLT1AntiBoost[3];
  double KLT1AntiBooststop[3];
  double KLT1AntiBoostKaonstop[3];
  int pdgcodeKLT1AntiBoost;

  double TimeKLT2AntiBoost;
  double EnergyDepKLT2AntiBoost;
  double kaonKinEKLT2AntiBoost;
  double gammaKinEKLT2AntiBoost;
  double lastkaonKinEKLT2AntiBoost;
  double XYZKLT2AntiBoost[3];
  double KLT2AntiBooststop[3];
  double KLT2AntiBoostKaonstop[3];
  int pdgcodeKLT2AntiBoost;

  double TimeKLBWAntiBoost;
  double EnergyDepKLBWAntiBoost;
  double kaonKinEKLBWAntiBoost;
  double lastkaonKinEKLBWAntiBoost;
  double XYZKLBWAntiBoost[3];
  double KLBWAntiBooststop[3];
  int pdgcodeKLBWAntiBoost;

  double TotalEnergyDepKLCZTAntiBoost;
  double gammaKinEKLCZTAntiBoost;
  double TimeKLCZTAntiBoost[MaxCZTHits];
  double EnergyDepKLCZTAntiBoost[MaxCZTHits];
  double kaonKinEKLCZTAntiBoost[MaxCZTHits];
  double XYZKLCZTAntiBoost[3][MaxCZTHits];
  int pdgcodeKLCZTAntiBoost[MaxCZTHits];

  double TimeKLDBoost;
  double EnergyDepKLDBoost;
  double kaonKinEKLDBoost;
  double lastkaonKinEKLDBoost;
  double XYZKLDBoost[3];
  double KLDBooststop[3];
  int pdgcodeKLDBoost;

  double TimeKLTBoost;
  double EnergyDepKLTBoost;
  double kaonKinEKLTBoost;
  double gammaKinEKLTBoost;
  double lastkaonKinEKLTBoost;
  double XYZKLTBoost[3];
  double KLTBooststop[3];
  double KLTBoostKaonstop[3];
  int pdgcodeKLTBoost;

  double TotalEnergyDepKLHPGeBoost;
  double gammaKinEKLHPGeBoost;
  double TimeKLHPGeBoost[MaxCZTHits];
  double EnergyDepKLHPGeBoost[MaxCZTHits];
  double kaonKinEKLHPGeBoost[MaxCZTHits];
  double XYZKLHPGeBoost[3][MaxCZTHits];
  int pdgcodeKLHPGeBoost[MaxCZTHits];

  double TimeGhost;
  double XYZGhost[3];
  double EnergyDepGhost;
  double KineticGhost;

  double TimeAnti1;
  double EnergyAnti1;
  double XYZAnti1[3];
  double XYZVertexAnti1[3];
  char particleNameAnti1[20];
  int parentIDAnti1;
  char parentNameAnti1[20];
  double TimeAnti2;
  double EnergyAnti2;
  double XYZAnti2[3];
  double XYZVertexAnti2[3];
  char particleNameAnti2[20];
  int parentIDAnti2;
  char parentNameAnti2[20];
  double TimeAnti3;
  double EnergyAnti3;
  double XYZAnti3[3];
  double XYZVertexAnti3[3];
  char particleNameAnti3[20];
  int parentIDAnti3;
  char parentNameAnti3[20];
  unsigned matVertIDFromMap_Veto1[20];
  double TimeAnti4;
  double EnergyAnti4;
  double XYZAnti4[3];
  double XYZVertexAnti4[3];
  char particleNameAnti4[20];
  int parentIDAnti4;
  char parentNameAnti4[20];
  double TimeAnti5;
  double EnergyAnti5;
  double XYZAnti5[3];
  double XYZVertexAnti5[3];
  char particleNameAnti5[20];
  int parentIDAnti5;
  char parentNameAnti5[20];

  int NbAnti;
  double TimeAnti[MaxAnti];
  double EnergyAnti[MaxAnti];
  double XYZAnti[MaxAnti][3];
  double XYZVertexAnti[MaxAnti][3];
  double KinVertexAnti[MaxAnti];
  char particleNameAnti[MaxAnti][20];
  int particleAntiPDGEncoding[MaxAnti];
  unsigned long int matvertAnti[MaxAnti];
  char materialVertexAnti[MaxAnti][20];
  int parentIDAnti[MaxAnti];
  char parentNameAnti[MaxAnti][20];
  char parentAntiPDGEncoding[MaxAnti];

  int nHitSciAnti;
  double TimeScintAnti[700];
  double EnergyScintAnti[700];
  double XYZScintAnti[700][3];
  int particleScintAntiPDGEncoding[700];
  unsigned matVertIDFromMap_Veto2[700];
  int copyScintAnti[700];
};

class SiddhartaHisto
{
public:
  SiddhartaHisto();
  ~SiddhartaHisto();

  void book();
  void save();
  void add1D(const G4String&, const G4String&, G4int nb=100, G4double x1=0.,
                                               G4double x2=1., G4double u=1.);
  void add3D(const G4String&, const G4String&, G4int nb1=100, G4double x11=0.,G4double x21=1., G4double u1=1.,
		                               G4int nb2=100, G4double x12=0.,G4double x22=1., G4double u2=1.,
		                               G4int nb3=100, G4double x13=0.,G4double x23=1., G4double u3=1.);

  void createHistogramWithAxes(TObject* object, TString xAxisName = "X Axis", TString yAxisName = "Y Axis", TString zAxisName = "Z Axis");
  void fillHistogram(const char* name, double xValue, doubleCheck yValue, doubleCheck zValue);
  void fillHistogramWithWeight(const char* name, double weight, double xValue, doubleCheck yValue, doubleCheck zValue);

  void setHisto1D(G4int, G4int, G4double, G4double, G4double);

  void fillHisto(G4int, G4double, G4double);
  void fillHisto(G4String, G4double, G4double);
  void fillHisto3(G4int, G4double, G4double, G4double, G4double);
  void fillHisto3(G4String, G4double, G4double, G4double, G4double);

  void scaleHisto(G4int, G4double);
  void addTuple(const G4String&, const G4String&, const G4String&);
  void fillTuple(G4int, const G4String&, G4double);
  void fillTuple(G4int, G4int, G4double);
  void fillTuple(G4int, const G4String&, G4String&);
  void fillTuple(G4int, const G4String&, G4bool);
  void fillTuple(const G4String&);

  void addRow(G4int);
  void setFileName(const G4String&);
  void setFileType(const G4String&);
  const G4String& FileType() const;

  std::map<G4String, G4int> hashhisto;
  std::map<G4String, G4int> hashhisto3;
  std::map<G4String, G4int> hashNtuple;
  std::map<G4String, G4double> hashunit;
  siddharta_t ntuData;

  template <typename T>
  T* getObject(const char* name)
  {
    TObject* tmp = fStats.FindObject(name);
    if (!tmp) {
      std::cerr << "getObject of " + std::string(name) + " returned nullptr" << std::endl;;
      return nullptr;
    }
    return dynamic_cast<T*>(tmp);
  }
  const THashTable* getStatsTable() const;

  void addMatIDFromString(std::string matName);
  int getMatIDFromString(std::string matName);

private:
  G4String histName;
  G4String histType;

  G4int nHisto;
  G4int nHisto3;
  G4int nTuple;
  G4int verbose;
  G4int defaultAct;

  TFile* hfileROOT;
  std::vector<TH1D*> ROOThisto;
  std::vector<TH3D*> ROOThisto3;
  std::vector<TTree*> ROOTntup;
  std::vector<std::vector<float>> Rarray;
  std::vector<G4int> Rcol;

  std::vector<G4int> active;
  std::vector<G4int> bins;
  std::vector<G4double> xmin;
  std::vector<G4double> xmax;
  std::vector<G4double> unit;
  std::vector<G4String> ids;
  std::vector<G4String> titles;
  std::vector<G4String> tupleName;
  std::vector<G4String> tupleId;
  std::vector<G4String> tupleList;
  std::vector<G4String> tupleListROOT;

  std::vector<G4int> active3;
  std::vector<G4int> bins31;
  std::vector<G4double> xmin31;
  std::vector<G4double> xmax31;
  std::vector<G4double> unit31;
  std::vector<G4int> bins32;
  std::vector<G4double> xmin32;
  std::vector<G4double> xmax32;
  std::vector<G4double> unit32;
  std::vector<G4int> bins33;
  std::vector<G4double> xmin33;
  std::vector<G4double> xmax33;
  std::vector<G4double> unit33;
  std::vector<G4String> ids3;
  std::vector<G4String> titles3;

  THashTable fStats;
  std::map<std::string, int> matID;
};

#endif
