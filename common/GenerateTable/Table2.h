#ifndef HYPERTABLE2_H
#define HYPERTABLE2_H

#include "Common.h"

#include <string>

#include <TLorentzVector.h>
#include <TTree.h>
#include <TF1.h>
#include <TFile.h>
#include <TVector3.h>

#include "AliAnalysisTaskHyperTriton2He3piML.h"

class Table2
{
public:
  Table2(std::string name, std::string title);
  void Fill(const RHyperTritonHe3pi &RHyperVec, const RCollision &RColl, bool fLambda);
  void Write() { tree->Write(); }

private:
  TTree *tree;
  TF1 *fHe3TPCcalib;
  float pt;
  float TPCnSigmaHe3;
  float ct;
  float m;
  float ArmenterosAlpha;
  float V0CosPA;
  float V0Chi2;
  float PiProngPt;
  float He3ProngPt;
  float ProngsDCA;
  float PiProngPvDCA;
  float He3ProngPvDCA;
  float NpidClustersHe3;
  float NitsClustersHe3;
  float NpidClustersPion;
  float TPCnSigmaPi;
  float Lrec;
  float centrality;
  float V0radius;
  float PiProngPvDCAXY;
  float He3ProngPvDCAXY;
  float Rapidity;
  float PseudoRapidityHe3;
  float PseudoRapidityPion;
  float Matter;
  float TOFnSigmaHe3;
  float TOFnSigmaPi;
  float decradius;
  float HasITSLayerHe[6];
  float HasITSLayerPi[6];
};

Table2::Table2(std::string name, std::string title)
{
  tree = new TTree(name.data(), title.data());

  string hypUtilsDir = getenv("HYPERML_UTILS");
  string calibFileArg = hypUtilsDir + "/He3TPCCalibration.root";

  TFile calibFile(calibFileArg.data(), "READ");
  fHe3TPCcalib = dynamic_cast<TF1 *>(calibFile.Get("He3TPCCalib")->Clone());
  calibFile.Close();

  tree->Branch("pt", &pt);
  tree->Branch("TPCnSigmaHe3", &TPCnSigmaHe3);
  tree->Branch("ct", &ct);
  tree->Branch("m", &m);
  tree->Branch("ArmenterosAlpha", &ArmenterosAlpha);
  tree->Branch("V0CosPA", &V0CosPA);
  tree->Branch("V0Chi2", &V0Chi2);
  tree->Branch("PiProngPt", &PiProngPt);
  tree->Branch("He3ProngPt", &He3ProngPt);
  tree->Branch("ProngsDCA", &ProngsDCA);
  tree->Branch("He3ProngPvDCA", &He3ProngPvDCA);
  tree->Branch("PiProngPvDCA", &PiProngPvDCA);
  tree->Branch("He3ProngPvDCAXY", &He3ProngPvDCAXY);
  tree->Branch("PiProngPvDCAXY", &PiProngPvDCAXY);
  tree->Branch("NpidClustersHe3", &NpidClustersHe3);
  tree->Branch("NpidClustersPion", &NpidClustersPion);
  tree->Branch("NitsClustersHe3", &NitsClustersHe3);  
  tree->Branch("TPCnSigmaPi", &TPCnSigmaPi);
  tree->Branch("Lrec", &Lrec);
  tree->Branch("centrality", &centrality);
  tree->Branch("V0radius", &V0radius);
  tree->Branch("Rapidity", &Rapidity);
  tree->Branch("PseudoRapidityHe3", &PseudoRapidityHe3);
  tree->Branch("PseudoRapidityPion", &PseudoRapidityPion);
  tree->Branch("Matter", &Matter);
  tree->Branch("TOFnSigmaHe3",&TOFnSigmaHe3);
  tree->Branch("TOFnSigmaPi",&TOFnSigmaPi);
  tree->Branch("ITS1He3",&HasITSLayerHe[0]);
  tree->Branch("ITS2He3",&HasITSLayerHe[1]);
  tree->Branch("ITS3He3",&HasITSLayerHe[2]);
  tree->Branch("ITS4He3",&HasITSLayerHe[3]);
  tree->Branch("ITS5He3",&HasITSLayerHe[4]);
  tree->Branch("ITS6He3",&HasITSLayerHe[5]);
  tree->Branch("ITS1pi",&HasITSLayerPi[0]);
  tree->Branch("ITS2pi",&HasITSLayerPi[1]);
  tree->Branch("ITS3pi",&HasITSLayerPi[2]);
  tree->Branch("ITS4pi",&HasITSLayerPi[3]);
  tree->Branch("ITS5pi",&HasITSLayerPi[4]);
  tree->Branch("ITS6pi",&HasITSLayerPi[5]);
};

void Table2::Fill(const RHyperTritonHe3pi &RHyper, const RCollision &RColl, bool fLambda = false)
{

  float mother_mass,fat_daughter_mass;
  if(!fLambda)
    fat_daughter_mass = AliPID::ParticleMass(AliPID::kHe3);
  else
    fat_daughter_mass = AliPID::ParticleMass(AliPID::kProton);

  centrality = RColl.fCent;
  double eHe3 = hypot(RHyper.fPxHe3, RHyper.fPyHe3, RHyper.fPzHe3, fat_daughter_mass);
  double ePi = hypot(RHyper.fPxPi, RHyper.fPyPi, RHyper.fPzPi, AliPID::ParticleMass(AliPID::kPion));

  TLorentzVector he3Vector, piVector, hyperVector;
  he3Vector.SetPxPyPzE(RHyper.fPxHe3, RHyper.fPyHe3, RHyper.fPzHe3, eHe3);
  piVector.SetPxPyPzE(RHyper.fPxPi, RHyper.fPyPi, RHyper.fPzPi, ePi);
  hyperVector = piVector + he3Vector;

  TVector3 v(RHyper.fDecayX, RHyper.fDecayY, RHyper.fDecayZ);
  float pointAngle = hyperVector.Angle(v);
  float CosPA = std::cos(pointAngle);

  float qP, qN, qT;
  if (RHyper.fMatter == true)
  {
    qP = SProd(hyperVector, he3Vector) / fabs(hyperVector.P());
    qN = SProd(hyperVector, piVector) / fabs(hyperVector.P());
    qT = VProd(hyperVector, he3Vector) / fabs(hyperVector.P());
  }
  else
  {
    qN = SProd(hyperVector, he3Vector) / fabs(hyperVector.P());
    qP = SProd(hyperVector, piVector) / fabs(hyperVector.P());
    qT = VProd(hyperVector, piVector) / fabs(hyperVector.P());
  }

  float alpha = (qP - qN) / (qP + qN);
  ct = kHypertritonMass * (hypot(RHyper.fDecayX, RHyper.fDecayY, RHyper.fDecayZ) / hyperVector.P());
  m = hyperVector.M();
  ArmenterosAlpha = alpha;
  V0CosPA = CosPA;
  V0Chi2 = RHyper.fChi2V0;
  PiProngPt = hypot(RHyper.fPxPi, RHyper.fPyPi);
  He3ProngPt = hypot(RHyper.fPxHe3, RHyper.fPyHe3);
  ProngsDCA = RHyper.fDcaV0daughters;
  PiProngPvDCA = RHyper.fDcaPi2PrimaryVertex;
  He3ProngPvDCA = RHyper.fDcaHe32PrimaryVertex;
  PiProngPvDCAXY = RHyper.fDcaPi2PrimaryVertexXY;
  He3ProngPvDCAXY = RHyper.fDcaHe32PrimaryVertexXY;
  Lrec = hypot(RHyper.fDecayX, RHyper.fDecayY, RHyper.fDecayZ);
  V0radius = hypot(RHyper.fDecayX, RHyper.fDecayY);
  NpidClustersHe3 = RHyper.fNpidClustersHe3;
  NitsClustersHe3 = RHyper.fITSclusHe3;
  NpidClustersPion = RHyper.fNpidClustersPi;
  TPCnSigmaPi = RHyper.fTPCnSigmaPi;
  TPCnSigmaHe3 = RHyper.fTPCnSigmaHe3;
  TOFnSigmaHe3 = RHyper.fTOFnSigmaHe3;
  TOFnSigmaPi = RHyper.fTOFnSigmaPi;
  pt = hyperVector.Pt();
  Rapidity = hyperVector.Rapidity();
  Matter = RHyper.fMatter;
  PseudoRapidityHe3 = he3Vector.PseudoRapidity();
  PseudoRapidityPion = piVector.PseudoRapidity();
  HasITSLayerHe[0] = RHyper.fITSclusHe3 & 1; // SPD interno
  HasITSLayerHe[1] = RHyper.fITSclusHe3 & (1 << 1); // SPD esterno
  HasITSLayerHe[2] = RHyper.fITSclusHe3 & (1 << 2); // SDD interno
  HasITSLayerHe[3] = RHyper.fITSclusHe3 & (1 << 3); // SPD interno
  HasITSLayerHe[4] = RHyper.fITSclusHe3 & (1 << 4); // SPD esterno
  HasITSLayerHe[5] = RHyper.fITSclusHe3 & (1 << 5); // SDD interno
  HasITSLayerPi[0] = RHyper.fITSclusPi & 1; // SPD interno
  HasITSLayerPi[1] = RHyper.fITSclusPi & (1 << 1); // SPD esterno
  HasITSLayerPi[2] = RHyper.fITSclusPi & (1 << 2); // SDD interno
  HasITSLayerPi[3] = RHyper.fITSclusPi & (1 << 3); // SPD interno
  HasITSLayerPi[4] = RHyper.fITSclusPi & (1 << 4); // SPD esterno
  HasITSLayerPi[5] = RHyper.fITSclusPi & (1 << 5); // SDD interno
  float ITSradii[6] = {3.9,7.6,15.0,23.9,38.0,43.0};//cm

  bool ITSreject = false;
  for(int itslayer=0; itslayer<6; itslayer++){
    if((HasITSLayerHe[itslayer] || HasITSLayerPi[itslayer]) && V0radius>ITSradii[itslayer]){
      ITSreject = true;
      break;
    }
  }

  if (He3ProngPt > 1.2 && ProngsDCA < 1.6 && NpidClustersHe3>30 && ITSreject)
    tree->Fill();
  else
  {
  }
}

#endif
