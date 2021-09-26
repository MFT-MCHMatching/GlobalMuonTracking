

#include "FairGenerator.h"

class MatcherBoxGen : public FairGenerator
{

 public:
  MatcherBoxGen(int nparticles, int pdgcode, float etamin, float etamax,
                float ptmin, float ptmax) : FairGenerator(),
                                            mPDGCode(pdgcode),
                                            mNParticles(nparticles),
                                            mEtaMin(etamin),
                                            mEtaMax(etamax),
                                            mPtMin(ptmin),
                                            mPtMax(ptmax){};
  ~MatcherBoxGen() = default;

  int mPDGCode;
  int mNParticles;
  float mEtaMin;
  float mEtaMax;
  float mPtMin;
  float mPtMax;
  bool mRandomizeCharge = true;
  void disableRandomCharge() { mRandomizeCharge = false; }

  Bool_t ReadEvent(FairPrimaryGenerator* primGen) override
  {

    int iPart = mNParticles;
    while (iPart) {
      float pt = gRandom->Uniform(mPtMin, mPtMax);
      float eta = gRandom->Uniform(mEtaMin, mEtaMax);
      float phi = gRandom->Uniform(0., TMath::TwoPi());
      float px = pt * TMath::Cos(phi);
      float py = pt * TMath::Sin(phi);
      float tanl = tan(TMath::Pi() / 2 - 2 * atan(exp(-eta)));
      float pz = tanl * pt;
      int charge = 1;

      if (mRandomizeCharge && (gRandom->Rndm() < 0.5)) {
        charge = -1;
      }
      primGen->AddTrack(charge * mPDGCode, px, py, pz, 0, 0, 0);
      printf("Add track %d %.2f %.2f %.2f  \n", charge * mPDGCode, px, py, pz);

      iPart--;
    }
    return kTRUE;
  }

 private:
};

FairGenerator* matcherMuBoxGen(int nParticles = 1, int pdgCode = 13, float etamin = -3.6f, float etamax = -2.4f, float ptmin = 0.01f, float ptmax = 20.f)
{

  auto gen = new MatcherBoxGen(nParticles, pdgCode, etamin, etamax, ptmin, ptmax);
  return gen;
}
