#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMOffline/Trigger/plugins/TriggerDQMBase.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
//---------------------------------------------------------
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TLorentzVector.h"
#include <cassert>

class ZplusJetsMonitor : public DQMEDAnalyzer, public TriggerDQMBase {
public:
  typedef dqm::reco::MonitorElement MonitorElement;
  typedef dqm::reco::DQMStore DQMStore;

  ZplusJetsMonitor(const edm::ParameterSet&);
  ~ZplusJetsMonitor() throw() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void dqmBeginRun(edm::Run const & iRun, edm::EventSetup const& iSetup) override;

  bool islowZPt(double ZPt);
  bool ismediumZPt(double ZPt);
  bool ishighZPt(double ZPt);

  bool isBarrel(double eta);
  bool isEndCap(double eta);
  bool isForward(double eta);

  void bookMESub(DQMStore::IBooker&,
                 ObjME* a_me,
                 const int len_,
                 const std::string& h_Name,
                 const std::string& h_Title,
                 const std::string& h_subOptName,
                 const std::string& h_subOptTitle,
                 const bool doDirectBalancevsZPt = true);

  void FillME(ObjME* a_me,
              const double ZJetB_,
              const double JetminusZ_,
	      const double Zpt_,
 	      const double Jetpt_,
              const bool doDirectBalancevsZPt = true);  

private:
  const std::string folderName_;

  const std::string processName_; // process name of (HLT) process for which to get HLT configuration
  /// The instance of the HLTConfigProvider as a data member
  HLTConfigProvider hltConfig_;

  edm::EDGetTokenT<trigger::TriggerEvent> triggerEventObject_;
  //edm::EDGetTokenT<trigger::TriggerEventWithRefs> triggerEventObject_; //for RAW?

  double ptcut_;
  bool isPFJetTrig;
  bool isCaloJetTrig;

  std::vector<double> directbalanceBinning_;
  std::vector<double> ZptBinning_;
  std::vector<double> jetptBinning_;
  //MEbinning ZJetB_binning_{50, 0.,3.0};
  MEbinning JoZ_binning_{50, -75.0, 75.0};

  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
 
  ObjME a_ME[5];
  ObjME a_ME_HB[5];
  ObjME a_ME_HE[5];
  ObjME a_ME_HF[5];
  ObjME a_ME_HB_lowZPt[5];
  ObjME a_ME_HE_lowZPt[5];
  ObjME a_ME_HF_lowZPt[5];
  ObjME a_ME_HB_mediumZPt[5];
  ObjME a_ME_HE_mediumZPt[5];
  ObjME a_ME_HF_mediumZPt[5];
  ObjME a_ME_HB_highZPt[5];
  ObjME a_ME_HE_highZPt[5];
  ObjME a_ME_HF_highZPt[5];
 
  ObjME mZMassME_;
  ObjME mDPhiZJetME_;
  ObjME mZJetAsymmetryME_;

  std::vector<double> v_jetpt;
  std::vector<double> v_jeteta;
  std::vector<double> v_jetphi;
  int numtrigobj;
  std::vector<float> trigobj_pt;
  std::vector<float> trigobj_eta;
  std::vector<float> trigobj_phi;
  double Met_px;
  double Met_py;
  TLorentzVector muon_1;
  TLorentzVector muon_2;
  TLorentzVector Zhltreco;

  int count_passes = 0;
  int passcut_muon = 0;
  int passcut_Z = 0;
  int passcut_Jet = 0;
  int passjets = 0;
  bool printbeginrun = false;
  bool print = true;

};

ZplusJetsMonitor::ZplusJetsMonitor(const edm::ParameterSet& iConfig)
    : folderName_(iConfig.getParameter<std::string>("FolderName")),
      processName_(iConfig.getParameter<std::string>("processName")),
      triggerEventObject_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerEventObject"))),
      ptcut_(iConfig.getParameter<double>("ptcut")),
      isPFJetTrig(iConfig.getParameter<bool>("ispfjettrg")),
      isCaloJetTrig(iConfig.getParameter<bool>("iscalojettrg")),
      directbalanceBinning_(
          iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<std::vector<double> >("directbalanceBinning")),
      ZptBinning_(
          iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<std::vector<double> >("ZptBinning")),
      jetptBinning_(
          iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<std::vector<double> >("jetptBinning")),
      triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResultsLabel")) ) {}
				    
ZplusJetsMonitor::~ZplusJetsMonitor() throw() {

}

void ZplusJetsMonitor::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const& iSetup) {
 
  std::string histname, histtitle;
  std::string hist_obtag = "";
  std::string histtitle_obtag = "";
  std::string currentFolder = folderName_;
  ibooker.setCurrentFolder(currentFolder);

  if (isPFJetTrig) {
    hist_obtag = "pfjet";
    histtitle_obtag = "PFJet";
  } else if (isCaloJetTrig) {
    hist_obtag = "calojet";
    histtitle_obtag = "CaloJet";
  } else {
    hist_obtag = "pfjet";
    histtitle_obtag = "PFJet";
  }  //default is pfjet

  histname = "DiMuonMass";
  histtitle = "DiMuonMass";
  bookME(ibooker, mZMassME_, histname, histtitle, 50, 71., 111.);
  histname = "DPhiZJ";
  histtitle = "DPhiZJ";
  bookME(ibooker, mDPhiZJetME_, histname, histtitle, 100, 0., acos(-1.));
  histname = "ZJetAsymmetry";
  histtitle = "ZJetAsymmetry";
  bookME(ibooker, mZJetAsymmetryME_, histname, histtitle, 100, -1., 1.);

  bookMESub(ibooker, a_ME, sizeof(a_ME) / sizeof(a_ME[0]), hist_obtag, histtitle_obtag, "", "");
  bookMESub(ibooker, a_ME_HB, sizeof(a_ME_HB) / sizeof(a_ME_HB[0]), hist_obtag, histtitle_obtag, "HB", "(HB)", true);
  bookMESub(ibooker, a_ME_HE, sizeof(a_ME_HE) / sizeof(a_ME_HE[0]), hist_obtag, histtitle_obtag, "HE", "(HE)", true);
  bookMESub(ibooker, a_ME_HF, sizeof(a_ME_HF) / sizeof(a_ME_HF[0]), hist_obtag, histtitle_obtag, "HF", "(HF)", true);
  bookMESub(ibooker, a_ME_HB_lowZPt, sizeof(a_ME_HB_lowZPt) / sizeof(a_ME_HB_lowZPt[0]), hist_obtag, histtitle_obtag, "HB_lowZPt", "(HB) lowZPt", true);
  bookMESub(ibooker, a_ME_HE_lowZPt, sizeof(a_ME_HE_lowZPt) / sizeof(a_ME_HE_lowZPt[0]), hist_obtag, histtitle_obtag, "HE_lowZPt", "(HE) lowZPt", true);
  bookMESub(ibooker, a_ME_HF_lowZPt, sizeof(a_ME_HF_lowZPt) / sizeof(a_ME_HF_lowZPt[0]), hist_obtag, histtitle_obtag, "HF_lowZPt", "(HF) lowZPt", true);
  bookMESub(ibooker, a_ME_HB_mediumZPt, sizeof(a_ME_HB_mediumZPt) / sizeof(a_ME_HB_mediumZPt[0]), hist_obtag, histtitle_obtag, "HB_mediumZPt", "(HB) mediumZPt", true);
  bookMESub(ibooker, a_ME_HE_mediumZPt, sizeof(a_ME_HE_mediumZPt) / sizeof(a_ME_HE_mediumZPt[0]), hist_obtag, histtitle_obtag, "HE_mediumZPt", "(HE) mediumZPt", true);
  bookMESub(ibooker, a_ME_HF_mediumZPt, sizeof(a_ME_HF_mediumZPt) / sizeof(a_ME_HF_mediumZPt[0]), hist_obtag, histtitle_obtag, "HF_mediumZPt", "(HF) mediumZPt", true);
  bookMESub(ibooker, a_ME_HB_highZPt, sizeof(a_ME_HB_highZPt) / sizeof(a_ME_HB_highZPt[0]), hist_obtag, histtitle_obtag, "HB_highZPt", "(HB) highZPt", true);
  bookMESub(ibooker, a_ME_HE_highZPt, sizeof(a_ME_HE_highZPt) / sizeof(a_ME_HE_highZPt[0]), hist_obtag, histtitle_obtag, "HE_highZPt", "(HE) highZPt", true);
  bookMESub(ibooker, a_ME_HF_highZPt, sizeof(a_ME_HF_highZPt) / sizeof(a_ME_HF_highZPt[0]), hist_obtag, histtitle_obtag, "HF_highZPt", "(HF) highZPt", true);

}

void ZplusJetsMonitor::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
 
  v_jetpt.clear();
  v_jeteta.clear();
  v_jetphi.clear();

  numtrigobj = 0;
  trigobj_pt.clear();
  trigobj_eta.clear();
  trigobj_phi.clear();
  
  Met_px = 0.;
  Met_py = 0.;

  //Get TriggerResults 
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken_, triggerResults);
  if (!triggerResults.isValid())
    return;

  edm::Handle<trigger::TriggerEvent> aodTriggerEvent;
  iEvent.getByToken(triggerEventObject_, aodTriggerEvent);
  if (!aodTriggerEvent.isValid())
    return;

  edm::TriggerNames triggerNames_ = iEvent.triggerNames(*triggerResults);//all trigger names available 

  //path & module
  std::string mu_pathName="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_PFJet40_v1"; 
  //module name for muons !!
  std::string mu_moduleName="hltDiMuon178Mass8Filtered";       
  
  bool passTrig= false;

 //===============
 //print all paths in each event so not so practical -- prefer to print from begin run if needed
 // for(unsigned i=0; i<triggerNames_.size(); i++){
	//  std::cout<<triggerNames_.triggerName(i)<<std::endl;
 // }
 //===============

  //unsigned hltcfgIndex = hltConfig_.triggerIndex(mu_pathName);
  unsigned index = triggerNames_.triggerIndex(mu_pathName);  // find in which index is that path
  //if(!(hltcfgIndex == index)){
//	std::cout<<"Error in trigger index"<<std::endl;
  //	return;
 // }
  
  if(index>0 && index < triggerNames_.size() && triggerResults->accept(index)){ // if trigger is accepted and INDEX IS VALID 
	if(print) std::cout<<"Trigger accepted"<<std::endl;
	passTrig = true;
	count_passes+=1;
	if(print) std::cout<<"# of trigger passed "<<count_passes<<std::endl;
 	}
  if(!passTrig){
	if(print) std::cout<<"Trigger did not pass"<<std::endl;  // skip event if trigger fails
	return;
   }

  //---------muons--
  //try to find module for this path
  
  const unsigned int module_size=hltConfig_.size(index);
  std::vector<std::string> module_names = hltConfig_.moduleLabels(index); // (pathname) works too 
//
  if(!(module_size == module_names.size())){
	std::cout<<"ERROR IN MODULES COUNTING"<<std::endl;
	return;
     }
  if(module_size == 0){
	std::cout<<"no modules in this path ?!?!"<<std::endl;
	return;
  } //// //forrehlt commended out
  

   // This works if you have the correct module name
   edm::InputTag moduleFilter;
   moduleFilter = edm::InputTag(mu_moduleName, "", processName_);
   if(print) std::cout<<moduleFilter<<"  "<<std::endl;

   // check whether the module is packed up in TriggerEvent product   
   trigger::size_type filterIndex_ = aodTriggerEvent->filterIndex(moduleFilter); 
   if(print) std::cout<< filterIndex_<<"  "<< aodTriggerEvent->sizeFilters() <<std::endl;
   if (filterIndex_>=aodTriggerEvent->sizeFilters()) {
     return;
    } 

   if(print) std:: cout<< " filter label/filter index" << mu_moduleName << "/" << filterIndex_<<std::endl;

   const trigger::Vids& VIDS_ = aodTriggerEvent->filterIds(filterIndex_);
   const trigger::Keys& KEYS_ = aodTriggerEvent->filterKeys(filterIndex_);
   const trigger::size_type nI_ = VIDS_.size();
   const trigger::size_type nK_ = KEYS_.size();
   assert(nI_ == nK_);
   const trigger::TriggerObjectCollection& TOC(aodTriggerEvent->getObjects());
   if(print){
      for (trigger::size_type idx = 0; idx < nI_; ++idx) {
        const trigger::TriggerObject& TO(TOC[KEYS_[idx]]);
        //VIDS :  muons-->83 / jets-->85  / met-->87
        std::cout << " idx:  " << idx << "| vid  " << VIDS_[idx] << "|" <<" keys  "<< KEYS_[idx] << "triggerobject: " << " obj_id " << TO.id() << " Pt " << TO.pt() << " eta " << TO.eta() << " phi " << TO.phi() << " mass " << TO.mass()<<std::endl;
        }
    }//if print 

  if(nK_<2){
	if(print) std:: cout<<" under 2 objects cant have a dimuon"<<std::endl;
	return;
     } 
     else {
        for (const auto& key : KEYS_) {
            double pt = TOC[key].pt();
            double eta = TOC[key].eta();
            double phi = TOC[key].phi();
	    double mass = TOC[key].mass();
            int id = TOC[key].id(); 
            unsigned int kCnt0 = 0;

            TLorentzVector v1;  //keep first muon 
            if (abs(id) == 13){  // check if it is a muon
              v1.SetPtEtaPhiM(pt, eta, phi, mass);
	    }
            else{
             v1.SetPtEtaPhiM(0., 0., 0., 0.);
	    }
            unsigned int kCnt1 = 0;
            for (const auto& key1 : KEYS_) {
              if (key != key1 && kCnt1 > kCnt0) {  // avoid double counting separate objs

                double pt2 = TOC[key1].pt();
                double eta2 = TOC[key1].eta();
                double phi2 = TOC[key1].phi();
		double mass2 = TOC[key1].mass();
                int id2 = TOC[key1].id();
     
                if ((id + id2) == 0) {  // check di-object system charge and flavor

                  TLorentzVector v2;
                  if (abs(id2) == 13){  // check if it is a muon
                    v2.SetPtEtaPhiM(pt2, eta2, phi2, mass2); //muon_m);
		  }
                  else{
                    v2.SetPtEtaPhiM(0., 0., 0., 0.0);
		  }

                  TLorentzVector dimuon = v1 + v2;

                  //selection cuts: these are from twiki but they may need some changes 
               //   if(v1.Pt()>20. && v2.Pt()>20. && abs(v1.Eta())<2.3  && abs(v2.Eta())<2.3){
                   //  if(fabs(dimuon.M()-91.2)<20.0 && dimuon.Pt()>30.0){
       			muon_1 = v1;
 			muon_2 = v2;
  			Zhltreco = muon_1+muon_2;
		//	std::cout<<"Z Pt "<<Zhltreco.Pt() <<std::endl;
		//	mZMassME_.numerator->Fill(Zhltreco.M());
 		 //   }//end Z cuts
		//  }//end muon cuts
           
                }//end check di-object
                else{ //if not di-object
		    return;
		}
             }// end avoid duplicate objects
	    kCnt1++;
 	   }// key1
           kCnt0++;
        }// key
     }// end else  
    //////}//endif 

  //=========== Jets ==============
  // index of last module executed in this Path
  const unsigned int moduleIndex= triggerResults->index(index); // here would be HLTBool at the end
  if(print) std::cout<<moduleIndex-1<<"  "<<module_names[moduleIndex-1]<<std::endl;  // the second to last would be the last module that is saved 
  assert(moduleIndex < module_size);

  // results from TriggerEvent product
  const std::string& ImoduleLabel=module_names[moduleIndex-1];   
  const std::string ImoduleType=hltConfig_.moduleType(ImoduleLabel);
  if(print) std::cout <<ImoduleLabel<< "   "<<ImoduleType<<std::endl;
  // check whether the module is packed up in TriggerEvent product
  const unsigned int filterIndex = aodTriggerEvent->filterIndex(edm::InputTag(ImoduleLabel, "", processName_));
  if (filterIndex < aodTriggerEvent->sizeFilters()) {
		if(print) std::cout<< " filter  - label/type " << ImoduleLabel << "/" << ImoduleType<<std::endl;
		const trigger::Vids& VIDS = aodTriggerEvent->filterIds(filterIndex);
		const trigger::Keys& KEYS = aodTriggerEvent->filterKeys(filterIndex);
		const trigger::size_type nI = VIDS.size();
		const trigger::size_type nK = KEYS.size();
		assert(nI == nK);
		const trigger::TriggerObjectCollection& objects(aodTriggerEvent->getObjects());
		if(print){
		      std::cout<<nI<<std::endl;
		      for (trigger::size_type idx = 0; idx < nI; ++idx) {
			const trigger::TriggerObject& TO_(objects[KEYS[idx]]);
			//VIDS :  muons-->83 / jets-->85  / met-->87
			std::cout << " idx  " << idx << " vid  " << VIDS[idx] << "/" <<" keys  "<< KEYS[idx] << ": " << " obj_id " << TO_.id() << " " << TO_.pt() << " " << TO_.eta() << " " << TO_.phi() << " " << TO_.mass()<<std::endl;
        		}
   		}

		for (const auto& key : KEYS) {
			//selection cuts: will change
			//if(objects[key].pt()>ptcut_ && fabs(objects[key].eta())<1.3){
			//	std::cout<<objects[key].pt()<<std::endl;
			v_jetpt.push_back(objects[key].pt());
			v_jeteta.push_back(objects[key].eta());
			v_jetphi.push_back(objects[key].phi());
			//}//se cuts
		}
	}
    
  double jetpt_ = -999999.;
  double jeteta_ = -99.;
  double jetphi_ = -99.;
  double Z_Pt = -99.;
  double dphi = -99.;
  double assymetry = -9999.;  
  double ZJetBalance_ = - 9999.;
  double J1Pt_minus_ZPt_ = -9999.;
  bool muon_pass = muon_1.Pt()>20. &&  muon_2.Pt()>20. && abs(muon_1.Eta())<2.3 && abs(muon_2.Eta())<2.3;
  bool Z_pass = fabs(Zhltreco.M()-91.2)<20.0 && Zhltreco.Pt()>15.0;
  bool Jet_pass = !v_jetpt.empty() && v_jetpt[0]>=ptcut_;
  if(!muon_pass){ // && Z_pass && Jet_pass){
	return;
   }

  passcut_muon+=1;
  std::cout<<"muon cut passed # "<<passcut_muon<<std::endl;
  if(!Z_pass){
	return;
  }
  passcut_Z+=1;
  std::cout<<"Z cut passed # "<<passcut_Z<<std::endl;
  if(!Jet_pass){
	return;
  }
  passcut_Jet+=1;
  std::cout<<"Jet cut passed # "<<passcut_Jet<<std::endl;

  jetpt_ = v_jetpt[0];
  jeteta_ = v_jeteta[0];
  jetphi_ = v_jetphi[0];
  Z_Pt = Zhltreco.Pt();
			
  dphi = fabs(jetphi_ - Zhltreco.Phi());
  if (dphi > acos(-1.)) {
	dphi = 2 * acos(-1.) - dphi;
  }
  if(dphi<2.5){
	  return;	 
  }
  mZMassME_.numerator->Fill(Zhltreco.M());
  mDPhiZJetME_.numerator->Fill(dphi);
  J1Pt_minus_ZPt_ = jetpt_ - Zhltreco.Pt();
  ZJetBalance_ = jetpt_ / Zhltreco.Pt();
  assymetry = (Z_Pt - jetpt_) / (Z_Pt + jetpt_);
  mZJetAsymmetryME_.numerator->Fill(assymetry);
	
  //=========================================================================================================================  
  FillME(a_ME, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
  if (isBarrel(jeteta_)) {
    FillME(a_ME_HB, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    if (islowZPt(Z_Pt)){
      FillME(a_ME_HB_lowZPt, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    }
    else if (ismediumZPt(Z_Pt)){
      FillME(a_ME_HB_mediumZPt, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    }
    else if (ishighZPt(Z_Pt)){
      FillME(a_ME_HB_highZPt, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    }
  }// end is barrel
  if (isEndCap(jeteta_)) {
    FillME(a_ME_HE, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    if (islowZPt(Z_Pt)){
      FillME(a_ME_HE_lowZPt, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    }
    else if (ismediumZPt(Z_Pt)){
      FillME(a_ME_HE_mediumZPt, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    }
    else if (ishighZPt(Z_Pt)){
      FillME(a_ME_HE_highZPt, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    }
  }// end is endcap
  if (isForward(jeteta_)) {
    FillME(a_ME_HF, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    if (islowZPt(Z_Pt)){
      FillME(a_ME_HF_lowZPt, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    }
    else if (ismediumZPt(Z_Pt)){
      FillME(a_ME_HF_mediumZPt, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    }
    else if (ishighZPt(Z_Pt)){
      FillME(a_ME_HF_highZPt, ZJetBalance_, J1Pt_minus_ZPt_, Z_Pt, jetpt_, true);
    }
  }// end is Forward
	 
}
// Z Pt categories: low, medium, high 
bool ZplusJetsMonitor::islowZPt(double ZPt) {      //(30<pTZ<90)
  bool output = false;
  if (ZPt >30. && ZPt < 90.)
    output = true;
  return output;
}

bool ZplusJetsMonitor::ismediumZPt(double ZPt) {   //(90<pTZ<140)
  bool output = false;
  if (ZPt > 90. && ZPt < 140.)
    output = true;
  return output;
}

bool ZplusJetsMonitor::ishighZPt(double ZPt) {   //(pTZ>140)
  bool output = false;
  if (ZPt > 140.)
    output = true;
  return output;
}
//detector areas
bool ZplusJetsMonitor::isBarrel(double eta) {  
  bool output = false;
  if (fabs(eta) <= 1.3)
    output = true;
  return output;
}

bool ZplusJetsMonitor::isEndCap(double eta) {
  bool output = false;
  if (fabs(eta) <= 3.0 && fabs(eta) > 1.3)
    output = true; 
  return output;
}

/// For Hcal Forward Area
bool ZplusJetsMonitor::isForward(double eta) {
  bool output = false;
  if (fabs(eta) > 3.0)
    output = true;
  return output;
}

void ZplusJetsMonitor::FillME(ObjME* a_me,
                        const double ZJetB_,
                        const double JetminusZ_,
                        const double Zpt_,
 	      		const double Jetpt_,
                        const bool doDirectBalancevsZPt) {

    a_me[0].numerator->Fill(ZJetB_);      //  index 0 = ZJetBalance
    a_me[1].numerator->Fill(JetminusZ_);  //  index 1 = J1Pt_minus_ZPt_
    a_me[2].numerator->Fill(Zpt_);	  //  index 2 = Pt Z
    a_me[3].numerator->Fill(Jetpt_);      //  index 3 = Jet Pt 
    if (doDirectBalancevsZPt){
      a_me[4].numerator->Fill(Zpt_,ZJetB_);    // index 4 = Balance vs Z Pt
     }
}

void ZplusJetsMonitor::bookMESub(DQMStore::IBooker& Ibooker,
                           ObjME* a_me,
                           const int len_,
                           const std::string& h_Name,
                           const std::string& h_Title,
                           const std::string& h_subOptName,
                           const std::string& hSubT,
                           const bool doDirectBalancevsZPt) {

  std::string hName = h_Name;
  std::string hTitle = h_Title;
  const std::string hSubN = h_subOptName.empty() ? "" : "_" + h_subOptName;

  int nbin_JoZ = JoZ_binning_.nbins;
  double maxbin_JoZ = JoZ_binning_.xmax;
  double minbin_JoZ = JoZ_binning_.xmin;

  hName = h_Name + "DirectBalance" + hSubN;
  hTitle = h_Title + " DirectBalance " + hSubT;
  bookME(Ibooker, a_me[0], hName, hTitle,  directbalanceBinning_); 
  setMETitle(a_me[0], h_Title + "", "events");

  hName = h_Name + "J1Pt_minus_ZPt" + hSubN;
  hTitle = h_Title + "J1Pt_minus_ZPt " + hSubT;
  bookME(Ibooker, a_me[1], hName, hTitle, nbin_JoZ, minbin_JoZ, maxbin_JoZ);
  setMETitle(a_me[1], h_Title, "events");

  hName = h_Name + "ZpT" + hSubN;
  hTitle = h_Title + " Z pT " + hSubT;
  bookME(Ibooker, a_me[2], hName, hTitle, ZptBinning_);
  setMETitle(a_me[2], h_Title + " pT [GeV]", "events / [GeV]");
 
  hName = h_Name + "JetpT" + hSubN;
  hTitle = h_Title + " Jet pT " + hSubT;
  bookME(Ibooker, a_me[3], hName, hTitle, jetptBinning_); 
  setMETitle(a_me[3], h_Title + " pT [GeV]", "events / [GeV]");

  if (doDirectBalancevsZPt) {
    hName = h_Name + "DirectBalancevsZPt" + hSubN;
    hTitle = h_Title + " Direct Balance vs ZPt " + hSubT;
    bookME(Ibooker, a_me[4], hName, hTitle, ZptBinning_,directbalanceBinning_); 
    setMETitle(a_me[4], h_Title +" Zpt", "Direct Balance");
  }

}

void ZplusJetsMonitor::dqmBeginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHighLevelTrigger#Access_to_the_HLT_configuration
  // "init" return value indicates whether intitialisation has succeeded
  // "changed" parameter indicates whether the config has actually changed
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    // if init returns TRUE, initialisation has succeeded!
    if (changed) {
     // The HLT config has actually changed wrt the previous Run, hence rebook your
     // histograms or do anything else dependent on the revised HLT config
     std::vector<std::string> triggerPaths = hltConfig_.triggerNames();
     for (const auto& pathName : triggerPaths) {
	if(printbeginrun) std::cout << "ZplusJetsMonitor::dqmBeginRun " << pathName << std::endl;	
     }
    }
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
	edm::LogError("ZplusJetsMonitor") << " HLT config extraction failure with process name " << processName_;
    // In this case, all access methods will return empty values!
  }
}

void ZplusJetsMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("FolderName", "HLT/JME/ZplusJets");
  desc.add<std::string>("processName","reHLT");
  desc.add<edm::InputTag>("triggerEventObject", edm::InputTag("hltTriggerSummaryAOD::reHLT"));

  desc.add<double>("ptcut", 40.);
  desc.add<bool>("ispfjettrg", true);
  desc.add<bool>("iscalojettrg", false);

  edm::ParameterSetDescription histoPSet;

  std::vector<double> bins = {-3.99,-3.97,-3.95,-3.93,-3.91,-3.89,-3.87,-3.85,-3.83,-3.81,-3.79,-3.77,-3.75,-3.73,-3.71,-3.69,-3.67,-3.65,-3.63,-3.61,-3.59,-3.57,-3.55,
-3.53,-3.51,-3.49,-3.47,-3.45,-3.43,-3.41,-3.39,-3.37,-3.35,-3.33,-3.31,-3.29,-3.27,-3.25,-3.23,-3.21,-3.19,-3.17,-3.15,-3.13,-3.11,-3.09,-3.07,
-3.05,-3.03,-3.01,-2.99,-2.97,-2.95,-2.93,-2.91,-2.89,-2.87,-2.85,-2.83,-2.81,-2.79,-2.77,-2.75,-2.73,-2.71,-2.69,-2.67,-2.65,-2.63,-2.61,-2.59,
-2.57,-2.55,-2.53,-2.51,-2.49,-2.47,-2.45,-2.43,-2.41,-2.39,-2.37,-2.35,-2.33,-2.31,-2.29,-2.27,-2.25,-2.23,-2.21,-2.19,-2.17,-2.15,-2.13,-2.11,
-2.09,-2.07,-2.05,-2.03,-2.01,-1.99,-1.97,-1.95,-1.93,-1.91,-1.89,-1.87,-1.85,-1.83,-1.81,-1.79,-1.77,-1.75,-1.73,-1.71,-1.69,-1.67,-1.65,-1.63,
-1.61,-1.59,-1.57,-1.55,-1.53,-1.51,-1.49,-1.47,-1.45,-1.43,-1.41,-1.39,-1.37,-1.35,-1.33,-1.31,-1.29,-1.27,-1.25,-1.23,-1.21,-1.19,-1.17,-1.15,
-1.13,-1.11,-1.09,-1.07,-1.05,-1.03,-1.01,-0.99,-0.97,-0.95,-0.93,-0.91,-0.89,-0.87,-0.85,-0.83,-0.81,-0.79,-0.77,-0.75,-0.73,-0.71,-0.69,-0.67,
-0.65,-0.63,-0.61,-0.59,-0.57,-0.55,-0.53,-0.51,-0.49,-0.47,-0.45,-0.43,-0.41,-0.39,-0.37,-0.35,-0.33,-0.31,-0.29,-0.27,-0.25,-0.23,-0.21,-0.19,
-0.17,-0.15,-0.13,-0.11,-0.09,-0.07,-0.05,-0.03,-0.01,0.01,0.03,0.05,0.07,0.09,0.11,0.13,0.15,0.17,0.19,0.21,0.23,0.25,0.27,0.29,0.31,0.33,0.35,
0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65,0.67,0.69,0.71,0.73,0.75,0.77,0.79,0.81,0.83,0.85,0.87,0.89,0.91,0.93,
0.95,0.97,0.99,1.01,1.03,1.05,1.07,1.09,1.11,1.13,1.15,1.17,1.19,1.21,1.23,1.25,1.27,1.29,1.31,1.33,1.35,1.37,1.39,1.41,1.43,1.45,1.47,1.49,1.51,
1.53,1.55,1.57,1.59,1.61,1.63,1.65,1.67,1.69,1.71,1.73,1.75,1.77,1.79,1.81,1.83,1.85,1.87,1.89,1.91,1.93,1.95,1.97,1.99,2.01,2.03,2.05,2.07,2.09,
2.11,2.13,2.15,2.17,2.19,2.21,2.23,2.25,2.27,2.29,2.31,2.33,2.35,2.37,2.39,2.41,2.43,2.45,2.47,2.49,2.51,2.53,2.55,2.57,2.59,2.61,2.63,2.65,2.67,
2.69,2.71,2.73,2.75,2.77,2.79,2.81,2.83,2.85,2.87,2.89,2.91,2.93,2.95,2.97,2.99,3.01,3.03,3.05,3.07,3.09,3.11,3.13,3.15,3.17,3.19,3.21,3.23,3.25,
3.27,3.29,3.31,3.33,3.35,3.37,3.39,3.41,3.43,3.45,3.47,3.49,3.51,3.53,3.55,3.57,3.59,3.61,3.63,3.65,3.67,3.69,3.71,3.73,3.75,3.77,3.79,3.81,3.83,
3.85,3.87,3.89,3.91,3.93,3.95,3.97,3.99,4.01,4.03,4.05,4.07,4.09,4.11,4.13,4.15,4.17,4.19,4.21,4.23,4.25,4.27,4.29,4.31,4.33,4.35,4.37,4.39,4.41,
4.43,4.45,4.47,4.49,4.51,4.53,4.55,4.57,4.59,4.61,4.63,4.65,4.67,4.69,4.71,4.73,4.75,4.77,4.79,4.81,4.83,4.85,4.87,4.89,4.91,4.93,4.95,4.97,4.99,
5.01,5.03,5.05,5.07,5.09,5.11,5.13,5.15,5.17,5.19,5.21,5.23,5.25,5.27,5.29,5.31,5.33,5.35,5.37,5.39,5.41,5.43,5.45,5.47,5.49,5.51,5.53,5.55,5.57,
5.59,5.61,5.63,5.65,5.67,5.69,5.71,5.73,5.75,5.77,5.79,5.81,5.83,5.85,5.87,5.89,5.91,5.93,5.95,5.97,5.99};

  histoPSet.add<std::vector<double> >("directbalanceBinning", bins);
  std::vector<double> bins_ = {12, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 1500};  // Z pT Binning
  histoPSet.add<std::vector<double> >("ZptBinning", bins_);
  std::vector<double> Jbins_ = {
      0.,   20.,  40.,  60.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160.,
      170., 180., 190., 200., 220., 240., 260., 280., 300., 350., 400., 450., 1000.};  // Jet pT Binning
  histoPSet.add<std::vector<double> >("jetptBinning", Jbins_);

  desc.add<edm::InputTag>("TriggerResultsLabel", edm::InputTag("TriggerResults::reHLT"));
  desc.add<edm::ParameterSetDescription>("histoPSet", histoPSet);

  descriptions.add("zjetsmonitoring", desc);
}

DEFINE_FWK_MODULE(ZplusJetsMonitor);
