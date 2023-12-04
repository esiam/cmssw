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
                 const bool doMPF = false);

  void FillME(ObjME* a_me,
              const double ZJetB_,
              const double JetminZ_,
              const double mpf_,
              const bool doMPF = false);  

private:
  const std::string folderName_;

  const std::string processName_; // process name of (HLT) process for which to get HLT configuration
  /// The instance of the HLTConfigProvider as a data member
  HLTConfigProvider hltConfig_;

  const bool requireValidHLTPaths_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerEventObject_;
  //edm::EDGetTokenT<trigger::TriggerEventWithRefs> triggerEventObject_; //for RAW?
  bool hltPathsAreValid_;

  double ptcut_;
  bool isPFJetTrig;
  bool isCaloJetTrig;

  const bool enableFullMonitoring_;

  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
 
  ObjME a_ME[3];
  ObjME a_ME_HB[3];
  ObjME a_ME_HE[3];
  ObjME a_ME_HF[3];
  ObjME a_ME_HB_lowZPt[3];
  ObjME a_ME_HE_lowZPt[3];
  ObjME a_ME_HF_lowZPt[3];
  ObjME a_ME_HB_mediumZPt[3];
  ObjME a_ME_HE_mediumZPt[3];
  ObjME a_ME_HF_mediumZPt[3];
  ObjME a_ME_HB_highZPt[3];
  ObjME a_ME_HE_highZPt[3];
  ObjME a_ME_HF_highZPt[3];
 
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

  MEbinning ZJetB_binning_{50, 0.,3.0};
  MEbinning JoZ_binning_{50, -75.0, 75.0};
  MEbinning MPF_binning_{50, 0., 2.0};

  bool printbeginrun = false;
  bool print = false;

};

ZplusJetsMonitor::ZplusJetsMonitor(const edm::ParameterSet& iConfig)
    : folderName_(iConfig.getParameter<std::string>("FolderName")),
      processName_(iConfig.getParameter<std::string>("processName")),
      requireValidHLTPaths_(iConfig.getParameter<bool>("requireValidHLTPaths")),
      triggerEventObject_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerEventObject"))),
      hltPathsAreValid_(false),
      ptcut_(iConfig.getParameter<double>("ptcut")),
      isPFJetTrig(iConfig.getParameter<bool>("ispfjettrg")),
      isCaloJetTrig(iConfig.getParameter<bool>("iscalojettrg")),
      enableFullMonitoring_(iConfig.getParameter<bool>("enableFullMonitoring")),
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

  //Jet paths
  std::vector<std::string> Jet_pathNames = {"HLT_PFJet80_v24", "HLT_PFJet140_v23", "HLT_PFJet200_v23", "HLT_PFJet260_v24", "HLT_PFJet320_v24", "HLT_PFJet400_v24", "HLT_PFJet450_v25", "HLT_PFJet500_v25"}; 
  
  //Jet modules
  std::vector<std::string> Jet_moduleNames = {"hltSinglePFJet80", "hltSinglePFJet140", "hltSinglePFJet200", "hltSinglePFJet260", "hltSinglePFJet320", "hltSinglePFJet400", "hltSinglePFJet450", "hltSinglePFJet500"};

  //muon path & module
  std::string mu_pathName="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v9"; 
  //Dont use module name hard coded -- it was incorrect !!
  std::string mu_moduleName="hltDiMuon178Mass8Filtered";      //"hltPreMu17TrkIsoVVLMu8TrkIsoVVLDZMass8";    
  
  bool mu_passTrig= false;
 
  unsigned hltcfgIndex = hltConfig_.triggerIndex(mu_pathName);
  unsigned index = triggerNames_.triggerIndex(mu_pathName);  // find in which index is that path
  if(!(hltcfgIndex == index)){
	std::cout<<"Error in trigger index"<<std::endl;
  	return;
  }
  if (index>0 && index < triggerNames_.size() && triggerResults->accept(index) ){ // if trigger is accepted and INDEX IS VALID --> do not forget
	if(print) std::cout<<"Moun trigger accepted"<<std::endl;
	mu_passTrig = true;
 	}
  if(!mu_passTrig){
	if(print) std::cout<<"Moun trigger did not pass"<<std::endl;  // skip event if muon trigger fails
	return;
 }

  //---------muons--
  //try to find module for this path
  const unsigned int module_size=hltConfig_.size(index);
  std::vector<std::string> module_names = hltConfig_.moduleLabels(index); // (pathname) works too 

  if(!(module_size == module_names.size())){
	std::cout<<"ERROR IN MODULES COUNTING"<<std::endl;
	return;
     }
  if(module_size == 0){
	std::cout<<"no modules in this path ?!?!"<<std::endl;
	return;
  } 
  

   // This works if you have the correct module name!!
   /*edm::InputTag moduleFilter;
   moduleFilter = edm::InputTag(mu_moduleName, "", processName_);
   if(print) std::cout<<moduleFilter<<"  "<<std::endl;

   // check whether the module is packed up in TriggerEvent product   
   trigger::size_type filterIndex_ = aodTriggerEvent->filterIndex(moduleFilter); 
   // std::cout<< filterIndex_<<"  "<< aodTriggerEvent->sizeFilters() <<std::endl;
   if (filterIndex_>=aodTriggerEvent->sizeFilters()) {
     return;
    }*/

  // index of last module executed in this Path
  const unsigned int moduleIndex= triggerResults->index(index); // here would be HLTBool at the end
  if(print) std::cout<<moduleIndex-1<<"  "<<module_names[moduleIndex-1]<<std::endl;  // the second to last would be the last module that is saved 
  assert(moduleIndex < module_size);

  // results from TriggerEvent product
  const std::string& ImoduleLabel=module_names[moduleIndex-1];   
  const std::string ImoduleType=hltConfig_.moduleType(ImoduleLabel);
  // std::cout <<ImoduleLabel<< "   "<<ImoduleType<<std::endl;
  // check whether the module is packed up in TriggerEvent product
  const unsigned int filterIndex = aodTriggerEvent->filterIndex(edm::InputTag(ImoduleLabel, "", processName_));
  if (filterIndex<aodTriggerEvent->sizeFilters()) {
	  
   if(print) std:: cout<< " filter label/type " << ImoduleLabel << "/" << ImoduleType<<std::endl;

   const trigger::Vids& VIDS_ = aodTriggerEvent->filterIds(filterIndex);
   const trigger::Keys& KEYS_ = aodTriggerEvent->filterKeys(filterIndex);
   const trigger::size_type nI_ = VIDS_.size();
   const trigger::size_type nK_ = KEYS_.size();
   assert(nI_ == nK_);
   const trigger::TriggerObjectCollection& TOC(aodTriggerEvent->getObjects());
   if(print){
      for (trigger::size_type idx = 0; idx < nI_; ++idx) {
      	const trigger::TriggerObject& TO(TOC[KEYS_[idx]]);
      	//VIDS :  muons-->83 / jets-->85  / met-->87
     	std::cout << " idx  " << idx << " vid  " << VIDS_[idx] << "/" <<" keys  "<< KEYS_[idx] << ": " << " obj_id " << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()<<std::endl;
    	}
   } 

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
                  if(v1.Pt()>20. && v2.Pt()>20. && abs(v1.Eta())<2.3  && abs(v2.Eta())<2.3){
                     if(fabs(dimuon.M()-91.2)<20.0 && dimuon.Pt()>30.0){
       			muon_1 = v1;
 			muon_2 = v2;
  			Zhltreco = muon_1+muon_2;
			//std::cout<<"Z Pt "<<Zhltreco.Pt() <<std::endl;
			mZMassME_.numerator->Fill(Zhltreco.M());
 		    }//end Z cuts
		  }//end muon cuts
           
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
    }//endif

 //------------
 //Jet triggers
 //------------
 bool passJetTrig = false;
 for(int k=0; k<8; k++){
   //std::cout<<k<<". "<<"Jet Path "<< Jet_pathNames[k] <<std::endl; 
   unsigned hltcfgIndex_J = hltConfig_.triggerIndex(Jet_pathNames[k]);
   unsigned index_J = triggerNames_.triggerIndex(Jet_pathNames[k]);
   if(hltcfgIndex_J != index_J){
      return;	
   }
   if (index_J>0 && index_J < triggerNames_.size() && triggerResults->accept(index_J)){    // to do this correct when a jet trigger pass, i need to stop searching
        passJetTrig = true;
	if(print) std::cout<<"Passed "<<Jet_pathNames[k] <<std::endl;
	const unsigned int module_size_J = hltConfig_.size(index_J);
	std::vector<std::string> module_names_J = hltConfig_.moduleLabels(Jet_pathNames[k]);
	if(module_size_J != module_names_J.size()){
		std::cout<<"ERROR IN MODULES COUNTING"<<std::endl;
        	return;
	}
	// index of last module executed in this Path
	const unsigned int moduleIndex_J= triggerResults->index(index_J);
	assert(moduleIndex_J < module_size_J);

	const std::string& ImoduleLabel=module_names_J[moduleIndex_J-1];
	const std::string ImoduleType=hltConfig_.moduleType(ImoduleLabel);
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
		for (const auto& key : KEYS) {
			int id = objects[key].id();
      			if(id == 1){  //is a jet
				//selection cuts: will change
				if(objects[key].pt()>ptcut_ && fabs(objects[key].eta())<1.3){
					v_jetpt.push_back(objects[key].pt());
           				v_jeteta.push_back(objects[key].eta());
           				v_jetphi.push_back(objects[key].phi());
          			}//se cuts
       			}
    		}
	}
	break;
    }// end if is accepted
   
   }//end run over Jet paths

   if(!passJetTrig){
        if(print) std::cout<<"DID NOT PASS JET TRIGGER"<<std::endl;
        return;  //skip event if not one Jet trigger at least pass the event
   }

   //Variables
   double jetpt_ = v_jetpt[0];
   double jeteta_ = v_jeteta[0];
   double jetphi_ = v_jetphi[0];
   double Z_Pt = Zhltreco.Pt();

   double assymetry = (Z_Pt - jetpt_) / (Z_Pt + jetpt_);
   double dphi = fabs(jetphi_ - Zhltreco.Phi());
   if (dphi > acos(-1.)) {
      dphi = 2 * acos(-1.) - dphi;
    }

  mDPhiZJetME_.numerator->Fill(dphi);
 // mZMassME_.numerator->Fill(Zhltreco.M());
  mZJetAsymmetryME_.numerator->Fill(assymetry);

  double ZJetBalance_ = jetpt_ / Zhltreco.Pt();
  double J1Pt_minus_ZPt_ = jetpt_ - Zhltreco.Pt();
  double MPF_ = 1. + (Met_px * Zhltreco.Px() + Met_py * Zhltreco.Py()) / (Zhltreco.Pt() * Zhltreco.Pt());
  
//=========================================================================================================================  
  FillME(a_ME, ZJetBalance_, J1Pt_minus_ZPt_, MPF_);
  if (isBarrel(jeteta_)) {
    FillME(a_ME_HB, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    if (islowZPt(Z_Pt)){
      FillME(a_ME_HB_lowZPt, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    }
    else if (ismediumZPt(Z_Pt)){
      FillME(a_ME_HB_mediumZPt, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    }
    else if (ishighZPt(Z_Pt)){
      FillME(a_ME_HB_highZPt, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    }
  }// end is barrel
  if (isEndCap(jeteta_)) {
    FillME(a_ME_HE, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    if (islowZPt(Z_Pt)){
      FillME(a_ME_HE_lowZPt, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    }
    else if (ismediumZPt(Z_Pt)){
      FillME(a_ME_HE_mediumZPt, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    }
    else if (ishighZPt(Z_Pt)){
      FillME(a_ME_HE_highZPt, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    }
  }// end is endcap
  if (isForward(jeteta_)) {
    FillME(a_ME_HF, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    if (islowZPt(Z_Pt)){
      FillME(a_ME_HF_lowZPt, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    }
    else if (ismediumZPt(Z_Pt)){
      FillME(a_ME_HF_mediumZPt, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    }
    else if (ishighZPt(Z_Pt)){
      FillME(a_ME_HF_highZPt, ZJetBalance_, J1Pt_minus_ZPt_, MPF_, true);
    }
  }// end is Forward  */
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
                        const double JetminZ_,
                        const double mpf_,
                        const bool doMPF) {

    a_me[0].numerator->Fill(ZJetB_);    //  index 0 = ZJetBalance
    a_me[1].numerator->Fill(JetminZ_);  //  index 1 = J1Pt_minus_ZPt_
    if (doMPF){
      a_me[2].numerator->Fill(mpf_);    //index 2 = MPF
     }
}

void ZplusJetsMonitor::bookMESub(DQMStore::IBooker& Ibooker,
                           ObjME* a_me,
                           const int len_,
                           const std::string& h_Name,
                           const std::string& h_Title,
                           const std::string& h_subOptName,
                           const std::string& hSubT,
                           const bool doMPF) {

  std::string hName = h_Name;
  std::string hTitle = h_Title;
  const std::string hSubN = h_subOptName.empty() ? "" : "_" + h_subOptName;

  int nbin_ZJetB = ZJetB_binning_.nbins;
  double maxbin_ZJetB = ZJetB_binning_.xmax;
  double minbin_ZJetB = ZJetB_binning_.xmin;

  int nbin_MPF = MPF_binning_.nbins;
  double maxbin_MPF = MPF_binning_.xmax;
  double minbin_MPF = MPF_binning_.xmin;

  int nbin_JoZ = JoZ_binning_.nbins;
  double maxbin_JoZ = JoZ_binning_.xmax;
  double minbin_JoZ = JoZ_binning_.xmin;

  hName = h_Name + "ZJetBalance" + hSubN;
  hTitle = h_Title + " ZJetBalance " + hSubT;
  bookME(Ibooker, a_me[0], hName, hTitle, nbin_ZJetB, minbin_ZJetB, maxbin_ZJetB);
  setMETitle(a_me[0], h_Title + " [GeV]", "events / [GeV]");

  hName = h_Name + "J1Pt_minus_ZPt" + hSubN;
  hTitle = h_Title + "J1Pt_minus_ZPt " + hSubT;
  bookME(Ibooker, a_me[1], hName, hTitle, nbin_JoZ, minbin_JoZ, maxbin_JoZ);
  setMETitle(a_me[1], h_Title, "events");


  if (doMPF) {
    hName = h_Name + "MPF" + hSubN;
    hTitle = h_Title + " MPF " + hSubT;
    bookME(Ibooker, a_me[2], hName, hTitle, nbin_MPF, minbin_MPF, maxbin_MPF);
    setMETitle(a_me[2], h_Title , "events");
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
  desc.add<std::string>("processName","HLT");
  desc.add<bool>("requireValidHLTPaths", true);
  desc.add<edm::InputTag>("triggerEventObject", edm::InputTag("hltTriggerSummaryAOD::HLT"));

  desc.add<double>("ptcut", 40.);
  desc.add<bool>("ispfjettrg", true);
  desc.add<bool>("iscalojettrg", false);

  desc.add<bool>("enableFullMonitoring", true);

  desc.add<edm::InputTag>("TriggerResultsLabel", edm::InputTag("TriggerResults::HLT"));

  descriptions.add("zjetsmonitoring", desc);
}

DEFINE_FWK_MODULE(ZplusJetsMonitor);
