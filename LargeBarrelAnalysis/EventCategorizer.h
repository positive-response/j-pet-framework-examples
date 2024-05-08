/**
 *  @copyright Copyright 2020 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file EventCategorizer.h
 */

#ifndef EVENTCATEGORIZER_H
#define EVENTCATEGORIZER_H

#include <JPetUserTask/JPetUserTask.h>
#include "EventCategorizerTools.h"
#include <JPetEvent/JPetEvent.h>
#include <JPetHit/JPetHit.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <string>
#include "TFile.h"


class JPetWriter;

/**
 * @brief User Task categorizing Events
 *
 * Task attempts to add types of events to each event. Each category/type
 * has separate method for checking, if current event fulfills set of conditions.
 * These methods are defined in tools class. More than one type can be added to an event.
 * Set of controll histograms are created, unless the user decides not to produce them.
 */
class EventCategorizer : public JPetUserTask{
public:
	EventCategorizer(const char * name);
	virtual ~EventCategorizer();
	virtual bool init() override;
	virtual bool exec() override;
	virtual bool terminate() override;

	static const double kUnknownEventType;

protected:
	const std::string kBack2BackSlotThetaDiffParamKey = "Back2Back_Categorizer_SlotThetaDiff_float";
	const std::string kScatterTOFTimeDiffParamKey = "Scatter_Categorizer_TOF_TimeDiff_float";
	const std::string kDeexTOTCutMinParamKey = "Deex_Categorizer_TOT_Cut_Min_float";
	const std::string kDeexTOTCutMaxParamKey = "Deex_Categorizer_TOT_Cut_Max_float";
	const std::string kMaxTimeDiffParamKey = "EventCategorizer_MaxTimeDiff_float";
	const std::string kSaveControlHistosParamKey = "Save_Control_Histograms_bool";
        const std::string kTOTCalculationType = "HitFinder_TOTCalculationType_std::string";
	const std::string kHitReqParamKey = "hitReq_int";
		
	void saveEvents(const std::vector<JPetEvent>& event);
	int fHitRequired = 5; 
	double fScatterTOFTimeDiff = 2000.0;
	double fB2BSlotThetaDiff = 3.0;
	double fDeexTOTCutMin = 30000.0;
	double fDeexTOTCutMax = 50000.0;
	double fMaxTimeDiff = 1000.;
	bool fSaveControlHistos = true;
        std::string fTOTCalculationType = "";
	void initialiseHistograms();

	int random = rand();
	std::string Name = "test"+std::to_string(random)+".root";
        //TFile *file = TFile::Open(Name.c_str(), "RECREATE");	
	TTree* pTree22 = nullptr;

	double flowEnergyCut = 30;
	double fAnnihilationEnergyCut = 650;
	

        int fTimeWindowNumber = 0;
        int fEventNumber = 0;
        int fNumberOfHits = 0;
        std::vector<float> fPosX {};
        std::vector<float> fPosY {};
        std::vector<float> fPosZ {};
        std::vector<float> fEnergy {};
        std::vector<float> fTime {};
	std::vector<unsigned int> fHitType {};
	std::vector<unsigned int> fVtxIndex {};
        float fRecoOrthoVtxPosX = 0;
        float fRecoOrthoVtxPosY = 0;
        float fRecoOrthoVtxPosZ = 0;
        bool fIsAcc = kFALSE;
        bool fIsOPs = kFALSE;
        bool fIsPickOff = kFALSE;
        bool fContainsPrompt = kFALSE;
        bool fIsScattered = kFALSE;
        bool fIsSecondary = kFALSE;
        float fTimeDiff = -100000.0;
};

#endif /* !EVENTCATEGORIZER_H */

