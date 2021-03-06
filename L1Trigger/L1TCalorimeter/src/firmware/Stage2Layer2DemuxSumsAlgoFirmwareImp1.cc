///
/// \class l1t::Stage2Layer2SumsAlgorithmFirmwareImp1
///
/// \author:
///
/// Description:

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2DemuxSumsAlgoFirmware.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"

#include <vector>
#include <algorithm>


l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::Stage2Layer2DemuxSumsAlgoFirmwareImp1(CaloParamsHelper* params) :
  params_(params), cordic_(Cordic(144*16,12,8))  // These are the settings in the hardware - should probably make this configurable
{
}


l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::~Stage2Layer2DemuxSumsAlgoFirmwareImp1() {


}


void l1t::Stage2Layer2DemuxSumsAlgoFirmwareImp1::processEvent(const std::vector<l1t::EtSum> & inputSums,
                                                              std::vector<l1t::EtSum> & outputSums) {

  int32_t et(0), metx(0), mety(0), metx2(0), mety2(0), ht(0), mhtx(0), mhty(0), metPhi(0), metPhi2(0), mhtPhi(0);
  uint32_t met(0), met2(0), mht(0);
  uint32_t mbp0(0), mbm0(0), mbp1(0), mbm1(0);

  // Add up the x, y and scalar components
  for (std::vector<l1t::EtSum>::const_iterator eSum = inputSums.begin() ; eSum != inputSums.end() ; ++eSum )
    {
      switch (eSum->getType()) {

      case l1t::EtSum::EtSumType::kTotalEt:
        et += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEtx:
        metx += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEty:
        mety += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalHt:
        ht += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalHtx:
        mhtx += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalHty:
        mhty += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEtx2:
        metx2 += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kTotalEty2:
        mety2 += eSum->hwPt();
        break;

      case l1t::EtSum::EtSumType::kMinBiasHFP0:
	mbp0 = eSum->hwPt();
	break;

      case l1t::EtSum::EtSumType::kMinBiasHFM0:
	mbm0 = eSum->hwPt();
	break;

      case l1t::EtSum::EtSumType::kMinBiasHFP1:
	mbp1 = eSum->hwPt();
	break;

      case l1t::EtSum::EtSumType::kMinBiasHFM1:
	mbm1 = eSum->hwPt();
	break;

      default:
        continue; // Should throw an exception or something?
      }
    }
  
  // leave out to preserve bitwise agreement with l1t-tsg-v6-cand:
  //if (et>0xFFF)   et   = 0xFFF;
  //if (metx>0xFFF) metx = 0xFFF;
  //if (mety>0xFFF) mety = 0xFFF;
  //if (ht>0xFFF)   ht   = 0xFFF;
  //if (mhtx>0xFFF) mhtx = 0xFFF;
  //if (mhty>0xFFF) mhty = 0xFFF;
  //if (metx2>0xFFF) metx2 = 0xFFF;
  //if (mety2>0xFFF) mety2 = 0xFFF;

  
  // Final MET calculation
  if (metx != 0 || mety != 0 ) cordic_( metx , mety , metPhi , met );
  // sets the met scale back to the original range for output into GT, this corresponds to
  // the previous scaling of sin/cos factors in calculation of metx and mety by 2^10 = 1024
  met >>= 10; 

  // Final MET2 calculation
  if (metx2 != 0 || mety2 != 0 ) cordic_( metx2 , mety2 , metPhi2 , met2 );

  // Final MHT calculation
  if (mhtx != 0 || mhty != 0 ) cordic_( mhtx , mhty , mhtPhi , mht );
  // sets the mht scale back to the original range for output into GT, the other 4
  // bits are brought back just before the accumulation of ring sum in MP jet sum algorithm
  mht >>= 6; 

  // Make final collection
  math::XYZTLorentzVector p4;

  l1t::EtSum etSumTotalEt(p4,l1t::EtSum::EtSumType::kTotalEt,et,0,0,0);
  l1t::EtSum etSumMissingEt(p4,l1t::EtSum::EtSumType::kMissingEt,met,0,metPhi>>4,0);
  l1t::EtSum etSumMissingEt2(p4,l1t::EtSum::EtSumType::kMissingEt2,met2,0,metPhi2>>4,0);
  l1t::EtSum htSumht(p4,l1t::EtSum::EtSumType::kTotalHt,ht,0,0,0);
  l1t::EtSum htSumMissingHt(p4,l1t::EtSum::EtSumType::kMissingHt,mht,0,mhtPhi>>4,0);
  l1t::EtSum etSumMinBiasHFP0(p4,l1t::EtSum::EtSumType::kMinBiasHFP0,mbp0,0,0,0);
  l1t::EtSum etSumMinBiasHFM0(p4,l1t::EtSum::EtSumType::kMinBiasHFM0,mbm0,0,0,0);
  l1t::EtSum etSumMinBiasHFP1(p4,l1t::EtSum::EtSumType::kMinBiasHFP1,mbp1,0,0,0);
  l1t::EtSum etSumMinBiasHFM1(p4,l1t::EtSum::EtSumType::kMinBiasHFM1,mbm1,0,0,0);

  outputSums.push_back(etSumTotalEt);
  outputSums.push_back(etSumMinBiasHFP0);
  outputSums.push_back(htSumht);
  outputSums.push_back(etSumMinBiasHFM0);
  outputSums.push_back(etSumMissingEt);
  outputSums.push_back(etSumMinBiasHFP1);
  outputSums.push_back(htSumMissingHt);
  outputSums.push_back(etSumMinBiasHFM1);
  outputSums.push_back(etSumMissingEt2);

}
