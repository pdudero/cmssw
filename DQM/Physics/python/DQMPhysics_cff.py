import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

from DQM.Physics.bphysicsOniaDQM_cfi import *
from DQM.Physics.ewkMuDQM_cfi import *
from DQM.Physics.ewkElecDQM_cfi import *
from DQM.Physics.ewkMuLumiMonitorDQM_cfi import *
from DQM.Physics.qcdPhotonsDQM_cfi import *
from DQM.Physics.topSingleLeptonDQM_cfi import *
from DQM.Physics.topDiLeptonOfflineDQM_cfi import *
from DQM.Physics.topSingleLeptonDQM_PU_cfi import *
from DQM.Physics.singleTopDQM_cfi import *
from DQM.Physics.ewkMuLumiMonitorDQM_cfi import *
from DQM.Physics.susyDQM_cfi import *
from DQM.Physics.HiggsDQM_cfi import *
from DQM.Physics.ExoticaDQM_cfi import *
from DQM.Physics.B2GDQM_cfi import *
from DQM.Physics.CentralityDQM_cfi import *
from DQM.Physics.topJetCorrectionHelper_cfi import *

dqmPhysics = cms.Sequence( bphysicsOniaDQM 
                           *ewkMuDQM
                           *ewkElecDQM
                           *ewkMuLumiMonitorDQM
                           *qcdPhotonsDQM
			   *topSingleMuonMediumDQM
                           *topSingleElectronMediumDQM	
                           *singleTopMuonMediumDQM
                           *singleTopElectronMediumDQM
                           *DiMuonDQM
			   *DiElectronDQM
			   *ElecMuonDQM
                           *susyDQM
                           *HiggsDQM
                           *ExoticaDQM
                           *B2GDQM
                           )
eras.phase1Pixel.toReplaceWith(dqmPhysics, dqmPhysics.copyAndExclude([ # FIXME
    ewkMuDQM,            # Excessive printouts because 2017 doesn't have HLT yet
    ewkElecDQM,          # Excessive printouts because 2017 doesn't have HLT yet
    ewkMuLumiMonitorDQM, # Excessive printouts because 2017 doesn't have HLT yet
]))

bphysicsOniaDQMHI = bphysicsOniaDQM.clone(vertex=cms.InputTag("hiSelectedVertex"))
dqmPhysicsHI = cms.Sequence(bphysicsOniaDQMHI+CentralityDQM)

from DQM.Physics.qcdPhotonsCosmicDQM_cff import *
dqmPhysicsCosmics = cms.Sequence(dqmPhysics)
dqmPhysicsCosmics.replace(qcdPhotonsDQM, qcdPhotonsCosmicDQM)
