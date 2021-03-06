import FWCore.ParameterSet.Config as cms

from RecoBTag.CTagging.candidateCombinedSecondaryVertexSoftLeptonCvsLComputer_cfi import *

candidateNegativeCombinedSecondaryVertexSoftLeptonCvsLComputer = candidateCombinedSecondaryVertexSoftLeptonCvsLComputer.clone(
   vertexFlip = True,
   trackFlip = True
)

candidateNegativeCombinedSecondaryVertexSoftLeptonCvsLComputer.trackSelection.sip3dSigMax = 0
candidateNegativeCombinedSecondaryVertexSoftLeptonCvsLComputer.trackPseudoSelection.sip3dSigMax = 0
candidateNegativeCombinedSecondaryVertexSoftLeptonCvsLComputer.trackPseudoSelection.sip2dSigMin = -99999.9
candidateNegativeCombinedSecondaryVertexSoftLeptonCvsLComputer.trackPseudoSelection.sip2dSigMax = -2.0
