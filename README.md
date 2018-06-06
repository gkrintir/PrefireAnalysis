PrefireAnalysis
===============
This package helps analyse the effect of prefiring (particularly from EE in late 2017) on an analysis.

# Installation
This should work in any `94X` CMSSW or newer
```bash
pushd $CMSSW_BASE/src
mkdir -p L1Trigger && cd L1Trigger
git clone git@github.com:nsmith-/PrefireAnalysis.git
cd ..
scram b
popd
```

# Usage
The idea is to select so-called 'unprefirable' events in whatever primary dataset is used for your analysis,
and then to check how often a (vetoed) L1A is present in the previous bunch crossing via the uGT record.

To select the events, the TCDS L1A history needs to be available.  It is found in AOD data as `L1AcceptBunchCrossingCollection`
with tag `scalersRawToDigi`, and an EDFilter is put in front of the analysis chain to read it: 
```python
process.prefireVetoFilter = cms.EDFilter("TriggerRulePrefireVetoFilter",
    l1AcceptRecordLabel = cms.InputTag("scalersRawToDigi"),
)

process.myPath = cms.Path(process.prefireVetoFilter + process.analysisChain)
```
If using MINIAOD, one can use the secondary file solution,
i.e. adding `secondaryFileNames` to `process.source`, or in the case of CRAB, `config.Data.useParent = True`.
A complete example that just saves the filtered MINIAOD events is in `test/testTriggerRulePrefireVeto_cfg.py`.

In your analysis chain, please save the uGT trigger decision (FinOR) of BX -1, accessed as follows:
```c++
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"

// class MyAna : public edm::one::EDAnalyzer<edm::one::SharedResources> {
// ...
    edm::EDGetTokenT<BXVector<GlobalAlgBlk>> l1GtToken_;

// constructor...
  l1GtToken_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("l1GtSrc")))

// analyze function...
  edm::Handle<BXVector<GlobalAlgBlk>> l1GtHandle;
  iEvent.getByToken(l1GtToken_, l1GtHandle);
  bool prefire = l1GtHandle->begin(-1)->getFinalOR();
```
Then check with final selection how often this is true.  If the fraction is less than about 1%, then this means
there is no prefiring, as there is always the possibility (~L1A rate / MinBias rate) that interesting physics is
in the previous bunch crossing.  If much larger, then this is telling you the fraction of data in your phase
space lost due to prefiring.
