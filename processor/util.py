import numpy as np

# jet id
def get_jetid_mask(jets, jet_type="CHS", use_lepton_veto=True):
    #year = events.metadata.get("year", year)
    #era = events.metadata.get("era", era)
    #if year is None:
    #    raise ValueError("year cannot be None")
    #if year[-1].isalpha():
    #    year = year[:-1]
    #    era = [-1]
    #if era is None:
    #    raise ValueError("era cannot be None")

    # https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV (v17)
    # Recommendations for the 13.6 TeV data analysis: runs 2022 re-reco CDE, 2022 FG and 2023 BCD
    #jets = events[jet_name]

    abs_eta = np.abs(jets.eta)

    # it is unclear whether to use nConstElecs or nElectrons fields...
    # from a quick test, it seems like jetId uses nConstElecs
    #nElectrons = jets["nConstElecs"] if "nConstElecs" in jets.fields else 
    # also for nConstMuons and nMuons...
    # from a quick test, it seems like jetId uses nConstMuons
    #nMuons = jets["nConstMuons"] if "nConstMuons" in jets.fields else jets["nMuons"]
    #nMuons = jets["nMuons"]
    #nPhotons = jets["nConstPhotons"] if "nConstPhotons" in jets.fields else jets["nPhotons"]
    #nCh = jets["nConstChHads"] if "nConstChHads" in jets.fields else jets["nCh"]
    #nNh = jets["nConstNeuHads"] if "nConstNeuHads" in jets.fields else jets["nNh"]

    nElectrons = jets["nElectrons"]
    nMuons = jets["nMuons"]
    nPhotons = jets["nPhotons"]
    nCh = jets["nCh"]
    nNh = jets["nNh"]

    #NumConst = jets["nConstituents"]
    NumConst = nElectrons+nMuons+nPhotons+nCh+nNh
    CHM = nCh + nElectrons + nMuons
    CEMF = jets["chEmEF"]
    CHF = jets["chHEF"]
    NumNeutralParticle = nNh + nPhotons
    NEMF = jets["neEmEF"]
    NHF = jets["neHEF"]
    #MUF = jets["muEF"] if "muEF" in jets.fields else jets["muEmEF"]
    MUF = jets["muEF"]

    # mask number corresponds to column in Jet ID table
    # mask0 |eta| in [0, 2.6]
    # mask1 |eta| in (2.6, 2.7]
    # mask2 |eta| in (2.7, 3.0]
    # mask3 |eta| in (3.0, 5.0]
    mask0 = (abs_eta<=2.6) & (CHM>0) & (CHF>0.01) & (NumConst>1) & (NEMF<0.9) & (NHF<0.99)
    if use_lepton_veto:
        mask0 = mask0 & (CEMF<0.8) & (MUF<0.8)
                    
    mask1 = (abs_eta>2.6) & (abs_eta<=2.7) & (NEMF<0.99) & (NHF<0.9) 
    if use_lepton_veto:
        mask1 = mask1 & (CEMF<0.8) & (MUF<0.8)
    if jet_type == "CHS":
        mask1 = mask1 & (CHM>0)
                    
    mask2 = (abs_eta>2.7) & (abs_eta<=3.0) & (NHF<0.99)
    if jet_type == "CHS":
        mask2 = mask2 & (NEMF<0.99) & (NumNeutralParticle>1)
                    
    mask3 = (abs_eta>3.0) & (NEMF<0.40)
    if jet_type == "PUPPI":
        mask3 = mask3 & (NumNeutralParticle>=2)
    elif jet_type == "CHS":
        mask3 = mask3 & (NumNeutralParticle>10)
                    
    mask = mask0 | mask1 | mask2 | mask3
    return mask

# trigger
def get_trigger_mask(events, trigger, prescale_dict=dict()):
    if isinstance(trigger, str):
        trigger = trigger.split()
    
    mask = events[tuple(trigger[0].split("_", 1))]
    if trigger[0] in prescale_dict:
        ps = prescale_dict[trigger[0]]
        mask = mask & np.random.choice([True, False], size=len(events), p=[1/ps, 1-1/ps])
        
    if len(trigger) == 1:
        return mask
    else:
        if trigger[1] in {"OR", "AND", "|", "&"}:
            for idx, trg in enumerate(trigger[1:]):
                if trg in {"OR", "AND", "|", "&"}:
                    continue
                else:
                    this_mask = events[tuple(trg.split("_", 1))]
                    if trg in prescale_dict:
                        ps = prescale_dict[trg]
                        this_mask = this_mask & np.random.choice([True, False], size=len(events), p=[1/ps, 1-1/ps])
                        
                    if trigger[idx] == "OR" or trigger[idx] == "|":
                        mask = mask | this_mask
                    else:
                        mask = mask & this_mask
        else:
            for trg in trigger[1:]:
                if trg in prescale_dict:
                    ps = prescale_dict[trg]
                    this_mask = this_mask & np.random.choice([True, False], size=len(events), p=[1/ps, 1-1/ps])
                mask = mask | this_mask
    return mask