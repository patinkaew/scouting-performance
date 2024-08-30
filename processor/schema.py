from coffea.nanoevents import NanoAODSchema

class ScoutingNanoAODSchema(NanoAODSchema):
    
    mixins = {
        **NanoAODSchema.mixins,
        "PFJet": "Jet",
        "ScoutingPFJet": "Jet",
        "ScoutingPFJetRecluster": "Jet",
        "ScoutingCHSJetRecluster": "Jet",
        "ScoutingFatCHSJetRecluster": "FatJet",
        "ScoutingRho": "Rho",
        "ScoutingMET": "MissingET",
        "ScoutingMuon": "Muon",
        "ScoutingMuonNoVtx" : "Muon",
        "ScoutingMuonVtx" : "Muon",
        # JMENano
        "JetCHS" : "Jet",
    }