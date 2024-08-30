import numpy as np
import awkward as ak
import hist
import hist.dask as hda
from coffea import processor as coffea_processor
from .util import *

pt_bins = np.array([0., 20., 40.,   50, 60.,  70, 80., 90, 100., 110, 120.,  130, 140., 150, 160.,
         170, 180.,  190, 200.,  210, 220., 230, 240., 250, 260., 280., 300.,  320., 340.,
         360., 380.,  400., 450., 500., 600., 800, 1000, 1200, 1400, 1800])

ht_bins = np.array([0., 20.,   40.,   60.,   80.,  100.,  120.,  140.,  160.,
         180.,  200.,  220.,  240.,  260.,  280.,  300.,  320.,  340.,
         360.,  380.,  400., 450., 500., 600., 700, 800, 900, 1000, 1200, 1400, 1800, 2000, 2500, 3000])

def get_necessary_trigger_names(trigger_tasks):
    necessary_trigger_names = set()
    for trigger_task in trigger_tasks:
        necessary_trigger_names.add(trigger_task["signal"]["name"])
        necessary_trigger_names.add(trigger_task["reference"]["name"])
        necessary_trigger_names.add(f"({trigger_task["signal"]["name"]}) AND ({trigger_task["reference"]["name"]})")
    return necessary_trigger_names

class TriggerProcessor(coffea_processor.ProcessorABC):
    def __init__(self, trigger_tasks):
        self._trigger_tasks = trigger_tasks
        
    def process(self, events):
        # selection
        events["Muon_Good"] = events["Muon"][
            (events.Muon.pt > 10)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.pfRelIso04_all < 0.25)
            & (events.Muon.looseId)
            ]

        muon_name = "ScoutingMuonNoVtx"
        events[muon_name+"_Good"] = events[muon_name][
            (events[muon_name].pt > 10)
            & (abs(events[muon_name].eta) < 2.4)
            & (abs(events[muon_name].trk_dxy) < 0.2)
            & (abs(events[muon_name].trackIso) < 0.15)
            & (abs(events[muon_name].trk_dz) < 0.5)
            #& (events[muon_name, "type"] == 2)
            & (events[muon_name].normchi2 < 10)
            & (events[muon_name].nValidRecoMuonHits > 0)
            & (events[muon_name].nRecoMuonMatchedStations > 1)
            & (events[muon_name].nValidPixelHits > 0)
            & (events[muon_name].nTrackerLayersWithMeasurement > 5)
        ]
        
        jet_types = ["Jet", "JetCHS", "PFJet", "ScoutingPFJet", "ScoutingPFJetRecluster", "ScoutingCHSJetRecluster"]
        fatjet_types = ["FatJet", "ScoutingFatCHSJetRecluster"]
        
        # undo correction
        events["Jet", "pt"] = (1-events["Jet", "rawFactor"])*events["Jet", "pt"]
        events["JetCHS", "pt"] = (1-events["JetCHS", "rawFactor"])*events["JetCHS", "pt"]
        events["FatJet", "pt"] = (1-events["FatJet", "rawFactor"])*events["FatJet", "pt"]
        
        for jet_type in jet_types:
            if jet_type in ["Jet", "PFJet", "JetCHS"]:
                good_muon_name = "Muon_Good"
            else:
                good_muon_name = "ScoutingMuonNoVtx_Good"
            events[jet_type+"_Good"] = events[jet_type][
                (events[jet_type].pt > 30)
                & (abs(events[jet_type].eta) < 2.5)
                & ak.all(events[jet_type].metric_table(events[good_muon_name]) > 0.4, axis=-1) # default metric: delta_r
            ]
            if jet_type in ["Jet", "JetCHS"]:
                jetid_mask = (events[jet_type+"_Good"].jetId >= 4)
                events[jet_type+"_Good"] = events[jet_type+"_Good"][jetid_mask]
            else:
                jetid_mask = get_jetid_mask(events[jet_type+"_Good"], jet_type="CHS", use_lepton_veto=True)
                events[jet_type+"_Good"] = events[jet_type+"_Good"][jetid_mask]
            
        for fatjet_type in fatjet_types:
            if fatjet_type == ["FatJet"]:
                good_muon_name = "Muon_Good"
            else:
                good_muon_name = "ScoutingMuonNoVtx_Good"
            events[fatjet_type+"_Good"] = events[fatjet_type][
                (events[fatjet_type].pt > 30)
                & (abs(events[fatjet_type].eta) < 2.5)
                & ak.all(events[fatjet_type].metric_table(events[good_muon_name]) > 0.8, axis=-1)  # default metric: delta_r
            ]

        # add HT
        for jet_type in jet_types+fatjet_types:
            events[jet_type+"_Good", "HT"] = ak.sum(events[jet_type+"_Good", "pt"], axis=-1)
        
        # define histograms
        dataset_axis = hist.axis.StrCategory([], name="dataset", label="Dataset", growth=True)
        jet_type_axis = hist.axis.StrCategory(jet_types+fatjet_types, name="jet_type", 
                                              label="Type of Jet", growth=False)
        trigger_axis = hist.axis.StrCategory(get_necessary_trigger_names(self._trigger_tasks), 
                                             name="trigger", label="Trigger name", growth=False)
        pt_axis = hist.axis.Variable(pt_bins, name="pt", label="jet $p_T$")
        ht_axis = hist.axis.Variable(ht_bins, name="ht", label="Event HT")
        
        h_jet = hda.hist.Hist(dataset_axis, jet_type_axis, trigger_axis, pt_axis, ht_axis)

        # fill histograms
        dataset = events.metadata.get("dataset", "untitled")
        
        filled_trigger_names = set()
        for trigger_task in self._trigger_tasks:
            signal_trigger_name = trigger_task["signal"]["name"]
            reference_trigger_name = trigger_task["reference"]["name"]
            intersect_trigger_name = f"({signal_trigger_name}) AND ({reference_trigger_name})"

            signal_trigger_mask = get_trigger_mask(events, signal_trigger_name)
            reference_trigger_mask = get_trigger_mask(events, reference_trigger_name)
            intersect_trigger_mask = signal_trigger_mask & reference_trigger_mask

            if signal_trigger_name not in filled_trigger_names:
                #print(signal_trigger_name)
                for jet_type in jet_types+fatjet_types:     
                    h_jet.fill(dataset=dataset, jet_type=jet_type,
                               trigger=signal_trigger_name,
                               pt=ak.flatten(events[signal_trigger_mask][jet_type+"_Good"][:, :1]["pt"]),
                               ht=ak.flatten(events[signal_trigger_mask][jet_type+"_Good"][:, :1]["HT"]),
                              )
            filled_trigger_names.add(signal_trigger_name)
            
            if reference_trigger_name not in filled_trigger_names:
                #print(reference_trigger_name)
                for jet_type in jet_types+fatjet_types:
                    h_jet.fill(dataset=dataset, jet_type=jet_type,
                               trigger=reference_trigger_name,
                               pt=ak.flatten(events[reference_trigger_mask][jet_type+"_Good"][:, :1]["pt"]),
                               ht=ak.flatten(events[reference_trigger_mask][jet_type+"_Good"][:, :1]["HT"]),
                              )
            filled_trigger_names.add(reference_trigger_name)
            
            if intersect_trigger_name not in filled_trigger_names:
                #print(intersect_trigger_name)
                for jet_type in jet_types+fatjet_types:
                    #print(jet_type)
                    h_jet.fill(dataset=dataset, jet_type=jet_type,
                               trigger=intersect_trigger_name,
                               pt=ak.flatten(events[intersect_trigger_mask][jet_type+"_Good"][:, :1]["pt"]),
                               ht=ak.flatten(events[intersect_trigger_mask][jet_type+"_Good"][:, :1]["HT"]),
                              )
            filled_trigger_names.add(intersect_trigger_name)
        
        return {
            "h_jet": h_jet,
        }
        
    def postprocess(self, accumulator):
        return accumulator