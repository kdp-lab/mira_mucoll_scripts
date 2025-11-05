from array import array
from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1D, TMath, TFile, TLorentzVector, TEfficiency
from math import *
from optparse import OptionParser
from itertools import combinations
import glob,os

occ_conversion = True

def getBin(detector, side, layer):

    """
    docstring
    """
    weight = 1.

    layer_area = [270.40, 270.40, 448.50, 448.50, 655.20, 655.20, 904.80, 904.80,                                    # VXD barrel
                  389.00, 389.00, 378.96, 378.96, 364.36, 364.36, 312.48, 312.48,                                    # VXD endcaps
                  389.00, 389.00, 378.96, 378.96, 364.36, 364.36, 312.48, 312.48, 8117.85, 22034.16, 51678.81,       # IT barrel
                  6639.65, 10611.59, 10078.04, 9900.19, 9307.37, 8595.98, 8299.56,                                   # IT endcaps
                  6639.65, 10611.59, 10078.04, 9900.19, 9307.37, 8595.98, 8299.56, 140032.91, 194828.39, 249623.88,  # OT barrel
                  69545.45, 69545.45, 69545.45, 69545.45,  # OT endcaps
                  69545.45, 69545.45, 69545.45, 69545.45]

    bin_n = 0

    if detector == 1:
        bin_n = layer
    elif detector == 2:
        if side > 0:
            bin_n = layer+8
        else:
            bin_n = layer+16
    elif detector == 3:
        bin_n = layer+24
    elif detector == 4:
        if side > 0:
            bin_n = layer+27
        else:
            bin_n = layer+34
    elif detector == 5:
        bin_n = layer+41
    elif detector == 6:
        if side > 0:
            bin_n = layer+44
        else:
            bin_n = layer+48

    weight = 1./layer_area[bin_n]

    return bin_n, weight


#########################
# parameters
parser = OptionParser()
parser.add_option('-i', '--inFiles', help='--inFile Output_REC.slcio',
                  type=str,default='*.slcio')
parser.add_option('-o', '--outFile', help='--outFile histos_occupancy.root',
                  type=str, default='histos_occupancy.root')
(options, args) = parser.parse_args()


vb_cut = 0.26
ib_cut = 0.50
ob_cut = 0.60

if ',' in options.inFiles:
    in_files = [f.strip() for f in options.inFiles.split(',')]
else:
    in_files = sorted(glob.glob(options.inFiles))
if len(in_files) == 0:
    raise RuntimeError(f'No files matched "{options.inFiles}"')
#########################
max_files = 100
h_nhits_nowei = TH1D('h_nhits_nowei', 'h_nhits_nowei', 52, 0., 52)
h_nhits = TH1D('h_nhits', 'h_nhits', 52, 0., 52)
h_ntimehits = TH1D('h_ntimehits', 'h_ntimehits', 52, 0., 52)
in_files = in_files[:max_files]
print(len(in_files))

#########################

# create a reader and open an LCIO file
# reader = IOIMPL.LCFactory.getInstance().createLCReader()
# reader.open(options.inFile)

allCollections = [
    "OuterTrackerEndcapCollection",
    "OuterTrackerBarrelCollection",
    "InnerTrackerEndcapCollection",
    "InnerTrackerBarrelCollection",
    "VertexEndcapCollection",
    "VertexBarrelCollection"]

timeCollections = [
    "OTBarrelHits",
    "OTEndcapHits",
    "ITBarrelHits",
    "ITEndcapHits",
    "VXDBarrelHits",
    "VXDEndcapHits"]

pixel_size = {1: 0.0625,
              2: 0.0625,
              3: 0.005,
              4: 0.005,
              5: 0.005,
              6: 0.005}

totEv = 0

# loop over all events in the file

for fname in in_files:
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(fname)
    print(fname)

    for ievt, event in enumerate(reader):

        totEv = totEv + 1

        if ievt % 10 == 0:
            print("Processing " + str(ievt))

        # do all first
        for coll in allCollections:

            #print("Collection: ", coll)

            # Check if collection exists by looking at the collection names
            collection_names = event.getCollectionNames()
            if coll in collection_names:
                hitsCollection = event.getCollection(coll)
                encoding = hitsCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                decoder = UTIL.BitField64(encoding)

                # looping over hit collections
                for ihit, hit in enumerate(hitsCollection):
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    layer = decoder['layer'].value()
                    detector = decoder["system"].value()
                    side = decoder["side"].value()

                    bin, wei = getBin(detector, side, layer)
                    h_nhits.Fill(bin, wei)
                    h_nhits_nowei.Fill(bin)

        # now time
        for coll in timeCollections:

            # Check if the collection exists by looking at the collection names
            collection_names = event.getCollectionNames()
            if coll in collection_names:
                hitsCollection = event.getCollection(coll)
                encoding = hitsCollection.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                decoder = UTIL.BitField64(encoding)

                # looping over hit collections
                for ihit, hit in enumerate(hitsCollection):
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    layer = decoder['layer'].value()
                    detector = decoder["system"].value()
                    side = decoder["side"].value()

                    bin, wei = getBin(detector, side, layer)
                    h_ntimehits.Fill(bin, wei)
    reader.close()


h_nhits.Scale(1./totEv)
h_ntimehits.Scale(1./totEv)

##################
# write histograms
output_file = TFile(options.outFile, 'RECREATE')
h_nhits.Write()
h_ntimehits.Write()
h_nhits_nowei.Write()
output_file.Close()