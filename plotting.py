import numpy as np
import ROOT
import argparse

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("textFile_repulsion", type=str, help="")
parser.add_argument("textFile_attraction", type=str, help="")
parser.add_argument("outputFile", type=str, help="")
args = parser.parse_args()


#get the histogram in the repulsive case:
data_repulsion = open(args.textFile_repulsion)
deltaR, weight = np.genfromtxt(data_repulsion, unpack=1, delimiter=",", usecols=(0,1))

#max_x = 3.141/180.
max_x = max(deltaR)
#division through max_weight is technical only...the constant quotient vanished during the normalisation
max_weight = max(weight)

hist_repulsion = ROOT.TH1F("hist_repulsion", "hist_repulsion", 100, 0.0, max_x)
for i in range(len(deltaR)):
  hist_repulsion.Fill(deltaR[i], weight[i]/max_weight)


#get the histogram in the attractive case:
data_attraction = open(args.textFile_attraction)
deltaR, weight = np.genfromtxt(data_attraction, unpack=1, delimiter=",", usecols=(0,1))

hist_attraction = ROOT.TH1F("hist_attraction", "hist_attraction", 100, 0.0, max_x)
for i in range(len(deltaR)):
  hist_attraction.Fill(deltaR[i], weight[i]/max_weight)


canvas = ROOT.TCanvas("canvas", "canvas", 800, 500)
hist_repulsion.SetStats(False)
hist_attraction.SetStats(False)
hist_repulsion.Scale(1./hist_repulsion.Integral())
hist_attraction.Scale(1./hist_attraction.Integral())
hist_repulsion.SetLineColor(1)
hist_attraction.SetLineColor(2)
legend = ROOT.TLegend(0.7, 0.5, 0.95, 0.7)
legend.AddEntry(hist_repulsion)
legend.AddEntry(hist_attraction)
hist_repulsion.GetXaxis().SetTitle("deltaR [rad]")
hist_repulsion.GetYaxis().SetTitle("N")
hist_repulsion.Draw("HIST")
hist_attraction.Draw("HISTSAME")
legend.Draw()
canvas.Print(args.outputFile)