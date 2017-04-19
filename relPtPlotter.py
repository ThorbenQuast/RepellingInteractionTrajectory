import numpy as np
import ROOT
import argparse

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser()
parser.add_argument("dataFile", type=str, help="")
parser.add_argument("outputFile", type=str, help="")
args = parser.parse_args()


column_names = ["DeltaR_initial","E_e1_in","E_e1_out","E_e2_in","E_e2_out","E_q1_final","E_q1_in","E_q1_out","E_q2_final","E_q2_in","E_q2_out","IP1_x","IP1_y","IP1_z","IP2_x","IP2_y","IP2_z","Q2_1","Q2_2","event_weight","mass_q1","mass_q2","px_e1_in","px_e1_out","px_e2_in","px_e2_out","px_q1_final","px_q1_in","px_q1_out","px_q2_final","px_q2_in","px_q2_out","py_e1_in","py_e1_out","py_e2_in","py_e2_out","py_q1_final","py_q1_in","py_q1_out","py_q2_final","py_q2_in","py_q2_out","pz_e1_in","pz_e1_out","pz_e2_in","pz_e2_out","pz_q1_final","pz_q1_in","pz_q1_out","pz_q2_final","pz_q2_in","pz_q2_out","relPt_final","relPt_initial"]

#get the histogram in the repulsive case:
data_file = open(args.dataFile)

data = np.genfromtxt(data_file, delimiter=",", skip_header=1, names = column_names)

max_x = max(data["relPt_final"]-data["relPt_initial"])
hist_relPtGain = ROOT.TH1F("gain_pT", "gain_pT", 100, 0.0, max_x)

for i in range(len(data["relPt_final"])):
  hist_relPtGain.Fill(data["relPt_final"][i]-data["relPt_initial"][i], data["event_weight"][i])

canvas = ROOT.TCanvas("canvas", "canvas", 800, 500)
hist_relPtGain.Scale(1./hist_relPtGain.Integral())
hist_relPtGain.GetXaxis().SetTitle("#Delta p_{t,rel}")
hist_relPtGain.GetYaxis().SetTitle("N [a.u.]")
hist_relPtGain.Draw("HIST")
canvas.Print(args.outputFile)

