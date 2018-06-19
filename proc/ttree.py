import ctypes
import logging
from array import array
from ctypes import POINTER
from random import random

logging.basicConfig(level=logging.DEBUG)

logging.info('importing ROOT')
import ROOT
from ROOT import TFile, TTree, AddressOf
from ROOT.Math import LorentzVector

ROOT.gROOT.Reset()

logging.debug("start")

t = TTree( 't1', 'tree with histos' )
 
maxn = 10
n = array( 'i', [ 0 ] )
d = array( 'f', maxn*[ 0. ] )
t.Branch( 'mynum', n, 'mynum/I' )
t.Branch( 'myval', d, 'myval[mynum]/F' )

# it creates
#root [3] t1->Print()
#******************************************************************************
#*Tree    :t1        : tree with histos                                       *
#*Entries :       25 : Total =            2448 bytes  File  Size =       1619 *
#*        :          : Tree compression factor =   1.00                       *
#******************************************************************************
#*Br    0 :mynum     : mynum/I                                                *
#*Entries :       25 : Total  Size=        644 bytes  File Size  =        170 *
#*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
#*............................................................................*
#*Br    1 :myval     : myval[mynum]/F                                         *
#*Entries :       25 : Total  Size=       1523 bytes  File Size  =        958 *
#*Baskets :        1 : Basket Size=      32000 bytes  Compression=   1.00     *
#*............................................................................*

# LorentzVector in ttree
#NTuple.Branch(#Name, #CLASS, &NT_##Name);
LorentzVector_Class = LorentzVector('ROOT::Math::PxPyPzE4D<double>')
lorvec = LorentzVector_Class(2., 3., 0., 0.)
t.Branch("lorvec", lorvec)

# vector of LorVecs
#define VECTOR_OBJECTs_in_NTuple(NTuple, VECTOR_CLASS, Name)   VECTOR_CLASS NT_##Name; VECTOR_CLASS* pt_NT_##Name ; NTuple.Branch(#Name, #VECTOR_CLASS, &pt_NT_##Name);
#VECTOR_OBJECTs_in_NTuple(OUTNTUPLE, std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >>, gen_t_w1_final_p4s)

ROOT.gROOT.ProcessLine("typedef std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >> LorentzVectorS;")
lorvecS = ROOT.LorentzVectorS()
t.Branch("lorvecS", lorvecS)


for i in range(25):
   n[0] = min(i,maxn)
   for j in range(n[0]):
      d[j] = i*0.1+j
   lorvecS.clear()
   some_lorvector = LorentzVector_Class(random(), 33., 0., 0.)
   lorvecS.push_back(some_lorvector)
   some_lorvector = LorentzVector_Class(random(), random(), 0., 0.)
   lorvecS.push_back(some_lorvector)

   t.Fill()

outfile = TFile("outfile_ttree.root", "RECREATE")
outfile.Write()
t.Write()
outfile.Close()




logging.debug("end")

