import matplotlib
import numpy
import Bio
import pandas
import logomaker 
from pysanger import *
from Bio import SeqIO

def visualizar_seq(archivo_ab1, template , region="all", archivo_salida="output.pdf"):
  abidata = abi_to_dict(archivo_ab1)
  fseq, rseq = generate_consensusseq(abidata)
  fig = visualize(abidata, template=template, region=region)
  return fseq, rseq
