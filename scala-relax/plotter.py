import matplotlib.pyplot as plt
import numpy as np
import sys

dcs3d = "energy,angle,dcs\n"
en_tcs = "energy,tcs\n"
ave_angle = "energy,angle,dE\n"


def plot_ave_angle(data):
  x = []; y = [];
  for d in data:
    now = d.split(',')
    x.append(now[0]); y.append(now[1]);
  plt.figure()
  plt.plot(x, y, 'bo')
  plt.xlabel("Energy [ev]")
  plt.ylabel("Average CM Scattering Angle [deg]")
  plt.savefig("./plots/ave_angle.png")

def plot_ave_dE(data):
  x = []; y = [];
  for d in data:
    now = d.split(',')
    x.append(now[0]); y.append(now[2]);
  plt.figure()
  plt.plot(x, y, 'ro')
  plt.xlabel("Energy [ev]")
  plt.ylabel("Average Energy Loss per Collision [eV]")
  plt.savefig("./plots/ave_dE.png")

def plot_tcs(data):
  x = []; y = [];
  for d in data:
    now = d.split(',')
    x.append(now[0]); y.append(now[1]);
  plt.figure()
  plt.plot(x, y, 'ro')
  plt.xlabel("Energy [ev]")
  plt.ylabel("Total Cross Section [a0^2]")
  plt.savefig("./plots/tcs.png")

def get_xy(list_xy):
  x = []; y = [];
  for l in list_xy:
    x.append(l[0])
    y.append(l[1])
  return (x,y)

def plot_dcs3d(data):
  cleaned_data = []
  cleaned_dict = {}
  for d in data:
    st = d.split(',')
    tt = []
    for s in st:
      tt.append(float(s.strip()))
    if tt[0] in cleaned_dict.keys():
      cleaned_dict[tt[0]].append((tt[1],tt[2]))
    else:
      cleaned_dict[tt[0]] = []
    cleaned_data.append(tt)
  en = cleaned_dict.keys()
  max_en = max(en); (x3,y3) = get_xy(cleaned_dict[max_en]);
  min_en = min(en); (x1,y1) = get_xy(cleaned_dict[min_en]);
  med_en = sorted(en)[len(en)//2]; (x2,y2) = get_xy(cleaned_dict[med_en]);
  plt.plot(x1, y1, 'bo', label=str(min_en)+" eV")
  plt.plot(x2, y2, 'ko', label=str(med_en)+" eV")
  plt.plot(x3, y3, 'ro', label=str(max_en)+" eV")
  plt.yscale('log')
  plt.xlabel("CM Theta [deg]")
  plt.ylabel("CM Differential Cross Section [a0^2]")
  plt.legend()
  plt.savefig("./plots/dcs.png")

def file_type(data):
  header = data[0]
  if header == dcs3d:
    print("plotting dcs3d")
    plot_dcs3d(data[1:])
  elif header == ave_angle:
    plot_ave_angle(data[1:])
    plot_ave_dE(data[1:])
  elif header == en_tcs:
    plot_tcs(data[1:])
  else:
    print("header -"+header+"- not recognized")

def read_file(filename):
  with open(filename) as f:
    data = f.readlines()
  return data

def plot_file(filename):
  data = read_file(filename)
  file_type(data)


def main():
  if len(sys.argv) != 2:
    print("!! ERROR\n!! Pass filename to plot\n")
    print("$ python plotter.py <filename>\n")
  else:
    filename = sys.argv[1]
    print("plotting "+filename)
    plot_file(filename)

if __name__ == "__main__":
  main()