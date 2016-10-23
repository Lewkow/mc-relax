import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

dcs3d = "energy,angle,dcs\n"
en_tcs = "energy,tcs\n"
ave_angle = "energy,angle,dE\n"

colors  = ['bo-','go-','ro-','mo-','co-']
colors += ['b*-','g*-','r*-','m*-','c*-']

def plot_ave_angle(data):
  plt.figure()
  ax = plt.subplot(111)
  i = 0
  for filename in data.keys():
    x = []; y = [];
    for d in data[filename][1:]:
      now = d.split(',')
      x.append(now[0]); y.append(now[1]);
    (proj, targ) = get_proj_targ(filename)
    print("projectile: -%s-\ttarget: -%s-" % (proj, targ))
    ax.plot(x, y, colors[i], label="%s-%s" % (proj, targ))
    i += 1

  plt.xlabel("Energy [ev]")
  plt.ylabel("Average CM Scattering Angle [deg]")
  box = ax.get_position()
  ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

  # Put a legend to the right of the current axis
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  plt.savefig("./plots/ave_angle.png")

def plot_ave_dE(data):
  plt.figure()
  ax = plt.subplot(111)
  i = 0
  for filename in data.keys():
    x = []; y = [];
    for d in data[filename][1:]:
      now = d.split(',')
      x.append(now[0]); y.append(now[2]);
    (proj, targ) = get_proj_targ(filename)
    print("projectile: -%s-\ttarget: -%s-" % (proj, targ))
    ax.plot(x, y, colors[i], label="%s-%s" % (proj, targ))
    i += 1

  plt.xlabel("Energy [ev]")
  plt.ylabel("Average Energy Loss per Collision [eV]")
  box = ax.get_position()
  ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

  # Put a legend to the right of the current axis
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  plt.savefig("./plots/ave_dE.png")

def plot_tcs(data):
  plt.figure()
  ax = plt.subplot(111)
  i = 0
  for filename in data.keys():
    x = []; y = [];
    for d in data[filename][1:]:
      now = d.split(',')
      x.append(now[0]); y.append(now[1]);
    (proj, targ) = get_proj_targ(filename)
    print("projectile: -%s-\ttarget: -%s-" % (proj, targ))
    ax.plot(x, y, colors[i], label="%s-%s" % (proj, targ))
    i += 1

  plt.xlabel("Energy [ev]")
  plt.ylabel("Total Cross Section [a0^2]")
  box = ax.get_position()
  ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

  # Put a legend to the right of the current axis
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  plt.savefig("./plots/tcs.png")

def plot_dcs3d(data):
  def get_xy(list_xy):
    x = []; y = [];
    for l in list_xy:
      x.append(l[0]); y.append(l[1]);
    return (x,y)
  return None
  # (proj, targ) = get_proj_targ(filename)
  # print("projectile: %s\ttarget: %s" % (proj, targ))
  # cleaned_data = []
  # cleaned_dict = {}
  # for d in data:
  #   st = d.split(',')
  #   tt = []
  #   for s in st:
  #     tt.append(float(s.strip()))
  #   if tt[0] in cleaned_dict.keys():
  #     cleaned_dict[tt[0]].append((tt[1],tt[2]))
  #   else:
  #     cleaned_dict[tt[0]] = []
  #   cleaned_data.append(tt)
  # en = cleaned_dict.keys()
  # max_en = max(en); (x3,y3) = get_xy(cleaned_dict[max_en]);
  # min_en = min(en); (x1,y1) = get_xy(cleaned_dict[min_en]);
  # med_en = sorted(en)[len(en)//2]; (x2,y2) = get_xy(cleaned_dict[med_en]);
  # plt.plot(x1, y1, 'bo', label=str(min_en)+" eV")
  # plt.plot(x2, y2, 'ko', label=str(med_en)+" eV")
  # plt.plot(x3, y3, 'ro', label=str(max_en)+" eV")
  # plt.yscale('log')
  # plt.xlabel("CM Theta [deg]")
  # plt.ylabel("CM Differential Cross Section [a0^2]")
  # plt.legend()
  # plt.savefig("./plots/dcs.png")

def make_plot(filenames):
  data = {}
  for filename in filenames:
    data[filename] = read_file(filename)
  header = data[filenames[0]][0]
  for filename in filenames:
    if data[filename][0] != header:
      print("!! ERROR\n!! Files passed to plotter must be of same type\n")
  if header == dcs3d:
    print("plotting dcs3d")
    plot_dcs3d(data)
  elif header == ave_angle:
    print("plotting average scattering angle")
    print("plotting average energy loss per collision")
    plot_ave_angle(data)
    plot_ave_dE(data)
  elif header == en_tcs:
    print("plotting total cross section")
    plot_tcs(data)
  else:
    print("header -"+header+"- not recognized")

def get_proj_targ(filename):
  proj = filename.split('_')[-1].split('-')[0];
  targ = filename.split('_')[-1].split('-')[1].split('.')[0];
  return (proj, targ)

def read_file(filename):
  with open(filename) as f:
    data = f.readlines()
  return data

def main():
  if len(sys.argv) < 2:
    print("!! ERROR\n!! Pass filename to plot\n")
    print("$ python plotter.py <filename>\n")
  else:
    filenames = sys.argv[1:]
    print("plotting: ")
    for filename in filenames:
        print filename
    make_plot(filenames)

if __name__ == "__main__":
  main()
