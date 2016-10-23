import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
import sys
from numpy.random import uniform
import math

dcs3d = "energy,angle,dcs\n"
en_tcs = "energy,tcs\n"
ave_angle = "energy,angle,dE\n"

colors  = ['bo-','go-','ro-','mo-','co-']
colors += ['b*-','g*-','r*-','m*-','c*-']
colors += ['bx-','gx-','rx-','mx-','cx-']

def plot_ave_angle(data):
  plt.figure()
  ax = plt.subplot(111)
  i = 0; savename = "./plots/ave_angle_";
  for filename in data.keys():
    x = []; y = [];
    for d in data[filename][1:]:
      now = d.split(',')
      x.append(now[0]); y.append(now[1]);
    (proj, targ) = get_proj_targ(filename)
    if i == len(data.keys())-1:
      savename += proj+targ+".png"
    else:
      savename += proj+targ+"-"
    print("projectile: -%s-\ttarget: -%s-" % (proj, targ))
    ax.plot(x, y, colors[i], label="%s-%s" % (proj, targ))
    i += 1

  plt.xlabel("Energy [ev]")
  plt.ylabel("Average CM Scattering Angle [deg]")
  box = ax.get_position()
  ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

  # Put a legend to the right of the current axis
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  plt.savefig(savename)

def plot_ave_dE(data):
  plt.figure()
  ax = plt.subplot(111)
  i = 0; savename = "./plots/ave_dE_";
  for filename in data.keys():
    x = []; y = [];
    for d in data[filename][1:]:
      now = d.split(',')
      x.append(now[0]); y.append(now[2]);
    (proj, targ) = get_proj_targ(filename)
    if i == len(data.keys())-1:
      savename += proj+targ+".png"
    else:
      savename += proj+targ+"-"
    print("projectile: -%s-\ttarget: -%s-" % (proj, targ))
    ax.plot(x, y, colors[i], label="%s-%s" % (proj, targ))
    i += 1

  plt.xlabel("Energy [ev]")
  plt.ylabel("Average Energy Loss per Collision [eV]")
  box = ax.get_position()
  ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

  # Put a legend to the right of the current axis
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  plt.savefig(savename)

def plot_tcs(data):
  plt.figure()
  ax = plt.subplot(111)
  i = 0; savename = "./plots/tcs_";
  for filename in data.keys():
    x = []; y = [];
    for d in data[filename][1:]:
      now = d.split(',')
      x.append(now[0]); y.append(now[1]);
    (proj, targ) = get_proj_targ(filename)
    if i == len(data.keys())-1:
      savename += proj+targ+".png"
    else:
      savename += proj+targ+"-"
    print("projectile: -%s-\ttarget: -%s-" % (proj, targ))
    ax.plot(x, y, colors[i], label="%s-%s" % (proj, targ))
    i += 1

  plt.xlabel("Energy [ev]")
  plt.ylabel("Total Cross Section [a0^2]")
  box = ax.get_position()
  ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

  # Put a legend to the right of the current axis
  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  plt.savefig(savename)


def plot_dcs3d_heat(data):

  for filename in data.keys():
    (proj, targ) = get_proj_targ(filename)
    print("projectile: %s\ttarget: %s" % (proj, targ))
    savename = "./plots/dcs3d_heat_"+proj+targ+".png"
    X_dat = []; Y_dat = []; Z_dat = [];

    for d in data[filename][1:]:
      line = d.split(',')
      X_dat.append(float(line[1])); Y_dat.append(float(line[0])); Z_dat.append(math.log10(float(line[2])));

    x = np.array(X_dat,); y = np.array(Y_dat,); z = np.array(Z_dat,);

    xi = np.linspace(min(x), max(x), 100)
    yi = np.linspace(min(y), max(y), 100)
    zi = griddata(x, y, z, xi, yi, interp='linear')
    CS = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow,vmax=abs(zi).max(), vmin=-abs(zi).max())
    plt.xlabel('CM Scattering Angle [deg]', fontsize=16)
    plt.ylabel('CM Collision Energy [eV]', fontsize=16)
    plt.title('%s - %s Collision' % (proj, targ))
    plt.colorbar()  # draw colorbar
    plt.savefig(savename)

def plot_dcs3d(data):
  def get_xy(list_xy):
    x = []; y = [];
    for l in list_xy:
      x.append(l[0]); y.append(l[1]);
    return (x,y)

  plt.figure()
  ax = plt.subplot(111)
  i = 0; savename = "./plots/dcs_";
  for filename in data.keys():
    (proj, targ) = get_proj_targ(filename)
    print("projectile: %s\ttarget: %s" % (proj, targ))
    if i == len(data.keys())-1:
      savename += proj+targ+".png"
    else:
      savename += proj+targ+"-"
    cleaned_data = []
    cleaned_dict = {}
    for d in data[filename][1:]:
      st = d.split(',')
      tt = []
      for s in st:
        tt.append(float(s.strip()))
      if tt[0] in cleaned_dict.keys():
        cleaned_dict[tt[0]].append((tt[1], tt[2]))
      else:
        cleaned_dict[tt[0]] = []
      cleaned_data.append(tt)
    en = cleaned_dict.keys()
    max_en = max(en); (x3,y3) = get_xy(cleaned_dict[max_en]);
    min_en = min(en); (x1,y1) = get_xy(cleaned_dict[min_en]);
    med_en = sorted(en)[len(en)//2]; (x2,y2) = get_xy(cleaned_dict[med_en]);
    plt.plot(x1, y1, colors[3*i+0], label=proj+"-"+targ+": "+str(min_en)+" eV")
    plt.plot(x2, y2, colors[3*i+1], label=proj+"-"+targ+": "+str(med_en)+" eV")
    plt.plot(x3, y3, colors[3*i+2], label=proj+"-"+targ+": "+str(max_en)+" eV")
    i+=1
  plt.yscale('log')
  plt.xlabel("CM Theta [deg]")
  plt.ylabel("CM Differential Cross Section [a0^2]")
  plt.legend()
  plt.savefig(savename)

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
    # plot_dcs3d(data)
    plot_dcs3d_heat(data)
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
