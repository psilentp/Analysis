# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>
from pylab import *
import psummary_chemocircuit as psum
groups = {
    'AVA_AVB_l1':[26,27,28,29,30,31,32,33,34,35,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59],
    'AVB_AVA_l1':[60,61,62,63,64,65,66,68,69,70,71,72,73,74,75],
    'AVB_AVA_ttx':[96,97,98,99,100,103,104,105,106,107],
    'AVB_AVA':[60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,96,97,98,99,100,103,104,105,106,107],
    'ASH_AVB_l1':[88,89,90,92,93,94,95]
}
# <codecell>


def hide_spines():
    """Hides the top and rightmost axis spines from view for all active
        figures and their respective axes."""
    
    # Retrieve a list of all current figures.
    figures = [x for x in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
    for figure in figures:
        # Get all Axis instances related to the figure.
        for ax in figure.canvas.figure.get_axes():
            # Disable spines.
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            # Disable ticks.
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

def rpot_conductance(key):
    groups = {
    'AVA_AVB_l1':[26,27,28,29,30,31,32,33,34,35,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59],
    'AVB_AVA_l1':[60,61,62,63,64,65,66,68,69,70,71,72,73,74,75],
    'AVB_AVA_ttx':[96,97,98,99,100,103,104,105,106,107],
    'AVB_AVA':[60,61,62,63,64,65,66,68,69,70,71,72,73,74,75,96,97,98,99,100,103,104,105,106,107],
    'ASH_AVB_l1':[88,89,90,92,93,94,95]
    }
    fi = open('output.txt','r')
    lns = fi.readlines()
    calcdict = dict()

    for l in lns:
        data = [float(x) for x in l.split()]
        calcdict.update({int(data[0]):{'rpot':data[1],'cond':data[2],'std_err':data[3],'p_val':data[4]}})
    fig = figure()
    group = groups[key]
    xs = [calcdict[cennum]['rpot'] for cennum in group]
    ys = [calcdict[cennum]['cond'] for cennum in group]
    yerr = [calcdict[cennum]['std_err'] for cennum in group]
    errorbar(xs,ys,yerr=yerr,marker='o',ls='None',color = 'k')
    ax = gca()
    psum.format_iv(ax)
    return calcdict

calcdict = rpot_conductance('AVA_AVB_l1')
ax = gca()
ax.set_xbound([-100,100])

# <codecell>
glist = []
for key in groups.keys():
    group = groups[key]
    xs = [calcdict[cennum]['rpot'] for cennum in group]
    ys = [calcdict[cennum]['cond'] for cennum in group]
    yerr = [calcdict[cennum]['std_err'] for cennum in group]
    glist.append(ys)
    from scipy import stats
    print key + ' \tp:%s'%(stats.wilcoxon(ys)[1]) + '\tn=%s'%(len(ys))

fig = figure(figsize = (15,7))
boxplot(glist)
ax = fig.axes[0]
ax.set_xticklabels(groups.keys())
for x in range(len(glist)):
    plot(ones(len(glist[x]))*(x+1),
         glist[x],'o', 
         markersize = 4,
         markerfacecolor = (0,0,0,1.0), 
         markeredgecolor = 'k',
         alpha = 0.5)
    print len(glist[x])
hide_spines()

# <codecell>



# <codecell>


# <codecell>



# <codecell>



# <markdowncell>

# calculate the hours spent recording CEN's 26-72

# <codecell>

f = open('workflow.csv','r')
import csv
worklog = list()
for l in csv.reader(f):
    worklog.append({'starttime':l[1],'endtime':l[2],'task':l[3],'duration':l[4]})
worklog = worklog[1:]

# <codecell>

ephys = filter(lambda i:'ephys' in i['task'].lower(),worklog[429:])

# <codecell>

sum(map(lambda item:float(item['duration']),ephys))/(293-249)

# <markdowncell>

# 30 cells in a group, 10hrs/day, 2.95 hrs/cell = 9 days per group.

# <codecell>

30/(10/3.0)

# <codecell>


