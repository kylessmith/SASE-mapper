from collections import OrderedDict
import numpy as np
from math import pi
    

def plot(out_long_fn, genome, prefix, plot_type="matplotlib"):
    '''
    plot pvalues in bokeh or matplotlib
    '''
    
    plot_array_x = []
    plot_array_y = []
    rect_x = []
    rect_w = []
    text_x = []
    chrom_names = []
    continuous_chrom_pos = OrderedDict()
    chrom_count = 0
    previous_end = 0
    total_chroms = genome.keys()
    
    for i in xrange(len(total_chroms)):
        chrom = total_chroms[i]
        continuous_chrom_pos[chrom] = chrom_count
        chrom_count += genome[chrom]
        w = genome[chrom] / 2.0
        text_x.append(continuous_chrom_pos[chrom] + w)
        if i%2 == 0:
            rect_x.append(continuous_chrom_pos[chrom] + w)
            rect_w.append(genome[chrom])
            
        
    for line in open(out_long_fn, 'r'):
        chrom, start, end, pval = line.strip().split('\t')
        pval = max(float(pval), 1e-200)

        if int(start)+continuous_chrom_pos[chrom] != previous_end:
            plot_array_y.append(0.0)
            plot_array_y.append(0.0)
            plot_array_x.append(int(previous_end)+1)
            plot_array_x.append(int(start)+continuous_chrom_pos[chrom]-1)
            
        plot_array_y.append(-np.log10(float(pval)))
        plot_array_y.append(-np.log10(float(pval)))
        plot_array_x.append(int(start)+continuous_chrom_pos[chrom])
        plot_array_x.append(int(end)+continuous_chrom_pos[chrom])
        previous_end = int(end)+continuous_chrom_pos[chrom]
        
    if plot_type == "bokeh":
        bokeh_plot(plot_array_x, plot_array_y, rect_x, rect_w, total_chroms, prefix, text_x)
    elif plot_type == "matplotlib":
        matplotlib_plot(plot_array_x, plot_array_y, continuous_chrom_pos, rect_w, total_chroms, prefix, text_x, genome)


def bokeh_plot(plot_array_x, plot_array_y, rect_x, rect_w, total_chroms, prefix, text_x):
    '''
    Create peaks plot in bokeh
    '''
    
    from bokeh.plotting import figure, output_file, show
    from bokeh.models import FixedTicker

    p = figure(title="Signatures of Accelerated Somatic Evolution", plot_width=700, plot_height=700)
    p.line(x=plot_array_x, y=plot_array_y)
    
    h = max(plot_array_y) + np.std(plot_array_y)
    rect_y = list(np.repeat(h/2.0, len(rect_x)))
    p.rect(x=rect_x, y=rect_y, width=rect_w, height=h,
            fill_color='grey', fill_alpha=0.3,
            line_color='grey', line_alpha=0.3)
    p.text(x=text_x, y = -0.05,
            text=total_chroms, angle=pi/2,
            text_baseline='middle', text_align='right',
            text_font_size='10pt', text_font='arial')
    p.xaxis[0].ticker = FixedTicker(ticks=rect_x)
    p.xaxis.axis_label='Chromosomes'
    p.xaxis.major_label_orientation = pi/4
    p.yaxis.axis_label = '-log10 p-values'
    
    output_file(prefix+'.html')
    show(p)
    
    
def matplotlib_plot(plot_array_x, plot_array_y, continuous_chrom_pos, rect_w, total_chroms, prefix, text_x, genome):
    '''
    Create peaks plot in matplotlib
    '''
    
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    
    chrom_starts = continuous_chrom_pos.values()
    chrom_names = continuous_chrom_pos.keys()
    rect_x = [chrom_starts[x] for x in xrange(len(chrom_starts)) if x%2 == 0]
    
    h = max(plot_array_y) + np.std(plot_array_y)
        
    fig1 = plt.figure(figsize=(11,5))
    ax1 = fig1.add_subplot(111)
    
    ax1.plot(plot_array_x, plot_array_y)
    
    for i in xrange(len(rect_x)):
        ax1.add_patch(patches.Rectangle((rect_x[i], 0), rect_w[i], h, facecolor="grey", alpha=0.3, edgecolor="none"))
        
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off'          # ticks along the top edge are off
        )
        
    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        right='off'      # ticks along the bottom edge are off
        )
        
    plt.xticks(text_x, chrom_names, rotation=90)
    plt.xlabel("Chromosomes")
    plt.title("Signatures of Accelerated Somatic Evolution")
    plt.ylabel("-log10 p-values")
    plt.xlim((0, chrom_starts[-1] + genome[chrom_names[-1]]))
    plt.tight_layout()
    
    fig1.savefig(prefix+".png")
    
