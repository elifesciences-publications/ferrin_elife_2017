'''Custom matplotlib functions'''

fontsize = {'paper': 8, 'presentation': 12}

def clean_axis(ax, spineoffset=['left', 'bottom']):
    '''Clean current matplotlib axis
    * Remove top and right lines
    * Offset left and bottom spines by 5pt
    ax: handle to axis
    keyword arguments:
    spineoffset: Either ['left', 'bottom'} (default) or 'left'
    '''
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    if type(spineoffset) != list:
        spineoffset = list(spineoffset)
    for loc, spine in ax.spines.items():
        if loc in spineoffset:
            spine.set_position(('outward', 5))
        if loc in ['top', 'right']:
            spine.set_color('none')


def set_fonts(mpl, format='paper'):
    '''Takes handle to matplotlib as input and sets params'''
    mpl.rcParams['font.sans-serif'] = 'Arial'
    mpl.rcParams['figure.subplot.hspace'] = 0.5
    mpl.rcParams['figure.subplot.wspace'] = 0.5
    mpl.rcParams['font.size'] = fontsize[format]
    mpl.rcParams['xtick.labelsize'] = fontsize[format]
    mpl.rcParams['ytick.labelsize'] = fontsize[format]
    mpl.rcParams['axes.labelsize'] = fontsize[format]
    mpl.rcParams['axes.titlesize'] = fontsize[format]
    mpl.rcParams['legend.fontsize'] = fontsize[format]
