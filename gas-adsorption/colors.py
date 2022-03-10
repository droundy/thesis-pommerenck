import matplotlib.pyplot as plt
import matplotlib, re
import numpy as np

from matplotlib.font_manager import FontProperties
small_font = FontProperties()
small_font.set_size('small')

kelly_colors = [# 'F2F3F4',
                '222222', 'F3C300', '875692', 'F38400', 'A1CAF1', 'BE0032', 'C2B280', '848482', '008856', 'E68FAC', '0067A5', 'F99379', '604E97', 'F6A600', 'B3446C', 'DCD300', '882D17', '8DB600', '654522', 'E25822', '2B3D26']
kelly_colors = ['#'+c for c in kelly_colors]

prop_cycle = plt.rcParams['axes.prop_cycle']
print(list(prop_cycle.by_key()['color']))
print(kelly_colors)
colors = kelly_colors # list(prop_cycle.by_key()['color'])

the_order = [
    'Ni2m-dobdc',
    'MOF519',
    'NJUBai43',
    'NU1000',
    'MOF5',
    'HKUST1',
    'HKUST1-mono',
    'PCN14',
    'AX21',
    'Mg2dobdc',
    'Ni2dobdc',
    'Co2dobdc',
    'NU125',
    'NOTT112',
    'rhtMOF7',
    'CuMOF74',
    'PCN250',
    'NU1000',
    'UiO67',
    'UiO68Ant',
    'CYCU3Al',
    'Zn2bdc2dabco',
    'COF102',
]

# cm = plt.get_cmap('viridis')
if len(colors) < len(the_order):
    for x in np.linspace(0, 1, len(the_order) - len(colors)):
        colors.append(matplotlib.colors.hsv_to_rgb((x,1,1-x*0.25)))
# colors = [cm(1.*i/(len(the_order)-1)) for i in range(len(the_order))]

def order(mofs):
    return [x for x in the_order if x in mofs]+[x for x in mofs if x not in the_order]

def color(mof):
    if mof in the_order:
        return colors[the_order.index(mof)]
    return 'k'

def symbol(gas):
    if gas == 'methane':
        return '.'
    else:
        return '+' 

latex_mof_names = {
    'Ni2m-dobdc': 'Ni$_2$($m$-dodbc)',
    'Ni2dobdc': 'Ni$_2$(dodbc)',
    'Co2dobdc': 'Co$_2$(dodbc)',
    'Mg2dobdc':'Mg$_2$(dodbc)',
    'Zn2bdc2dabco': 'Zn$_2$(bdc)$_2$(dabco)',
    'CuMOF74': 'Cu-MOF-74',
    'CYCU3Al': 'CYCU-3-Al',
    'UiO68Ant': 'UiO-68-Ant',
    'rhtMOF7': 'rht-MOF-7',
    'HKUST1-mono': 'HKUST-1-mono',
    'NJUBai43': 'NJU-Bai-43',
}
def latex_mof(mof):
    if mof in latex_mof_names:
        return latex_mof_names[mof]
    print(mof, re.split('(\d+)', mof))
    letters, numbers,empty = re.split('(\d+)', mof)
    if empty == '':
        return letters + '-' + numbers
    return mof

#print(latex_mof('AX231'))
#print(latex_mof('COF102'))
