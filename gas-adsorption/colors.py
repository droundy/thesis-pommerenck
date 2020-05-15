import matplotlib.pyplot as plt
import matplotlib, re
import numpy as np

from matplotlib.font_manager import FontProperties
small_font = FontProperties()
small_font.set_size('small')

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = list(prop_cycle.by_key()['color'])

the_order = [
    'MOF5',
    'HKUST1',
    'PCN14',
    'Ni2dobdc',
    'Co2dobdc',
    'Mg2dobdc',
    'AX21',
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
    'Ni2dobdc': 'Ni$_2$(dodbc)',
    'Co2dobdc': 'Co$_2$(dodbc)',
    'Mg2dobdc':'Mg$_2$(dodbc)',
    'Zn2bdc2dabco': 'Zn$_2$(bdc)$_2$(dabco)',
    'CYCU3Al': 'CYCU$_3$Al',
    'UiO68Ant': 'UiO68Ant',
    'rhtMOF7': 'rht-MOF-7',
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
