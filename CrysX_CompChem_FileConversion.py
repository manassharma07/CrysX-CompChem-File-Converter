import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import subprocess
import sys
import time
from io import StringIO

try:
    # from openbabel import OBMol, OBConversion
    import openbabel
except ModuleNotFoundError as e:
    subprocess.Popen([f'{sys.executable} -m pip install --global-option=build_ext --global-option="-I/usr/include/openbabel3" --global-option="-L/usr/lib/openbabel" openbabel'], shell=True)
    subprocess.Popen([f'{sys.executable} -m pip install --global-option=build_ext --global-option="-I/home/appuser/include/openbabel3" --global-option="-L/home/appuser/lib/openbabel" openbabel'], shell=True)
    subprocess.Popen([f'{sys.executable} -m pip install --global-option=build_ext --global-option="-I/home/appuser/usr/include/openbabel3" --global-option="-L/home/appuser/usr/lib/openbabel" openbabel'], shell=True)
    # wait for subprocess to install package before running your actual code below
    time.sleep(90)
    
import os
from openbabel import pybel


if os.path.exists('viz.html'):
    os.remove('viz.html')

# Set page config
st.set_page_config(page_title='CrysX - CompChem File Converter', layout='wide', page_icon="🧊",
menu_items={
         'About': "# This online tool allows you to enter a molecular geometry file in various formats, and convert it to another format that you desire."
     })


# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write('### Made By [Manas Sharma](https://www.bragitoff.com/about/)')
st.sidebar.write('### *Powered by*')
st.sidebar.write('* [Py3Dmol]() for Visualization')
st.sidebar.write('* [Open Babel](http://openbabel.org/) for Format Conversion')
st.sidebar.write('## Brought to you by [CrysX](https://www.bragitoff.com/crysx/)')
st.sidebar.write('## Cite us:')
st.sidebar.write('[Sharma, M. & Mishra, D. (2019). J. Appl. Cryst. 52, 1449-1454.](http://scripts.iucr.org/cgi-bin/paper?S1600576719013682)')

# Main app
st.write('# CrysX - CompChem File Converter')
st.write('This online tool allows you to enter a molecular geometry file in various formats, and convert it to another format that you desire.')

placeholder_xyz_str = '''12
Comment
C         -0.65914       -1.21034        3.98683
C          0.73798       -1.21034        4.02059
C         -1.35771       -0.00006        3.96990
C          1.43653       -0.00004        4.03741
C         -0.65915        1.21024        3.98685
C          0.73797        1.21024        4.02061
H         -1.20447       -2.15520        3.97369
H          1.28332       -2.15517        4.03382
H         -2.44839       -0.00006        3.94342
H          2.52722       -0.00004        4.06369
H         -1.20448        2.15509        3.97373
H          1.28330        2.15508        4.03386
'''



### CONVERSION ###

## INPUT ##
col1, col2 = st.columns(2)
col1.write('## INPUT')
input_format = col1.selectbox('Select the input file format',
     ( 'xyz', 'tmol', 'sdf', 'cif', 'poscar','cub','cube','fhiaims','mcif','mmcif','mdl', 'mol', 'mol2', 'outmol', 'pwscf', 'smi', 'pdb', 'smiles','txt','txyz','text'))
input_text_area = col1.empty()
input_geom_str = input_text_area.text_area(label='Enter the contents of the source file here', value = placeholder_xyz_str, placeholder = 'Put your text here', height=400, key = 'input_text_area')
# Get rid of empty lines
input_geom_str = os.linesep.join([s for s in input_geom_str.splitlines() if s])
uploaded_file = col1.file_uploader("You can also choose a file on your system")
if uploaded_file is not None:
    # To read file as bytes:
    bytes_data = uploaded_file.getvalue()

    # To convert to a string based IO:
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))

    # To read file as string:
    string_data = stringio.read()
    placeholder_xyz_str = string_data
    input_geom_str = input_geom_str = input_text_area.text_area(label='Enter the contents of the source file here', value = placeholder_xyz_str, placeholder = 'Put your text here', height=400)
    

## OUTPUT ##
col2.write('## OUTPUT')
output_format = col2.selectbox('Select the output file format',
     ( 'tmol', 'xyz', 'sdf', 'cif','cub','cube','fh','fhiaims','mcif','mmcif','mdl','mol','mol2','outmol','smi','pdb','smiles','txt','txyz','text'))


# References:
# https://openbabel.org/docs/dev/UseTheLibrary/PythonDoc.html
# https://github.com/kbsezginel/chem-tools-tutorials/blob/master/openbabel/openbabel.ipynb
# http://www.biotech.fyicenter.com/1000088_List_of_File_Formats_Supported_by_Open_Babel.html
# https://open-babel.readthedocs.io/en/latest/UseTheLibrary/PythonInstall.html
# https://open-babel.readthedocs.io/en/latest/UseTheLibrary/Python_Pybel.html
# https://open-babel.readthedocs.io/en/latest/UseTheLibrary/Python_PybelAPI.html#pybel.descs
# https://open-babel.readthedocs.io/en/latest/UseTheLibrary/Python_PybelAPI.html#pybel.descs

# obconversion = OBConversion()
# obconversion.SetInAndOutFormats(input_format, output_format)  # Set input and output formats
# obmol = OBMol()                                # Create openbabel molecule instance
# obconversion.ReadString(obmol, input_geom_str)         # Read file (file is read into obmol object)
# output_geom_str = obconversion.WriteString(obmol )   # Convert file to output format and save
output_geom_str = 'None'
try:
    mol = pybel.readstring(input_format, input_geom_str)
    # mol.make3D()
    output_geom_str = mol.write(output_format)
except Exception as e:
    print('There was a problem with the conversion', e)

col2.text_area(label='Converted geometry file in the format selected by you',value=output_geom_str, height=400)
col2.download_button(
     label="Download the converted file",
     data=output_geom_str,
     file_name='converted_output.'+output_format,
     mime='text/csv',
 )

### VISUALIZATION ####
style = st.selectbox('Visualization style',['ball-stick','line','cross','stick','sphere','cartoon','clicksphere'])
col1, col2 = st.columns(2)
spin = col1.checkbox('Spin', value = False)
showLabels = col2.checkbox('Show Labels', value = False)
# style='stick'
# style='cartoon'
# style='sphere'
view = py3Dmol.view(width=500, height=300)
structure_for_visualization = ''
try:
    mol = pybel.readstring(input_format, input_geom_str)
    # mol.make3D()
    if style=='cartoon':
        structure_for_visualization = mol.write('pdb')
    else:
        structure_for_visualization = mol.write('xyz')
except Exception as e:
    print('There was a problem with the conversion', e)
if style=='cartoon':
    view.addModel(structure_for_visualization, 'pdb')
else:
    view.addModel(structure_for_visualization, 'xyz')
if style=='ball-stick': # my own custom style
    view.setStyle({'sphere':{'colorscheme':'Jmol','scale':0.3},
                       'stick':{'colorscheme':'Jmol', 'radius':0.}})
else:
    view.setStyle({style:{'colorscheme':'Jmol'}})
# Label addition template
# view.addLabel('Aromatic', {'position': {'x':-6.89, 'y':0.75, 'z':0.35}, 
#             'backgroundColor': 'white', 'backgroundOpacity': 0.5,'fontSize':18,'fontColor':'black',
#                 'fontOpacity':1,'borderThickness':0.0,'inFront':'true','showBackground':'false'})
if showLabels:
    for atom in mol:
        view.addLabel(str(atom.idx), {'position': {'x':atom.coords[0], 'y':atom.coords[1], 'z':atom.coords[2]}, 
            'backgroundColor': 'white', 'backgroundOpacity': 0.5,'fontSize':18,'fontColor':'black',
                'fontOpacity':1,'borderThickness':0.0,'inFront':'true','showBackground':'false'})
view.zoomTo()
view.spin(spin)
view.setClickable({'clickable':'true'});
view.enableContextMenu({'contextMenuEnabled':'true'})
view.show()
view.render()
# view.png()
t = view.js()
f = open('viz.html', 'w')
f.write(t.startjs)
f.write(t.endjs)
f.close()

HtmlFile = open("viz.html", 'r', encoding='utf-8')
source_code = HtmlFile.read() 
components.html(source_code, height = 300, width=900)
HtmlFile.close()
st.write('## Properties of the given chemical system')
st.write('### Weight ')
st.write(str(mol.molwt))
st.write('### Formula ')
st.write(str(mol.formula))
st.write('### Exact mass ')
st.write(str(mol.exactmass))
st.write('### Spin multiplicity')
st.write(str(mol.spin))
st.write('### Charge')
st.write(str(mol.charge))
# st.write('### Conformers')
# st.write((mol.conformers))
# st.write(mol.data)
# for atom in mol:
#     st.write(atom.coords[1])