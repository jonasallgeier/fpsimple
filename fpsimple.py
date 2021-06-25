from io import DEFAULT_BUFFER_SIZE
import numpy as np
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CheckboxButtonGroup, Panel, Tabs, CustomJS, Div, Slider, WheelZoomTool, RadioButtonGroup, FuncTickFormatter
from bokeh.plotting import Figure
from bokeh.resources import CDN
from bokeh.embed import components
from colorcet import coolwarm
from bokeh.transform import linear_cmap
from datetime import date

# set up two figures (flow net and travel time distribution)
figFN = Figure(tools="pan,reset,save,wheel_zoom",
            match_aspect=True,
            x_axis_label='x',
            y_axis_label='y')
figFN.xgrid.grid_line_color = None
figFN.ygrid.grid_line_color = None

figTTD = Figure(tools="pan,reset,save,wheel_zoom",
            y_axis_label='F(t) in %',
            x_axis_label='t in % of tMax',
            x_range=(0,100), y_range=(0, 100))

for fig in [figFN,figTTD] :
    fig.yaxis.major_label_text_font_size    = "25px"
    fig.xaxis.major_label_text_font_size    = "25px"
    fig.xaxis.axis_label_text_font_size     = "25px"
    fig.yaxis.axis_label_text_font_size     = "25px"
    fig.toolbar.logo                        = 'grey' # 'grey'/'normal'/None
    fig.toolbar.active_scroll               = figFN.select_one(WheelZoomTool)
    fig.toolbar.active_drag                 = None # disable pan, because it's bad on mobile
    fig.margin                              = (5,0,0,0)
    fig.min_width                           = 250
    fig.height_policy                       = "fit"
    fig.width_policy                        = "max"

# these settings define the sliders (min, max, ini, step)
bL      = [100, 3000, 1000, 10]
bI      = [np.log(0.0003), np.log(0.03), np.log(0.003), (np.log(0.03)-np.log(0.0003))/100]
bX      = [0.1, 0.5, 0.3, 0.2/100]
bN      = [0.4, 1.0, 1.0, 0.6/100]
bT      = [np.log(0.000001), np.log(0.005), (np.log(0.005)+np.log(1e-6))/2, (np.log(0.005)-np.log(0.000001))/100]
bA      = [-0.5, 0.5, 0, 2/100]
bQ      = [0.0, 3.0, 0.0, 3/100]
bP      = [0.1, 5.0, 1.0, 0.1]
shapes  = ["cosinusoidal", "bump", "composite"]

# define initial data set
iniL        = bL[2]
iniWmin     = bL[2]*bX[2]
x           = np.linspace(0,iniL,50)
yTop        = 0*x+iniWmin
yBot        = 0*x
xI1         = np.linspace(0,iniL,11)
yI1         = np.linspace(0,iniWmin,11)
xI1         = xI1[1:-1]
yI1         = yI1[1:-1]

XSH, YSH, XSP, YSP, CSH = [], [], [], [], []
for xT in xI1:
    XSH.append([xT,xT])
    YSH.append([0,iniWmin])
    CSH.append([1-xT/iniL])
for yT in yI1:
    XSP.append([0,iniL])
    YSP.append([yT,yT])

# define CDSs that the plots can use
srcESW      = ColumnDataSource(data=dict(xW=np.array([0,0]), yW=np.array([0,iniWmin]),
                                            xE=np.array([iniL,iniL]), yE=np.array([0,iniWmin]),
                                            xS=np.array([0,iniL]), yS=np.array([0,0])))
srcNorth    = ColumnDataSource(data=dict(x=x, y=yTop))
srcHE       = ColumnDataSource(data=dict(x=x, yTop=(np.NaN+yTop), yBot=(np.NaN+yBot)))
srcHS       = ColumnDataSource(data=dict(x=x, yTop=(np.NaN+yTop), yBot=(np.NaN+yBot)))
srcTTD      = ColumnDataSource(data=dict(t=[0,50,100],f=[0,np.NaN,100]))
srcIsoH     = ColumnDataSource(data=dict(xs=XSH,ys=YSH,clr=CSH))
srcIsoP     = ColumnDataSource(data=dict(xs=XSP,ys=YSP))
mapper      = linear_cmap(field_name='clr', palette=coolwarm ,low=0,high=1)

# initialize plotting glyphs; "h" --> handle
    # isopotential lines of stream function
hIsoP   = figFN.multi_line(xs='xs',ys='ys', line_color='rgba(230,232,237, 1.0)',  line_width=5, source=srcIsoP,line_alpha=0.7)
    # hillslope area
hHS     = figFN.varea(x='x', y1='yTop',y2='yBot', source=srcHS, fill_alpha=0.4,fill_color=(201,174,105),visible=True)
    # hyporheic exchange area
hHE     = figFN.varea(x='x', y1='yTop',y2='yBot', source=srcHE, fill_alpha=0.4,fill_color=(105,186,201),visible=True)
    # isopotential lines of hydraulic head
hIsoH   = figFN.multi_line(xs='xs',ys='ys', line_color=mapper,  line_width=5, source=srcIsoH,line_alpha=1.0)
    # western boundary
hW      = figFN.line(x='xW', y='yW', source=srcESW, line_width=5, color='rgba(192,2,6,1.0)', line_cap="round",visible=True)
    # eastern boundary
hE      = figFN.line(x='xE', y='yE', source=srcESW, line_width=5, color='rgba(32, 81, 219, 1.0)', line_cap="round",visible=True)
    # southern boundary
hS      = figFN.line(x='xS', y='yS', source=srcESW, line_width=5, color='rgba(105,186,201,1.0)', line_cap="round",visible=True)
    # northern boundary
hN      = figFN.line(x='x', y='y', source=srcNorth, line_width=5, color='rgba(201,174,105,1.0)', line_cap="round",visible=True)
    # travel time distribution
ttD     = figTTD.line(x='t',y='f',source=srcTTD, line_width=5, line_cap="round",color='rgba(219,33,82,1.0)')

# define sliders for input
sliL            = Slider(start=bL[0], end=bL[1], value=bL[2], step=bL[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(0)+' [L]'"""))
sliI            = Slider(start=bI[0], end=bI[1], value=bI[2], step=bI[3],
                    format=FuncTickFormatter(code="""return Math.exp(tick).toExponential(1).toString()+' [L/L]'"""))
sliLrat         = Slider(start=bX[0], end=bX[1], value=bX[2], step=bX[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(2)+' [L/L]'"""))
sliWrat         = Slider(start=bN[0], end=bN[1], value=bN[2], step=bN[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(2)+' [L/L]'"""))
sliT            = Slider(start=bT[0], end=bT[1], value=bT[2], step=bT[3],
                    format=FuncTickFormatter(code="return Math.exp(tick).toExponential(1).toString() + ' [L^2/T]'"))
sliAni          = Slider(start=bA[0], end=bA[1], value=bA[2], step=bA[3],
                    format=FuncTickFormatter(code="""return (10**tick).toFixed(2)+' [-]'"""))
sliQnor         = Slider(start=bQ[0], end=bQ[1], value=bQ[2], step=bQ[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(2)+' [-]'"""))
sliPor          = Slider(start=bP[0], end=bP[1], value=bP[2], step=bP[3],
                    format=FuncTickFormatter(code="""return tick.toFixed(1)+' [L]'"""))

for sli in [sliL, sliI, sliLrat, sliWrat, sliT, sliAni, sliQnor, sliPor] :
    sli.bar_color   = 'rgba(60, 60, 60, 0.6)'

# define slider labels, which are placed left to sliders
LblL         = Div(text="""<i>L</i>""")
LblWrat      = Div(text="""<i>w<sub>max</sub></i>/<i>w<sub>min</sub></i>""")
LblI         = Div(text="""<i>Δh/L</i>""")
LblT         = Div(text="""(<i>T<sub>x</sub></i>·<i>T<sub>y</sub></i>)<sup>1/2</sup>""")
LblAni       = Div(text="""<i>T<sub>x</sub></i>/<i>T<sub>y</sub></i>""")
LblQnor      = Div(text="""<i>Q</i><sup>~</sup><sub>north</sub>""")
LblLrat      = Div(text="""<i>w<sub>max</sub></i>/<i>L</i>""")
LblPor       = Div(text="""<i>Φ</i>""")

for lbl in [LblL,LblWrat,LblI,LblT,LblAni,LblQnor,LblLrat,LblPor]:
    lbl.style={'font-size': '115%'}

#  define radio button group for shape selection
radioGroup  = RadioButtonGroup(labels=shapes, active=0,margin=(10,3,10,3))

# define two check box button groups for visibility toggles
toggles = CheckboxButtonGroup(labels=["west","south","east","north"], active=[0,1,2,3],margin=(10,3,10,3))
toggles.js_on_click(CustomJS(args=dict(hIsoP=hIsoP,hIsoH=hIsoH,hE=hE,hW=hW,hS=hS,hN=hN),code="""
    var act         = this.active
    hW.visible      = act.includes(0)
    hS.visible      = act.includes(1)
    hE.visible      = act.includes(2)
    hN.visible      = act.includes(3)
"""))

toggles2 = CheckboxButtonGroup(labels=["heads","stream","zone south", "zone north"], active=[0,1,2,3],margin=(10,3,10,3))
toggles2.js_on_click(CustomJS(args=dict(hHE=hHE,hHS=hHS,hIsoP=hIsoP,hIsoH=hIsoH),code="""
    var act         = this.active
    hIsoH.visible   = act.includes(0)
    hIsoP.visible   = act.includes(1)
    hHE.visible     = act.includes(2)
    hHS.visible     = act.includes(3)
"""))

# define div that displays results as text
divResult   = Div(text="""There is no exchange zone.""",sizing_mode="stretch_width",visible=True,css_classes=["divResult"],style={'font-size': '135%'})

# first callback; adapts changes to geometry
with open('callbackGeometry.js','r') as file:
    cbCode = file.read()
callbackG = CustomJS(args=dict( srcNorth=srcNorth,
                                srcLin=srcESW,
                                sliL=sliL,
                                sliLrat=sliLrat,
                                sliWrat=sliWrat,
                                rG=radioGroup     ), code=cbCode)

#  second callback; updates the flow field and results
with open('callbackFlowSolution.js','r') as file:
    cbCode = file.read()
callbackF = CustomJS(args=dict( sliI=sliI,
                                sliT=sliT,
                                sliPor=sliPor,
                                sliAni=sliAni,
                                sliQnor=sliQnor,
                                sliL=sliL,
                                sliLrat=sliLrat,
                                sliWrat=sliWrat,
                                srcIsoH=srcIsoH,
                                srcHE=srcHE,
                                srcHS=srcHS,
                                srcIsoP=srcIsoP,
                                srcTTD=srcTTD,
                                divResult=divResult,
                                rG=radioGroup     ), code=cbCode)

# attach geometry callback to geometry properties
for sli in [sliL,sliLrat,sliWrat] :
    sli.js_on_change('value',callbackG)
radioGroup.js_on_click(callbackG)

# attach flow field callback to all properties
for sli in [sliL,sliLrat,sliWrat,sliI,sliT,sliAni,sliQnor,sliPor] :
    sli.js_on_change('value',callbackF)
radioGroup.js_on_click(callbackF)

# define the layout
layout = column(row(LblL,sliL,sizing_mode="stretch_width"),
                row(LblLrat,sliLrat,sizing_mode="stretch_width"),
                row(LblWrat,sliWrat,sizing_mode="stretch_width"),
                row(LblI,sliI,sizing_mode="stretch_width"),
                row(LblT,sliT,sizing_mode="stretch_width"),
                row(LblAni,sliAni,sizing_mode="stretch_width"),
                row(LblQnor,sliQnor,sizing_mode="stretch_width"),
                row(LblPor,sliPor,sizing_mode="stretch_width"),
                radioGroup,toggles,toggles2,divResult,sizing_mode="stretch_width")


tab1 = Panel(child=figFN, title="Flow Net")
tab2 = Panel(child=figTTD, title="Travel Times")

plots = Tabs(tabs=[tab1,tab2])

# print the script to file
script, (div1, div2) = components((layout,plots))
f = open("./themodel.js", "w")
script = "\n".join(script.split("\n")[2:-1])
f.write(script)
f.close()

# read in the template file
with open('template', 'r') as file :
  filedata = file.read()

# replace the target strings
filedata = filedata.replace('+LEFT+', div1)
filedata = filedata.replace('+RIGHT+', div2)
filedata = filedata.replace('+DATE+',date.today().strftime("%Y-%m-%d"))

# write to html file
with open('fpsimple.html', 'w') as file:
  file.write(filedata)