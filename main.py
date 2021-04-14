import numpy as np
from colorcet import coolwarm
from bokeh.io import curdoc
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Div, Slider, TextInput, WheelZoomTool, Button, RadioButtonGroup, FuncTickFormatter
from bokeh.plotting import figure
from bokeh.transform import linear_cmap
from bokeh.palettes import Spectral6
from skimage import measure
from bokeh.server.server import Server
from random import randint

# Set up the figure
fig = figure(tools="pan,reset,save,wheel_zoom",min_width=250,min_height=400,height_policy="fit",width_policy="max",match_aspect=True,name="rightPanel",css_classes=["main"])
fig.toolbar.autohide                    = False
fig.toolbar.logo                        = 'grey' # 'grey'/'normal'/None
fig.toolbar.active_scroll               = fig.select_one(WheelZoomTool)
fig.toolbar.active_drag                 = None # disable pan, because it's bad on mobile
fig.background_fill_color               = "#f2f2f2"
fig.border_fill_color                   = "#f2f2f2"
fig.xgrid.grid_line_color               = "#bfbfbf"
fig.ygrid.grid_line_color               = "#bfbfbf"
fig.yaxis.major_label_text_font_size    = "1.5em"
fig.xaxis.major_label_text_font_size    = "1.5em"

# Set up data containers (with dummy data)
srcOL       = ColumnDataSource(data=dict(x=[0], y=[0]))
srcH        = ColumnDataSource(data=dict(xs=[[0], [0]],ys=[[0], [0]],clr=[[0], [0]]))
srcP        = ColumnDataSource(data=dict(xs=[[0], [0]],ys=[[0], [0]],clr=[[0], [0]]))

# Initialize plotting glyphs
mapper      = linear_cmap(field_name='clr', palette=coolwarm ,low=0,high=1)
mapperP     = linear_cmap(field_name='clr', palette=('#67757E','#2D3439') ,low=0,high=1)
fig.patch(x='x', y='y', source=srcOL, fill_color=(160, 202, 230), line_color=(33, 97, 140), line_width=5)
isoH        = fig.multi_line(xs='xs',ys='ys', line_color=mapper,  line_width=4, source=srcH)
isoP        = fig.multi_line(xs='xs',ys='ys', line_color=mapperP, line_width=4, source=srcP)
fig.patch(x='x', y='y', source=srcOL, fill_color=(160, 202, 230), fill_alpha=0, line_color=(33, 97, 140), line_width=5)

bL = [100, 3000, 1550, 50]
bI = [np.log(0.0003), np.log(0.03), np.log(0.003), (np.log(0.03)-np.log(0.0003))/100]
bX = [0.1, 0.5, 0.3, 0.2/100]
bN = [0.4, 1.0, 0.7, 0.6/100]
bT = [np.log(0.000001), np.log(0.005), (np.log(0.005)+np.log(1e-6))/2, (np.log(0.005)-np.log(0.000001))/100]
bA = [-0.5, 0.5, 0, 2/100]
bQ = [0.0, 3.0, 0.0, 3/100]

# Set up widgets
sliL        = Slider(start=bL[0], end=bL[1], value=bL[2], step=bL[3], tooltips=True, bar_color=(160, 202, 230), title="length")
sliI        = Slider(start=bI[0], end=bI[1], value=bI[2], step=bI[3], tooltips=True, bar_color=(160, 202, 230), title="gradient",
                        format=FuncTickFormatter(code="return Math.exp(tick).toExponential(1)"))
sliLrat     = Slider(start=bX[0], end=bX[1], value=bX[2], step=bX[3], tooltips=True, bar_color=(160, 202, 230), title="aspect ratio")
sliWrat     = Slider(start=bN[0], end=bN[1], value=bN[2], step=bN[3], tooltips=True, bar_color=(160, 202, 230), title="width ratio")
sliT        = Slider(start=bT[0], end=bT[1], value=bT[2], step=bT[3], tooltips=True, bar_color=(160, 202, 230), title="transmissivity",
                format=FuncTickFormatter(code="return Math.exp(tick).toExponential(1)"))
sliAni      = Slider(start=bA[0], end=bA[1], value=bA[2], step=bA[3], tooltips=True, bar_color=(160, 202, 230), title="anisotropy")
sliQnor     = Slider(start=bQ[0], end=bQ[1], value=bQ[2], step=bQ[3], tooltips=True, bar_color=(160, 202, 230), title="northern influx")

divResult   = Div(text="""Initial text.""",sizing_mode="stretch_width",visible=True,css_classes=["divResult"],style={'font-size': '135%'})
shapes      = ["cosinusoidal", "bump", "composite"]
radioGroup  = RadioButtonGroup(labels=shapes, active=0)
btnRandom   = Button(label="Random", button_type="primary")
btnReset    = Button(label="Reset", button_type="primary")

# Initial values
fNorth  = lambda x: x+1
N       = 10
M       = np.ceil((0.5*N)**2)
L       = 1000
Tx      = 1e-3
Ty      = 1e-3
wMax    = 360
wMin    = 250
h1      = 10
h2      = 7
q       = 0
A       = np.NaN

# helper function
def F(x,y):
    return -q*x/Tx+(h1-h2)/L*y # opposite sign to paper!

# function to determine coefficients A
def getA():
    # x     = 0.5*L+0.5*L*np.cos((2*(np.arange(1,M+1,1))-1)/(2*M)*np.pi)
    x     = np.linspace(0,L,int(M))
    y     = fNorth(x)

    Umtx  = U(np.arange(0,N+1,1),x,y)
    Fvec  = F(x,y)
    rhs   = np.matmul(Umtx,Fvec[:,np.newaxis])
    lhs   = np.matmul(Umtx,Umtx.transpose())

    return np.linalg.lstsq(lhs, rhs, rcond=None)[0]

# helper function
def U(n,x,y):
    c = n[:,np.newaxis]*np.pi/L
    kap = np.sqrt(Tx/Ty)
    b1 = np.matmul(c,kap*(np.expand_dims(y,axis=y.ndim-1)-wMax))
    b2 = np.matmul(c,kap*(-np.expand_dims(y,axis=y.ndim-1)-wMax))
    b3 = np.matmul(c,np.expand_dims(x,axis=x.ndim-1))

    v = (np.exp(b1)+np.exp(b2))/(1+np.exp(-2*kap*c*wMax)) * np.cos(b3)
    return v/kap

# helper function
def V(n,x,y):
    c   = n*np.pi/L
    kap = np.sqrt(Tx/Ty)
    b1 = c * (kap*(y-wMax))
    b2 = c * (kap*(-y-wMax))
    return (np.exp(b1) - np.exp(b2))/( 1 + np.exp(-2*kap*c*wMax) )*np.sin(c*x)

# function to retrieve hydraulic head
def head(x,y):
    h = h1 + (h2-h1)/L*x
    for n in range(1,N+1):
        h = h + A[n]*V(n,x,y)
    return h

# function to retrieve stream function values
def psi(x,y):
    ps = A[0] + (h2-h1)/L*y
    An = A[1:]
    Ut  = U(np.arange(1,N+1,1),x,y)
    ps = ps + np.squeeze(np.matmul(An.transpose(),Ut))
    return -ps*Tx

# function to obtain exchange flux
def getQex():
    xT  = np.linspace(0, L, 100)
    yT  = 0*xT
    ps  = psi(xT,yT)
    return np.amin([ps[0], ps[-1]])-np.amin(ps)

# Set up callbacks
def drawOutline(attr,old,new):
    # Generate the new curve
    xN = np.linspace(L, 0, 50)
    yN = fNorth(xN)

    # Set up data
    x = np.concatenate(([0, L, L], xN, [0]))
    y = np.concatenate(([0, 0, wMin], yN, [wMin]))

    # update outline data
    srcOL.data = dict(x=x, y=y)

    resetSolution(0,0,0)

def drawFlowNet(attr,old,new):
    global A
    # solve the problem
    A = getA()

    # update scatter data
    val     = 25*np.min([wMax,L])
    xSc     = np.linspace(0, L, (np.ceil(val/wMax)).astype(np.int))
    ySc     = np.linspace(0,1, (np.ceil(val/L)).astype(np.int))
    XSc, YSc = np.meshgrid(xSc,ySc)
    YSc     = np.multiply(YSc,fNorth(XSc))

    # get contour lines
    for case in ["head","psi"]:
        # get the field
        if (case=="head"):
            field = head(XSc,YSc)
            mapper['transform'].low = np.min(field)
            mapper['transform'].high = np.max(field)
        elif (case =="psi"):
            field = psi(XSc,YSc)

        # initialize
        lvls    = np.linspace(np.min(field),np.max(field),11)
        lvls    = lvls[1:-1]
        if (case =="psi"):
            lvls = np.append(lvls, psi(0*xSc[0:1],0*ySc[0:1]))
        xs      = []
        ys      = []
        sclX    = np.size(field,1)-1
        sclY    = np.size(field,0)-1

        # iterate through contour levels
        try:
            for lvl in lvls:
                contours = measure.find_contours(field, lvl)
                xC = np.empty([0])
                yC = np.empty([0])
                for contour in contours:
                    xC = np.append(xC,contour[:,1]*L/sclX)
                    yC = np.append(yC,contour[:,0]*fNorth(xC)/sclY)
                    xC = np.append(xC,np.NaN)
                    yC = np.append(yC,np.NaN)
                xs.append(xC.tolist())
                ys.append(yC.tolist())
            if (case=="head"):
                clrs = lvls.flatten()
                srcH.data = dict(xs=xs,ys=ys,clr =clrs)
                isoH.visible = True
            elif (case =="psi"):
                clrs = lvls.flatten()
                clrs[0:-2] = 0
                clrs[-1] = 1
                srcP.data = dict(xs=xs,ys=ys,clr =clrs)
                isoP.visible = True
        except:
            isoP.visible = False
            isoH.visible = False

    Qex = getQex()
    if (Qex==0):
        divResult.text = "The exchange flux is Q<sub>ex</sub>  = 0.00 m<sup>3</sup>/s."
    else:
        post = np.floor(np.log10(Qex))
        pre = Qex/(10**post)
        divResult.text = "The exchange flux is Q<sub>ex</sub>  = {:.2f}·10<sup>{}</sup> m<sup>3</sup>/s.".format(pre,np.int(post))
    
    divResult.visible = True

def resetSolution(attrname,old,new):
    global A
    A = np.NaN
    isoH.visible = False
    isoP.visible = False
    hideText()

def hideText():
    divResult.visible = False

# more callbacks
def changedL(attrname,old,new):
    global L, wMax, wMin
    L = new
    wMax = sliLrat.value*L
    wMin = sliWrat.value*wMax

def changedWmin(attrname,old,new):
    global wMin,wMax
    wMin = wMax*new

def changedWmax(attrname,old,new):
    global wMax, L, wMin
    wMax    = L*new
    wMin    = sliWrat.value*wMax

def changedShape(attrname,old,new):
    global fNorth
    shape = radioGroup.labels[new]

    if shape=="bump":
        fNorth = lambda x: wMin+(wMax-wMin)*np.exp(1-1/(1-(2*x/L-1)**2+1e-9))    
    elif shape=="cosinusoidal":
        fNorth = lambda x: wMin+0.5*(wMax-wMin)*(1-np.cos(2*np.pi*x/L))
    else:
        fNorth = lambda x: wMin+(wMax-wMin)*0.5*( (1-np.cos(np.pi/(14/40)*(x/L-1/40)) )*((x/L>=1/40) & (x/L<15/40)) + (1+np.cos(np.pi/(14/40)*(x/L-25/40)) )*((x/L>=25/40) & (x/L <= 39/40)) + 2*((x/L >=15/40) & (x/L < 25/40)) )
    hideText()

def changedAni(attrname,old,new):
    global Tx,Ty
    Tm = np.sqrt(Tx*Ty)
    Tx = np.sqrt(Tm**2*(10**new))
    Ty = np.sqrt(Tm**2/(10**new))

def changedI(attrname,old,new):
    global h1, L
    h1 = h2 + L*np.exp(new)
    hideText()

def changedT(attrname,old,new):
    global Tx, Ty
    Tm = np.exp(new)
    Tx = np.sqrt(Tm**2*(10**sliAni.value))
    Ty = np.sqrt(Tm**2/(10**sliAni.value))
    hideText()

def changedQnor(attrname,old,new):
    global q, h1, h2, L, Tx, wMax, wMin
    Q0 = (h1-h2)/L*Tx*(wMax-wMin)
    q = new*Q0/L
    hideText()

# reset functionality
def resetSettings(attrname):
    setSettings(np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.0]))

# randomize functionality
def randSettings(attrname):
    setSettings(np.random.rand(7))

# functions to turn callbacks off and on
def callbacksOff():
    sliL.remove_on_change('value', changedL,drawOutline)
    sliL.remove_on_change('value_throttled', drawFlowNet)
    sliI.remove_on_change('value',changedI)
    sliI.remove_on_change('value_throttled',drawFlowNet)
    sliLrat.remove_on_change('value', changedWmax,drawOutline)
    sliLrat.remove_on_change('value_throttled', drawFlowNet)
    sliWrat.remove_on_change('value', changedWmin,drawOutline)
    sliWrat.remove_on_change('value_throttled', drawFlowNet)
    sliT.remove_on_change('value',changedT)
    sliT.remove_on_change('value_throttled',drawFlowNet)
    sliAni.remove_on_change('value',changedAni,drawOutline)
    sliAni.remove_on_change('value_throttled',drawFlowNet)
    sliQnor.remove_on_change('value',changedQnor)
    sliQnor.remove_on_change('value_throttled',drawFlowNet)
    
def callbacksOn():
    sliL.on_change('value', changedL,drawOutline)
    sliL.on_change('value_throttled', drawFlowNet)
    sliI.on_change('value',changedI)
    sliI.on_change('value_throttled',drawFlowNet)
    sliLrat.on_change('value', changedWmax,drawOutline)
    sliLrat.on_change('value_throttled', drawFlowNet)
    sliWrat.on_change('value', changedWmin,drawOutline)
    sliWrat.on_change('value_throttled', drawFlowNet)
    sliT.on_change('value',changedT)
    sliT.on_change('value_throttled',drawFlowNet)
    sliAni.on_change('value',changedAni,drawOutline)
    sliAni.on_change('value_throttled',drawFlowNet)
    sliQnor.on_change('value',changedQnor)
    sliQnor.on_change('value_throttled',drawFlowNet)

# setting function for randomized cases
def setSettings(pIn):
    global L,Tx,Ty,wMax,wMin,h1,h2,q
    callbacksOff()
    hideText()

    sliL.value      = bL[0]+ (bL[1]-bL[0])*pIn[0]
    sliLrat.value   = bX[0]+ (bX[1]-bX[0])*pIn[1]
    sliWrat.value   = bN[0]+ (bN[1]-bN[0])*pIn[2]
    sliI.value      = bI[0]+ (bI[1]-bI[0])*pIn[3]
    sliAni.value    = bA[0]+ (bA[1]-bA[0])*pIn[4]
    sliT.value      = bT[0]+ (bT[1]-bT[0])*pIn[5] 
    sliQnor.value   = bQ[0]+ (bQ[1]-bQ[0])*pIn[6]

    L       = sliL.value
    h1      = h2 + L*np.exp(sliI.value)
    wMax    = L*sliLrat.value
    wMin    = wMax *sliWrat.value
    Tm      = np.exp(sliT.value)
    An      = sliAni.value
    Tx      = np.sqrt(Tm**2*(10**An))
    Ty      = np.sqrt(Tm**2/(10**An))
    Q0      = (h1-h2)/L*Tx*(wMax-wMin)
    q       = sliQnor.value*Q0/L

    callbacksOn()
    drawOutline(0,0,0)
    drawFlowNet(0,0,0)

# main code
callbacksOn()
btnReset.on_click(resetSettings)
btnRandom.on_click(randSettings)

radioGroup.on_change('active',changedShape,drawOutline,drawFlowNet)

# Set up layouts and add to document
sliderLayout = column(sliL,sliI,sliLrat,sliWrat,sliT,sliAni,sliQnor,radioGroup,row(btnRandom,btnReset,sizing_mode="stretch_width"),
                        divResult,width=400,sizing_mode="stretch_width",name="leftPanel",css_classes=["side"])

curdoc().add_root(sliderLayout)
curdoc().add_root(fig)

# Update plot once initially to get rid of dummy data
changedShape(0,0,0)
resetSettings(0)
curdoc().title = 'fpSimple'