// retrieve current settings
const L     = sliL.value
const Lrat  = sliLrat.value
const Wrat  = sliWrat.value
const wMax  = Lrat*L
const wMin  = wMax*Wrat
const Tmean = Math.exp(sliT.value)
const aniso = sliAni.value
const Tx    = Math.sqrt(Tmean**2 * (10**aniso))
const Ty    = Math.sqrt(Tmean**2 / (10**aniso))
const h2    = 0
const h1    = h2 + L*Math.exp(sliI.value)
const Q0    = (h1-h2)/L*Tx*(wMax-wMin)
const q     = sliQnor.value*Q0/L
const shape = rG.active
const por   = sliPor.value

// determine coefficients A
const N = 10
const M = 25
var A   = getA(N,M,L,wMin,wMax,shape,Tx,Ty,h1,h2,q)

// get fields of h and psi
const [hField, pField, pMin, pMax] = getFields(L,wMin,wMax,shape,h1,h2,A,Tx,Ty);

// extract contour line sets
var lvlsH = Bokeh.LinAlg.linspace(h2,h1,11)
lvlsH.shift()
lvlsH.pop()

const [XsH,YsH] = getIsolines(hField,lvlsH,L,wMin,wMax,shape)
srcIsoH.data['xs']  = XsH
srcIsoH.data['ys']  = YsH
for (var i = 0; i < lvlsH.length; i++) {
    lvlsH[i] = (lvlsH[i]-h2)/(h1-h2) 
} 
srcIsoH.data['clr'] = lvlsH

var lvlsP = Bokeh.LinAlg.linspace(pMin,pMax,11)
lvlsP.shift()
lvlsP.pop()
const [XsP,YsP] = getIsolines(pField,lvlsP,L,wMin,wMax,shape)
srcIsoP.data['xs'] = XsP
srcIsoP.data['ys'] = YsP

// extract single isoline for exchange zone (he) 
const [heX,heY] = getIsoline(pField,pField[0][0],L,wMin,wMax,shape)
    
if (heX.length==0) {
    srcHE.data['x'] = [0,L]
    srcHE.data['yBot']    = [0,0]
    srcHE.data['yTop']    = [0,0]
} else {
    srcHE.data['x']       = heX
    var heYbot = []
    for (var i = 0; i < heX.length; i++) {
       heYbot.push(0) 
    } 
    srcHE.data['yBot']    = heYbot
    srcHE.data['yTop']    = heY
}

// extract single isoline for hillslope zone (hs)
const pVal = pField[pField.length-1][0]
const [hsX,hsY] = getIsoline(pField,pVal,L,wMin,wMax,shape)

if (hsX.length==0) {
    srcHS.data['x'] = [0,L]
    srcHS.data['yBot']    = [0,0]
    srcHS.data['yTop']    = [0,0]
} else {
    const xMax = Math.max(...hsX)
    if ((xMax < L) && (Math.min(...hsY) < 0.001*wMin))  {
        for (var i = Math.ceil(20*xMax/L); i < 21; i++) {
            hsX.push(L*i/20)
            hsY.push(0)
        } 
    }
    srcHS.data['x']       = hsX
    var hsYtop = []
    for (var i = 0; i < hsX.length; i++) {
       hsYtop.push(fNorth(hsX[i],L,wMin,wMax,shape)) 
    } 
    srcHS.data['yBot']    = hsY
    srcHS.data['yTop']    = hsYtop
}

// determine exchange zone area
var xtemp1      = Array.from(srcHE.data['x'])
var ytemp1      = Array.from(srcHE.data['yBot'])
const xtemp2    = xtemp1.slice()
const ytemp2    = Array.from(srcHE.data['yTop'])
const xPoly     = xtemp1.concat(xtemp2.reverse())
const yPoly     = ytemp1.concat(ytemp2.reverse())
var Aex         = polyarea(xPoly,yPoly)

// define boolean to tell if zone exists
var Qex = getQex(h1,h2,L,wMin,wMax,Tx,Ty,A)
var zoneExists =  (Qex >= 1e-19)

// if zone exists: get TTD and adjust results text
if (zoneExists) {
    const [t,f] = getTTD(pField,L,wMin,wMax,h1,h2,Tx,Ty,A,shape,por);

    // normalize ttD
    const Tmax = Math.max(...t)
    for (var i = 0; i < t.length; i++) {
        t[i] = 100*t[i]/Tmax
        f[i] = 100*f[i]
    } 
    srcTTD.data.t = t
    srcTTD.data.f = f
    
    const postQ = Math.floor(Math.log10(Qex))
    const preQ = Qex/(10**postQ)
    const postA = Math.floor(Math.log10(Aex))
    const preA = Aex/(10**postA)
    const postT = Math.floor(Math.log10(Tmax))
    const preT = Tmax/(10**postT)
    divResult.text = 'The exchange flux is Q<sub>ex</sub>  = '
                         + preQ.toFixed(2) 
                         + '·10<sup>' 
                         + postQ.toFixed(0) 
                         + ' </sup> L<sup>3</sup>/T. '
                         + 'The exchange zone area is A<sub>ex</sub> = '
                         + preA.toFixed(2)
                         + "·10<sup>" 
                         + postA.toFixed(0) 
                         +  "</sup>  L<sup>2</sup>."
                         + 'The longest exchange travel time is t<sub>max</sub> = '
                         + preT.toFixed(2)
                         + "·10<sup>" 
                         + postT.toFixed(0) 
                         +  "</sup>  T."
} else {
    divResult.text = 'There is no exchange zone.'
    srcTTD.data.t = [0,50,100]
    srcTTD.data.f = [0,NaN,100]
}

// communicate updates
srcIsoH.change.emit();
srcIsoP.change.emit();
srcTTD.change.emit();
srcHE.change.emit();
srcHS.change.emit();