var yNorth = srcNorth.data['y']
var xNorth = srcNorth.data['x']
const L = sliL.value
const Lrat = sliLrat.value
const Wrat = sliWrat.value
const wMax = Lrat*L
const wMin = wMax*Wrat
const val = rG.active

for (var i = 0; i < xNorth.length; i++) {
    xNorth[i] = i/(xNorth.length-1)*L
    yNorth[i] = fNorth(xNorth[i],L,wMin,wMax,val)
}

srcLin.data['yW'][1] = wMin
srcLin.data['yE'][1] = wMin
srcLin.data['xE'][0] = L
srcLin.data['xE'][1] = L
srcLin.data['xS'][1] = L

srcLin.change.emit();
srcNorth.change.emit();