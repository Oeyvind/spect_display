<Cabbage>
form caption("SVG FFT") size(600, 240), colour(30, 35, 40),  guiMode("queue") pluginId("def1")
;nslider bounds(0,5,20,10), channel("dBmin"), range(-96,-6,-60,1)
label bounds(8, 5, 50, 10), text("scale"), align("left")
nslider bounds(5, 18, 60, 22), channel("dBmin"), range(6,96,60,1), valuePostfix(" dB")
label bounds(73, 5, 65, 10), text("resolution"), align("left")
nslider bounds(70, 18, 65, 22), channel("resolution"), range(0.2,4.0,0.5,1, 0.1), valuePostfix(" semi")
label bounds(142, 5, 45, 10), text("skew"), align("left")
nslider bounds(140, 18, 45, 22), channel("skew"), range(0.3,1.5,1.0,1,0.01)
label bounds(192, 5, 70, 10), text("upd.rate"), align("left")
nslider bounds(190, 18, 50, 22), channel("updaterate"), range(5,50,15,1)

label bounds(245, 5, 180, 10), text("-  -  -  -  -  colours  -  -  -  -  -"), align("centre")
combobox bounds(245, 19, 90,17), channel("colour1"), channelType("string"), items( "red", "orange", "chocolate", "navajowhite", "khaki", "gold", "goldenrod", "yellow", "lime", "springgreen",  "aquamarine", "aqua", "deepskyblue", "pink", "hotpink", "fuchsia", "violet", "mediumvioletred"), value("chocolate")
combobox bounds(340, 19, 90,17), channel("colour2"), channelType("string"),  items( "red", "orange", "chocolate", "navajowhite", "khaki", "gold", "goldenrod", "yellow", "lime", "springgreen",  "aquamarine", "aqua", "deepskyblue", "pink", "hotpink", "fuchsia", "violet", "mediumvioletred"), value("navajowhite")

label bounds(438, 5, 60, 10), text("linewidth"), align("left")
nslider bounds(435, 18, 50, 22), channel("linewidth"), range(0.0,2.0,0.8)
label bounds(493, 5, 45, 10), text("fill"), align("left")
nslider bounds(490, 18, 45, 22), channel("fill"), range(0.0,1.0,0.55, 1, 0.01)
label bounds(543, 5, 60, 10), text("gradient"), align("left")
nslider bounds(540, 18, 50, 22), channel("gradient"), range(0.0,1.0,0.4, 1, 0.01)


image bounds(0, 40, 600, 200), channel("displayGrid"), colour(0,0,0,0) 
image bounds(0, 40, 600, 200), channel("displaySig1"), colour(0,0,0,0) 
image bounds(0, 40, 600, 200), channel("displaySig2"), colour(0,0,0,0) 
image bounds(0, 40, 600, 200), channel("displaySig3"), colour(0,0,0,0) 
;csoundoutput bounds(0, 240, 600, 200)
</Cabbage>
<CsoundSynthesizer>
<CsOptions>
-n -d -+rtmidi=NULL -M0 -m0d 

</CsOptions>
<CsInstruments>
; Initialize the global variables. 
ksmps = 32
nchnls = 2
0dbfs = 1

opcode PowArray, i[], i[]i
  ; normalize array, then shape, and rescale back to original 
  iArr[], ipow xin
  imax maxarray iArr
  iArr /= imax
  iArr2[] pow iArr, ipow
  iArr2 *= imax
  xout iArr2
endop

opcode MakeSplitpoints, i[]i[], iii
  ; make splitpoints for fft display so that we reduce the data array size according to the desired frequency resolution (in semitones)
  ; input is fft size, resolution (in semitones)
  ; output is
  ; array with the indices for the splitpoints, and
  ; array with the frequency for each display bin
  ifftsize, ifft_binsize, iresolution xin
  iFftfreqs[] genarray 0, sr/2, ifft_binsize
  iSplitindices[] init ifftsize/2 ; long enough, we will trim it later
  iDispfreqs[] init ifftsize/2 ; long enough, we will trim it later
  inindex = 0
  ioutindex = 0
  while iFftfreqs[inindex] <= sr/2-ifft_binsize do
    iDispfreqs[ioutindex] = iFftfreqs[inindex]
    iSplitindices[ioutindex] = inindex
    inindex += 1
    while iFftfreqs[inindex] <= min(iDispfreqs[ioutindex]*semitone(iresolution),sr/2-ifft_binsize) do
      inindex += 1
    od
    ioutindex += 1
  od
  trim iDispfreqs, ioutindex
  trim iSplitindices, ioutindex
  xout iSplitindices, iDispfreqs
endop

opcode FindFreqBins, i[], i[]i[]
  ; given a set of frequencies (vertical lines in the frequency grid), find the corresponding display bin indices
  iGridfreqs[], iDispfreqs[] xin
  iFreqBinindices[] init lenarray(iGridfreqs)+1
  index = 0
  index2 = 0
  while index < lenarray(iGridfreqs) do
    while iGridfreqs[index] > iDispfreqs[index2] do
      iFreqBinindices[index] = index2
      index2 += 1
    od
    index += 1
  od
  iFreqBinindices[lenarray(iFreqBinindices)-1] = lenarray(iDispfreqs); write last grid position
  iFreqBinindices += 1
  ; since the display grid will most often not fall on a display bin, we need to find the fractional bin indices
  iFreqBinindices_frac[] init lenarray(iFreqBinindices)
  index = 0
  while index < lenarray(iFreqBinindices)-1 do
    ifrac = (iDispfreqs[iFreqBinindices[index]]-iGridfreqs[index]) / (iDispfreqs[iFreqBinindices[index]]-iDispfreqs[iFreqBinindices[index]-1])
    iFreqBinindices_frac[index] = iFreqBinindices[index]-(ifrac^(1/12))
    index += 1
  od
  iFreqBinindices_frac[index] = iFreqBinindices[index]
  xout iFreqBinindices_frac
endop

opcode ArrayReduce, k[], k[]i[]k
  ; reduce an array :
  ; the split-indices divides the input array into regions,
  ; then we take the maximum value in each region and write to the output array
  kInamps[], iSplitindices[], kupdate xin
  insize = lenarray(kInamps)
  ioutsize = lenarray(iSplitindices)
  
  kOutamps[] init ioutsize
  if kupdate > 0 then
    kindex = 0
    kstartindex = 0
    while kindex < ioutsize do
      ; iterate through the input array, find the max amp value within each region defined by splitindices  and write that to the output bin
      kendindex = kindex > ioutsize-2 ? insize-1 : iSplitindices[kindex+1]
      kmax = 0
      while kstartindex <= kendindex-1 do
        kamp = kInamps[kstartindex]
        kmax max kmax, kamp
        kOutamps[kindex] = kmax
        kstartindex += 1
      od
      kindex += 1
    od
  endif
  xout kOutamps
endop

opcode svgGrid, 0, Si[]i[]iii
  ; make the grid for spectral display
  ; frequency lines vertically, at the frequencies set by iFreqs and the positions set by iPos
  ; amplitude lines at dB locations vertically
  SChannel, iFreqs[], iPos[],ilen, idBmin, inum_amplines xin
  iBounds[] cabbageGet SChannel, "bounds" 
  itopspace = 10 ; space at top for labels
  iheight = iBounds[3]-itopspace

  SPath init ""
  Slabels init ""
  
  ; amp lines
  index init 0
  while index < inum_amplines do
    SPath strcat SPath, sprintf({{
    <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:rgb(100,100,100);stroke-width:1" />
    }}, 0, itopspace+(iheight/inum_amplines)*index, iBounds[2], itopspace+(iheight/inum_amplines)*index)
    Samp sprintf "-%d dB", (idBmin/inum_amplines)*index
    Slabels strcat Slabels, sprintf({{<text x="%d" y="%d" fill="grey" font-size="10px" text-anchor="right">%s </text>}},
      1, 11+itopspace+(iheight/inum_amplines)*index, Samp)
    index +=1
  od
  SPath strcat SPath, sprintf({{
    <line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:rgb(100,100,100);stroke-width:1" />
    }}, 0, iBounds[3], iBounds[2], iBounds[3]); add bottom line
    
  ; frequency lines
  index = 0
  while index < lenarray(iFreqs) do
    SPath strcat SPath, sprintf({{<line x1="%d" y1="%d" x2="%d" y2="%d" style="stroke:rgb(100,100,100);stroke-width:1" />}}, iPos[index], itopspace, iPos[index], iBounds[3])
    if iFreqs[index] < 1000 then
      Sfreq sprintf "%d",  iFreqs[index]
    else
      Sfreq sprintf "%dk",  iFreqs[index]/1000
    endif
    Slabels strcat Slabels, sprintf({{<text x="%d" y="%d" fill="grey" font-size="10px" text-anchor="middle">%s </text>}}, iPos[index], 8, Sfreq)
    index +=1
  od
  SPath strcat SPath, Slabels
  cabbageSet SChannel, "svgElement", SPath
endop

opcode svgSpect, 0,Sk[]ii[]kSiii
  SChannel, kAmps[], ilen, iPos[], kupdate, Scolor, ilinewidth, ifillopacity, igradient xin
  iBounds[] cabbageGet SChannel, "bounds"
  kindex = 0 
  irange = sr/2
    
  if kupdate == 1 then
    SPath sprintfk {{
      <defs>
        <linearGradient id="Gradient1" x1="0" y1="0" x2="0" y2="1">
        <stop offset="0" stop-color="%s" stop-opacity="%f" />
        <stop offset="%f" stop-color="black" stop-opacity="0" />
        </linearGradient>
      </defs>}}, Scolor, ifillopacity, igradient
    SPath strcatk SPath, sprintfk("<path d=\"M0 %d", iBounds[3])  
    while kindex < ilen do
      kamp = kAmps[kindex]*iBounds[3]
      knextindex = kindex==lenarray(kAmps)-1?kindex:kindex+1
      kamp1 = (kAmps[kindex]-(kAmps[kindex]-kAmps[knextindex])*0.5)*iBounds[3]
      kx = iPos[kindex]
      kx1 = ((iPos[kindex]+iPos[knextindex])*0.5)
      ky = iBounds[3]-kamp
      ky1 = iBounds[3]-kamp1
      if kindex == 0 then
        SPath strcatk SPath, sprintfk(" %d %d ", kx, ky)      ; line
      else
        SPath strcatk SPath, sprintfk("Q %d %d %d %d ", kx, ky, kx1, ky1); curve
      endif
      kindex += 1
    od
    if igradient > 0 then
      SPath strcatk SPath, sprintfk("L%d %d \" style=\"fill:url(#Gradient1);fill-opacity:%f;stroke:%s;stroke-width:%f;\"/>", iBounds[2], iBounds[3], ifillopacity, Scolor, ilinewidth)
    else
      SPath strcatk SPath, sprintfk("L%d %d \" style=\"fill:%s;fill-opacity:%f;stroke:%s;stroke-width:%f;\"/>", iBounds[2], iBounds[3], Scolor, ifillopacity, Scolor, ilinewidth)
    endif
    cabbageSet 1, SChannel, "svgElement", SPath 
  endif
endop

opcode svgSignalDisplay, 0, SfkSiiii[]i[]iii
  Schannel, fsig, kupdate, Scolor, ifftsize, idisplen, idBmin, iSplitindices[], iPos[], ilinewidth, ifillopacity, igradient xin
  fsmooth pvsmooth fsig, 0.2, 0.2
  kAmps[] init ifftsize/2+1
  kPvsfreqs[] init ifftsize/2+1
  kframe  pvs2array kAmps, kPvsfreqs, fsmooth
  kDispAmps[] ArrayReduce kAmps, iSplitindices, kupdate
  kDispAmps += ampdbfs(-140) ; avoid overflow when input is silent
  ; convert to dB
  if kupdate > 0 then
  iscale = 0.8/(pow(idBmin/96,0.1)*96)
  
    kindex = 0
      while kindex < lenarray(kDispAmps)-1 do
      kDispAmps[kindex] = (dbfsamp(kDispAmps[kindex]))*iscale+1
      kindex += 1
  od
  endif
  svgSpect Schannel, kDispAmps, idisplen, iPos, kupdate, Scolor, ilinewidth, ifillopacity, igradient
endop

instr 1

  ; test med gui reset didplay params
  ; reinit grid etc
  SWidgetChannels[] cabbageGetWidgetChannels
  SChannel, ktrig cabbageChanged SWidgetChannels
  if ktrig > 0 then
    reinit settings
  endif
  settings:
  
  ifftsize = 2048
  ifft_binsize = (sr/ifftsize)
  iFftfreqs[] genarray ifft_binsize, sr/2, ifft_binsize
  iupdaterate chnget "updaterate"
  idBmin chnget "dBmin"
  Scolour1 chnget "colour1"
  Scolour2 chnget "colour2"
  inum_amplines = 5 ; number of amplitude lines in grid
  iresolution chnget "resolution" ; semitones
  iskew chnget "skew" ; frequency skew *after* we have remapped the linear fft array into semitone-spaced bins
  iskew = iskew*pow(tanh(iresolution*0.6),0.35)*0.9
  ilinewidth chnget "linewidth"
  ifillopacity chnget "fill"
  igradient chnget "gradient"
  iGridfreqs[] fillarray 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000

  iSplitindices[], iDispfreqs[] MakeSplitpoints ifftsize, ifft_binsize, iresolution 
  iBounds[] cabbageGet "displayGrid", "bounds"
  idispsize = iBounds[2]
  idisplen = lenarray(iSplitindices)
  iGridindices[] FindFreqBins iGridfreqs, iDispfreqs
  iGridpos[] = iGridindices*(idispsize/idisplen)
  iGridpos PowArray iGridpos, iskew
  svgGrid "displayGrid", iGridfreqs, iGridpos, idisplen, idBmin, inum_amplines
  ; signal display positions array
  iPos[] genarray 0, idisplen, 1
  iPos = iPos*(idispsize/idisplen)
  iPos PowArray iPos, iskew

  kupdate metro iupdaterate

    a1, a2 inch 1, 2
    ;anois rnd31 ampdbfs(-96), 1
    ;a1 = a1+anois
    ;a2 = a2+anois
  krms1 rms a1
  fsig1   pvsanal a1, ifftsize,ifftsize/2,ifftsize, 1
  svgSignalDisplay "displaySig1", fsig1, kupdate, Scolour1, ifftsize, idisplen, idBmin, iSplitindices, iPos, ilinewidth, ifillopacity, igradient
  krms2 rms a2
  fsig2   pvsanal a2, ifftsize,ifftsize/2,ifftsize, 1
  svgSignalDisplay "displaySig2", fsig2, kupdate, Scolour2, ifftsize, idisplen, idBmin, iSplitindices, iPos, ilinewidth, ifillopacity, igradient
      
endin

</CsInstruments>
<CsScore>
i1 0 [60*60*24*7] 
</CsScore>
</CsoundSynthesizer>
