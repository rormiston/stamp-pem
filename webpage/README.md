# Plan for webpage

## Homepage:

### Calendar:
* Points to specific days

### Description of Analysis:
* Coherence calculation

```latex
$\rho(f) = \frac{|s_1(f)*\times s_2(f)|^2}{P_1(f)P_2(f)}$
$P_I(f) = |s_I(f)|^2.$
```

* Coherence SNR calculation

```latex
SNR = $N\rho(f)$
$N$ = number of segments averaged
Note: $\frac{1}{N} is theoretical $\rho(f)$
```
* Descriptions of individual plots

```
Coherence Matrix: 
    Colors are coherence, x-axis is frequency, 
    y-axis is channels. Nice for looking at time-integrated
    coherence for many channels at once. 
    Grouped by subsystems.

Coherence Spectrogram:
    Colors are coherence, x-axis is time, 
    y-axis is frequency, made for each individual 
    channel. Nice for looking at time-progression 
    of coherence for single channel. Good for 
    checking for wandering lines.
```
* Where to find lists of channels
* Which subsystems where chosen

### Links
* Code for stamp_pem
* Code for gwpy
* CIS
* Summary pages

## Daily Main Page:

### Plots:
* Dynamic opening for coherence matrices for subsystems

### Links & Misc.:
* Link to each channel's page for coherence spectrogram for all locks from that day
    * Organized by subsystems on page
* Segment information for that day
* Links to coherence matrices for each lock.

## Individual Daily Channel Page:
### Plots:
* Individual coherence spectrogram for all locks from that day.
* (?) PSDs for that channel for that day (see if there's interest)

### Links:
* CIS link for that channel

## Navigation Bar:

### Links:
* Homepage/Top level calendar (may not be needed if there's a date-picker)
* Date-picker (see bootstrap page) that takes you to similar pages for next day
* Link to coherence matrices
* Link to code?
* Daily landing page
