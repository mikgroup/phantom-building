M104 S200 ;Print temperatures
M109 S200

G21        ;metric values
G90        ;absolute positioning
M82        ;set extruder to absolute mode
M107       ;start with the fan off
G32		; Define bed plane
G92 Z3.95  ; Measured value of nozzle offset
G1 Z5.0 F1200 ;move up

G1 X0 Y0 Z0.3 F3000;move the platform to purge extrusion
G92 E0 ;zero the extruded length
G1 F200 X100 E30 ;extrude 30mm of feed stock
G92 E0 ;zero the extruded length again
G0 X0 Y0 F3000
G1 F200 Y100 E30 ;draw 0,0 cross
G92 E0

G1 X100 Y100 Z5 F1200 ; Move to arbitrary position

M600 ; pause print, attach acrylic
G1 X10 Y10 Z3.22 F3000; acrylic thickness + 0.3 (CHANGE!!)
G92 E0 ;zero the extruded length
G1 F200 X50 E15 ;extrude 30mm of feed stock
G92 E0 ;zero again

M117 Printing...
