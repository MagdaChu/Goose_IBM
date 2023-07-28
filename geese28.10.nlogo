; Author:      Magda Chudzinska, AU
; Date:        Created (started) 2013-09-24
; Description: Geese movement model.The step size is 1h, the patch size is 50x50m (0.05kmX0.05km)
; The real world is 28.5 by 33.5 km so 571 by 675 patches
; Model starts 06.04 at 01:00 and finishes no later than 19.05 24:00

; INITIALLY MODEL WAS SUPPOSE TO START 11.04 SINCE THIS IS THE DATA I HAVE HABITAT MAPS FROM BUT IN REALITY AT 11.04 THERE SHOULD ALREADY BE
; 15220 AT THE STUDY AREA SO I DECIDED TO START THE MODEL 6.04 WHEN THERE SHOULD BE 5130 GEESE AT THE AREA AND I ASSUME THAT THERE IS NOT MUCH SNOW
; AND FIELD HABITAT AS AS I STARTED MAPPING 11.04
extensions [gis]

globals[
  water
  bcg_fields
  roost-list
  fw1-list
  fw2-list
  fw3-list
  fw4-list
  xllcorner
  yllcorner
]

breed [fieldsw1 fieldw1] ; must be declared before geese so geese will be on top of fields
breed [fieldsw2 fieldw2]
breed [fieldsw3 fieldw3]
breed [fieldsw4 fieldw4]
breed [roosts roost]
breed [geese goose]
breed [suns sun]

geese-own[
  energy
]

suns-own [day?]

roosts-own [
  rx-coord
  ry-coord
  id
  capacity
  ]

fieldsw1-own[
  fw1x-coord
  fw1y-coord
  areaw1
  habw1
  initial_hab_values1
  id1]

fieldsw2-own[
  fw2x-coord
  fw2y-coord
  areaw2
  habw2
  id2
  initial_hab_values2]

fieldsw3-own[
  fw3x-coord
  fw3y-coord
  areaw3
  habw3
  id3
  initial_hab_values3]

fieldsw4-own[
  fw4x-coord
  fw4y-coord
  areaw4
  habw4
  id4
  initial_hab_values4]

patches-own [background]

to load-background
  ;; this layer is just for an overview, it is not going to influence geese movement
  ;; (for this model to work with NetLogo's new plotting features,
  ;; __clear-all-and-reset-ticks should be replaced with clear-all at
  ;; the beginning of your setup procedure and reset-ticks at the end
  ;; of the procedure.)
  clear-all
  ; Load datasets
  set water gis:load-dataset "ascii/water_ascii.asc" ; I chose to use ascii files because I can define lower left coordinated of the world
  set bcg_fields gis:load-dataset "ascii/fields_ascii.asc"
  ; Set the world envelope to the union of dataset's envelopes
  gis:set-world-envelope (gis:envelope-union-of
                                                (gis:envelope-of water)
                                               (gis:envelope-of bcg_fields)
                                               )
  set xllcorner 605230.32995747 ; define the lower left corner (from ascii files)
  set yllcorner 7071665.4296497


    ask patches
  [
    set pcolor white
    set background "nothing"
    let value_w gis:raster-sample water patch pxcor pycor
    if (value_w > -9999)
    [set pcolor cyan
    set background "water"]
    let value_f gis:raster-sample bcg_fields patch pxcor pycor
    if (value_f > -9999)
    [set pcolor gray
    set background "patches"]
  ]



 reset-ticks
end



to load-roosts-from-txt

let roost_pos "input/roost_utm.txt"
let line-txt ""
let line-lst (list " " " " " " " ") ;as many as columns in my text file
set roost-list (list line-lst)
let line-lgt 0
file-open roost_pos
let i 0 ;characters in the text file
let j 0 ; words in the text file
let k 0 ; line number in the text file
let one-char ""
let one-word ""
let in-word true
    while [ not file-at-end? ] [
      set line-txt file-read-line
      set line-txt word line-txt " "
      set line-lgt length line-txt
      set i 0
      set j 0
      set in-word true
      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
        set one-char substring line-txt i (i + 1)
        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word) [  ; that is, prev. char was part of a word, current isn't; if a character after is n/ delimited or tab delimited or space delimited than is not a word
          set in-word false
          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
          set line-lst replace-item j line-lst one-word
          set j j + 1
        ]
        if ((one-char = "\t" or one-char = " ") and (not in-word)) [
          set one-word ""
        ]
        if (one-char != "\t" and one-char != " ") [
          set in-word true
          set one-word word one-word one-char
        ]
        set i i + 1
      ]
      if (k > 0) [ set roost-list lput line-lst roost-list ]
      set k k + 1
    ]
    file-close
end

to import_roosts

  let i 0
  if (not (length roost-list < 2)) [ ; First line contains header
    set roost-list remove-item 0 roost-list ; remove header
    let nb-of-roosts length roost-list
    create-roosts nb-of-roosts

      ask roosts
      [
      set id (item 0 (item i roost-list))
      set rx-coord ((item 1 (item i roost-list)) - xllcorner) * 0.020   ; comes from: the real environment is 28500m and there is 571 patches (571/28500 = 0.02)
      set ry-coord ((item 2 (item i roost-list)) - yllcorner) * 0.020   ; 675/33500 = 0.02
      set capacity item 3 (item i roost-list)
      setxy rx-coord ry-coord
      set i i + 1
      set color 14
      set shape "square"
      set size 8
      ]
       ]
  print count roosts
end

to load_fields_w1 ; from 11.04 00:00 till 18.04 23:59

let fields_w1_pos "input/habitat_eol_mac_16.12.txt" ; This is a new complete habitat map created 28.10 based on R script (details in making new habitat map folder)
; to the file I added initial value of the habitats as done in R script 'initial habitat values'
; than file has to be saved as Mac-type txt file as done in the same R scipt
let line-txt ""
let line-lst (list " " " " " " " " " " " " " " " " " " " " " " " " " ") ;as many as columns in my text file
set fw1-list (list line-lst)
let line-lgt 0
file-open fields_w1_pos
let i 0 ;characters in the text file
let j 0 ; words in the text file
let k 0 ; line number in the text file
let one-char ""
let one-word ""
let in-word true
    while [ not file-at-end? ] [
      set line-txt file-read-line
      set line-txt word line-txt " "
      set line-lgt length line-txt
      set i 0
      set j 0
      set in-word true
      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
        set one-char substring line-txt i (i + 1)
        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word) [  ; that is, prev. char was part of a word, current isn't; if a character after is n/ delimited or tab delimited or space delimited than is not a word
          set in-word false
          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
          set line-lst replace-item j line-lst one-word
          set j j + 1
        ]
        if ((one-char = "\t" or one-char = " ") and (not in-word)) [
          set one-word ""
        ]
        if (one-char != "\t" and one-char != " ") [
          set in-word true
          set one-word word one-word one-char
        ]
        set i i + 1
      ]
      if (k > 0) [ set fw1-list lput line-lst fw1-list ]
      set k k + 1
    ]
    file-close
end

to import_fields_w1

  let i 0
  if (not (length fw1-list < 2)) [ ; First line contains header
    set fw1-list remove-item 0 fw1-list ; remove header
    let nb-of-fields-w1 length fw1-list
    create-fieldsw1 nb-of-fields-w1

      ask fieldsw1
      [
      set habw1 (item 2 (item i fw1-list))
      set id1 (item 0 (item i fw1-list))
      set fw1x-coord ((item 6 (item i fw1-list)) - xllcorner) * 0.020   ; comes from: the real environment is 28500m and there is 571 patches (571/28500 = 0.02)
      set fw1y-coord ((item 7 (item i fw1-list)) - yllcorner) * 0.020   ; 675/33500 = 0.02
      set areaw1 item 1 (item i fw1-list) ; area is in m2
      set initial_hab_values1 (item 9 (item i fw1-list))
      setxy fw1x-coord fw1y-coord
      set i i + 1
      if habw1 = "grass" [set color green]
      if habw1 = "grain" [set color orange]
      if habw1 = "other" [set color 3]
      if habw1 = "plough" [set color black]
      if habw1 = "potato" [set color brown]
      if habw1 = "stubble" [set color yellow]
      set shape "circle"
      set size 4
      ]
       ]
  ;;;;;; assigning initial crop values for each habitat  ;;;;;
  ;;;;;; all detail in ODD ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Initial values must be in the input file so there is one less source of stochasticity ;;


  print count fieldsw1
end

to load_fields_w2 ;from 19.04 00:00 till 26.04 23:59

let fields_w2_pos "input/habitat_eol_mac_16.12.txt"
let line-txt ""
let line-lst (list " " " " " " " " " " " " " " " " " " " " " " " " " ") ;as many as columns in my text file
set fw2-list (list line-lst)
let line-lgt 0
file-open fields_w2_pos
let i 0 ;characters in the text file
let j 0 ; words in the text file
let k 0 ; line number in the text file
let one-char ""
let one-word ""
let in-word true
    while [ not file-at-end? ] [
      set line-txt file-read-line
      set line-txt word line-txt " "
      set line-lgt length line-txt
      set i 0
      set j 0
      set in-word true
      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
        set one-char substring line-txt i (i + 1)
        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word) [  ; that is, prev. char was part of a word, current isn't; if a character after is n/ delimited or tab delimited or space delimited than is not a word
          set in-word false
          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
          set line-lst replace-item j line-lst one-word
          set j j + 1
        ]
        if ((one-char = "\t" or one-char = " ") and (not in-word)) [
          set one-word ""
        ]
        if (one-char != "\t" and one-char != " ") [
          set in-word true
          set one-word word one-word one-char
        ]
        set i i + 1
      ]
      if (k > 0) [ set fw2-list lput line-lst fw2-list ]
      set k k + 1
    ]
    file-close
end

to import_fields_w2

  let i 0
  if (not (length fw2-list < 2)) [ ; First line contains header
    set fw2-list remove-item 0 fw2-list ; remove header
    let nb-of-fields-w2 length fw2-list
    create-fieldsw2 nb-of-fields-w2

      ask fieldsw2
      [
      set habw2 (item 3 (item i fw2-list))
      set fw2x-coord ((item 6 (item i fw2-list)) - xllcorner) * 0.020   ; comes from: the real environment is 28500m and there is 571 patches (571/28500 = 0.02)
      set fw2y-coord ((item 7 (item i fw2-list)) - yllcorner) * 0.020   ; 675/33500 = 0.02
      set areaw2 item 1 (item i fw2-list)
      set id2 (item 0 (item i fw2-list))
      set initial_hab_values2 (item 10 (item i fw2-list))
      setxy fw2x-coord fw2y-coord
      set i i + 1
      if habw2 = "grass" [set color green]
      if habw2 = "grain" [set color orange]
      if habw2 = "other" [set color 3]
      if habw2 = "plough" [set color black]
      if habw2 = "potato" [set color brown]
      if habw2 = "stubble" [set color yellow]
      set shape "circle"
      set size 4
      ]
       ]
  print count fieldsw2
end

to load_fields_w3 ;from 27.04 00:00 till 04.05 23:59

let fields_w3_pos "input/habitat_eol_mac_16.12.txt"
let line-txt ""
let line-lst (list " " " " " " " " " " " " " " " " " " " " " " " " " ") ;as many as columns in my text file
set fw3-list (list line-lst)
let line-lgt 0
file-open fields_w3_pos
let i 0 ;characters in the text file
let j 0 ; words in the text file
let k 0 ; line number in the text file
let one-char ""
let one-word ""
let in-word true
    while [ not file-at-end? ] [
      set line-txt file-read-line
      set line-txt word line-txt " "
      set line-lgt length line-txt
      set i 0
      set j 0
      set in-word true
      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
        set one-char substring line-txt i (i + 1)
        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word) [  ; that is, prev. char was part of a word, current isn't; if a character after is n/ delimited or tab delimited or space delimited than is not a word
          set in-word false
          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
          set line-lst replace-item j line-lst one-word
          set j j + 1
        ]
        if ((one-char = "\t" or one-char = " ") and (not in-word)) [
          set one-word ""
        ]
        if (one-char != "\t" and one-char != " ") [
          set in-word true
          set one-word word one-word one-char
        ]
        set i i + 1
      ]
      if (k > 0) [ set fw3-list lput line-lst fw3-list ]
      set k k + 1
    ]
    file-close
end

to import_fields_w3

  let i 0
  if (not (length fw3-list < 2)) [ ; First line contains header
    set fw3-list remove-item 0 fw3-list ; remove header
    let nb-of-fields-w3 length fw3-list
    create-fieldsw3 nb-of-fields-w3

      ask fieldsw3
      [
      set habw3 (item 4 (item i fw3-list))
      set fw3x-coord ((item 6 (item i fw3-list)) - xllcorner) * 0.020   ; comes from: the real environment is 28500m and there is 571 patches (571/28500 = 0.02)
      set fw3y-coord ((item 7 (item i fw3-list)) - yllcorner) * 0.020   ; 675/33500 = 0.02
      set areaw3 item 1 (item i fw3-list)
      set id3 (item 0 (item i fw3-list))
       set initial_hab_values3 (item 11 (item i fw3-list))
      setxy fw3x-coord fw3y-coord
      set i i + 1
      if habw3 = "grass" [set color green]
      if habw3 = "grain" [set color orange]
      if habw3 = "other" [set color 3]
      if habw3 = "plough" [set color black]
      if habw3 = "potato" [set color brown]
      if habw3 = "stubble" [set color yellow]
      set shape "circle"
      set size 4
      ]
       ]
  print count fieldsw3
end

to load_fields_w4 ;from 05.05 00:00 till 14.05 23:59

let fields_w4_pos "input/habitat_eol_mac_16.12.txt"
let line-txt ""
let line-lst (list " " " " " " " " " " " " " " " " " " " " " " " " " ") ;as many as columns in my text file
set fw4-list (list line-lst)
let line-lgt 0
file-open fields_w4_pos
let i 0 ;characters in the text file
let j 0 ; words in the text file
let k 0 ; line number in the text file
let one-char ""
let one-word ""
let in-word true
    while [ not file-at-end? ] [
      set line-txt file-read-line
      set line-txt word line-txt " "
      set line-lgt length line-txt
      set i 0
      set j 0
      set in-word true
      while [ i < line-lgt and k > 0 ] [  ; first line (k = 0) contains the header, must be omitted
        set one-char substring line-txt i (i + 1)
        if ((one-char = "\n" or one-char = "\t" or one-char = " ") and in-word) [  ; that is, prev. char was part of a word, current isn't; if a character after is n/ delimited or tab delimited or space delimited than is not a word
          set in-word false
          if (j > 0) [ set one-word read-from-string one-word ] ; cast to numeric
          set line-lst replace-item j line-lst one-word
          set j j + 1
        ]
        if ((one-char = "\t" or one-char = " ") and (not in-word)) [
          set one-word ""
        ]
        if (one-char != "\t" and one-char != " ") [
          set in-word true
          set one-word word one-word one-char
        ]
        set i i + 1
      ]
      if (k > 0) [ set fw4-list lput line-lst fw4-list ]
      set k k + 1
    ]
    file-close
end

to import_fields_w4

  let i 0
  if (not (length fw4-list < 2)) [ ; First line contains header
    set fw4-list remove-item 0 fw4-list ; remove header
    let nb-of-fields-w4 length fw4-list
    create-fieldsw4 nb-of-fields-w4

      ask fieldsw4
      [
      set habw4 (item 5 (item i fw4-list))
      set fw4x-coord ((item 6 (item i fw4-list)) - xllcorner) * 0.020   ; comes from: the real environment is 28500m and there is 571 patches (571/28500 = 0.02)
      set fw4y-coord ((item 7 (item i fw4-list)) - yllcorner) * 0.020   ; 675/33500 = 0.02
      set areaw4 item 1 (item i fw4-list)
      set id4 (item 0 (item i fw4-list))
       set initial_hab_values4 (item 12 (item i fw4-list))
      setxy fw4x-coord fw4y-coord
      set i i + 1
      if habw4 = "grass" [set color green]
      if habw4 = "grain" [set color orange]
      if habw4 = "other" [set color 3]
      if habw4 = "plough" [set color black]
      if habw4 = "potato" [set color brown]
      if habw4 = "stubble" [set color yellow]
      set shape "circle"
      set size 4
      ]
       ]
  print count fieldsw4
end

to starting_setups

  load-background
  load-roosts-from-txt
  import_roosts
  load_fields_w1
  import_fields_w1
  create_sun



end

to create_sun
  create-suns 1
  [
    set color black
    set size 40
    set shape "circle"
    setxy  530 651
    ;change appearance based on whether it is day or night
  ]
end

to create_geese ; for the moment a fixed amount of geese arrives to the model and no new geese arrive during the model.The proper arrival dynamic I will model later
   create-geese initial-number-of-geese [
    set color pink
    set energy (random-float 6214) + 9321
    ; an arriving goose should have an API asign between 2 +/- 0.5 API which is 12428 +/- 3107 kJ
    ;(that does not take fat synthesis efficiency into account and it shouldnt)
    ; the number should be than from 9321 and 15535
    ; from now on I will only keep energy in kJ in the model
    set size 5
    setxy  0 0 ; starting point does not matter cause I tell geese to move to a random roost immidietaly
    move-to one-of roosts
    ]

end

;to geese-ariving ; geese can arrive evry 3-4 days
  ; random-normal 55403 / 20 13606 / 20 ; mean and sd for 2012
;end

to show-day
  ; time of civil twilights are in 'civil twilight.xlsx'
  if ticks > 4 [ask turtles with [breed = suns]
  [set color yellow
  set day? TRUE]]

  if ticks > 20 [ask turtles with [breed = suns]
  [set color black
  set day? FALSE]]


end

to go

  tick ; otherwise there is no ticks
  grass_growth
  if ticks = 20 * 24 [ask turtles with [breed = fieldsw1] [die]]
  if ticks = 20 * 24 + 1 [load_fields_w2 import_fields_w2] ; the new fields appear 26.04 at 01:00
  if ticks = 28 * 24 [ask turtles with [breed = fieldsw2] [die]]
  if ticks = 28 * 24 + 1 [load_fields_w3 import_fields_w3] ;the new fields appear 03.05 at 01:00
  if ticks = 36 * 24 [ask turtles with [breed = fieldsw3] [die]]
  if ticks = 36 * 24 + 1 [load_fields_w4 import_fields_w4] ;the new fields appear 11.05 at 01:00
  if ticks = 44 * 24 [stop] ; the model runs from 06.04 and ends 19.05 the latest
  show-day

end

to grass_growth ; grass is growing everyday irrespective of the grazing preasssure
  ; grass is growing 0.06 cm every day during periods 1-3 and 0.47 cm/day in period 4. Conversion into biomass is done based on equetion provided in ODD


  ask fieldsw1 [
    if habw1 = "grass" [set initial_hab_values1 initial_hab_values1 + (((0.06 / 24) * 16.4) * areaw1)]

    ]

  ask fieldsw2 [
    ; set prev-crop-value pres-crop-value
    ]

;it would never work because you intend to compare id1 of this turtle with id2 of EVERY OTHER TURTLE...
;this is not what you are doing at all... and it is not what you would want to do.

;Maybe you need to use HATCH when you create the second set of fields, in the first place.
;HATCH causes the turtle to create an exact copy of itself (except for WHO number).
;So you use HATCH to create the copy, then reset any properties that should NOT be exactly the same.
end
@#$#@#$#@
GRAPHICS-WINDOW
290
10
862
695
-1
-1
0.7
1
10
1
1
1
0
1
1
1
0
571
0
675
0
0
1
ticks
30.0

BUTTON
906
11
988
44
clear w1
ask turtles with [breed = fieldsw1] [die]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
66
49
121
82
p2
load_fields_w2\nimport_fields_w2
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
907
54
985
87
clear w2
ask turtles with [breed = fieldsw2] [die]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
126
49
181
82
p3
load_fields_w3\nimport_fields_w3
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
186
49
241
82
p4
load_fields_w4\nimport_fields_w4
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
907
95
985
128
clear w3
ask turtles with [breed = fieldsw3] [die]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
907
136
985
169
clear w4
ask turtles with [breed = fieldsw4] [die]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
12
100
75
133
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
6
10
173
43
background, roosts & p1
starting_setups
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
92
103
200
136
NIL
create_geese
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
15
146
193
179
initial-number-of-geese
initial-number-of-geese
10
100
30.0
20
1
NIL
HORIZONTAL

BUTTON
7
49
62
82
p1
load_fields_w1\nimport_fields_w1
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
