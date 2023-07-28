; Author:      Magda Chudzinska, AU
; Date:        Created (started) 2013-09-24
; Description: Geese movement model.The step size is 1h, the patch size is 50x50m (0.05kmX0.05km)
; The real world is 28.5 by 33.5 km so 571 by 675 patches
; Model starts 06.04 at 00:00 and finishes no later than 19.05 24:00

; INITIALLY MODEL WAS SUPPOSE TO START 11.04 SINCE THIS IS THE DATA I HAVE HABITAT MAPS FROM BUT IN REALITY AT 11.04 THERE SHOULD ALREADY BE 
; 17380 AT THE STUDY AREA SO I DECIDED TO START THE MODEL 6.04 WHEN THERE SHOULD BE 7040 GEESE AT THE AREA AND I ASSUME THAT THERE IS NOT MUCH SNOW 
; AND FIELD HABITAT AS AS I STARTED MAPPING 11.04
; calculations in 'arrival dates.r'

;; PERIODS

; 1 - 6.04 - 25.04
; 2 - 26.04 - 03.05
; 3 - 04.05 - 11.05
; 4 - 12.05 - 19.05

extensions [gis array matrix]

globals[
  water                      ; this variable is just a background, it does not influence the goose movement/behaviour
  bcg_fields                 ; this variable is just a background, it does not influence the goose movement/behaviour
  roost-list                 ; roost list created together with Jesper during MIGRAPOP meeting in 2013
  fw-list                    ; list of all fields
  xllcorner                  ; x coordinates of lower left corner of the GIS file
  yllcorner                  ; y coordinates of lower left corner of the GIS file
  daynight-list              ; list of ticks assigned with day or night. Day and night is set as a civil twilight found at http://aa.usno.navy.mil/data/docs/RS_OneYear.php and set for year 2012, Steinkjer (11 29' and 64 0') with UTM + 2
  day                        ; a list as long as number of ticks, for each tick is says 'D' if it is a day and 'N' if it is a night
  patch_size                 ; I need to know what is a patch size in order to be able to calculate the distance between two points
  tick_list                  ; list of ticks when new geese should arrive to the model
  arrival_list               ; number of goose agents arriving to the model according to tick_list
  c                          ; counter
  d                          ; counter
  n-of-starved               ; number of geese which left the model due to starvation
  n-of-vesteralen            ; number of geese which left the model becuase they reached a desired API to move north
  gl_ticks_since_morning     ; 1= first step after civil twilight
  gl_max_num_geese_per_roost ; one list showing the maximum number of geese visiting each roost during the entire model and the id of that roost
  gl_mean_num_geese_per_roost ; one list showing the mean number of geese visiting each roost during the entire model and the id of that roost. this is mean of daily maximum number of visiting geese 
  eg_reduction_dist          ; the amount of metabolised energy intake which is substructed from a goose after disturbance. this is to mimic a situation when geese are disturbed within a time step and are not able to forage throughout the entire step
  n_geese_large_fields       ; number of geese visiting large fields (defined by field_area_threshold)
  n_geese_small_fields       ; number of geese visiting small fields (defined by field_area_threshold)
  n_disturbance_event        ; number of time geese were disturbed
  
  ;;; these two below parameters are only used for the analysis of the results, they do not take part or influence any processes in the model
  period                     ; defines which period of the stopover season it is (ticks 0-479 - period 1; 480-671 - period 2; 672-863 period 3; 864-1057 period 4)
  time_of_day                ; defines which time of day it is based on the 'gl_ticks_since_morning'. If 'gl_ticks_since_morning' is equal to a number as in the below table then time of the day is either morning, midday or evening (everything else is night):
  
  ;;;            morning   midday  evening
  ;;; period 1    1-6       7-11     12-16
  ;;; period 2    2-7       8-12     13-17
  ;;; period 3    3-8       9-13     14-18
  ;;; period 4    4-9       10-14    15-19
  
  ;;; global values which may be parameterised/or in sensitivity analyses
  ;;; they are all moved to the interface
  
;  n_geese_in_super_goose ; number of individuals consisting one super-goose
;  big_radius             ; radius within which geese move after leaving a night roost, when leave a field due to rule 2 and 4 and some aspects of rule 3
;  small_radius           ; radius within which geese move when leave a field due to rule 1 and some aspects of rule 3
;  roost_dist_radius      ; a radius defining how close to a roost geese should be in order to fly to that roost after disturbance 
;  cum_time_threshold     ; time after which geese move to a roost
;  alpha                  ; memory factor (equation 5)
;  field_area_threshold   ; a field size [m2] over which goose behaviour differes with disturbance
;  av_ee                  ; average energy expenditure (kJ/h) no matter if a goose is on a field or on a roost
;  ee_fly                 ; average energy expenditure of one hour of flying (kJ/h)
;  starvation_energy_stores         ; goose energy stores (kJ) below which geese starve and leave the model
;  prob-of-disturbance-small-fields ; probability a goose agent will be disturbed on small fields
;  prob-of-disturbance-large-fields ; probability a goose agent will be disturbed on small fields 
;  extra_initial_energy   ; energy value (kJ) which will increase or decrease the intial energy stores when geese arrive, used for sensitivity analysis
;  extra_leaving_energy   ; energy value (kJ) which will increase or decrease the energy stores threshold which defines when geese leave the model to migrate north
;  max_steps-to_stay_on_roost ; maximum number of steps a goose can stay on a day roost 8minimum is always set to 2) 
;  temp_disturbance       ; a fixed value added to the probability of disturbance on small and large fields - used in sensitivity analysis 
 
  ; initialisation
  
;  initial_grass_length          ; intial grass length of grass fields which  haven't been rewosn this year (cm)
;  initial_grain_density_grain   ; intial grain density on grain fields (g/m2)
;  grain_mean_stubble            ; mean value of grain on stubble fields (g/m2)
;  grain_sd_stubble              ; stadard deviation of number of grain on stubble field (g/m2)
  initial_grain_density_stubble ; intial grain density on stubble fields (g/m2)
  
  ; metabolised energy intake - equation 5 in ODD. These values won't be parameterised, they are based on empirical data
  ; however if one day a new, better data appear i can change it. These parameters will be hard coded and not put in the interface of the model
  
  gf_grass             ; ash-free energy content of grass (kJ of (g org DW))
  gf_grain             ; ash-free energy content of grain (kJ of (g org DW))
  gf_stubble           ; ash-free energy content of stubble (kJ of (g org DW))
  dr_grass             ; dropping rate of grass (per hour)
  dr_grain             ; dropping rate of grain (per hour)
  dr_stubble           ; dropping rate of stubble (per hour)
  gd_grass             ; ash-free energy content of grass droppings (kJ of (g org DW))
  gd_grain             ; ash-free energy content of grain droppings (kJ of (g org DW))
  gd_stubble           ; ash-free energy content of stubble droppings (kJ of (g org DW))
  m_grass              ; average ash-free mas of one grass dropping (g)
  m_grain              ; average ash-free mas of one grain dropping (g)
  m_stubble            ; average ash-free mas of one stubble dropping (g)
  mei_potato           ; metabolisable enrgy intake on potato fields (kJ/h)
  
  ; habitat depletion - equation 3 and 4 in ODD. data are based on the studies by Baveco and Nolet (unpublished data) where they empricially measured how geese firage in response to food biomass
  ; some information are in Baveco et al (2001), Ecological modelling 222
  ; hard coded as well as explained above
  
  b1                   ; b1-b3 are regression coefficients
  b2
  b3
  tc                   ; cropping time (s)
  R_max                ; maximum rate of chewing (g/s)
  a_grain              ; attack rate on grain (m2/h)  
  a_stubble            ; attack rate on stubble (m2/h)               
  h_stubble            ; handling time on stubble (s/g)
  h_grain              ; handling time on grain (s/g)  
  
  
  ;; Matrix just for calibration
  calibration-data
]

breed [norths north]    ; geese which leave the model because they reached a desired mass go there, used for DEBUGGING
breed [starvs starv]    ; geese which leave the model because they dropped to API 0 go there, used for DEBUGGING
breed [fields field]    ; must be declared before geese so geese will be on top of fields     
breed [roosts roost]
breed [suns sun]        ; this is just to show day and night, this breed does not influence goose behaviour/movement
breed [geese goose]


geese-own [tot_nei              ; total net energy intake over the entire model
           dnei                 ; daily net energy intake calculated from morning civil twilight till the end of the night
           dee                  ; daily energy expenditure 
           deg                  ; daily metabolised energy intake rate 
           step_ee              ; energy expenditure per step EXCLUDING energy for flying  
           step_eg              ; metabolised energy intake per step
           step_nei             ; net energy intake per step EXCLUDING flight costs, I need it after scaring events, than the total nei is reduced by a certain amount if geese are disturbed
           tot_ee               ; total energy expenditure over the entire model
           tot_eg               ; total energy gain during over the entire model
           ir                   ; intake rate on a given field
           pr_feed              ; proportion of time spent on feeding per stime step
           time_fly             ; time necessary for flying from the previous location (distance flown (km) / 50)
           en_fly               ; energy required for flying (calculated as 416.35 (kJ/h) * distance flown (km) / 50 since we assume average speed of 50 km/h)
           ticks_since_morning  ; 1= first step after civil twilight, same as the global variable gl_ticks_since_morning but for debugging reason it is easier for me to keep is as goose variable as well
           cum_time_feeding     ; cumulative time spent on feeding since leaving a morning roost and going back to another roost
           cum_tot_time_feeding ; cumulative time spent on feeding since leaving a morning roost, not set to zero when a goose gets to a day roost
           previous_loc         ; a previously visited field or a roost
           night_roost_location ; a location of the night roost, for DEBUGGING
           tot_dist_flown       ; the total distance a goose flown during the model (km)
           tot_en_fly           ; the total energy expenditures spent on flying during the entire model
           in_energy_stores     ; energy stores just after arrival to Trondelag (API 2 +/- 0.5; 12428 ± 3107 kJ)
           steps_to_stay_on_roost ; number of time steps a goose has to stay on roost after going back to a day roost (this number differs depending what caused geese to move to such a roost)
           day_roost_500m       ; reporter saying whether a goose is on a day roost because is there after being scared with 300 m rules or other
           left_model           ; if goose left the model it has this value updated
           ;moved?               ; indicates whether was goose was aready dragged by another goose or not
           mei_grass            ; variable used to calculate intitial gamma (initial expected gain rate) in calculate_initial_gamma procedure
           mei_stubble          ; variable used to calculate intitial gamma (initial expected gain rate) in calculate_initial_gamma procedure
           mei_ploughed         ; variable used to calculate intitial gamma (initial expected gain rate) in calculate_initial_gamma procedure
           mei_grain            ; variable used to calculate intitial gamma (initial expected gain rate) in calculate_initial_gamma procedure
           mei_pot              ; variable used to calculate intitial gamma (initial expected gain rate) in calculate_initial_gamma procedure
           gamma                ; expected gain rate
           moving_north_threshold ; since a mass threshold after which geese move to North-Norway is not a fixed value but a range (3122 +/- 90g), each goose has this value defined at the beginning of the model
           departure_date       ; day of departure based on GDD
           steps_spent_on-day_roost ; sum of steps a goose spent on a day roost. it is set to zero at the beginning of each day. It is used to calculte how much time during daylight, a goose spend on roost
           dist_closest_roost   ; distance to the closest roost, used for plotting output (plot_distance_to_roost)
           disturbed?           ; after a disturbance event is set to true, otherwise false
           on_roost?            ; if geese are on roost it is set to true, neccessary for calculating distance to the closest roost (calculate dist_to_roost procedure)
           id_field_with_max_grass; a who number of a potential grass field to which a goose may move if this is a field with highest mei within small_radius radius (applies to DMR2)
           id_field_with_max_stubble; a who number of a potential stubble field to which a goose may move if this is a field with highest mei within small_radius radius (applies to DMR2)
           id_field_with_max_grain; a who number of a potential grain field to which a goose may move if this is a field with highest mei within small_radius radius (applies to DMR2)
           id_field_with_max_potato; a who number of a potential potato field to which a goose may move if this is a field with highest mei within small_radius radius (applies to DMR2)
           mei_list             ; a list to which each goose stores a mei obtained at a given field and a who of that field. List is reset at the end of each day. The list is used in DMR3

]

roosts-own 
[
  rx-coord              ; x coordinate of roost
  ry-coord              ; x coordinate of roost
  rid                   ; id of a roost as set during MIGRAPOP meeting
  capacity              ; the maximum number of geese ever observed on a roost during counst on years 2005-2013, carrying capacity of roosts is not modelled in that model but I keep it for a future
  hab                   ; habitat on roost will always be roost
  num_geese_step        ; a list showing number of goose agents which visit a given roost per step (only if at least one goose is on that roost so there is no 0 in that list). the list is reset at the beginning of each day
  max_geese_day         ; a list showing the maximum number of geese (not cumulative) which visited a given roost during day
  ]

fields-own[
  fx-coord           ; x coordinate of a field, fields are represented as the points located in the middle of each field
  fy-coord           ; y coordinate of a field, fields are represented as the points located in the middle of each field
  area               ; the field area of a field, does not change from period to period [m2]
  hab                ; habitat type on a field in a givn moment
  hab1               ; habitat type on a field in period 1
  hab2               ; habitat type on a field in period 2
  hab3               ; habitat type on a field in period 3
  hab4               ; habitat type on a field in period 4
  hab_biomass         ; current value (in a given tick) of habitat biomass  [g/field]
  id_gis             ; id of the field as in GIS file, id does not change between the periods
  id                 ; id just as a consecutive numbers 
  gr_length          ; shows the grass length in a current moment in cm, for debugging only
  grain_dens         ; shows the grain density [number of grains/m2], for debugging only
  dist_prob          ; in each step each field has this value assign as random-float 1. if this value is lower than the probability of disturbance, geese on this field leave this field
  occupied?          ; yes if there is at least one goose agent on that field, no if a patch is empty
]


patches-own [background]


;;;;;;;;;;;;;;; THESE ARE INITIAL MODEL VALUES WHEN THE MODEL STARTS (HABITAT BIOMASS; GEESE NUMBER AND POSITIONS ;;;;;;;;;;;;;;;;;;;;;;;;;;


to setups
  ; random-seed 3 ;it means that every time I run the simUlAtion, it will choose the same random number
  ; it is sufficient to set random-seed once, just here and if I want certain procedures to have different seed or no seed specified I have to write it that procedure
  ca ; clear all so pressing setup clears everything
  load_background
  read_roosts_from_txt
  load_roosts
  read_fields_from_txt
  load_fields
  read_daynight_from_txt


  set gl_ticks_since_morning 20 
  set period 1
  set time_of_day "night"
  
  ; these are the lists showing at which tick, what number of goose agents has to arrive to the model. Geese only arrive until 30.04 as desribed in ODD
  set c 1 ; neccesary to run a list on geese-arrival
  set d 1 ; neccessary to update gdd and days since Denmark
  set tick_list [0 25 49 73 97 121 145 169 193 217 241 265 289 313 337 361 385 409 433 457 481 505 529 553 577 ] ; I had to add an extra number at the end of each list, otherwise model would stop
  set arrival_list [7040 1954 2028 2084 2124 2146 2151 2142 2117 2077 2023 1955 1873 1780 1672 1554 1424 1283 1132  971  800  620  432  235   30 ]
  ;set arrival_list [352 98 101 104 106 108 107 107 106 104 101 98 94 89 83 78 71 64 57 48 40 31 22 12 1 0]
  
  set patch_size 0.05 ; patch size in km (length of one edge of a patch)
  
  ; setting global values which may be parameterised and used for sensitivity analysis
  ; all these values are moved to the interface
  
;  set n_geese_in_super_goose 20
;  set big_radius 5
;  set small_radius 1
;  set alpha 0.001
;  set av_ee 54.28
;  set ee_fly 416.35
;  set starvation_energy_stores 9620
;  set prob-of-disturbance-small-fields 0.3
;  set prob-of-disturbance-large-fields 0.2
;  set cum_time_threshold 3
;  set field_area_threshold 60000 ; Jesper showed (unpublished data) that geese flee from on average 120 m in Mid-Norway as a reaction to disturbance so a field has to be at least 240*240m in order for geese to forage relatively happy
;  set steps_on_500m_roost 1
;  set max_steps-to_stay_on_roost 3
;  set percentage_rule 0.8
;  set roost_dist_radius 0.3
;  set extra_initial_energy 0
;  set extra_leaving_energy 0

  
  create-norths 1 
  [
    set color green
    set size 20
    set shape "square"
    setxy  115 656
    set label "N"
  ]
  
  create-starvs 1 
  [
    set color 63
    set size 20
    set shape "square"
    setxy  174 655
    set label "S"
  ]
  
  
  
  create-suns 1 
  [
    set color black
    set size 40
    set shape "moon"
    setxy  25 612
    ;change appearance based on whether it is day or night
  ]
  
 

  ; metabolised energy intake - equation 5 in ODD
  
         
  set gf_grass 16.18           
  set gf_grain 14.55           
  set gf_stubble 17.24             
  set dr_grass 9.8            
  set dr_grain 4.5            
  set dr_stubble 5.6          
  set gd_grass 12.82           
  set gd_grain 11.47           
  set gd_stubble 13.97          
  set m_grass 0.88             
  set m_grain 0.97            
  set m_stubble 1.30  
  set mei_potato 879.5         
  
  ; habitat depletion - equation 3 and 4 in ODD
  
  set b1 0.28                   
  set b2 9.6
  set b3 2.8
  set tc 0.42                
  set R_max 0.0085
  set a_grain 0.00163                
  set a_stubble 0.00325                         
  set h_stubble 80         
  set h_grain 180  
  
  ask fields [
    
    set initial_grain_density_stubble ((random grain_sd_stubble + 1) + grain_mean_stubble)
    set hab hab1 
    ; all initial values are in grams per field
    if hab = "grass"   [set hab_biomass (16.4 * initial_grass_length * area)  set gr_length (hab_biomass / area / 16.4) ] ;the initial grass value for all the fields is 2cm and change into biomass (g/field) as described in the ODD. 
    if hab = "plough"  [set hab_biomass 0]
    if hab = "potato"  [set hab_biomass 0] 
    if hab = "stubble" [set hab_biomass (initial_grain_density_stubble * area) set grain_dens initial_grain_density_stubble] ; average grain density as measured by me and Bart in 2013
    if hab = "grain"   [set hab_biomass (initial_grain_density_grain * area) set grain_dens initial_grain_density_grain] ; average grain density as measured by me and Bart in 2013
    if hab = "grass"   [set color green ] ; if grass length is defined after hab_biomass, it can already appear is settings, otherwise after the first step
    if hab = "grain"   [set color orange]
    if hab = "other"   [set color 127]
    if hab = "plough"  [set color gray]
    if hab = "potato"  [set color brown]
    if hab = "stubble" [set color yellow] 

    set shape "circle"
    set size 3
    set occupied? FALSE
    ] 
  
  create_geese
  ask geese [calculate_initial_gamma]
 
  ask roosts [set num_geese_step [] set max_geese_day [] ]
  
  set gl_max_num_geese_per_roost []
  set gl_mean_num_geese_per_roost []
end
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;; ALL LOADING PROCEDURES: BACKGROUND; ROOSTS; FIELDS; DAY AND NIGHT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to load_background
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
  gis:set-world-envelope (gis:envelope-union-of (gis:envelope-of water) (gis:envelope-of bcg_fields))
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
    [set pcolor 9
    set background "patches"]
  ]   
  

                                       
 reset-ticks
end



to read_roosts_from_txt

; code taken from Jacob's porpoise model

let roost_pos "input/roost_utm_24.02.txt" ; this list has roost with 34 and 7 deleted becuase they were outside my study area
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

to load_roosts
  
  let i 0
  if (not (length roost-list < 2)) [ ; First line contains header
    set roost-list remove-item 0 roost-list ; remove header
    let nb-of-roosts length roost-list
    create-roosts nb-of-roosts

      ask roosts 
      [
      set rid (item 0 (item i roost-list))
      set rx-coord ((item 1 (item i roost-list)) - xllcorner) * 0.020   ; comes from: the real environment is 28500m and there is 571 patches (571/28500 = 0.02)
      set ry-coord ((item 2 (item i roost-list)) - yllcorner) * 0.020   ; 675/33500 = 0.02 
      set capacity item 3 (item i roost-list) 
      setxy rx-coord ry-coord
      set hab "roost"
      set i i + 1
      set color 14
      set shape "square"
      set size 8
      ]
       ]
  print count roosts
end

to read_fields_from_txt
  
let fields_pos "input/habitat_eol_mac_06.02.txt" ; This is a new complete habitat map created 28.10 based on R script (details in making new habitat map folder)
; than file has to be saved as Mac-type txt file as done in the same R scipt 
let line-txt ""
let line-lst (list " " " " " " " " " " " " " " " " " " " " " " " " " ") ;as many as columns in my text file
set fw-list (list line-lst)
let line-lgt 0
file-open fields_pos
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
      if (k > 0) [ set fw-list lput line-lst fw-list ]
      set k k + 1
    ]
    file-close    
end

to load_fields
  
  let i 0
  if (not (length fw-list < 2)) [ ; First line contains header
    set fw-list remove-item 0 fw-list ; remove header
    let nb-of-fields length fw-list
    create-fields nb-of-fields

      ask fields 
      [
      set hab1 (item 2 (item i fw-list))
      set hab2 (item 3 (item i fw-list))
      set hab3 (item 4 (item i fw-list))
      set hab4 (item 5 (item i fw-list))
      set id (item 0 (item i fw-list))
      set id_gis (item 8 (item i fw-list))
      set fx-coord ((item 6 (item i fw-list)) - xllcorner) * 0.020   ; comes from: the real environment is 28500m and there is 571 patches (571/28500 = 0.02)
      set fy-coord ((item 7 (item i fw-list)) - yllcorner) * 0.020   ; 675/33500 = 0.02 
      set area (item 1 (item i fw-list))  ; area is in m2 i
      setxy fx-coord fy-coord
      set i i + 1

      ]
       ]
  
  print count fields
end  



to read_daynight_from_txt
; this file contains as many rows and ticks in the model. every row says whether it is a day (D) or night (N)
; day is set as the time between two civil twilights (source: aa.usno.navy.mil)
  file-open "input/daynight.txt" 
  if file-read = "day" [set daynight-list file-read set day item 0 daynight-list]
  file-close  
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; GEESE POPULATION SETUP ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; creating geese 


to create_geese 
   create-geese round ((item 0 arrival_list) / n_geese_in_super_goose) [ ;at the start of the model there should be 352 agents
    set color black
    setxy  0 0 ; starting point does not matter cause I tell geese to move to a random roost immidietaly
    move-to one-of roosts
    set tot_nei 0
    set tot_ee 0 ; geese haven't spent any energy when the model starts
    set tot_eg 0 ; geese haven't gained any energy when the model starts
    set step_ee 0
    set step_eg 0
    set step_nei 0
    set en_fly 0 ; geese haven't spent any energy on flying when the model starts
    set time_fly 0 ; geese haven't spent any time on flying when the model starts
    set tot_en_fly 0 
    set in_energy_stores ((random-float 6214.1) + 18941 + extra_initial_energy) ;The initial goose energy stores (kJ) is taken as the quivalent of API 2 +/- 0.5 (18941 - 25155 kJ)
    set ticks_since_morning 20 ; geese arrive to the model at miDnight
    set cum_time_feeding 0
    set cum_tot_time_feeding 0
    set previous_loc one-of roosts-here ; a goose always starts the model from a roost so it has to remember which roost it was so in the next step it can calculate the distance to this roost. there must be one-of, otherwise it reports agenset
    set night_roost_location one-of roosts-here
    set tot_dist_flown 0 ; geese haven't moved from the beginning of the model
    set day_roost_500m false
    set left_model false
    set moving_north_threshold ((random-float 4256.1) + 42908 + extra_leaving_energy) ; the API of geese leaving Trondelag is between 3.75-4.25 API (42908 - 47164 kJ, fat efficiency taken into account (0.80))
    set departure_date ((random 376) + 624) ; if geese obtain a defined energy, they should leave between 1st and 10th of May what is day 26-35 of the model and therefore tick 624-840
    set steps_spent_on-day_roost 0
    set size 4
    set disturbed? false
    ;set moved? false
    set on_roost? true
    set mei_list []
    
    ]
   
  
end


to calculate_initial_gamma
 
; the initial gamma (expected gain rate) is calculated as the possible average gain rate of fields within 5km from the first roost geese fly to when they arrive to Mid-Norway 
; this is according to 'classical' Charnov marginal value theorem and is based on an idea that geese have some expectation based on  experience from previous years / other stopover sites
;we assume energy expenditure equal to average energy expenditure of one hour
  

  if any? fields with [hab = "grass"] with [distance myself < (5 / patch_size)]
  [
    let av_gr_length (mean [gr_length] of fields with [hab = "grass"] with [distance myself < (5 / patch_size)])
    let ir_grass 3600 * (1 / ((1 + b2 * av_gr_length / 100) / (b1 * av_gr_length / 100)*(tc + b3 * av_gr_length / 100) + 1 / R_max))
    set mei_grass ((ir_grass * gf_grass) - (dr_grass *  m_grass * gd_grass))
  ]

  if any? fields with [hab = "stubble"] with [distance myself < (5 / patch_size)]
  [
    let av_stubble (mean [grain_dens] of fields with [hab = "stubble"] with [distance myself < (5 / patch_size)])
    let ir_stubble 3600 * ((a_stubble * av_stubble) / (1 + a_stubble * h_stubble * av_stubble))
    set mei_stubble ((ir_stubble * gf_stubble) - (dr_stubble * m_stubble * gd_stubble))
  ]

 
  if any? fields with [hab = "grain"] with [distance myself < (5 / patch_size)]
  [
    let av_grain (mean [grain_dens] of fields with [hab = "grain"] with [distance myself < (5 / patch_size)])
    let ir_grain 3600  * ((a_grain * av_grain) / (1 + a_grain * h_grain * av_grain)) 
    set mei_grain ((ir_grain * gf_grain) - (dr_grain * m_grain * gd_grain))
  ]

  if any? fields with [hab = "potato"] with [distance myself < (5 / patch_size)]
  [
    set mei_pot mei_potato
  ]
  
    if any? fields with [hab = "plough"] with [distance myself < (5 / patch_size)]
  [
    set mei_ploughed 0
  ]
                                    

   if decision_making_rule = "random" or decision_making_rule = "asocial learning" or decision_making_rule = "social learning" or decision_making_rule = "asocial and social" [set gamma (mean (list mei_grass mei_stubble mei_ploughed mei_grain mei_pot) - av_ee)]
   if decision_making_rule = "max energy" [set gamma (mean (list mei_grass mei_stubble mei_grain mei_pot) - av_ee)]


;;;;;;;;;;;;  PARAMETERS TO CALIBRATE ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  file-open "experiment-plan-calibration.txt"      ;; Import matrix for calibration 
  let calibration-matrix file-read 
  set calibration-data matrix:from-row-list calibration-matrix
  file-close
  
  set cum_time_threshold   (matrix:get calibration-data expindex 0)
  set roost_dist_radius    (matrix:get calibration-data expindex 1)
  set alpha                (matrix:get calibration-data expindex 2)


 
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;; PROCEDURES IN 'GO' ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to go

;if ticks = 1 [random-seed 3] ; I also have to specify it here, otherwise in later ticks it does not work
if ticks = 792 [output-number-of-geese-on-roost]
;if ticks > 0 [output-debugging-geese-movement-roost-rules]
;if ticks > 0 [output-debugging-energetics] ; if this is after comment 'tick', it starts printing tick 1 as tick 2 but IT STILL DOES NOT PRINT THE LAST TICK
;  if ticks > 0 [output-debugging-hab-depletion]
;  if ticks = 480 or ticks = 481 [output-hab-biomass] ; this only prints the output at two ticks
;  let N n-values 15 [?] foreach N [if ticks = 100 * (item ? N) [output]] ; in the output table I should have a value if hab_biomass every 100 ticks, [?] is a reporter-task
  
  tick ; otherwise there is no ticks
  
;  ifelse count geese = 0
;    [stop] ; model either stops when or geese migrated north or starved
;    [
      if ticks = 1057 [stop]
;    ] ; or when is no later then 19.05 
  
;;  profiler:start ; this shows how much time each procedure takes
  update_daynight
  show-day
  update_day_ticks
  habitat_change
  goose_arrival
  ;print gl_ticks_since_morning
  ask fields with [hab = "grass"] [Habitat_growth]
  
;order does matter here. Geese first fly to a field, then something happens on this field
;because it is all in one bracket, everything is done by one agent, then by another etc. it matters in depletion process

  
  
   
      set eg_reduction_dist random-float 1
      ask fields [set dist_prob random-float 1]
         
      ask geese 
      
      [ 
        ifelse day = "D" 
        [          
          movement_during_day
          calculate_dist_to_roost
        ]
        [
          movement_during_night
        ]
        energy_expenditure
        metabolised_energy_intake
        net_energy_intake
        if any? fields-here [hab_depletion] ; geese cannot deplete roosts
        leaving_the_model
      ]
      
      ask fields ; everytime there a goose on a patch, it is marked as occupied. Necessary for DMR 4 and 5
      [
        ifelse any? geese-here
        [set occupied? TRUE]
        [set occupied? FALSE]
      ]
      
      ; making a list with number of goose agents visiting each roost (for pattern 5 in POM). 
      ;In real world number of geese roosting at each roost site was only monitored during daylight so same should apply in the model
      ; I also need information if number of visiting geese = 0
      
      
      ask roosts
      [
        if day = "D" 
        [
          let number_of-roosting_geese count geese-here
          set num_geese_step lput number_of-roosting_geese num_geese_step
        ]
        
      ]
      
      ; every morning the list is reset
      if gl_ticks_since_morning = 1
        [
          ask roosts [set num_geese_step []]
        ]
      
      ; for the partial analysis I only need a list showing the maximum number of roosting geese per day per roost
      ; so I need a list showing the max number of geese visiting each roost per day
      
      if gl_ticks_since_morning = 23 and ticks > 3 ; (during the first night of the model the list is empty)
        [
          ask roosts [set max_geese_day lput max num_geese_step max_geese_day]
        ]
        
      ; for the final analysis I only need a list showing the maximum and the mean of daily maxima number of roosting geese during the entire model per roost
      ; this is a global variable which will be calculated 08.05 becuase this is the latest date from which I have roost counts 
      
      
      if ticks = 792 ; 08.05
      [
        ask roosts 
        [
          let max_number_geese_per_roost max max_geese_day
          ;let mean_number_geese_per_roost mean max_geese_day
          
          set gl_max_num_geese_per_roost  lput list rid max_number_geese_per_roost gl_max_num_geese_per_roost
          ;set gl_mean_num_geese_per_roost lput list rid mean_number_geese_per_roost  gl_mean_num_geese_per_roost
          ;show gl_max_num_geese_per_roost
        ]
      ]
      
;; setting these two below global parameters (period and time_of_day) are only used for the analysis of the results, they do not take part or influence any processes in the model
  
  ifelse ticks < 479 [set period 1]
  [ifelse ticks > 479 and ticks <= 671 [set period 2]
    [ifelse ticks > 671 and ticks <= 863 [set period 3]
      [set period 4]]]
  

  ifelse period = 1 and gl_ticks_since_morning >= 1 and gl_ticks_since_morning <= 6 [set time_of_day "morning"]
  [ifelse period = 1 and gl_ticks_since_morning >= 7 and gl_ticks_since_morning <= 11 [set time_of_day "midday"]
    [ifelse period = 1 and gl_ticks_since_morning >= 12 and gl_ticks_since_morning <= 16 [set time_of_day "evening"]
      [ifelse period = 2 and gl_ticks_since_morning >= 2 and gl_ticks_since_morning <= 7 [set time_of_day "morning"]
        [ifelse period = 2 and gl_ticks_since_morning >= 8 and gl_ticks_since_morning <= 12 [set time_of_day "midday"]
          [ifelse period = 2 and gl_ticks_since_morning >= 13 and gl_ticks_since_morning <= 17 [set time_of_day "evening"]
            [ifelse period = 3 and gl_ticks_since_morning >= 3 and gl_ticks_since_morning <= 8 [set time_of_day "morning"]
              [ifelse period = 3 and gl_ticks_since_morning >= 9 and gl_ticks_since_morning <= 13 [set time_of_day "midday"]
                [ifelse period = 3 and gl_ticks_since_morning >= 14 and gl_ticks_since_morning <= 18 [set time_of_day "evening"]
                  [ifelse period = 4 and gl_ticks_since_morning >= 4 and gl_ticks_since_morning <= 9 [set time_of_day "morning"]
                    [ifelse period = 4 and gl_ticks_since_morning >= 10 and gl_ticks_since_morning <= 14 [set time_of_day "midday"]
                      [ifelse period = 4 and gl_ticks_since_morning >= 15 and gl_ticks_since_morning <= 19 [set time_of_day "evening"]
                        [set time_of_day "night"]]]]]]]]]]]]
  
  print time_of_day
  
;;    plot_mean_dnei ; at the end of each day shows mean daily net energy intake of all geese
;;    plot_time_on_day_roosts ; at the end of each day shows mean number of steps geese spent on day roosts as proportion of daylight hours
;    plot_distance_to_roost
;;    plot_phenology ; at the end of each day plots number of geese present in the model

;;    plot_grass_growth
  ;  plot_n_grass
  ;  if ticks > 480 [plot_grain_m2]
  ;  plot_stubble_m2
  ;  if ticks > 481 [plot_grain_m2] ; there is no grain fields before so it will complain
  ;  plot_grain_m2
  ;  plot_plough
  
;;  profiler:stop
  ;let N n-values 1056 [?] foreach N [if ticks = 1 * (item ? N) [print profiler:report]] ; every 1 ticks
;;  profiler:reset
   
end

to update_daynight
  set day item (ticks - 1) daynight-list ; every tick is taken a consecutive value from the daynight file
end


to show-day  ; this is for sun to change colour depending whether it is day or night
ifelse day = "D" [ask suns [set color yellow set shape "sun"]] 
                 [ask suns [set color black set shape "moon"]]
end

to update_day_ticks
  
; ticks are restarted when civil twilight starts to show how many ticks has passes since morning so daily energetics can be calculated as period from sunrise to sunrise
  
  let day_at_previous_tick "N" 
  if ticks > 1 [set day_at_previous_tick item (ticks  - 2) daynight-list]
  ;print (word "pervious" day_at_previous_tick "current " day)
  ifelse day = "D" and day_at_previous_tick = "N" 
   [set gl_ticks_since_morning 1 ]
   [set gl_ticks_since_morning gl_ticks_since_morning + 1]
   
  ask geese [set ticks_since_morning gl_ticks_since_morning]
end


to habitat_change ; changes habitat on the periodical basis
  
  if ticks = 20 * 24 + 1 [ask fields 
                                    [ 
                                      set initial_grain_density_stubble ((random grain_sd_stubble + 1) + grain_mean_stubble)
                                      set hab hab2
                                      ifelse hab = "grass"   [set color green] 
                                      [ifelse hab = "grain"   [set color orange]
                                      [ifelse hab = "other"   [set color 127] ; plum colour
                                      [ifelse hab = "plough"  [set color gray]
                                      [ifelse hab = "potato"  [set color brown]
                                      [set color yellow]]]]] ; In case of "stubble"
                                      
                                      if hab = "grass" and hab1 = "grass"  [set gr_length gr_length set hab_biomass hab_biomass ifelse gr_length = 0 [set color red] [set color green] stop] ; the grass gets the values of grass biomass from the previous step. 
                                      if hab = "grass" and hab1 != "grass" [set gr_length 1 set hab_biomass (16.4 * gr_length * area) set grain_dens 0 stop] ;if it wasn't grass before (could have been other or stubble or ploughed), it should have value 1 because this is new grass which just started growing
                                      if hab = "grain" and hab1 = "grass"  [set gr_length 0 set grain_dens initial_grain_density_grain set hab_biomass (initial_grain_density_grain * area)  stop]; fields which were grass and later were ploughed, need to have grass length reset to zero
                                      ;if hab = "plough" and hab1 = "grass" [set gr_length 0 set hab_biomass 0 set grain_dens 0 stop]
                                      if hab = "grain"   [set grain_dens initial_grain_density_grain set hab_biomass (initial_grain_density_grain * area) stop] 
                                      if hab = "plough"  [set hab_biomass 0 set grain_dens 0 set gr_length 0 stop]
                                      if hab = "other"   [set hab_biomass 0 set grain_dens 0 set gr_length 0 stop]                                      
                                      if hab = "potato"  [set hab_biomass hab_biomass stop]
                                      if hab = "stubble" and hab1 = "stubble" [set hab_biomass hab_biomass set grain_dens grain_dens stop] 
                                      if hab = "stubble" and hab1 != "stubble"[set grain_dens initial_grain_density_stubble set hab_biomass (initial_grain_density_stubble * area) stop]

                                      ]
                                    ]
  
  if ticks = 28 * 24 + 1 [ask fields 
                                    [
                                      set initial_grain_density_stubble ((random grain_sd_stubble + 1) + grain_mean_stubble)
                                      set hab hab3 ;the new fields appear 03.05 at 01:00
                                      ifelse hab = "grass"   [set color green] 
                                      [ifelse hab = "grain"   [set color orange]
                                      [ifelse hab = "other"   [set color 127] ; plum colour
                                      [ifelse hab = "plough"  [set color gray]
                                      [ifelse hab = "potato"  [set color brown]
                                      [set color yellow]]]]] ; In case of "stubble"
                                      
                                      if hab = "grass" and hab2 = "grass"    [set gr_length gr_length set hab_biomass hab_biomass ifelse gr_length = 0 [set color red] [set color green] stop] ; the grass gets the values of grass biomass from the previous step. 
                                      if hab = "grass" and hab2 != "grass"   [set gr_length 1 set hab_biomass (16.4 * gr_length * area) set grain_dens 0 stop] ;if it wasn't grass before (could have been other or stubble or ploughed), it should have value 1 because this is new grass which just started growing, it should have value 1 because this is new grass which just started growing
                                      if hab = "grain" and hab2 = "grass"    [set gr_length 0 set grain_dens initial_grain_density_grain set hab_biomass (initial_grain_density_grain * area)  stop]; fields which were grass and later were ploughed, need to have grass length reset to zero
                                      if hab = "plough" and hab2 = "grass"   [set gr_length 0 set hab_biomass 0 set grain_dens 0 stop] 
                                      if hab = "plough" and hab2 = "stubble" [set gr_length 0 set hab_biomass 0 set grain_dens 0 stop] 
                                      if hab = "grain" and hab2 = "grain"    [set hab_biomass hab_biomass set grain_dens grain_dens stop] ; it has to get the grain values from the previous step
                                      if hab = "grain" and hab2 != "grain"   [set grain_dens initial_grain_density_grain set hab_biomass (initial_grain_density_grain * area) set gr_length 0 stop]
                                      if hab = "other"   [set hab_biomass 0 set grain_dens 0 set gr_length 0 stop]
                                      if hab = "plough"  [set hab_biomass 0 set grain_dens 0 set gr_length 0 stop]
                                      if hab = "potato"  [set hab_biomass hab_biomass set grain_dens 0 set gr_length 0 stop]
                                      if hab = "stubble" and hab2 = "stubble" [set hab_biomass hab_biomass set grain_dens grain_dens stop] 
                                      if hab = "stubble" and hab2 != "stubble"[set grain_dens initial_grain_density_stubble set hab_biomass (initial_grain_density_stubble * area)  stop]
 
                                      ]
                                    ]
  
   if ticks = 36 * 24 + 1 [ask fields 
                                     [
                                      set initial_grain_density_stubble ((random grain_sd_stubble + 1) + grain_mean_stubble)
                                      set hab hab4;the new fields appear 11.05 at 01:00
                                      ifelse hab = "grass"   [set color green] 
                                      [ifelse hab = "grain"   [set color orange]
                                      [ifelse hab = "other"   [set color 127] ; plum colour
                                      [ifelse hab = "plough"  [set color gray]
                                      [ifelse hab = "potato"  [set color brown]
                                      [set color yellow]]]]] ; In case of "stubble"
                                      
                                      if hab = "grass" and hab3 = "grass"    [set gr_length gr_length set hab_biomass hab_biomass ifelse gr_length = 0 [set color red] [set color green] stop] ; the grass gets the values of grass biomass from the previous step. 
                                      if hab = "grass" and hab3 != "grass"   [set gr_length 1 set hab_biomass (16.4 * gr_length * area) set grain_dens 0 stop] ;if it wasn't grass before (could have been other or stubble or ploughed), it should have value 1 because this is new grass which just started growing, it should have value 1 because this is new grass which just started growing
                                      if hab = "grain" and hab3 = "grass"    [set gr_length 0 set grain_dens initial_grain_density_grain set hab_biomass (initial_grain_density_grain * area)  stop]; fields which were grass and later were ploughed, need to have grass length reset to zero
                                      if hab = "plough" and hab3 = "grass"   [set gr_length 0 set hab_biomass 0 set grain_dens 0 stop] 
                                      if hab = "grain" and hab3 = "grain"    [set hab_biomass hab_biomass set grain_dens grain_dens stop] ; it has to get the grain values from the previous step
                                      if hab = "grain" and hab3 != "grain"   [set grain_dens initial_grain_density_grain set hab_biomass (initial_grain_density_grain * area) set gr_length 0 stop]
                                      if hab = "other"   [set hab_biomass 0 set grain_dens 0 set gr_length 0 stop]
                                      if hab = "plough"  [set hab_biomass 0 set grain_dens 0 set gr_length 0 stop]
                                      if hab = "potato"  [set hab_biomass hab_biomass set grain_dens 0 set gr_length 0 stop]
                                      if hab = "stubble" and hab3 = "stubble" [set hab_biomass hab_biomass set grain_dens grain_dens stop] 
                                      if hab = "stubble" and hab3 != "stubble"[set grain_dens initial_grain_density_stubble set hab_biomass (initial_grain_density_stubble * area) stop]
                                     ]
                                   ]
end

 to goose_arrival
   
;geese arrive to the model according to field observations and polynomial relationship given by Bart and Hand. Geese arrive until 30.04 when around 43000 geese (2150 agents)
;should be in my study area. Details about calculations are in 'arrival dates.r'

  while [ticks = item c tick_list] 
  [create-geese round((item c arrival_list) /  n_geese_in_super_goose)
  [ set color black
    setxy  0 0 ; starting point does not matter cause I tell geese to move to a random roost immidietaly
    move-to  one-of roosts
    set tot_nei 0
    set tot_ee 0 ; geese haven't spent any energy when the model starts
    set tot_eg 0 ; geese haven't gained any energy when the model starts
    set en_fly 0 ; geese haven't spent any energy on flying when the model starts
    set time_fly 0 ; geese haven't spent any time on flying when the model starts
    set tot_en_fly 0 
    set step_ee 0
    set step_eg 0
    set step_nei 0
    set in_energy_stores ((random-float 6214.1) + 18941 + extra_initial_energy); The initial goose energy stores (kJ) is taken as the quivalent of API 2 +/- 0.5
    set ticks_since_morning 0 ; geese arrive to the model at mignight
    set cum_time_feeding 0
    set cum_tot_time_feeding 0
    set previous_loc one-of roosts-here ; a goose always starts the model from a roost so it has to remember which roost it was so in the next step it can calculate the distance to this roost. there must be one-of, otherwise it reports agenset
    set night_roost_location one-of roosts-here
    set tot_dist_flown 0 ; geese haven't moved from the beginning of the model
    set day_roost_500m false
    set left_model false
    set moving_north_threshold ((random-float 4256.1) + 42908 + extra_leaving_energy) ; the API of geese leaving Trondelag is between 3.75-4.25 API (42908 - 47164 kJ, fat efficiency taken into account (0.80))
    set departure_date ((random 376) + 624) ; if geese obtain a defined energy, they should leave between 1st and 10th of May what is day 26-35 of the model and therefore tick 624-840
    set size 8
    set steps_spent_on-day_roost 0
    calculate_initial_gamma
    set disturbed? false
    ;set moved? false
    set on_roost? true
    set mei_list []
    
   ]
    set c c + 1 if c > 24 [set c 0] ; by having this I don't have to have this extra 0 in 'arrival_list'
    ;print c
   ]
  
  
  ;print count geese
 end 


to Habitat_growth 
  
  ; grass is growing everyday irrespective of the grazing preasssure
  ; grass is growing 0.038 cm/day during periods 1, 0.126 cm/day in period 2, 0.279 cm/day in period 3 and 0.459 cm/day in period 4
  ; Conversion into biomass is done based on equation provided in ODD
    

     if ticks < 20 * 24 + 1                             
     [set gr_length gr_length + (0.038 / 24) 
      ifelse gr_length = 0 [set color red] [set color green] 
      set hab_biomass (16.4 * gr_length * area)
     ]
          
     if ticks >= 20 * 24 + 1 and ticks < 28 * 24 + 1    
     [set gr_length gr_length + (0.126 / 24) 
      ifelse gr_length = 0 [set color red] [set color green] 
      set hab_biomass (16.4 * gr_length * area)
      ]
     
     if ticks >= 28 * 24 + 1 and ticks < 36 * 24 + 1    
     [set gr_length gr_length + (0.279 / 24)
      ifelse gr_length = 0 [set color red] [set color green] 
      set hab_biomass (16.4 * gr_length * area)
      ]
     
     if ticks >= 36 * 24 + 1                            
     [set gr_length gr_length + (0.459 / 24)
      ifelse gr_length = 0 [set color red] [set color green] 
      set hab_biomass (16.4 * gr_length * area)
      ]
end

to calculate_dist_to_roost
  
  ; in period 1 gl_ticks_since_morning = 1 corresponds to 5 AM so gl_ticks_since_morning = 2 would be 6AM, =8 would be 12AM and =14 would be 6PM
  ; in period 2 gl_ticks_since_morning = 1 corresponds to 4 AM so gl_ticks_since_morning = 3 would be 6AM, =9 would be 12AM and =15 would be 6PM
  ; in period 3 gl_ticks_since_morning = 1 corresponds to 3 AM so gl_ticks_since_morning = 4 would be 6AM, =10 would be 12AM and =16 would be 6PM
  ; in period 4 gl_ticks_since_morning = 1 corresponds to 2 AM so gl_ticks_since_morning = 5 would be 6AM, =11 would be 12AM and =17 would be 6PM
  ; civil twilights are summerised in 'civil twilight.xlsx' file
  

        let closest_roost min-one-of roosts [distance myself]
        let dist (((distance closest_roost) * patch_size))
        set dist_closest_roost dist

  
end

to movement_during_day
  
  ; this procedure was first debugged in the little model (maly_model-disturbance) and then debugged in this model by following the reported comments
  
  ; in the lower right corner of the model there are few fields which are further away than 5km from any roost
  ; (they are up 6.5 km) so I had to make an exception and if geese are on any of these fields, they go to the closest roost within 7 km
  
  ; using distance myself is much faster than in-radius
  
  let rand random-float 1
  
  ifelse left_model = false ; it only applies to geese which haven't left the model
  [ 
    
    ifelse ticks_since_morning > 1
    [
      ifelse any? fields-here
      [let closest_field min-one-of fields-here [distance myself]
       set on_roost? false
        ; what happens on larger fields:
        ifelse [area] of closest_field > field_area_threshold 
        [
          
          ; everytime a goose lands on a large field it is counted. Used for plotting number of geese foraging on large and small fields
          ;set n_geese_large_fields n_geese_large_fields + 1
          
          ; what happens if there is disturbance on larger fields:
         ; if rand > prob-of-disturbance-large-fields [print word rand " so no dist"]
          ifelse [dist_prob] of closest_field < prob-of-disturbance-large-fields + temp_disturbance
          [
            set disturbed? true
            set n_disturbance_event n_disturbance_event + 1
            movement_disturbance
 
          ]
          
          ; what happens if there is no disturbance on larger fields:
          [
            set disturbed? false
            movement_no_disturbance
            
          ]
          
        ]
        
        ; what happens on smaller fields:
        [
           ;everytime a goose lands on a small field it is counted. Used for plotting number of geese foraging on large and small fields
           ;set n_geese_small_fields n_geese_small_fields + 1
           
          ; what happens if there is disturbance on smaller fields, the same as if there was a disturbance on larger fields:
         ; if rand > prob-of-disturbance-small-fields [print word rand " so no dist"]
          ifelse [dist_prob] of closest_field < prob-of-disturbance-small-fields + temp_disturbance
          [
            set disturbed? true
            set n_disturbance_event n_disturbance_event + 1
            movement_disturbance
          ]
          
          ; what happens if there is no disturbance on small fields, the same as if there is no disturbance on larger fields:
          [
           set disturbed? false
           movement_no_disturbance
          ]
        ]
      ]
      ; what happens if geese get to the roost during the day:
      
      [
        set on_roost? true
        set disturbed? false
        
        ; if they get to the roost becuase they need time to digest
        
        ifelse day_roost_500m != true
        [
          ifelse steps_to_stay_on_roost > 0
          [
           set steps_to_stay_on_roost steps_to_stay_on_roost - 1 fd 0
           ;show (word "Im on a day roost and leave in " steps_to_stay_on_roost " steps") ; if steps are larger then 0, then decrease one step per tick and stay on this roost
          ]
          
          ; and this is what happens when geese leave the day roost
          
          [
            if decision_making_rule = "random" 
            [
              let field5km one-of fields with [distance myself < (big_radius / patch_size)]
              move-to field5km ; if step is 0, geese should move to a random field within 5 km
              ;show "I left a day roost to a random field within 5km" 
            ]
           
            if decision_making_rule = "max energy"
            [
              
              let field1km  fields with [distance myself < (small_radius / patch_size)]
              ifelse any? field1km 
              
              [ 
                ;show "I left a day roost to a max field within 1km" 
                field-with_highest_nei_within_1km
              ] 
              
              ;; in cases when there is no field within 1km, goose moves to a field with max nei within 5km
              
              [
                ;show "I left a day roost to a max field within 5km becuase there was no field within 1km" 
                field-with_highest_nei_within_5km
                
                
              ]
              ;show "I left a day roost to a max field within 1km" 
              ;field-with_highest_nei_within_1km 
            ]
            
            if decision_making_rule = "asocial learning"
            [
              let field1km  fields with [distance myself < (small_radius / patch_size)]
              ifelse any? field1km 
              
              [ 
                move-to one-of field1km ; if step is 0, geese should move to a random field within 1 km
                ;show "I left a day roost to a random field within 1km"
              ] 
              
              ;; in cases when there is no field within 1km, goose moves to a random field within 5km
              
              [
                let potentialfields5km fields with [distance myself < (big_radius / patch_size)]
                move-to one-of potentialfields5km
                ;show "I left a day roost to a random field within 5km because there was no field within 1km"
              ]
;              let max_field with-max-first mei_list 
;              move-to one-of fields with [who = max_field]
           
  
            ]
            
            if decision_making_rule = "social learning" or decision_making_rule = "asocial and social"
            [
              ;geese should move to a random field already occupied by geese within 5km 
              let occupied_field5km  fields with [distance myself < (big_radius / patch_size) and occupied? = TRUE]
              ifelse any? occupied_field5km 
              [
                move-to one-of occupied_field5km
                ;show "left day roost to random occupied field within 5km"
              ]
              [
                ; if none of the fields is occupied within 5km, they should move to a random field within 5km
                let random_field5km  fields with [distance myself < (big_radius / patch_size)]
                
                move-to one-of random_field5km
                ;show "left day roost to random field within 5km becuase none was occupied within 5km"
              ]
              
            ]
            set on_roost? false
          ]
        ]
        
        ; if the get to the roost after disturbance becuase they were within 500m from a roost
        
        [
          
          ifelse steps_to_stay_on_roost > 0
          [set steps_to_stay_on_roost steps_to_stay_on_roost - 1 fd 0
          ;show (word "Im on 500m day roost and leave in " steps_to_stay_on_roost " steps") ; if steps are larger then 0, then decrease one step per tick and stay on this roost
          ]
          [
            
            if decision_making_rule = "random" or decision_making_rule = "asocial learning"
            [
              
              let field1km fields with [distance myself < (small_radius / patch_size)]
              ifelse any? field1km 
              
              [ 
                move-to one-of field1km ; if step is 0, geese should move to a random field within 1 km
                set on_roost? false
                set day_roost_500m false
                ;show "I left a 500m day roost to a random field within 1km"
              ] 
              
              ;; in cases when there is no field within 1km, goose moves to a random field within 5km
              
              [
                let potentialfields5km fields with [distance myself < (big_radius / patch_size)]
                move-to one-of potentialfields5km
                set on_roost? false
                set day_roost_500m false
                ;show "I left a 500m day roost to a random field within 5km because there was no field within 1km"
              ]
              
            ]
            
            if decision_making_rule = "max energy"
            [

                let field1km  fields with [distance myself < (small_radius / patch_size)]
                ifelse any? field1km 
                
                [ 
                  ;show "I left a day roost to a max field within 1km" 
                  field-with_highest_nei_within_1km
                  set on_roost? false
                  set day_roost_500m false
                ] 
                
                ;; in cases when there is no field within 1km, goose moves to a field with max nei within 5km
                
                [
                  ;show "I left a day roost to a max field within 5km becuase there was no field within 1km" 
                  field-with_highest_nei_within_5km
                  set on_roost? false
                  set day_roost_500m false
                    
                ]

              ]
              ;show "I left a day roost to a max field within 1km" 
;              field-with_highest_nei_within_1km 
;              set on_roost? false
;              set day_roost_500m false 

            if decision_making_rule = "social learning" or decision_making_rule = "asocial and social"
            [
               ; first geese should move to a random field already occupied by geese within 1km 
              let occupied_field1km  fields with [distance myself < (small_radius / patch_size) and occupied? = TRUE]
              ifelse any? occupied_field1km 
              [
                move-to one-of occupied_field1km
                set on_roost? false
                set day_roost_500m false
                ;show "left day roost to random occupied field within 1km"
              ]
              [
                ; if none of the fields is occupied within 1km, they should move to a random, unoccupied field within 1km
                let random_field1km  fields with [distance myself < (small_radius / patch_size) and occupied? = FALSE]
                ifelse any? random_field1km 
                [
                  move-to one-of random_field1km
                  set on_roost? false
                  set day_roost_500m false
                  ;show "left day roost to random field within 1km becuase none was occupied within 1km"
                ]
                [
                  ; if there is no field of any type within 1km, geese should move to a random, occupied field within 5km 
                  let occupied_field5km  fields with [distance myself < (big_radius / patch_size) and occupied? = TRUE]
                  ifelse any? occupied_field5km 
                  [
                    move-to one-of occupied_field5km 
                    set on_roost? false
                    set day_roost_500m false
                    ;show "left day roost to random field within 5km becuase there was no field within 1km"
                  ]
                  [
                    ; if none of the fields within 5km is occupied and there is no fields of any type within 1km, geese move to a random (unoccupied) field within 5km
                    let potentialfields5km fields with [distance myself < (big_radius / patch_size) and occupied? = FALSE]
                    move-to one-of potentialfields5km
                    set on_roost? false
                    set day_roost_500m false
                    ;show "left day roost to random field within 5km becuase there was no occupied field within 1 and 5km"
                  ]
                ]
              ]             
              
            ]      

          ]
        ]
        
        
      ]
      
    ]
    
    ; what happens if there is a first tick of the day and geese have to leave a night roost
    
    
    [
      ; this should be here becuase if I move to to the night procedure it will update it every tick until night is gone and it will eat a lot of time and energy of the model
      ; but it has to be before move procedure so geese do not change the previous_loc yet
      set cum_tot_time_feeding 0
      set steps_spent_on-day_roost 0
      
      if decision_making_rule = "random"
      [
        let field5km one-of fields with [distance myself < (big_radius / patch_size)] ; patch_size = 0.05 so 100 will be 5 km 
        move-to field5km
        ;show "I left night roost to a random field within 5km"
      ]
      
      if decision_making_rule = "max energy"
      [
        ;show "I left night roost to a max field within 5km"
        field-with_highest_nei_within_5km      
      ]
      
      if decision_making_rule = "asocial learning" or decision_making_rule = "asocial and social"
      [

; if this is a first day in Mid_norway for a goose, it goes to a random field within 5km from the roost to which it arrived to the model
; such goose cannot go a field which yesterday resulted in highest MEI becuase there was no yesterday for that goose

        ifelse empty? mei_list
        [
        let field5km one-of fields with [distance myself < (big_radius / patch_size)] ; patch_size = 0.05 so 100 will be 5 km 
        move-to field5km
        ;show "I left n.r. to a random 5km becuase I just arrived to Tr"
        ]
        
        [
        ;show "I left night roost to a max field from yesterday"           
        let max_field with-max-first mei_list 
        ;show max_field
        move-to one-of fields with [who = max_field]
        set mei_list [] ; at the beginning of each day, the list is set to zero so geese only rememebr the visited fields from the previous day

        ]

   
      ]
      
      if decision_making_rule = "social learning" 
        [
          ;geese should move to a random field already occupied by geese within 5km
          let occupied_field5km  fields with [distance myself < (big_radius / patch_size) and occupied? = TRUE]
          ifelse any? occupied_field5km 
            [
              move-to one-of occupied_field5km
              ;show "left night roost to random occupied field within 5km"
            ]
            [
              ; if none of the fields is occupied within 5km, they should move to a random field within 5km
              let random_field5km  fields with [distance myself < (big_radius / patch_size)]
              
              move-to one-of random_field5km
              ;show "left night roost to random field within 5km becuase none was occupied within 5km"
            ]
          
        ]
      
      set on_roost? false 
         
    ]
  ]
  
  ; this is what happens if geese leave the model and move either to starved or move north field
  
  [fd 0
    ;show "I left the model"
  ] 

end

; the reporter below show the second item of sublist in mei_list (the second item always shows who of a field) for 
; which a first item of that sublist is max of the entire mei_list
; this reported was copied from http://netlogo-users.18673.x6.nabble.com/How-to-select-items-with-max-values-in-a-list-td5001876.html

to-report with-max-first [lists] 
  let highest-first max map first lists 
  let winning? task [first ? = highest-first] 
  report first butfirst one-of filter winning? lists 
end 

to movement_disturbance
  
  ifelse any? roosts with [distance myself < (roost_dist_radius / patch_size)] ; in a situation if (a or/and b) and c, c must on a seperate if statement because the procedure will run already of a or b is true, reagrdless c
    [ 
      
      move-to one-of roosts with [distance myself < (roost_dist_radius / patch_size)] ; 6*0.05 (patch_size) is 0.3 km, 
      set day_roost_500m true
      set on_roost? true
      set steps_to_stay_on_roost steps_on_500m_roost ; here cumulative time is not set to zero becuase geese go there for a short time
      set steps_spent_on-day_roost steps_spent_on-day_roost + steps_to_stay_on_roost
      ;show "500m rule"
    ]
  
  ; if geese are not close to a roost they go to a random field within 1 km
    [ 
      
      if decision_making_rule = "random" or decision_making_rule = "asocial learning" 
      [
        move-to one-of fields with [distance myself < (small_radius / patch_size)]
        ;show "dist, no cum, no 500m I moved to a random field within 1km"
      ] ; within 1 km, but geese CAN also come back to the same field 
      
      if decision_making_rule = "max energy" 
      [
        
        let field1km  fields with [distance myself < (small_radius / patch_size)]
        ifelse any? field1km 
        
          [ 
            ;show "I left to a max field within 1km after disturbance" 
            field-with_highest_nei_within_1km
          ] 
        
        ;; in cases when there is no field within 1km, goose moves to a field with max nei within 5km
        
          [
           ;show "I left to a max field within 5km after disturbance becuase there was no field within 1km"
            field-with_highest_nei_within_5km          
          ]

        ;field-with_highest_nei_within_1km
      ]
      
      if decision_making_rule = "social learning" or decision_making_rule = "asocial and social"
        [
          ; first geese should move to a random field already occupied by geese within 1km but not a field goose is on
          let closest_field min-one-of fields-here [distance myself]
          let id_of_field_im_on [id] of closest_field
          let occupied_field1km  fields with [distance myself < (small_radius / patch_size) and occupied? = TRUE and id != id_of_field_im_on]
          ifelse any? occupied_field1km 
            [
              move-to one-of occupied_field1km
              ;show "dist,moved to random occupied field within 1km"
            ]
            [
              ; if none of the fields is occupied within 1km, they should move to a random, unoccupied field within 1km
              let random_field1km  fields with [distance myself < (small_radius / patch_size) and occupied? = FALSE]
              ifelse any? random_field1km 
                [
                  move-to one-of random_field1km
                  ;show "dist, moved to random field within 1km becuase none was occupied within 1km"
                ]
                [
                  ; if there is no field of any type within 1km, geese should move to a random, occupied field within 5km but not a field goose is on
                  let occupied_field5km  fields with [distance myself < (big_radius / patch_size) and occupied? = TRUE and id != id_of_field_im_on]
                  ifelse any? occupied_field5km 
                  [
                    move-to one-of occupied_field5km 
                    ;show "dist, moved to random field within 5km becuase there was no field within 1km"
                  ]
                  [
                    ; if none of the fields within 5km is occupied and there is no fields of any type within 1km, geese move to a random (unoccupiep) field within 5km
                    let potentialfields5km fields with [distance myself < (big_radius / patch_size) and occupied? = FALSE]
                    move-to one-of potentialfields5km
                    ;show "dist, moved to random field within 5km becuase there was no occupied field within 1 and 5km"
                  ]
                ]
            ]             
          
        ]
      
      set on_roost? false
      set steps_to_stay_on_roost 0
      
    ]
     

end



to movement_no_disturbance
  
  if decision_making_rule = "random"
  [
    
    ifelse cum_time_feeding >= cum_time_threshold ; ifelse must be first. geese first obtain energy > x and then flies to a roost. this is the mean of 2683-4096 estimated by Hammond and Diamond
    [ let roost5km one-of roosts with [distance myself < (big_radius / patch_size)]
      
      ; there are very few fields which are further away than 5 km from any roost. if geese are on these field, than they move to a roost within 7 km
      ifelse roost5km = nobody
      
        [let roost7km one-of roosts with [distance myself < (7 / patch_size)] ; patch_size = 0.05 so 100 will be 5 km 
          move-to roost7km
          set on_roost? true]
        [move-to roost5km
          set on_roost? true]
      
      set cum_time_feeding 0
      set steps_to_stay_on_roost (random (max_steps-to_stay_on_roost - 1)) + 2
      set steps_spent_on-day_roost steps_spent_on-day_roost + steps_to_stay_on_roost
      ;show "no disturbance but cum >=3"
    ]
    
    
    
    [ 
      giving-up_procedure
    ]
  ]
  
  
    if decision_making_rule = "max energy" or decision_making_rule = "asocial learning" or decision_making_rule = "social learning" or decision_making_rule = "asocial and social"; geese minimise energy so they always go to the closest roost
  [
    
    ifelse cum_time_feeding >= cum_time_threshold ; ifelse must be first. geese first obtain energy > x and then flies to a roost. this is the mean of 2683-4096 estimated by Hammond and Diamond
    [ 
      let closest_roost min-one-of roosts [distance myself]
      move-to closest_roost set on_roost? true
      set on_roost? true
      set cum_time_feeding 0
      set steps_to_stay_on_roost (random (max_steps-to_stay_on_roost - 1)) + 2
      set steps_spent_on-day_roost steps_spent_on-day_roost + steps_to_stay_on_roost
      ;show "no disturbance but cum >=3-closest roost"
    ]
            
    [ 
      giving-up_procedure
    ]
  ]
  

end



to giving-up_procedure
  
; calculation of expected gain rate, equation 5

  let exp_gamma (alpha * step_nei + (1 - alpha) * gamma)
  
  ifelse step_nei < exp_gamma
  
  [
    if decision_making_rule = "random" or decision_making_rule = "asocial learning"  
    [
      move-to one-of fields with [distance myself < (small_radius / patch_size)]
      ;show "I left this field to a random in 1km because of GU"
    ]
    
    if decision_making_rule = "max energy" 
    [
      ;show "I left this field to a max because of GU"
      field-with_highest_nei_within_1km    
    ]
    
    if decision_making_rule = "social learning" or decision_making_rule = "asocial and social"
      [
        ; first geese should move to a random field already occupied by geese within 1km but not field goose is on
        let closest_field min-one-of fields-here [distance myself]
        let id_of_field_im_on [id] of closest_field
        let occupied_field1km  fields with [distance myself < (small_radius / patch_size) and occupied? = TRUE and id != id_of_field_im_on]
        ifelse any? occupied_field1km 
          [
            move-to one-of occupied_field1km
            ;show "GU, left this field to random occupied field within 1km"
          ]
          [
            ; if none of the fields is occupied within 1km, they should move to a random, unoccupied field within 1km
            let random_field1km  fields with [distance myself < (small_radius / patch_size) and occupied? = FALSE]
            ifelse any? random_field1km 
              [
                move-to one-of random_field1km
                ;show "GU, left to random field within 1km becuase none was occupied within 1km"
              ]
              [
                ; if there is no field of any type within 1km, geese should move to a random, occupied field within 5km but field a goose is on
                let occupied_field5km  fields with [distance myself < (big_radius / patch_size) and occupied? = TRUE and id != id_of_field_im_on]
                ifelse any? occupied_field5km 
                  [
                    move-to one-of occupied_field5km 
                    ;show "GU, left to random field within 5km becuase there was no field within 1km"
                  ]
                  [
                    ; if none of the fields within 5km is occupied and there is no fields of any type within 1km, geese move to a random (unoccupiep) field within 5km
                    let potentialfields5km fields with [distance myself < (big_radius / patch_size) and occupied? = FALSE]
                    move-to one-of potentialfields5km
                    ;show "GU, left to random field within 5km becuase there was no occupied field within 1 and 5km"
                  ]
              ]
          ]             
        
      ]
    
    set on_roost? false
    set gamma exp_gamma
   
  ]
  
  [
    ;show "I am happy on that field, i stay"
  ]

 
end

to field-with_highest_nei_within_1km
  
    ifelse any? fields with [hab = "grass"] with [distance myself < (small_radius / patch_size)]
  [
    let max_gr_length (max [gr_length] of fields with [hab = "grass"] with [distance myself < (small_radius / patch_size)])
    set id_field_with_max_grass ([who] of one-of fields with [gr_length = max_gr_length]  with [distance myself < (small_radius / patch_size)])
    let ir_grass 3600 * (1 / ((1 + b2 * max_gr_length / 100) / (b1 * max_gr_length / 100)*(tc + b3 * max_gr_length / 100) + 1 / R_max))
    set mei_grass ((ir_grass * gf_grass) - (dr_grass *  m_grass * gd_grass))
  ]
  [
    set mei_grass 0
  ]

  ifelse any? fields with [hab = "stubble"] with [distance myself < (small_radius / patch_size)]
  [
    let max_stubble (max [grain_dens] of fields with [hab = "stubble"] with [distance myself < (small_radius / patch_size)])
    set id_field_with_max_stubble ([who] of one-of fields with [grain_dens = max_stubble] with [distance myself < (small_radius / patch_size)])
    let ir_stubble 3600 * ((a_stubble * max_stubble) / (1 + a_stubble * h_stubble * max_stubble))
    set mei_stubble ((ir_stubble * gf_stubble) - (dr_stubble * m_stubble * gd_stubble))
  ]
  [
    set mei_stubble 0 
  ]

 
  ifelse any? fields with [hab = "grain"] with [distance myself < (small_radius / patch_size)]
  [
    let max_grain (max [grain_dens] of fields with [hab = "grain"] with [distance myself < (small_radius / patch_size)])
    set id_field_with_max_grain ([who] of one-of fields with [grain_dens = max_grain] with [distance myself < (small_radius / patch_size)])
    let ir_grain 3600  * ((a_grain * max_grain) / (1 + a_grain * h_grain * max_grain)) 
    set mei_grain ((ir_grain * gf_grain) - (dr_grain * m_grain * gd_grain))
  ]
  [
   set mei_grain 0 
  ]

  ifelse any? fields with [hab = "potato"] with [distance myself < (small_radius / patch_size)]
  [ 
    let potential_potato one-of fields with [hab = "potato"] with [distance myself < (small_radius / patch_size)]
    set mei_pot mei_potato
    set id_field_with_max_potato ([who] of potential_potato)
  ]
  [
   set mei_pot 0 
  ]
  
   let max_mei (max (list mei_grass mei_stubble mei_grain mei_pot)) 
     
   ifelse max_mei = mei_grass 
   [
     move-to one-of fields with [who = [id_field_with_max_grass] of myself]
     ;show "Im on grass with highest mei"
   ]
   [
     ifelse max_mei = mei_stubble 
     [
       move-to one-of fields with [who = [id_field_with_max_stubble] of myself]
       ;show "Im on stubble with highest mei"
     ]
     [
       ifelse max_mei = mei_grain 
       [
         move-to one-of fields with [who = [id_field_with_max_grain] of myself]
         ;show "Im on grain with highest mei"
       ]
       [
         move-to one-of fields with [who = [id_field_with_max_potato] of myself] 
         ;show "Im on potato with highest mei"
         
       ]
     ]
   ]                            


end


to field-with_highest_nei_within_5km
  
    ifelse any? fields with [hab = "grass"] with [distance myself < (big_radius / patch_size)]
  [
    let max_gr_length (max [gr_length] of fields with [hab = "grass"] with [distance myself < (big_radius / patch_size)])
    set id_field_with_max_grass ([who] of one-of fields with [gr_length = max_gr_length]  with [distance myself < (big_radius / patch_size)])
    let ir_grass 3600 * (1 / ((1 + b2 * max_gr_length / 100) / (b1 * max_gr_length / 100)*(tc + b3 * max_gr_length / 100) + 1 / R_max))
    set mei_grass ((ir_grass * gf_grass) - (dr_grass *  m_grass * gd_grass))
  ]
  [
    set mei_grass 0
  ]

  ifelse any? fields with [hab = "stubble"] with [distance myself < (big_radius / patch_size)]
  [
    let max_stubble (max [grain_dens] of fields with [hab = "stubble"] with [distance myself < (big_radius / patch_size)])
    set id_field_with_max_stubble ([who] of one-of fields with [grain_dens = max_stubble] with [distance myself < (big_radius / patch_size)])
    let ir_stubble 3600 * ((a_stubble * max_stubble) / (1 + a_stubble * h_stubble * max_stubble))
    set mei_stubble ((ir_stubble * gf_stubble) - (dr_stubble * m_stubble * gd_stubble))
  ]
  [
    set mei_stubble 0 
  ]

 
  ifelse any? fields with [hab = "grain"] with [distance myself < (big_radius / patch_size)]
  [
    let max_grain (max [grain_dens] of fields with [hab = "grain"] with [distance myself < (big_radius / patch_size)])
    set id_field_with_max_grain ([who] of one-of fields with [grain_dens = max_grain] with [distance myself < (big_radius / patch_size)])
    let ir_grain 3600  * ((a_grain * max_grain) / (1 + a_grain * h_grain * max_grain)) 
    set mei_grain ((ir_grain * gf_grain) - (dr_grain * m_grain * gd_grain))
  ]
  [
   set mei_grain 0 
  ]

  ifelse any? fields with [hab = "potato"] with [distance myself < (big_radius / patch_size)]
  [ 
    let potential_potato one-of fields with [hab = "potato"] with [distance myself < (big_radius / patch_size)]
    set mei_pot mei_potato
    set id_field_with_max_potato ([who] of potential_potato)
  ]
  [
   set mei_pot 0 
  ]
  
   let max_mei (max (list mei_grass mei_stubble mei_grain mei_pot)) 
     
   ifelse max_mei = mei_grass 
   [
     move-to one-of fields with [who = [id_field_with_max_grass] of myself]
     ;show "Im on grass with highest mei"
   ]
   [
     ifelse max_mei = mei_stubble 
     [
       move-to one-of fields with [who = [id_field_with_max_stubble] of myself]
       ;show "Im on stubble with highest mei"
     ]
     [
       ifelse max_mei = mei_grain 
       [
         move-to one-of fields with [who = [id_field_with_max_grain] of myself]
         ;show "Im on grain with highest mei"
       ]
       [
         move-to one-of fields with [who = [id_field_with_max_potato] of myself] 
         ;show "Im on potato with highest mei"
         
       ]
     ]
   ]                            


end



to movement_during_night ; patch leaving rule 4
  
; this procedure was first debugged in the little model (maly_model-disturbance) and then debugged in this model by following the reported comments
  
; if the nights comes and geese are not on a roost yet they move to a roost within a radius. 
;if they are already on a roost they stay there until day comes
  
ifelse left_model = false ; it only applies to geese which haven't left the model
[
  ifelse any? roosts-here 
    [ 
      fd 0
      set on_roost? true
      ;show "it is night, im on the roost"
      let my_roost min-one-of roosts-here [distance myself] ; it cannot be a list
      set cum_time_feeding 0
      set cum_tot_time_feeding 0
      set night_roost_location one-of roosts-here
    ]
    [
      
      if decision_making_rule = "random"
      [
        let potentialRoosts5km roosts with [distance myself < (big_radius / patch_size)]
        ; there are very few fields which are further away than 5 km from any roost. if geese are on these field, than they move to a roost within 7 km
        ifelse any? potentialRoosts5km       
        [move-to one-of potentialRoosts5km set on_roost? true]
        [move-to one-of roosts with [distance myself < (7 / patch_size)] set on_roost? true]  
        ;show "I just got to the night roost"
      ]
      
      if decision_making_rule = "max energy" or decision_making_rule = "asocial learning" or decision_making_rule = "social learning" or decision_making_rule = "asocial and social"
      [
        let closest_roost min-one-of roosts [distance myself]
        move-to closest_roost set on_roost? true
        ;show "I just got to the closest roost-night"
        
      ]
    ]
]
   
; this is what happens if geese leave the model and move either to starved or move north field

[fd 0
 ;show "I left the model" 
]
end
 
     
;;;;;;;;;;;;;;;;;;;;; WHAT HAPPEND TO A GOOSE AGENT IF IS ON A FIELD ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;; UPDATING ENERGETOCS STATEMENT ;;;;;;;;;;;;;

to energy_expenditure
  

  
;calculting the time and energy required for flying from the previous location
; we assume av. speed of a goose 50 km /h. Distance is the number of patches between two points 
;(two middle points of patches on which my agents are) so in order to change into km, I need to multiply it by the size of one patch (km)
; DISTANCE DEPENDS WHETHER THE WORLD IS WRAPPED OR NOT; MY WORLD CANNOT BE WRAPPED!!!
;The energy expenditure does not differ between habitats and fields and roosts


if any? fields-here or any? roosts-here
[
  set time_fly (((distance previous_loc) * patch_size)) / 50
  set en_fly (time_fly * ee_fly)
  set tot_dist_flown tot_dist_flown + ((distance previous_loc) * patch_size)
  set tot_ee tot_ee + en_fly +  ((1 - time_fly) * av_ee)
  set step_ee ((1 - time_fly) * av_ee) ; this is WITHOUT energy exp of flying. No matter if geese are disturbed or not they still have to spend energy on flying from the previous place
  ifelse ticks_since_morning > 1 
    [set dee dee + en_fly +  ((1 - time_fly) * av_ee)]
    [set dee step_ee + en_fly] ; at first tick of the day, dee must only be equal to step_dee plus energy spent on flying
  set tot_en_fly tot_en_fly + en_fly
; energy expenditure is proportional to time of one step (1h) minus time required for flying from the previous location and also includes the enrgy required for flying from the previous location

]

; after the entire procedure this current field is set as the previous field
ifelse any? fields-here 
[set previous_loc one-of fields-here] 
[set previous_loc one-of roosts-here] 

end 

to metabolised_energy_intake
  


; energy intake is proportional to time of one step (1h) minus time required for flying from the previous location 
; and proportional to time spent feeding during one step (details in ODD)

; I have to first calcuate the distance to the closest field-here. If 2 fields are close to each other (are on the same patch), a goose will report two fields-here
; and if I used any? field here, it would gain energy from both of them at the same time (if they were from two different habitat type)

; if geese happen to forage on habitat were habitat biomass is close to 0, I wrote a code to be sure that habitat biomass canot go below 0
if any? fields-here
[
 let closest_field min-one-of fields-here [distance myself]

 if [hab] of closest_field = "grass" 
   [
     ifelse [gr_length] of closest_field > 0
     [
       set pr_feed ((random-float 0.241) + 0.54)
       set ir 3600 *  (1 / ((1 + b2 * [gr_length] of closest_field / 100) / (b1 * [gr_length] of closest_field / 100)*(tc + b3 * [gr_length] of closest_field / 100) + (1 / R_max))) 
       set step_eg ((ir * gf_grass) - (dr_grass * m_grass * gd_grass))
       
       ifelse disturbed? = true
       [
         set step_eg (step_eg - (eg_reduction_dist * step_eg)); 
       ]
       [
         set step_eg step_eg
       ]
       if step_eg < 0 [set step_eg 0]
       
       if decision_making_rule = "asocial learning" or decision_making_rule = "asocial and social"
       [
         ; at each step, on patches where step_eg>0, each goose stores in the below list the step_eg and a who number of a field on which a goose is currently
         set mei_list lput list step_eg [who] of closest_field mei_list
       ]
     ]
     [
       set pr_feed 0 
       set step_eg 0
     ]
   ]
 
 if [hab] of closest_field = "grain"
   [
     ifelse [grain_dens] of closest_field > 0
     [
       set pr_feed ((random-float 0.161) + 0.70)
       set ir 3600 * ((a_grain * [grain_dens] of closest_field) / (1 + a_grain * h_grain * [grain_dens] of closest_field)) 
       set step_eg ((ir * gf_grain) - (dr_grain * m_grain * gd_grain)) 
       ifelse disturbed? = true
       [
         set step_eg (step_eg - (eg_reduction_dist * step_eg))
       ]
       [
         set step_eg step_eg
       ] 
       if step_eg < 0 [set step_eg 0]
       
       if decision_making_rule = "asocial learning" or decision_making_rule = "asocial and social"
         [
           set mei_list lput list step_eg [who] of closest_field mei_list
         ]
     ]
     [
       set pr_feed 0 
       set step_eg 0
     ]
   ]
 
 if [hab] of closest_field = "stubble"
   [
     ifelse [grain_dens] of closest_field > 0
     [
       set pr_feed ((random-float 0.241) + 0.54)
       set ir 3600 * ((a_stubble * [grain_dens] of closest_field) / (1 + a_stubble * h_stubble * [grain_dens] of closest_field)) 
       set step_eg ((ir * gf_stubble) - (dr_stubble * m_stubble * gd_stubble))
       ifelse disturbed? = true
       [
         set step_eg (step_eg - (eg_reduction_dist * step_eg))
       ]
       [
         set step_eg step_eg
       ]
       if step_eg < 0 [set step_eg 0]
       
       if decision_making_rule = "asocial learning" or decision_making_rule = "asocial and social"
       [
         set mei_list lput list step_eg [who] of closest_field mei_list
        ]
     ]
     [
       set pr_feed 0 
       set step_eg 0
     ]
   ]

 if [hab] of closest_field = "potato"
   [
    set pr_feed ((random-float 0.241) + 0.54) 
    set step_eg (mei_potato)
   ifelse disturbed? = true
      [
        set step_eg (step_eg - (eg_reduction_dist * step_eg))
             ]
      
      [
        set step_eg step_eg
      ]
      
   if decision_making_rule = "asocial learning" or decision_making_rule = "asocial and social"
     [
       set mei_list lput list step_eg [who] of closest_field mei_list
     ]
   ]

 if [hab] of closest_field = "other" or [hab] of closest_field = "plough"
   [
    set pr_feed 0 
    set step_eg 0
   ]
]

if any? roosts-here or any? norths-here or any? starvs-here
[
 set pr_feed 0
 set step_eg 0
]

set tot_eg tot_eg + step_eg

ifelse ticks_since_morning > 1 
  [
   set deg deg + step_eg   
   ifelse disturbed? = true  
   [      
     set cum_time_feeding (cum_time_feeding + (pr_feed - pr_feed * eg_reduction_dist))
     set cum_tot_time_feeding (cum_tot_time_feeding + (pr_feed - pr_feed * eg_reduction_dist))
   ]
   [
     set cum_time_feeding cum_time_feeding + pr_feed
     set cum_tot_time_feeding cum_tot_time_feeding + pr_feed
   ]
  ]
  [
   set deg step_eg
   set cum_time_feeding pr_feed
   set cum_tot_time_feeding pr_feed
  ]

end



to net_energy_intake 
 set tot_nei (tot_eg - tot_ee)
 set dnei (deg - dee)
 set step_nei (step_eg - step_ee) ; excluding energy for flying 
end

;;;;;; END OF UPDATING ENERGETICS ;;;;;;;;



to hab_depletion 
  
  let closest_field min-one-of fields-here [distance myself] 
  if [hab] of closest_field = "grass"
    [
      
      ask closest_field 
      [
        set hab_biomass (hab_biomass - (n_geese_in_super_goose * [ir] of myself))
        if hab_biomass <= 0 [set hab_biomass 0]
        set gr_length (hab_biomass / area / 16.4)
        
      ]
    ]
  
  if [hab] of closest_field = "grain"
    [
      
      ask closest_field 
      [
        set hab_biomass (hab_biomass - (n_geese_in_super_goose * [ir] of myself))
        if hab_biomass <= 0 [set hab_biomass 0]
        set grain_dens (hab_biomass / area)
        
      ]
    ]
  
  if [hab] of closest_field = "stubble"
    [
      
      ask closest_field 
      [
        set hab_biomass (hab_biomass - (n_geese_in_super_goose * [ir] of myself))        
        if hab_biomass <= 0 [set hab_biomass 0]
        set grain_dens (hab_biomass / area)
        
      ]
    ]
 
end



to leaving_the_model
  
; geese in the same step will first go to a field and fly north/starv so the energy intake and expenditure from this field is added to general energetics but this
; field is not going to be written in the output
  
  
  if (tot_nei + in_energy_stores) <= starvation_energy_stores 
  [
    move-to one-of starvs 
    set n-of-starved n-of-starved + 1
    set left_model "starved"
    die
  ] ; geese move to the starvation square and stay there



  if (tot_nei + in_energy_stores) >= moving_north_threshold and ticks >= departure_date   

  [
    move-to one-of norths 
    set n-of-vesteralen n-of-vesteralen + 1
    set left_model "moved north"
    die
 
  ]  

end

;;;;;;;;;;;;;;;;;;; patern 5 in POM - output file ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


to output-number-of-geese-on-roost
  file-open "n_geese_on_roost.csv"
  file-type ticks       file-type ","
  ;file-type gl_mean_num_geese_per_roost file-type ","
  file-print gl_max_num_geese_per_roost
  file-close
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; MODEL OUTPUT GRAPHS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

to plot_mean_dnei
  set-current-plot "Mean DNEI"
  
  let day_at_next_tick "N" 
  if ticks > 1 and ticks < 1055 [set day_at_next_tick item (ticks) daynight-list]
  ;print (word "next" day_at_next_tick "current " day)
  
  if day = "N" and day_at_next_tick = "D" 
  [
    let mean_dnei (mean [dnei] of geese)
    plot mean_dnei
  ]
   
end



to plot_time_on_day_roosts
  
  ; calculates proportion of time available for foraging (daylight) spent on day roost
  ; period 1 - average day length 17.28h
  ; period 2 - 18.80
  ; period 3 - 20.04
  ; period 4 - 20.97
  
  set-current-plot "Time on day roost"
  let day_at_next_tick "N" 
  if ticks > 1 and ticks < 1055 [set day_at_next_tick item (ticks) daynight-list]
  if day = "N" and day_at_next_tick = "D" 
  [
    if ticks <= 20 * 24
      [
        let time_on_day_roost (mean [steps_spent_on-day_roost] of geese)
        plot time_on_day_roost / 17.28
      ]
    
    if ticks <= 28 * 24 and ticks > 20 * 24
      [
        let time_on_day_roost (mean [steps_spent_on-day_roost] of geese)
        plot time_on_day_roost / 18.80
      ]
    
    if ticks <= 36 * 24 and ticks > 28 * 24
      [
        let time_on_day_roost (mean [steps_spent_on-day_roost] of geese)
        plot time_on_day_roost / 20.04
      ]
    
    if ticks > 36 * 24
      [
        let time_on_day_roost (mean [steps_spent_on-day_roost] of geese)
        plot time_on_day_roost / 20.97
      ]
  ]
end



;to plot_distance_to_roost
;  
;  ; in period 1 gl_ticks_since_morning = 1 corresponds to 5 AM so gl_ticks_since_morning = 2 would be 6AM, =8 would be 12AM and =14 would be 6PM
;  ; in period 2 gl_ticks_since_morning = 1 corresponds to 4 AM so gl_ticks_since_morning = 3 would be 6AM, =9 would be 12AM and =15 would be 6PM
;  ; in period 3 gl_ticks_since_morning = 1 corresponds to 3 AM so gl_ticks_since_morning = 4 would be 6AM, =10 would be 12AM and =16 would be 6PM
;  ; in period 4 gl_ticks_since_morning = 1 corresponds to 2 AM so gl_ticks_since_morning = 5 would be 6AM, =11 would be 12AM and =17 would be 6PM
;  ; civil twilights are summerised in 'civil twilight.xlsx' file
;  
;  set-current-plot "Distance to roost"
;  
;  
;  if ticks <= 20 * 24
;  [
;    
;    
;    if gl_ticks_since_morning = 2
;      [
;        let mean_dist_roost_6  (mean [dist_closest_roost] of geese)
;        set-current-plot-pen "morning"
;        plot mean_dist_roost_6
;      ]
;      
;    if gl_ticks_since_morning = 8
;      [
;        let mean_dist_roost_12 (mean [dist_closest_roost] of geese)
;         set-current-plot-pen "midday"
;        plot mean_dist_roost_12
;      ]
;      
;    if gl_ticks_since_morning = 14
;      [
;        let mean_dist_roost_18 (mean [dist_closest_roost] of geese)
;        set-current-plot-pen "afternoon"
;        plot mean_dist_roost_18
;      ]
;    
;  ]
;  
;  if ticks <= 28 * 24 and ticks > 20 * 24
;  [
;    
;    
;    if gl_ticks_since_morning = 3
;      [
;        let mean_dist_roost_6  (mean [dist_closest_roost] of geese)
;        set-current-plot-pen "morning"
;        plot mean_dist_roost_6
;      ]
;      
;    if gl_ticks_since_morning = 9
;      [
;        let mean_dist_roost_12 (mean [dist_closest_roost] of geese)
;         set-current-plot-pen "midday"
;        plot mean_dist_roost_12
;      ]
;      
;    if gl_ticks_since_morning = 15
;      [
;        let mean_dist_roost_18 (mean [dist_closest_roost] of geese)
;        set-current-plot-pen "afternoon"
;        plot mean_dist_roost_18
;      ]
;    
;  ]
;  
;  if ticks <= 36 * 24 and ticks > 28 * 24
;  [
;    
;    if gl_ticks_since_morning = 4
;      [
;        let mean_dist_roost_6  (mean [dist_closest_roost] of geese)
;        set-current-plot-pen "morning"
;        plot mean_dist_roost_6
;      ]
;      
;    if gl_ticks_since_morning = 10
;      [
;        let mean_dist_roost_12 (mean [dist_closest_roost] of geese)
;         set-current-plot-pen "midday"
;        plot mean_dist_roost_12
;      ]
;      
;    if gl_ticks_since_morning = 16
;      [
;        let mean_dist_roost_18 (mean [dist_closest_roost] of geese)
;        set-current-plot-pen "afternoon"
;        plot mean_dist_roost_18
;      ]
;  ]
;  
;  if ticks > 36 * 24 
;  [
;    
;    if gl_ticks_since_morning = 5
;      [
;        let mean_dist_roost_6  (mean [dist_closest_roost] of geese)
;        set-current-plot-pen "morning"
;        plot mean_dist_roost_6
;      ]
;      
;    if gl_ticks_since_morning = 11
;      [
;        let mean_dist_roost_12 (mean [dist_closest_roost] of geese)
;         set-current-plot-pen "midday"
;        plot mean_dist_roost_12
;      ]
;      
;    if gl_ticks_since_morning = 17
;      [
;        let mean_dist_roost_18 (mean [dist_closest_roost] of geese)
;        set-current-plot-pen "afternoon"
;        plot mean_dist_roost_18
;      ]
;    
;  ]
;  
;  
;  
;end


to plot_phenology
  set-current-plot "Count super-geese"
  
  let day_at_next_tick "N" 
  if ticks > 1 and ticks < 1055 [set day_at_next_tick item (ticks) daynight-list]
  if day = "N" and day_at_next_tick = "D" 
  [
    let n_super_geese count geese
    plot n_super_geese
  ]

end







;;;;;;;;;;;;; DEBUGING   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;to output-hab-biomass
;  
;;; Printing the output table for testing
;
;file-open "hab_biomass_output.csv"
; ask fields[
;     file-type ticks         file-type ","
;     file-type id            file-type ","
;     file-type hab           file-type ","
;     file-type area          file-type ","
;     file-print hab_biomass]     
;file-close
;end
;
;to output-debugging-energetics
;  file-open "debugging_energetics.csv"
; 
; ask geese
; [ if any? fields-here
;   [let closest_field min-one-of fields-here [distance myself]
;  
;     ask closest_field
;     [
;       file-type ticks         file-type ","
;       file-type day           file-type ","
;       file-type hab           file-type ","
;       file-type gr_length     file-type ","
;       file-type grain_dens    file-type ","
;       ]  
; file-type pr_feed            file-type ","
; file-type time_fly            file-type ","
; file-type step_ee             file-type ","
; file-type step_eg             file-type ","
; file-type en_fly              file-type ","
; file-type dnei                file-type ","
; file-type dee                 file-type ","
; file-type tot_ee              file-type ","
; file-type tot_eg              file-type ","
; file-type tot_nei             file-type ","
; file-type tot_mass            file-type ","
; file-print deg
; ]
;   
; if any? roosts-here
;   [let closest_roost min-one-of roosts-here [distance myself]
;     
;     ask closest_roost
;     [
;       file-type ticks             file-type ","
;       file-type day               file-type ","
;       file-type hab               file-type ","
;      ]  
;     file-type pr_feed             file-type ","
;     file-type time_fly            file-type ","
;     file-type step_ee             file-type ","
;     file-type step_eg             file-type ","
;     file-type en_fly              file-type ","
;     file-type dnei                file-type ","
;     file-type dee                 file-type ","
;     file-type tot_ee              file-type ","
;     file-type tot_eg              file-type ","
;     file-type tot_nei             file-type ","
;     file-type tot_mass            file-type ","
;     file-print deg
;   ]
;   
; if any? norths-here or any? starvs-here
;   [
;     file-type ticks               file-type ","
;     file-type day                 file-type ","
;     file-type breed               file-type ","
;     file-type 0                   file-type ","
;     file-type 0                   file-type ","
;     file-type  pr_feed            file-type ","
;     file-type time_fly            file-type ","
;     file-type step_ee             file-type ","
;     file-type step_eg             file-type ","
;     file-type en_fly              file-type ","
;     file-type dnei                file-type ","
;     file-type dee                 file-type ","
;     file-type tot_ee              file-type ","
;     file-type tot_eg              file-type ","
;     file-type tot_nei             file-type ","
;     file-type tot_mass            file-type ","
;     file-print deg
;   ]
; ]
;file-close
;end
;
;to output-debugging-hab-depletion
;  file-open "debugging_hab_depletion.csv"
;  ask geese
;  [
;   let closest_field min-one-of fields-here [distance myself]
;  
;     ask closest_field
;     [
;       file-type ticks         file-type ","
;       file-type hab           file-type ","
;       file-type area          file-type ","
;       file-type hab_biomass   file-type ","
;       file-type gr_length     file-type ","
;       ]  
;  file-type time_fly           file-type ","
;  file-print pr_feed]
;
;file-close
;end
;
;to output-debugging-geese-movement-roost-rules
;  file-open "debugging_roost_rules.csv"
;  
;  
;  ask geese
;  [
;    if any? fields-here 
;   [
;   let closest_field min-one-of fields-here [distance myself]
;  
;  
;     ask closest_field
;     [
;       file-type ticks         file-type ","
;       file-type hab           file-type ","
;       file-type day           file-type ","
;       ]
; 
; file-type  pr_feed            file-type ","
; file-type time_fly            file-type ","
; file-type cum_time_feeding    file-type ","
; file-type dnei        file-type ","
; file-type dee         file-type ","
; file-print deg         ;file-type ","
;  ]
;   
;  if any? roosts-here
;    [let closest_roost min-one-of roosts-here [distance myself]
;         ask closest_roost
;     [
;       file-type ticks         file-type ","
;       file-type hab           file-type ","
;       file-type day           file-type ","
;       ]
; 
; file-type  pr_feed            file-type ","
; file-type time_fly            file-type ","
; file-type cum_time_feeding    file-type ","
; file-type dnei        file-type ","
; file-type dee         file-type ","
; file-print deg         ;file-type ","
;  ]]
;
;file-close
;end


 ; 1. I first checked whether the number of grass fields for each period is the same as in the input file and it is

;to plot_n_grass
;  set-current-plot "N of grass fields"
;  plot count fields with [hab = "grass"]
;end

; 2. below I check the average grass length per field. Some grass fields are resown and 'appear' from the stubble fields (undersown stubble) and these new grass fields starts from grass length 1cm and not 2 cm as the grass fields from the previous season
; that is why there are these 'steps' in the graph 

to plot_grass_growth ; if I plot grass length without grass-growth procedure and setting grass value to 2 instead of 1cm, than the values is constant
  set-current-plot "Grass length"
  let grass_length []
  set grass_length ((sum [hab_biomass / area / 16.4] of fields with [hab = "grass"])  / count fields with [hab = "grass"])
  plot grass_length
end
;
to plot_stubble_m2 ; it should be constant if there is no depletion and it is
  set-current-plot "Stubble m2"
  let grain_dens_stubble ((sum [grain_dens] of fields with [hab = "stubble"]) / (count fields with [hab = "stubble"]))
  plot grain_dens_stubble
  ;print grain_count
end

to plot_grain_m2 ; it should be constant if there is no depletion and it is
  set-current-plot "Grain m2"
  let grain_dens_grain 0
  if any? fields with [hab = "grain"]
   [set grain_dens_grain ((sum [grain_dens] of fields with [hab = "grain"]) / (count fields with [hab = "grain"]))]
  plot grain_dens_grain
end


;
;to plot_plough ; it should be constant and it is
;  set-current-plot "Plough"
;  let plough []
;  set plough ((sum [hab_biomass / area ] of fields with [hab = "plough"]) / (count fields with [hab = "plough"]))
;  plot plough
;end
@#$#@#$#@
GRAPHICS-WINDOW
185
10
767
717
-1
-1
0.4
1
10
1
1
1
0
0
0
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
72
10
135
43
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
68
43
setups
setups
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
4
48
59
81
step
go
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
61
48
116
81
10 steps
repeat 10 [go]
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
118
48
181
81
480 step
repeat 480 [go]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
892
546
1065
669
Grass length
NIL
NIL
0.0
1056.0
0.0
10.0
true
false
"" ""
PENS
"pen-0" 1.0 0 -7500403 true "" ""

CHOOSER
5
87
148
132
decision_making_rule
decision_making_rule
"random" "max energy" "asocial learning" "social learning" "asocial and social"
0

MONITOR
18
139
75
184
starved
n-of-starved
17
1
11

MONITOR
85
140
142
185
north
n-of-vesteralen
17
1
11

PLOT
1072
547
1246
673
Stubble m2
NIL
NIL
0.0
1056.0
0.0
10.0
true
false
"" ""
PENS
"pen-0" 1.0 0 -7500403 true "" ""

PLOT
890
683
1064
820
Grain m2
NIL
NIL
0.0
1056.0
0.0
10.0
true
false
"" ""
PENS
"pen-0" 1.0 0 -7500403 true "" ""

BUTTON
34
197
120
230
pen down
ask geese [pd]
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
9
271
179
407
Mean DNEI
days
kJ
0.0
50.0
0.0
2000.0
true
false
"" ""
PENS
"pen-0" 1.0 0 -7500403 true "" ""

PLOT
9
556
179
694
Count super-geese
days
N of geese
0.0
50.0
0.0
2300.0
true
false
"" ""
PENS
"pen-0" 1.0 0 -7500403 true "" ""

PLOT
9
414
179
549
Time on day roost
days
prop of daylight
0.0
50.0
0.0
1.0
true
false
"" ""
PENS
"pen-0" 1.0 0 -7500403 true "" ""

TEXTBOX
895
521
1045
542
Plots for depugging
11
0.0
1

PLOT
9
725
331
875
Distance to roost
days
km
0.0
50.0
0.0
5.0
true
true
"" ""
PENS
"morning" 1.0 0 -2674135 true "" ""
"midday" 1.0 0 -11221820 true "" ""
"afternoon" 1.0 0 -16777216 true "" ""

TEXTBOX
895
10
1241
38
Parameteers used in sensitivity analysis and parameterisation, see description in ODD and globals at the top of the code
11
0.0
1

INPUTBOX
1182
371
1321
431
n_geese_in_super_goose
20
1
0
Number

INPUTBOX
1184
436
1322
496
big_radius
5
1
0
Number

INPUTBOX
1035
242
1170
302
small_radius
1
1
0
Number

INPUTBOX
1033
309
1170
369
alpha
0.036
1
0
Number

INPUTBOX
1036
114
1172
174
av_ee
54.28
1
0
Number

INPUTBOX
1036
178
1172
238
ee_fly
416.35
1
0
Number

INPUTBOX
1181
179
1323
239
starvation_energy_stores
9620
1
0
Number

INPUTBOX
1033
372
1171
432
prob-of-disturbance-small-fields
0.3
1
0
Number

INPUTBOX
1181
50
1322
110
prob-of-disturbance-large-fields
0.2
1
0
Number

INPUTBOX
1038
51
1174
111
cum_time_threshold
3
1
0
Number

INPUTBOX
891
309
1029
369
field_area_threshold
60000
1
0
Number

INPUTBOX
1181
305
1323
365
max_steps-to_stay_on_roost
3
1
0
Number

INPUTBOX
893
375
1029
435
extra_initial_energy
0
1
0
Number

INPUTBOX
1178
241
1321
301
extra_leaving_energy
0
1
0
Number

MONITOR
1072
685
1162
730
# disturbance
n_disturbance_event
0
1
11

MONITOR
1178
685
1235
730
geese
count geese
17
1
11

INPUTBOX
892
51
1032
111
initial_grass_length
3
1
0
Number

INPUTBOX
892
115
1030
175
initial_grain_density_grain
17
1
0
Number

INPUTBOX
891
181
1031
241
grain_mean_stubble
22
1
0
Number

INPUTBOX
892
245
1030
305
grain_sd_stubble
22
1
0
Number

INPUTBOX
1181
114
1321
174
temp_disturbance
0
1
0
Number

INPUTBOX
893
442
1029
502
roost_dist_radius
0.1
1
0
Number

INPUTBOX
1038
441
1172
501
steps_on_500m_roost
1
1
0
Number

SLIDER
425
783
597
816
expindex
expindex
0
1000
0
1
1
NIL
HORIZONTAL

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

moon
false
0
Polygon -7500403 true true 175 7 83 36 25 108 27 186 79 250 134 271 205 274 281 239 207 233 152 216 113 185 104 132 110 77 132 51

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

sun
false
0
Circle -7500403 true true 75 75 150
Polygon -7500403 true true 300 150 240 120 240 180
Polygon -7500403 true true 150 0 120 60 180 60
Polygon -7500403 true true 150 300 120 240 180 240
Polygon -7500403 true true 0 150 60 120 60 180
Polygon -7500403 true true 60 195 105 240 45 255
Polygon -7500403 true true 60 105 105 60 45 45
Polygon -7500403 true true 195 60 240 105 255 45
Polygon -7500403 true true 240 195 195 240 255 255

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
NetLogo 5.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="calibration-random" repetitions="10" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <timeLimit steps="1057"/>
    <exitCondition>count geese = 0</exitCondition>
    <metric>count geese</metric>
    <metric>day</metric>
    <metric>gl_ticks_since_morning</metric>
    <metric>ticks</metric>
    <metric>mean [dnei] of geese</metric>
    <metric>mean [steps_spent_on-day_roost] of geese</metric>
    <metric>mean [dist_closest_roost] of geese with [on_roost? = false]</metric>
    <metric>n-of-vesteralen</metric>
    <metric>n-of-starved</metric>
    <metric>period</metric>
    <metric>time_of_day</metric>
    <steppedValueSet variable="expindex" first="0" step="1" last="999"/>
    <enumeratedValueSet variable="decision_making_rule">
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="calibration-maxenergy" repetitions="10" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <timeLimit steps="1057"/>
    <exitCondition>count geese = 0</exitCondition>
    <metric>count geese</metric>
    <metric>day</metric>
    <metric>gl_ticks_since_morning</metric>
    <metric>ticks</metric>
    <metric>mean [dnei] of geese</metric>
    <metric>mean [steps_spent_on-day_roost] of geese</metric>
    <metric>mean [dist_closest_roost] of geese with [on_roost? = false]</metric>
    <metric>n-of-vesteralen</metric>
    <metric>n-of-starved</metric>
    <metric>period</metric>
    <metric>time_of_day</metric>
    <steppedValueSet variable="expindex" first="0" step="1" last="999"/>
    <enumeratedValueSet variable="decision_making_rule">
      <value value="&quot;max energy&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="calibration-asociallearning" repetitions="10" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <timeLimit steps="1057"/>
    <exitCondition>count geese = 0</exitCondition>
    <metric>count geese</metric>
    <metric>day</metric>
    <metric>gl_ticks_since_morning</metric>
    <metric>ticks</metric>
    <metric>mean [dnei] of geese</metric>
    <metric>mean [steps_spent_on-day_roost] of geese</metric>
    <metric>mean [dist_closest_roost] of geese with [on_roost? = false]</metric>
    <metric>n-of-vesteralen</metric>
    <metric>n-of-starved</metric>
    <metric>period</metric>
    <metric>time_of_day</metric>
    <steppedValueSet variable="expindex" first="0" step="1" last="999"/>
    <enumeratedValueSet variable="decision_making_rule">
      <value value="&quot;asocial learning&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="calibration-sociallearning" repetitions="10" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <timeLimit steps="1057"/>
    <exitCondition>count geese = 0</exitCondition>
    <metric>count geese</metric>
    <metric>day</metric>
    <metric>gl_ticks_since_morning</metric>
    <metric>ticks</metric>
    <metric>mean [dnei] of geese</metric>
    <metric>mean [steps_spent_on-day_roost] of geese</metric>
    <metric>mean [dist_closest_roost] of geese with [on_roost? = false]</metric>
    <metric>n-of-vesteralen</metric>
    <metric>n-of-starved</metric>
    <metric>period</metric>
    <metric>time_of_day</metric>
    <steppedValueSet variable="expindex" first="0" step="1" last="999"/>
    <enumeratedValueSet variable="decision_making_rule">
      <value value="&quot;social learning&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="calibration-asociaandsocial" repetitions="10" runMetricsEveryStep="true">
    <setup>setups</setup>
    <go>go</go>
    <timeLimit steps="1057"/>
    <exitCondition>count geese = 0</exitCondition>
    <metric>count geese</metric>
    <metric>day</metric>
    <metric>gl_ticks_since_morning</metric>
    <metric>ticks</metric>
    <metric>mean [dnei] of geese</metric>
    <metric>mean [steps_spent_on-day_roost] of geese</metric>
    <metric>mean [dist_closest_roost] of geese with [on_roost? = false]</metric>
    <metric>n-of-vesteralen</metric>
    <metric>n-of-starved</metric>
    <metric>period</metric>
    <metric>time_of_day</metric>
    <steppedValueSet variable="expindex" first="0" step="1" last="999"/>
    <enumeratedValueSet variable="decision_making_rule">
      <value value="&quot;asocial and social&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 1.0 0.0
0.0 1 1.0 0.0
0.2 0 1.0 0.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
