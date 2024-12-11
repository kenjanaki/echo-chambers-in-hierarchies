;; 2-D representation of a belief propagation network written by Toby D. Pilditch and Jan-Philipp Fränken (2019)


;;##################################################
; EXTENSIONS
;;##################################################


extensions [
  stats
]


;;##################################################
; VARIABLES
;;##################################################


globals [
  linked-neighbours       ;; how many agents an agent is linked to
  number-agents           ;; total number of agents
  neutral-agents          ;; agents that have not committed to an opinion (keeping their prior value)
  opinion-A-agents        ;; agents holding opinion A
  opinion-B-agents        ;; agents holding opinion B
  match-counter           ;; checking whether number of  A agents = number of B agents and counts how many times there was no change (i.e. totals have been the same)
  tick-tot-odd            ;; counters used to check when no further learning takes place (i.e. if tick-tot-odd and tick-tot-even are equal)
  tick-tot-even
  peak-spread
  cl-prop-same            ;; mean of proportion of like-minded out of total neighbours
  cl-prop-diff            ;; mean of proportion of anti-minded out of total neighbours
  dist
  help
  generation-count
  num
  initial-turts
  central-node
]

turtles-own [
  p-H                     ;; prior belief of agent (i.e. starting point for learning)
  communication-direction ;; variable that stores the direction an agent has been 'pushed' to (i.e. what it was told by its neighbor before updating); it dictates the flipping in the BMSC function
  p-E-mean-red            ;; stores mean p-E of neighbors holding opinion 1
  p-E-mean-blue           ;; stores mean p-E of neighbors holding opinion 2
  p-T-mean-red            ;; stores mean p-T of neighbors holding opinion 1
  p-T-mean-blue           ;; stores mean p-T of neighbors holding opinion 2
  SCval1                  ;; count of neighbors holding same opinion as communicating agent or target agent dependent on the conditional used
  SCval2                  ;; count of neighbors holding different opinion as communicating agent or target agent dependent on the conditional used
  p-H-R                   ;; posterior belief of agent dictating posterior-opinion
  posterior-opinion       ;; posterior opinion after learning process
  prior-con               ;; confidence in prior-opinion
  p-E                     ;; p-Expertise of persuader / other agent that communicates opinion
  p-T                     ;; p-Trustworthiness of persuader / other agent that communicates opinion
  p-E-weighted            ;; conformity weighted expertise value used in the BSCM
  p-T-weighted            ;; conformity weighted trustworthiness value used in the BSCM

  network-size            ;; size of each agent's network (i.e. number of links)
  opinion-threshold       ;; threshold that defines individual levels of scepticism
  recieve-belief          ;; checks whether the learner has received an opinion from the last tick
  belief-state            ;; whether the agent holds an opinion (0 = neutral agent, 1 = opinion-A-agent, 2 = opinion-B-agent)
  com-num                 ;; counts how often a learner has communicated

  clust-num-same          ;; value for number of "like-minded" linked neighbours
  clust-num-diff          ;; value for number of "anti-minded" linked neighbours
  clust-prop-same         ;; proportion of like-minded out of total neighbours
  clust-prop-diff         ;; proportion of anti-minded out of total neighbours
  nearest-learners

  center-node?
  generation
  model?

  parent
  copy-number

]

patches-own []

breed [learners learner]
breed [ex-learners ex-learner]



;;##################################################
;; SETUP
;;##################################################


to setup
  clear-all
  if(use-seed?)[random-seed seed]
  ;setup-patches          ;; function setting up patches
  setup-turtles          ;; function setting up turtles
  reset-ticks
end

to setup-priors
  ;; setting up P(H), p(E), and P(T)
  ask turtles[
    set p-H random-normal prior-mean prior-sd
    set p-E random-normal prior-mean prior-sd
    set p-T random-normal prior-mean prior-sd
    while [p-H > 1 or p-H < 0][
      set p-H random-normal prior-mean prior-sd]
    while [p-E > 1 or p-E < 0][
      set p-E random-normal prior-mean prior-sd]
    while [p-T > 1 or p-T < 0][
      set p-T random-normal prior-mean prior-sd]]
end


to setup-links
  if (network != "hierarchical1")[
    let link-counter 0
    ask turtles [
      ifelse prox-YN [
        set nearest-learners min-n-of network-size (other turtles) [ distance myself ]
      ][
        set nearest-learners n-of network-size (other turtles)]
      let near-learners turtle-set (nearest-learners)
      create-links-with near-learners
    ]
  ]
  ask links
  [set color grey - 1
    set thickness 0.1]
end


to setup-turtles
  if network = "random" [                        ;; random network
    create-learners n_agents
    [ setxy ((random (2 * max-pxcor)) + min-pxcor)
      ((random (2 * max-pycor)) + min-pycor)
      set-characteristics
    ]
  ]

  if network = "scale-free" [                     ;; Scale free network setup make the initial network of two turtles and an edge
    make-node nobody                              ;; first node, unattached
    make-node turtle 0                            ;; second node, attached to first node
    let SF-create-count 2
    while [SF-create-count < n_agents]            ;;loop for generating scale-free network
    [ask links [ set color gray ]
      make-node find-partner
      layout
      set SF-create-count (SF-create-count + 1)
    ]
    ask turtles
        [set-characteristics
    ]
    ask links                                   ;; Severing "builder" links
    [die]
  ]

  if network = "hierarchical1" [
    setup-hierarchy
    repeat max_links[create-generation1]
    ask turtles
    [set-characteristics
    ]
  ]

  if network = "hierarchical2" [
    setup-hierarchy
    repeat max_links[create-generation2]
    ask turtles
    [
      set-characteristics
    ]

  ]

  ;
  ifelse(network = "hierarchical1" or network = "hierarchical2")[
    ask central-node[ ifelse neut-event-YN
      [set color green
        set belief-state 3]
      [set color red
        set belief-state 1]
      ifelse n_init_believers < 2
      [setxy 0 0]
      [setxy ((random (2 * max-pxcor)) + min-pxcor)
        ((random (2 * max-pycor)) + min-pycor)]
      set shape "circle"
      set size 1
      set network-size max_links                ;;idea: (random maxLinks) + 1 *** Could use random-gamma alpha lambda (where alpha is mean/variance, and lambda is 1/(variance/mean))... reflects social media connections (few with lots...most with some)
      set recieve-belief 3
      set com-Num 0
    ]
  ]
  [
    create-ex-learners n_init_believers[;;create single ex-learner believer in closest to center (or random pos?)
      ifelse neut-event-YN
      [set color green
        set belief-state 3]
      [set color red
        set belief-state 1]
      ifelse n_init_believers < 2
      [setxy 0 0]
      [setxy ((random (2 * max-pxcor)) + min-pxcor)
        ((random (2 * max-pycor)) + min-pycor)]
      set shape "circle"
      set size 1
      set network-size max_links                ;;idea: (random maxLinks) + 1 *** Could use random-gamma alpha lambda (where alpha is mean/variance, and lambda is 1/(variance/mean))... reflects social media connections (few with lots...most with some)
      set recieve-belief 3
      set com-Num 0
    ]
  ]

  setup-priors
  setup-links
  global-initialization
end


to set-characteristics
  set color grey
  set shape "circle"
  set size .5
  set belief-state 0
  set network-size -1
  while [( network-size  < 1) or ( network-size  > 999)][
    set network-size random-normal max_links round(max_links / 2)];(random maxLinks) + 1]                                   ;;connectivity density dependent
  set recieve-belief 0
  set com-Num 0
  set opinion-threshold 0.5
  set posterior-opinion 0 ;; placeholder for posterior regarding belief
end


to global-initialization
  set opinion-A-agents count turtles with [color = red]
  set opinion-B-agents count turtles with [color = blue]
  set number-agents count turtles
  set neutral-agents count turtles with [color = grey]
  set tick-tot-odd 1
  set tick-tot-even 0
end



;;##################################################
;; GO
;;##################################################



to go
  ifelse tick-tot-odd = tick-tot-even ;; if there are no more agents getting updated
    [set match-counter match-counter + 1]
  [set match-counter 0]
  if match-counter > 1
    [ stop ]
  ask turtles with [recieve-belief = 1 or recieve-belief = 2]
    [ learn ]
  ask turtles with [belief-state > 0]
  [ propagate ]
  ask turtles
    [
      if belief-state = 1
      [ set clust-num-same count link-neighbors with [color = red]
        set clust-num-diff count link-neighbors with [color = blue]
        if clust-num-same > 0
        [ set clust-prop-same ((clust-num-same / (count link-neighbors)) * 100)]
        if clust-num-diff > 0
        [ set clust-prop-diff ((clust-num-diff / (count link-neighbors)) * 100)]
      ]
      if belief-state = 2
      [ set clust-num-same count link-neighbors with [color = blue]
        set clust-num-diff count link-neighbors with [color = red]
        if clust-num-same > 0
        [ set clust-prop-same ((clust-num-same / (count link-neighbors)) * 100)]
        if clust-num-diff > 0
        [ set clust-prop-diff ((clust-num-diff / (count link-neighbors)) * 100)]
      ]
  ]
  set opinion-A-agents count turtles with [color = red]
  set opinion-B-agents count turtles with [color = blue]
  set neutral-agents count turtles with [color = grey]
  set cl-prop-same mean [clust-prop-same] of turtles
  set cl-prop-diff mean [clust-prop-diff] of turtles
  ;;; Peak Spread
  if ((abs(tick-tot-odd - tick-tot-even)) / n_agents) * 100 > peak-spread
    [set peak-spread ((abs(tick-tot-odd - tick-tot-even)) / n_agents) * 100]
  ;;;
  ifelse (remainder(ticks + 2) 2) = 1
  [set tick-tot-odd opinion-A-agents + opinion-B-agents]
  [set tick-tot-even opinion-A-agents + opinion-B-agents]
  tick
end



;;##################################################
;; LEARNING AND PROPAGATING
;;##################################################


to learn
    set posterior-opinion BSCm-integration p-H p-E-weighted p-T-weighted recieve-belief

  ifelse (random 2 = 1)                                      ;; given enforced coding of >=/=, this prevents one-sided declaration level  ;; evaluation-process resulting in a belief-state congruent with the posterior-opinion
    [ifelse posterior-opinion < 0.5
      [  ;; believers
        set belief-state 2
        set recieve-belief 3
        set color blue
      ]
      [ ;; refuters
        set belief-state 1
        set recieve-belief 3
        set color red
      ]
  ]
  [
    ifelse posterior-opinion <= 0.5
    [  ;; believers
      set belief-state 2
      set recieve-belief 3
      set color blue
    ]
    [ ;; refuters
      set belief-state 1
      set recieve-belief 3
      set color red
    ]
  ]

end


to-report BSCm-integration [p-Hyp p-Exp p-Trust flip]   ;; function that is used for contrasting the prior-opinion with evidence to form posterior using modified version of Bayes theorem

  let expertise_influence_spec expertise_influence

  let p-nH (1 - p-Hyp)     ;; anteriors
  let p-nE (1 - p-Exp)
  let p-nT (1 - p-Trust)


  let p-H-E-T 0.5
  let p-H-nE-T 0.5
  let p-H-nE-nT 0.5
  let p-H-E-nT 0.5

  let p-nH-E-T 0.5
  let p-nH-nE-T 0.5
  let p-nH-nE-nT 0.5
  let p-nH-E-nT 0.5


  if (flip = 2) [
    set p-H-E-T 0.5 - (expertise_influence_spec * 2) ;; specifying CPT
    set p-H-nE-T 0.5 - expertise_influence_spec
    set p-H-nE-nT 0.5 + expertise_influence_spec
    set p-H-E-nT 0.5 + (expertise_influence_spec * 2)

    set p-nH-E-T 1 - (0.5 - (expertise_influence_spec * 2))
    set p-nH-nE-T 1 - (0.5 - expertise_influence_spec)
    set p-nH-nE-nT 1 - (0.5 + expertise_influence_spec)
    set p-nH-E-nT 1 - (0.5 + (expertise_influence_spec * 2))

  ]

  if (flip = 1) [
    set p-H-E-T 0.5 + (expertise_influence_spec * 2) ;; specifying CPT
    set p-H-nE-T 0.5 + expertise_influence_spec
    set p-H-nE-nT 0.5 - expertise_influence_spec
    set p-H-E-nT 0.5 - (expertise_influence_spec * 2)

    set p-nH-E-T 1 - (0.5 + (expertise_influence_spec * 2))
    set p-nH-nE-T 1 - (0.5 + expertise_influence_spec)
    set p-nH-nE-nT 1 - (0.5 - expertise_influence_spec)
    set p-nH-E-nT 1 - (0.5 - (expertise_influence_spec * 2))

  ]


  ;; BAYESIAN SOURCE CREDIBILITY MODEL INTEGRATION


  ;; Calculating P(Rep|H):
  let p-R-H ((p-H-E-T * p-Exp * p-Trust) + (p-H-nE-T * p-nE * p-Trust) + (p-H-nE-nT * p-nE * p-nT) + (p-H-E-nT * p-Exp * p-nT))
  ;; Calculating P(Rep|¬H):
  let p-R-nH ((p-nH-E-T * p-Exp * p-Trust) + (p-nH-nE-T * p-nE * p-Trust) + (p-nH-nE-nT * p-nE * p-nT) + (p-nH-E-nT * p-Exp * p-nT))

  ;; plugging into BSCm Central Equation to calculate P(H|Rep):
  set p-H-R (p-H * p-R-H)/(p-H * p-R-H + p-nH * p-R-nH)


  report p-H-R

end


to propagate                  ;; reaching out to other agents and communicating opinion

  if any? link-neighbors with [color = grey] [
    if belief-state = 1 and (random 100 < (prop-likelihood * 100))[                   ;;P(declaration)
      pick-target (1) (red) (blue) (1)
    ]
    if belief-state = 2 and (random 100 < (prop-likelihood * 100)) [
      pick-target (2) (blue) (red) (-1)
    ]
    ;; neutral starting point - how does a neutral agent communicate

    if belief-state = 3 [
      let com-target link-neighbors with [ color = grey ]
      if com-target != nobody [
        let own-No [who] of self
        let targ-No [who] of com-target
        ask com-target [
          set recieve-belief (random 2 + 1)
          ask link ([who] of self) own-No
          [set color green]
        ]
      ]
      set com-Num com-Num + 1
    ]
  ]
end


to pick-target [ bel-state opin-col-A opin-col-B pos-neg ]
  let com-target link-neighbors with [ color = grey ]
  if com-target != nobody [
    let own-No [who] of self
    let targ-No [who] of com-target
    ask com-target [
      set recieve-belief bel-state



      ;; counting neighbors congruent what has been communicated and computing their weighted mean P(E) and P(T) values



      if (recieve-belief = 2) [     ;; if they communicated that they believe blue (i.e. zero), here encoded as receive-belief = 2

        set SCval1 count link-neighbors with [color = blue]
        ifelse (SCval1 >= 1)[
          set p-E-mean-blue mean [p-E] of link-neighbors with [color = blue]
          set p-T-mean-blue mean [p-T] of link-neighbors with [color = blue]]
        [set p-E-mean-blue 0.5
          set p-T-mean-blue 0.5]

        set SCval2 count link-neighbors with [color = red]
        ifelse (SCval2 >= 1) [
          set p-E-mean-red mean [p-E] of link-neighbors with [color = red]
          set p-T-mean-red mean [p-T] of link-neighbors with [color = red]]
        [set p-E-mean-red 0.5
          set p-T-mean-red 0.5]

        set p-E-weighted ((SCval1) / (SCval1 + SCval2))
        set p-T-weighted ((SCval1) / (SCval1 + SCval2))
        if (modulate-weight-by-mean = true)[
          set p-E-weighted ((SCval1 * p-E-mean-blue) / (SCval1 * p-E-mean-blue + SCval2 * p-E-mean-red))
          set p-T-weighted ((SCval1 * p-T-mean-blue) / (SCval1 * p-T-mean-blue + SCval2 * p-T-mean-red))] ]

      if (recieve-belief = 1) [
        set SCval1 count link-neighbors with [color = red]    ;; if they communicated that they believe red (i.e. 1), here encoded as receive-belief = 1
        ifelse (SCval1 >= 1)[
          set p-E-mean-red mean [p-E] of link-neighbors with [color = red]
          set p-T-mean-red mean [p-T] of link-neighbors with [color = red]]
        [set p-E-mean-red 0.5
          set p-T-mean-red 0.5]

        set SCval2 count link-neighbors with [color = blue]
        ifelse (SCval2 >= 1) [
          set p-E-mean-blue mean [p-E] of link-neighbors with [color = blue]
          set p-T-mean-blue mean [p-T] of link-neighbors with [color = blue]]
        [set p-E-mean-blue 0.5
          set p-T-mean-blue 0.5]


        set p-E-weighted ((SCval1) / (SCval1 + SCval2))
        set p-T-weighted ((SCval1) / (SCval1 + SCval2))
        if (modulate-weight-by-mean = true)[
          set p-E-weighted ((SCval1 * p-E-mean-red) / (SCval1 * p-E-mean-red + SCval2 * p-E-mean-blue))
          set p-T-weighted ((SCval1 * p-T-mean-red) / (SCval1 * p-T-mean-red + SCval2 * p-T-mean-blue))]]

      if (no-social-influence = true)[
        set p-E-weighted random-float 1
        set p-T-weighted random-float 1]
      ask link ([who] of self) own-No
      [set color opin-col-A]
    ]

  ]

end

;;##################################################
;; FUNCTIONS SOURCED FROM OTHER MODELS
;;##################################################

;;## 1 ## Wilensky, U. (2005).  NetLogo Preferential Attachment model.  http://ccl.northwestern.edu/netlogo/models/PreferentialAttachment. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

to make-node [old-node]
  create-turtles 1
  [
    set color red
    if old-node != nobody
      [ create-link-with old-node [ set color green ]
        ;; position the new node near its partner
        move-to old-node
        fd 8
    ]
  ]
end

;; This code is the heart of the "preferential attachment" mechanism
to-report find-partner
  report [one-of both-ends] of one-of links
end


;;#####################################
;; HIERARCHY V3
;;#####################################

to setup-hierarchy
  ca
  set dist 5
  set generation-count 1
  set num n_agents
  crt num[
    create-links-with other turtles
    set shape "dot"
    set size 1.5
    set copy-number 0
  ]
  let c 0
  repeat (num - 1)[
    ask one-of turtles with[xcor = 0 and ycor = 0][
      set heading 360 / (num - 1) * c fd dist
    ]
    set c c + 1
  ]
  ask turtles[
    ifelse(xcor = 0 and ycor = 0)[set center-node? true set central-node self][set center-node? false]
    set color red
    set generation 1
    set model? true
  ]
end

to layout
  ;; the number 3 here is arbitrary; more repetitions slows down the
  ;; model, but too few gives poor layouts
  repeat 10 [
    ;; Refactoring the link lengths (MODIFY DENSITY HERE?)
    ;; the more turtles we have to fit into the same amount of space,
    ;; the smaller the inputs to layout-spring we'll need to use
    let factor ((sqrt count turtles) / 3) ;; Here SF-density-mod influences the distance factor across the network - will impact search function...
                                          ;; numbers here are arbitrarily chosen for pleasing appearance
    layout-spring turtles links (1 / factor) (7 / factor) (1 / factor)
    display  ;; for smooth animation
  ]
  ;;; Centering network ;;;
  ;; don't bump the edges of the world
  let x-offset max [xcor] of turtles + min [xcor] of turtles
  let y-offset max [ycor] of turtles + min [ycor] of turtles
  ;; big jumps look funny, so only adjust a little each time
  set x-offset limit-magnitude x-offset 0.1
  set y-offset limit-magnitude y-offset 0.1
  ask turtles [ setxy (xcor - x-offset / 2) (ycor - y-offset / 2) ]
end

to-report limit-magnitude [number limit]
  if number > limit [ report limit ]
  if number < (- limit) [ report (- limit) ]
  report number
end

to create-generation1
  our-hatch
  replicate-links
  ask turtles [set heading (72 * (copy-number - 1)) fd 2 * generation * copy-number]
  let c 0
  ask turtles with[center-node? = false and generation = generation-count][if(c < (num - 1) ^ (generation))[create-link-with central-node set c c + 1]]
  layout
end

to create-generation2
  our-hatch
  replicate-links
  ask turtles [set heading (72 * (copy-number - 1)) fd 2 * generation * copy-number]
  let c 0
  ask turtles with[generation = generation-count][if(c <= (count turtles with [generation = generation-count]) * (p ^ (generation - 1)))[create-link-with one-of turtles with [generation = generation-count - 1] set c c + 1]]
  layout
end

to our-hatch
  set initial-turts count turtles + 1
  set initial-turts initial-turts - 1
  set generation-count generation-count + 1
  let copy 1
  while[count turtles < num * initial-turts][
    let c 0
    while[c < initial-turts][
      ask turtle c [hatch 1[rt 20 fd 2 set parent (who mod (initial-turts)) set generation generation-count set copy-number copy]]
      ;ask turtle c [hatch 1[set parent (who mod (initial-turts)) set generation generation-count set copy-number copy]]
      set c c + 1
    ]
    set copy copy + 1
  ]
end

to replicate-links
  ask turtles with[generation = generation-count][
    let me self
    let c copy-number
    ask turtle parent[
      ask link-neighbors[
        ask turtle (who + (c * initial-turts))[create-link-with me]
      ]
    ]
  ]
end

;;#####################################
@#$#@#$#@
GRAPHICS-WINDOW
380
10
826
457
-1
-1
3.62
1
1
1
1
1
0
1
1
1
-60
60
-60
60
0
0
1
ticks
30.0

BUTTON
235
340
290
409
NIL
setup
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
292
340
347
409
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

SLIDER
115
40
225
73
n_agents
n_agents
0
10000
6.0
1
1
NIL
HORIZONTAL

PLOT
380
460
540
610
P(H)
NIL
NIL
0.0
1.0
0.0
100.0
true
false
"" "set-plot-y-range 0 round(n_agents / 20)"
PENS
"P(H)" 0.01 1 -16777216 true "" "histogram [p-H] of turtles"

PLOT
10
460
370
610
P(T), P(E)
NIL
NIL
0.0
1.0
0.0
2.0
true
true
"" "set-plot-y-range 0 round(n_agents / 20)"
PENS
"p-expertise" 0.01 1 -7858858 true "" "histogram [p-E] of turtles"
"p-trust" 0.01 1 -13840069 true "" "histogram [p-T] of turtles"

TEXTBOX
259
319
337
337
Setup and Go
11
0.0
1

TEXTBOX
120
20
220
40
number of agents or nodes in a module
8
0.0
1

TEXTBOX
125
80
215
106
random or scale-free or hierarchical
8
0.0
1

TEXTBOX
1108
77
1284
103
NIL
10
0.0
1

SLIDER
245
40
355
73
max_links
max_links
0
500
3.0
1
1
NIL
HORIZONTAL

TEXTBOX
250
20
355
46
maximum agent links or iterations
8
0.0
1

SLIDER
600
515
692
548
prior-mean
prior-mean
0
1
0.5
0.01
1
NIL
HORIZONTAL

SLIDER
695
515
787
548
prior-sd
prior-sd
0
1
0.2
0.01
1
NIL
HORIZONTAL

SWITCH
245
170
355
203
neut-event-YN
neut-event-YN
1
1
-1000

SLIDER
115
240
230
273
n_init_believers
n_init_believers
0
100
1.0
1
1
NIL
HORIZONTAL

SWITCH
115
175
230
208
prox-YN
prox-YN
0
1
-1000

SLIDER
245
105
355
138
prop-likelihood
prop-likelihood
0
1
1.0
0.01
1
NIL
HORIZONTAL

PLOT
1035
10
1235
160
Propagation rate
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot (abs(tick-tot-odd - tick-tot-even) / n_agents) * 100"

PLOT
1240
10
1441
160
Proportion Time Plot
Time
Proportion
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"Opinion A" 1.0 0 -2674135 true "" "plot opinion-A-agents "
"Opinion B" 1.0 0 -13345367 true "" "plot opinion-B-agents"

TEXTBOX
263
159
413
177
first event neutral?
8
0.0
1

TEXTBOX
140
155
236
181
close links vs. links everywhere
8
0.0
1

TEXTBOX
140
220
206
246
number of initial opinion holders 
8
0.0
1

TEXTBOX
270
90
420
108
P(declaration)
8
0.0
1

TEXTBOX
1629
16
1779
34
Parameters\n
11
0.0
1

PLOT
830
10
1030
160
Degree of Clustering Post-Cascade
Time
Percentage of Neighbours 
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"Similar " 1.0 0 -16777216 true "" "plot cl-prop-same"
"Different" 1.0 0 -13840069 true "" "plot cl-prop-diff"

TEXTBOX
1180
190
1355
218
CPT for the source credibility model
11
0.0
1

TEXTBOX
1095
240
1149
258
H-true
11
0.0
1

TEXTBOX
1085
280
1140
298
H-n-true
11
0.0
1

TEXTBOX
1145
215
1172
233
T-E
11
0.0
1

TEXTBOX
1205
215
1245
233
T-nE
11
0.0
1

TEXTBOX
1270
215
1313
233
nT-E
11
0.0
1

TEXTBOX
1335
215
1393
233
nT-nE
11
0.0
1

MONITOR
1128
230
1188
267
p-H-E-T
0.5 + (expertise_influence  * 2)
2
1
9

TEXTBOX
1190
345
1340
378
determines the magnitude to which expertise influences opinion formation
9
0.0
1

MONITOR
1128
266
1187
303
p-nH-T-E
1 - (0.5 + (2 * expertise_influence))
2
1
9

MONITOR
1185
230
1250
267
p-H-T-nE
0.5 + expertise_influence
2
1
9

MONITOR
1188
266
1251
303
p-nH-T-nE
1 - (0.5 + expertise_influence)
2
1
9

MONITOR
1250
230
1315
267
p-H-nT-E
0.5 - (expertise_influence * 2)
2
1
9

MONITOR
1250
265
1316
302
p-nH-nT-E
1 - (0.5 - (expertise_influence * 2))
2
1
9

MONITOR
1315
230
1387
267
p-H-nT-nE
0.5 - expertise_influence
2
1
9

MONITOR
1315
265
1387
302
p-nH-nT-nE
1 - (0.5 - expertise_influence)
2
1
9

PLOT
870
205
1042
341
 P(H|rep)
NIL
NIL
0.0
1.0
0.0
100.0
false
false
"" ";set-plot-y-range 0 round(n_agents / 20)"
PENS
"default" 0.01 1 -16777216 true "" "histogram [posterior-opinion] of turtles"

SLIDER
1187
306
1337
339
expertise_influence
expertise_influence
0
1
1.0
0.01
1
NIL
HORIZONTAL

SWITCH
599
555
784
588
no-social-influence
no-social-influence
1
1
-1000

INPUTBOX
140
315
190
375
seed
3.0
1
0
Number

SWITCH
120
380
210
413
use-seed?
use-seed?
1
1
-1000

SWITCH
599
475
784
508
modulate-weight-by-mean
modulate-weight-by-mean
0
1
-1000

CHOOSER
115
100
225
145
network
network
"random" "scale-free" "hierarchical1" "hierarchical2"
3

SLIDER
245
240
355
273
p
p
0
1
1.0
0.1
1
NIL
HORIZONTAL

TEXTBOX
250
220
335
240
fraction for stochastic hierarchy
8
0.0
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
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Replication_2017" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <metric>cl-prop-same</metric>
    <metric>ticks</metric>
    <metric>peak-spread</metric>
    <steppedValueSet variable="max_links" first="1" step="10" last="500"/>
    <enumeratedValueSet variable="min_con">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.58"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="gaussian-prior">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <metric>cl-prop-same</metric>
    <metric>ticks</metric>
    <metric>peak-spread</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0.04"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="4.6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="3.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="33"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (1)" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0.42"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="-0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="social-conformity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (2)" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="169"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="-0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="social-conformity">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (3)" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="169"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.49"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="-0.09"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (4)" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="116"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="-0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (5)" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="134"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (6)" repetitions="29" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="134"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="-0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (7)" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="134"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="-0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (8)" repetitions="30" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="134"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="-0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (9)" repetitions="20" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="29"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="-0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (10)" repetitions="100" sequentialRunOrder="false" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <metric>cl-prop-same</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="-0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (11)" repetitions="200" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (12)" repetitions="500" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (13)" repetitions="500" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (14)" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (15)" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (16)" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (17)" repetitions="500" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (18)" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (19)" repetitions="50" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="-0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (20)" repetitions="20" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (21)" repetitions="20" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (22)" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="23"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (23)" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="99"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (24)" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="99"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="700"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (25)" repetitions="30" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <metric>cl-prop-same</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment (26)" repetitions="500" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <enumeratedValueSet variable="max_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_links">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_con">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-sd">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-learning">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quadratic-con-prior?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth">
      <value value="0.69"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-threshold">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="con-in-prior">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-h-given-c">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_agents">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-Bayes">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ground-truth-impact?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning-through-RL">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-variance">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="scale-free">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="change_shape?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-sd">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-con-mean">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sc-bel-prop">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_con">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-E-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="learning_rate">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="nb-prop-YN">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_con">
      <value value="2.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p-T-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_con">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="beta_T">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="evidence">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max_E">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="min_T">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="Hierarchical Networks" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>opinion-A-agents</metric>
    <metric>opinion-B-agents</metric>
    <metric>cl-prop-same</metric>
    <metric>ticks</metric>
    <metric>peak-spread</metric>
    <steppedValueSet variable="max_links" first="1" step="1" last="3"/>
    <enumeratedValueSet variable="network">
      <value value="&quot;hierarchical1&quot;"/>
      <value value="&quot;hierarchical2&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prop-likelihood">
      <value value="0.1"/>
      <value value="0.5"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n_init_believers">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="neut-event-YN">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prox-YN">
      <value value="true"/>
    </enumeratedValueSet>
    <steppedValueSet variable="n_agents" first="1" step="1" last="5"/>
    <enumeratedValueSet variable="prior-mean">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prior-sd">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="no-social-influence">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <steppedValueSet variable="p" first="0" step="0.5" last="1"/>
    <enumeratedValueSet variable="modulate-weight-by-mean">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="expertise_influence">
      <value value="0"/>
      <value value="0.2"/>
      <value value="0.4"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
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
1
@#$#@#$#@
