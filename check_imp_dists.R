#check distribution imputed variables

completedData <- read.csv(paste0(cas_dir,'projects/W1_W2_pubertal_timing/MItests.csv'))

#all missing hormones

hormonesi1 <- completedData[ which(completedData$tagid=='TAG006' & completedData$wave=='2'), ]
hormonesi2 <- completedData[ which(completedData$tagid=='TAG009' & completedData$wave=='2'), ]
hormonesi3 <- completedData[ which(completedData$tagid=='TAG012' & completedData$wave=='1'), ]
hormonesi4 <- completedData[ which(completedData$tagid=='TAG012' & completedData$wave=='2'), ]
hormonesi5 <- completedData[ which(completedData$tagid=='TAG013' & completedData$wave=='2'), ]
hormonesi6 <- completedData[ which(completedData$tagid=='TAG015' & completedData$wave=='2'), ]
hormonesi7 <- completedData[ which(completedData$tagid=='TAG017' & completedData$wave=='1'), ]
hormonesi8 <- completedData[ which(completedData$tagid=='TAG023' & completedData$wave=='2'), ]
hormonesi9 <- completedData[ which(completedData$tagid=='TAG030' & completedData$wave=='2'), ]
hormonesi10 <- completedData[ which(completedData$tagid=='TAG040' & completedData$wave=='2'), ]
hormonesi11 <- completedData[ which(completedData$tagid=='TAG045' & completedData$wave=='1'), ]
hormonesi12 <- completedData[ which(completedData$tagid=='TAG050' & completedData$wave=='1'), ]
hormonesi13 <- completedData[ which(completedData$tagid=='TAG052' & completedData$wave=='1'), ]
hormonesi14 <- completedData[ which(completedData$tagid=='TAG058' & completedData$wave=='2'), ]
hormonesi15 <- completedData[ which(completedData$tagid=='TAG064' & completedData$wave=='2'), ]
hormonesi16 <- completedData[ which(completedData$tagid=='TAG083' & completedData$wave=='2'), ]
hormonesi17 <- completedData[ which(completedData$tagid=='TAG091' & completedData$wave=='2'), ]
hormonesi18 <- completedData[ which(completedData$tagid=='TAG098' & completedData$wave=='2'), ]
hormonesi19 <- completedData[ which(completedData$tagid=='TAG102' & completedData$wave=='2'), ]
hormonesi20 <- completedData[ which(completedData$tagid=='TAG105' & completedData$wave=='2'), ]
hormonesi21 <- completedData[ which(completedData$tagid=='TAG109' & completedData$wave=='2'), ]
hormonesi22 <- completedData[ which(completedData$tagid=='TAG117' & completedData$wave=='1'), ]
hormonesi23 <- completedData[ which(completedData$tagid=='TAG127' & completedData$wave=='2'), ]
hormonesi24 <- completedData[ which(completedData$tagid=='TAG137' & completedData$wave=='2'), ]
hormonesi25 <- completedData[ which(completedData$tagid=='TAG143' & completedData$wave=='2'), ]
hormonesi26 <- completedData[ which(completedData$tagid=='TAG145' & completedData$wave=='2'), ]
hormonesi27 <- completedData[ which(completedData$tagid=='TAG165' & completedData$wave=='2'), ]
hormonesi28 <- completedData[ which(completedData$tagid=='TAG169' & completedData$wave=='2'), ]
hormonesi29 <- completedData[ which(completedData$tagid=='TAG170' & completedData$wave=='2'), ]
hormonesi30 <- completedData[ which(completedData$tagid=='TAG174' & completedData$wave=='1'), ]
hormonesi31 <- completedData[ which(completedData$tagid=='TAG174' & completedData$wave=='2'), ]

shapiro.test(hormonesi1$DHEA_cor)
shapiro.test(hormonesi2$DHEA_cor)
shapiro.test(hormonesi3$DHEA_cor)
shapiro.test(hormonesi4$DHEA_cor)
shapiro.test(hormonesi5$DHEA_cor)
shapiro.test(hormonesi6$DHEA_cor)
shapiro.test(hormonesi7$DHEA_cor)
shapiro.test(hormonesi8$DHEA_cor)
shapiro.test(hormonesi9$DHEA_cor)
shapiro.test(hormonesi10$DHEA_cor)
shapiro.test(hormonesi11$DHEA_cor)
shapiro.test(hormonesi12$DHEA_cor)
shapiro.test(hormonesi13$DHEA_cor)
shapiro.test(hormonesi14$DHEA_cor)
shapiro.test(hormonesi15$DHEA_cor)
shapiro.test(hormonesi16$DHEA_cor)
shapiro.test(hormonesi17$DHEA_cor)
shapiro.test(hormonesi18$DHEA_cor)
shapiro.test(hormonesi19$DHEA_cor)
shapiro.test(hormonesi20$DHEA_cor)
shapiro.test(hormonesi21$DHEA_cor)
shapiro.test(hormonesi22$DHEA_cor)
shapiro.test(hormonesi23$DHEA_cor)
shapiro.test(hormonesi24$DHEA_cor)
shapiro.test(hormonesi25$DHEA_cor)
shapiro.test(hormonesi26$DHEA_cor)
shapiro.test(hormonesi27$DHEA_cor)
shapiro.test(hormonesi28$DHEA_cor)
shapiro.test(hormonesi29$DHEA_cor)
shapiro.test(hormonesi30$DHEA_cor)
shapiro.test(hormonesi31$DHEA_cor)


hist(hormonesi18$DHEA_cor)
hist(hormonesi27$DHEA_cor)

#qqplots

shapiro.test(hormonesi1$TEST_cor)
shapiro.test(hormonesi2$TEST_cor)
shapiro.test(hormonesi3$TEST_cor)
shapiro.test(hormonesi4$TEST_cor)
shapiro.test(hormonesi5$TEST_cor)
shapiro.test(hormonesi6$TEST_cor)
shapiro.test(hormonesi7$TEST_cor)
shapiro.test(hormonesi8$TEST_cor)
shapiro.test(hormonesi9$TEST_cor)
shapiro.test(hormonesi10$TEST_cor)
shapiro.test(hormonesi11$TEST_cor)
shapiro.test(hormonesi12$TEST_cor)
shapiro.test(hormonesi13$TEST_cor)
shapiro.test(hormonesi14$TEST_cor)
shapiro.test(hormonesi15$TEST_cor)
shapiro.test(hormonesi16$TEST_cor)
shapiro.test(hormonesi17$TEST_cor)
shapiro.test(hormonesi18$TEST_cor)
shapiro.test(hormonesi19$TEST_cor)
shapiro.test(hormonesi20$TEST_cor)
shapiro.test(hormonesi21$TEST_cor)
shapiro.test(hormonesi22$TEST_cor)
shapiro.test(hormonesi23$TEST_cor)
shapiro.test(hormonesi24$TEST_cor)
shapiro.test(hormonesi25$TEST_cor)
shapiro.test(hormonesi26$TEST_cor)
shapiro.test(hormonesi27$TEST_cor)
shapiro.test(hormonesi28$TEST_cor)
shapiro.test(hormonesi29$TEST_cor)
shapiro.test(hormonesi30$TEST_cor)
shapiro.test(hormonesi31$TEST_cor)

hist(hormonesi7$DHEA_cor)
hist(hormonesi10$TEST_cor)

#qqplots

shapiro.test(hormonesi1$EST_cor)
shapiro.test(hormonesi2$EST_cor)
shapiro.test(hormonesi3$EST_cor)
shapiro.test(hormonesi4$EST_cor)
shapiro.test(hormonesi5$EST_cor)
shapiro.test(hormonesi6$EST_cor)
shapiro.test(hormonesi7$EST_cor)
shapiro.test(hormonesi8$EST_cor)
shapiro.test(hormonesi9$EST_cor)
shapiro.test(hormonesi10$EST_cor)
shapiro.test(hormonesi11$EST_cor)
shapiro.test(hormonesi12$EST_cor)
shapiro.test(hormonesi13$EST_cor)
shapiro.test(hormonesi14$EST_cor)
shapiro.test(hormonesi15$EST_cor)
shapiro.test(hormonesi16$EST_cor)
shapiro.test(hormonesi17$EST_cor)
shapiro.test(hormonesi18$EST_cor)
shapiro.test(hormonesi19$EST_cor)
shapiro.test(hormonesi20$EST_cor)
shapiro.test(hormonesi21$EST_cor)
shapiro.test(hormonesi22$EST_cor)
shapiro.test(hormonesi23$EST_cor)
shapiro.test(hormonesi24$EST_cor)
shapiro.test(hormonesi25$EST_cor)
shapiro.test(hormonesi26$EST_cor)
shapiro.test(hormonesi27$EST_cor)
shapiro.test(hormonesi28$EST_cor)
shapiro.test(hormonesi29$EST_cor)
shapiro.test(hormonesi30$EST_cor)
shapiro.test(hormonesi31$EST_cor)


#all missing CESDC

cesdi1 <- completedData[ which(completedData$tagid=='TAG004' & completedData$wave=='1'), ]
cesdi2 <- completedData[ which(completedData$tagid=='TAG006' & completedData$wave=='2'), ]
cesdi3 <- completedData[ which(completedData$tagid=='TAG012' & completedData$wave=='2'), ]
cesdi4 <- completedData[ which(completedData$tagid=='TAG013' & completedData$wave=='2'), ]
cesdi5 <- completedData[ which(completedData$tagid=='TAG022' & completedData$wave=='1'), ]
cesdi6 <- completedData[ which(completedData$tagid=='TAG023' & completedData$wave=='2'), ]
cesdi7 <- completedData[ which(completedData$tagid=='TAG024' & completedData$wave=='1'), ]
cesdi8 <- completedData[ which(completedData$tagid=='TAG038' & completedData$wave=='1'), ]
cesdi9 <- completedData[ which(completedData$tagid=='TAG044' & completedData$wave=='1'), ]
cesdi10 <- completedData[ which(completedData$tagid=='TAG077' & completedData$wave=='1'), ]
cesdi11 <- completedData[ which(completedData$tagid=='TAG078' & completedData$wave=='1'), ]
cesdi12 <- completedData[ which(completedData$tagid=='TAG083' & completedData$wave=='2'), ]
cesdi13 <- completedData[ which(completedData$tagid=='TAG095' & completedData$wave=='1'), ]
cesdi14 <- completedData[ which(completedData$tagid=='TAG102' & completedData$wave=='2'), ]
cesdi15 <- completedData[ which(completedData$tagid=='TAG109' & completedData$wave=='2'), ]
cesdi16 <- completedData[ which(completedData$tagid=='TAG114' & completedData$wave=='1'), ]
cesdi17 <- completedData[ which(completedData$tagid=='TAG124' & completedData$wave=='1'), ]
cesdi18 <- completedData[ which(completedData$tagid=='TAG127' & completedData$wave=='2'), ]
cesdi19 <- completedData[ which(completedData$tagid=='TAG137' & completedData$wave=='2'), ]
cesdi20 <- completedData[ which(completedData$tagid=='TAG145' & completedData$wave=='2'), ]
cesdi21 <- completedData[ which(completedData$tagid=='TAG165' & completedData$wave=='2'), ]
cesdi22 <- completedData[ which(completedData$tagid=='TAG170' & completedData$wave=='2'), ]
cesdi23 <- completedData[ which(completedData$tagid=='TAG175' & completedData$wave=='2'), ]
cesdi24 <- completedData[ which(completedData$tagid=='TAG209' & completedData$wave=='2'), ]
cesdi25 <- completedData[ which(completedData$tagid=='TAG215' & completedData$wave=='2'), ]

shapiro.test(cesdi1$CESDC_total)
shapiro.test(cesdi2$CESDC_total)
shapiro.test(cesdi3$CESDC_total)
shapiro.test(cesdi4$CESDC_total)
shapiro.test(cesdi5$CESDC_total)
shapiro.test(cesdi6$CESDC_total)
shapiro.test(cesdi7$CESDC_total)
shapiro.test(cesdi8$CESDC_total)
shapiro.test(cesdi9$CESDC_total)
shapiro.test(cesdi10$CESDC_total)
shapiro.test(cesdi11$CESDC_total)
shapiro.test(cesdi12$CESDC_total)
shapiro.test(cesdi13$CESDC_total)
shapiro.test(cesdi14$CESDC_total)
shapiro.test(cesdi15$CESDC_total)
shapiro.test(cesdi16$CESDC_total)
shapiro.test(cesdi17$CESDC_total)
shapiro.test(cesdi18$CESDC_total)
shapiro.test(cesdi19$CESDC_total)
shapiro.test(cesdi20$CESDC_total)
shapiro.test(cesdi21$CESDC_total)
shapiro.test(cesdi22$CESDC_total)
shapiro.test(cesdi23$CESDC_total)
shapiro.test(cesdi24$CESDC_total)
shapiro.test(cesdi25$CESDC_total)

#all missing SCARED anxiety mean

scaredi1 <- completedData[ which(completedData$tagid=='TAG006' & completedData$wave=='2'), ]
scaredi2 <- completedData[ which(completedData$tagid=='TAG012' & completedData$wave=='2'), ]
scaredi3 <- completedData[ which(completedData$tagid=='TAG013' & completedData$wave=='2'), ]
scaredi4 <- completedData[ which(completedData$tagid=='TAG019' & completedData$wave=='1'), ]
scaredi5 <- completedData[ which(completedData$tagid=='TAG020' & completedData$wave=='1'), ]
scaredi6 <- completedData[ which(completedData$tagid=='TAG022' & completedData$wave=='1'), ]
scaredi7 <- completedData[ which(completedData$tagid=='TAG023' & completedData$wave=='2'), ]
scaredi8 <- completedData[ which(completedData$tagid=='TAG023' & completedData$wave=='1'), ]
scaredi9 <- completedData[ which(completedData$tagid=='TAG024' & completedData$wave=='1'), ]
scaredi10 <- completedData[ which(completedData$tagid=='TAG028' & completedData$wave=='1'), ]
scaredi11 <- completedData[ which(completedData$tagid=='TAG034' & completedData$wave=='1'), ]
scaredi12 <- completedData[ which(completedData$tagid=='TAG038' & completedData$wave=='1'), ]
scaredi13 <- completedData[ which(completedData$tagid=='TAG044' & completedData$wave=='1'), ]
scaredi14 <- completedData[ which(completedData$tagid=='TAG077' & completedData$wave=='1'), ]
scaredi15 <- completedData[ which(completedData$tagid=='TAG083' & completedData$wave=='2'), ]
scaredi16 <- completedData[ which(completedData$tagid=='TAG095' & completedData$wave=='1'), ]
scaredi17 <- completedData[ which(completedData$tagid=='TAG102' & completedData$wave=='2'), ]
scaredi18 <- completedData[ which(completedData$tagid=='TAG109' & completedData$wave=='2'), ]
scaredi19 <- completedData[ which(completedData$tagid=='TAG114' & completedData$wave=='1'), ]
scaredi20 <- completedData[ which(completedData$tagid=='TAG124' & completedData$wave=='1'), ]
scaredi21 <- completedData[ which(completedData$tagid=='TAG127' & completedData$wave=='2'), ]
scaredi22 <- completedData[ which(completedData$tagid=='TAG137' & completedData$wave=='2'), ]
scaredi23 <- completedData[ which(completedData$tagid=='TAG145' & completedData$wave=='2'), ]
scaredi24 <- completedData[ which(completedData$tagid=='TAG165' & completedData$wave=='2'), ]
scaredi25 <- completedData[ which(completedData$tagid=='TAG170' & completedData$wave=='2'), ]
scaredi26 <- completedData[ which(completedData$tagid=='TAG175' & completedData$wave=='2'), ]
scaredi27 <- completedData[ which(completedData$tagid=='TAG209' & completedData$wave=='2'), ]
scaredi28 <- completedData[ which(completedData$tagid=='TAG215' & completedData$wave=='2'), ]
scaredi29 <- completedData[ which(completedData$tagid=='TAG233' & completedData$wave=='1'), ]

shapiro.test(scaredi1$SCARED_anxiety_mean)
shapiro.test(scaredi2$SCARED_anxiety_mean)
shapiro.test(scaredi3$SCARED_anxiety_mean)
shapiro.test(scaredi4$SCARED_anxiety_mean)
shapiro.test(scaredi5$SCARED_anxiety_mean)
shapiro.test(scaredi6$SCARED_anxiety_mean)
shapiro.test(scaredi7$SCARED_anxiety_mean)
shapiro.test(scaredi8$SCARED_anxiety_mean)
shapiro.test(scaredi9$SCARED_anxiety_mean)
shapiro.test(scaredi10$SCARED_anxiety_mean)
shapiro.test(scaredi11$SCARED_anxiety_mean)
shapiro.test(scaredi12$SCARED_anxiety_mean)
shapiro.test(scaredi13$SCARED_anxiety_mean)
shapiro.test(scaredi14$SCARED_anxiety_mean)
shapiro.test(scaredi15$SCARED_anxiety_mean)
shapiro.test(scaredi16$SCARED_anxiety_mean)
shapiro.test(scaredi17$SCARED_anxiety_mean)
shapiro.test(scaredi18$SCARED_anxiety_mean)
shapiro.test(scaredi19$SCARED_anxiety_mean)
shapiro.test(scaredi20$SCARED_anxiety_mean)
shapiro.test(scaredi21$SCARED_anxiety_mean)
shapiro.test(scaredi22$SCARED_anxiety_mean)
shapiro.test(scaredi23$SCARED_anxiety_mean)
shapiro.test(scaredi24$SCARED_anxiety_mean)
shapiro.test(scaredi25$SCARED_anxiety_mean)
shapiro.test(scaredi26$SCARED_anxiety_mean)
shapiro.test(scaredi27$SCARED_anxiety_mean)
shapiro.test(scaredi28$SCARED_anxiety_mean)
shapiro.test(scaredi29$SCARED_anxiety_mean)

hist(scaredi13$SCARED_anxiety_mean)
