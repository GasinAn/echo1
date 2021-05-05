from matplotlib import pyplot as plt
import numpy as np
t = np.arange(37159,37311)

data=[
 13.9197580409484,
 14.8093372788811,
 16.7892457650880,
 18.3745368851591,
 19.7602048353211,
 21.3462303575682,
 23.3427502931939,
 25.4197152197698,
 28.8071227500749,
 31.1952029991782,
 33.8437925247457,
 36.6430028442664,
 39.5728979478026,
 43.1434201420243,
 46.0746795461115,
 49.1667479645348,
 52.2594864828177,
 55.4130298675788,
 58.6273738440170,
 61.7826064015181,
 65.0285861262920,
 68.2054205472399,
 71.1231665673394,
 74.3916957848086,
 77.5311070534904,
 80.6013913468933,
 82.9026033434793,
 86.2046457005281,
 89.3975263156077,
 92.4912964100244,
 95.5958993480767,
 98.3013766190052,
101.527706424318 ,
104.534665102181 ,
106.942795249203 ,
110.101350650887 ,
113.230608645674 ,
116.260594600240 ,
119.271405407830 ,
121.712863288191 ,
124.704865312062 ,
127.787352070415 ,
130.620505372436 ,
133.544084645130 ,
136.658241326599 ,
140.382475940538 ,
143.287225845901 ,
146.512435698997 ,
149.617831997488 ,
152.453536316297 ,
155.769444163112 ,
159.295549856609 ,
162.661668095187 ,
166.357928345044 ,
169.374111852345 ,
172.920315699876 ,
176.376402380892 ,
179.832350925759 ,
183.308077219899 ,
187.023651487651 ,
190.418956776771 ,
193.903784046862 ,
197.398333276686 ,
201.112501907487 ,
204.776304371775 ,
208.319621429165 ,
211.982278300090 ,
215.604362792407 ,
219.316026362973 ,
223.226047937048 ,
227.135977429509 ,
229.866665440091 ,
233.775393176431 ,
237.563416541651 ,
239.880649043775 ,
243.897036368139 ,
247.732580865648 ,
251.897222947583 ,
255.091095821502 ,
258.694230806576 ,
262.216416510311 ,
265.737675893533 ,
269.238007005060 ,
272.627404300658 ,
276.285937154439 ,
279.833705583403 ,
283.360512916240 ,
286.766555004244 ,
290.281651454181 ,
294.135943763271 ,
297.439256794666 ,
300.912779543975 ,
304.463820519494 ,
307.834846267191 ,
311.375254414133 ,
314.984978699556 ,
318.483979512478 ,
321.722262950377 ,
325.540047834453 ,
329.077158684756 ,
332.513797677133 ,
336.029959112887 ,
339.525687553291 ,
342.910983146271 ,
346.475938336645 ,
350.160644109364 ,
353.714876228774 ,
357.218917975792 ,
360.662969621192 ,
364.156763327036 ,
367.660400453776 ,
371.164085299536 ,
374.967805971049 ,
378.361548676001 ,
381.905336214895 ,
385.239340088991 ,
388.633866032459 ,
392.248095848575 ,
395.492865751421 ,
398.667887487350 ,
402.393426979944 ,
405.709444893018 ,
409.325683063783 ,
412.392895335104 ,
415.740470662444 ,
418.938584903161 ,
422.337454721344 ,
425.656938914212 ,
428.997117124920 ,
432.278044287171 ,
435.579792681108 ,
438.862376180859 ,
442.135702686251 ,
445.489882688543 ,
448.664888888397 ,
452.480808258946 ,
455.307617908541 ,
458.305345246073 ,
461.723957234209 ,
464.913524050328 ,
468.384032210095 ,
471.705482158961 ,
475.047800980551 ,
478.461031684861 ,
481.865138398190 ,
485.300144391986 ,
488.716004075668 ,
492.092731333652 ,
495.590262425602 ,
498.988567946745 ,
502.447701147772 ,
506.017587486750 ,
]
srp0=[
 13.9197580409484,
 19.4456212713386,
 18.7367104009782,
 24.6595191630456,
 27.2950453330366,
 29.7679067045389,
 33.5126230342986,
 32.7040650218556,
 39.2908564907440,
 40.7064149420468,
 45.6039157618233,
 47.0752273503336,
 47.4302839071346,
 53.0268303756885,
 54.3222015596315,
 60.9888071604701,
 60.4852743439405,
 62.5879680108816,
 65.9927997949776,
 68.3997809846283,
 75.5132599272821,
 74.1604941155391,
 77.7381313752657,
 78.6565563944048,
 83.0099170738295,
 89.0751313645609,
 88.3858012325178,
 92.4709324589623,
 91.5689799534181,
 97.8665463137835,
101.920362054905 ,
103.209637378117 ,
106.568873680613 ,
105.211257741134 ,
112.441256845944 ,
114.576156301052 ,
118.399848701159 ,
120.180736159201 ,
119.834147115248 ,
126.438481117105 ,
127.650698443648 ,
133.527412111586 ,
133.640140091920 ,
135.380704787312 ,
139.949648282942 ,
141.523801153820 ,
148.181778089710 ,
147.315337483390 ,
151.347798329634 ,
153.381493331180 ,
156.401869774398 ,
162.071617969427 ,
161.567294972796 ,
167.208993626675 ,
167.268203147517 ,
172.032392531873 ,
175.479106725512 ,
176.727332497633 ,
182.675722085544 ,
181.952424921937 ,
187.863320143949 ,
188.825419237456 ,
192.685456803552 ,
197.539878936760 ,
197.684431476719 ,
203.388728039808 ,
202.620948469954 ,
208.870888185205 ,
211.935971016626 ,
214.494337812564 ,
218.498291349293 ,
217.467712049142 ,
224.726727393921 ,
226.388006082604 ,
231.889307053622 ,
233.291404681780 ,
233.533577445274 ,
239.828226697763 ,
241.376178106466 ,
249.004671880460 ,
248.192560872427 ,
250.567709230948 ,
254.226561594445 ,
257.206707600580 ,
265.176791345368 ,
263.620879239843 ,
267.764610990491 ,
268.495799414624 ,
273.804936213238 ,
280.268400639161 ,
279.740554328480 ,
284.325789563163 ,
283.235052235242 ,
290.540173616554 ,
294.504531268777 ,
296.407169556666 ,
300.037224898895 ,
298.866134283076 ,
306.733500787998 ,
308.540898045611 ,
313.259574472711 ,
315.028619845380 ,
315.337708854239 ,
321.998994333627 ,
322.922501816188 ,
329.530100784353 ,
329.551160601951 ,
332.327817352660 ,
336.398107575606 ,
337.953049637829 ,
344.699749544882 ,
344.000831701906 ,
349.190912786119 ,
350.467160670018 ,
353.659612289638 ,
358.701389467054 ,
358.816682749230 ,
365.230047701167 ,
364.676527538299 ,
369.469802955889 ,
371.850216730491 ,
374.040138442720 ,
380.053637651109 ,
379.356234248454 ,
384.732486136086 ,
384.755785956192 ,
389.438509200506 ,
393.877895397036 ,
394.611749911733 ,
399.235081232281 ,
398.090135196291 ,
404.518024488970 ,
407.157995327807 ,
410.221802634032 ,
413.027901182495 ,
412.174108655473 ,
418.750395833368 ,
420.293499904931 ,
425.758548751200 ,
426.392688906972 ,
426.903298708203 ,
432.117420141028 ,
433.776226056569 ,
440.681808920833 ,
439.716731351499 ,
442.041152172890 ,
444.984086706297 ,
447.904490692886 ,
454.796354353863 ,
453.413490115175 ,
457.176220755767 ,
457.751039556058 ,
]
srp1=[
 13.9197580409484,
 17.5861778754492,
 15.5447771650777,
 19.8114262960309,
 21.4422323026945,
 22.6743621827926,
 25.9711834404298,
 25.2661874833450,
 30.7684214715663,
 32.2655422103656,
 36.1290714546046,
 38.3407067724402,
 39.4110705882526,
 44.4760347368938,
 46.3207423672362,
 51.6589782629327,
 52.8363477522831,
 55.4306905559075,
 59.0587788942672,
 61.7454743487382,
 67.2975646661356,
 68.1186661474960,
 71.5369080542359,
 73.8199585239674,
 77.4384557227941,
 82.3993332615754,
 83.4734702036463,
 87.0760910385819,
 88.6075713428320,
 92.8672415464688,
 96.9941750133163,
 98.5899060585367,
102.060277886833 ,
103.401753546107 ,
107.902613026989 ,
111.344542850304 ,
113.468585147000 ,
116.780176964247 ,
118.246102153998 ,
122.749603220760 ,
125.777480940234 ,
128.331624416534 ,
131.607332594833 ,
133.339092352450 ,
137.767094302397 ,
140.603231685446 ,
143.525238051840 ,
146.927972557443 ,
148.986106399986 ,
153.368658172314 ,
156.146377745000 ,
159.424642547306 ,
163.089437145961 ,
165.484598617056 ,
169.862832204061 ,
172.639369175616 ,
176.241861087867 ,
180.163976626844 ,
182.781749567742 ,
187.104527972503 ,
189.856165745013 ,
193.697452111866 ,
197.844245306025 ,
200.556085003378 ,
204.778801974902 ,
207.515919464233 ,
211.487278127840 ,
215.832980200046 ,
218.532870354185 ,
222.619661009771 ,
225.394649674529 ,
229.371389698472 ,
233.861506744545 ,
236.471568957383 ,
240.374170198953 ,
243.284152678782 ,
247.134509049223 ,
251.677826354880 ,
254.200828778709 ,
257.854830844536 ,
261.051893236357 ,
264.669206109389 ,
269.128250060949 ,
271.716463867364 ,
275.050535275590 ,
278.695528272995 ,
282.067251050932 ,
286.217960136519 ,
289.123218511326 ,
292.093251080875 ,
296.201241517319 ,
299.488199189332 ,
303.079373972799 ,
306.490541319463 ,
309.207682576877 ,
313.498071397692 ,
317.012066161298 ,
319.994496963470 ,
323.802657700586 ,
326.620899146581 ,
330.567701044447 ,
334.505621512296 ,
337.264405389788 ,
340.986097420077 ,
344.320636173095 ,
347.568263367316 ,
351.643155050499 ,
354.855680026922 ,
358.020959381652 ,
361.901705867595 ,
364.822929058919 ,
368.414962601162 ,
372.335516226931 ,
375.116473192358 ,
378.845150172302 ,
382.283709013403 ,
385.360938303183 ,
389.170392318322 ,
392.290475729619 ,
395.153451781251 ,
399.085589337691 ,
402.559753932110 ,
405.445267439364 ,
408.993530437513 ,
411.740402563243 ,
415.168860507565 ,
419.098476000005 ,
421.824736682188 ,
424.818217704809 ,
428.530708166553 ,
431.823407698946 ,
434.856804583690 ,
438.184366782059 ,
440.725041784867 ,
444.359450127424 ,
448.203183857816 ,
450.752499904492 ,
453.826913007806 ,
457.280810679200 ,
460.493153289939 ,
463.751623112278 ,
467.002288939771 ,
469.600478617316 ,
473.431198457046 ,
477.191229304540 ,
479.793840664434 ,
483.150798731082 ,
486.524121145052 ,
489.780184244514 ,
493.422076950250 ,
496.723828861594 ,
499.511049858173 ,
]
srp2=[
 13.9197580409484,
 17.5861778754492,
 15.5447771650777,
 19.8114262960309,
 21.4422323026945,
 22.6743621827926,
 25.9711834404298,
 25.2661874833450,
 30.7684214715663,
 32.2655422103656,
 36.1290714546046,
 38.3407067724402,
 39.4110705882526,
 44.4760347368938,
 46.3207423672362,
 51.6589782629327,
 52.8363477522831,
 55.4306905559075,
 59.0587788942672,
 61.7454743487382,
 67.2975646661356,
 68.1186661474960,
 71.5369080542359,
 73.8199585239674,
 77.4384557227941,
 82.3993332615754,
 83.4734702036463,
 87.0760910385819,
 88.6075713428320,
 92.8672415464688,
 96.9941750133163,
 98.5899060585367,
102.060277886833 ,
103.401753546107 ,
107.902613026989 ,
111.344542850304 ,
113.468585147000 ,
116.780176964247 ,
118.246102153998 ,
122.749603220760 ,
125.780432614172 ,
128.329836752088 ,
131.609580226559 ,
133.339430349642 ,
137.766738231072 ,
140.614640150888 ,
143.516752720923 ,
146.933289592708 ,
148.984732197508 ,
153.368163731786 ,
156.164524978732 ,
159.410602634755 ,
163.098824625752 ,
165.481739657368 ,
169.862422012522 ,
172.662570872952 ,
176.223504371122 ,
180.178024272105 ,
182.777401563610 ,
187.103054612599 ,
189.883800275561 ,
193.676537527448 ,
197.862670064841 ,
200.550850741431 ,
204.773704909201 ,
207.548827695024 ,
211.465850930844 ,
215.853546801317 ,
218.528925553425 ,
222.607690708895 ,
225.434885209716 ,
229.351903580838 ,
233.878593173543 ,
236.472610984618 ,
240.351823484905 ,
243.334302357136 ,
247.121730490432 ,
251.682586610972 ,
254.211747911882 ,
257.820515966836 ,
261.111025159963 ,
264.672490426681 ,
269.111644071125 ,
271.740155364566 ,
275.008462739877 ,
278.752948470568 ,
282.098278970441 ,
286.178615351522 ,
289.154582299623 ,
292.057858487797 ,
296.235071796372 ,
299.549730504951 ,
303.034246867310 ,
306.512810718991 ,
309.202425861203 ,
313.487706123422 ,
317.079831613547 ,
319.977178913265 ,
323.795658849405 ,
326.660054522154 ,
330.525720132202 ,
334.530688522233 ,
337.295413183056 ,
340.949085731234 ,
344.382390151767 ,
347.557101760385 ,
351.603173361116 ,
354.908121648997 ,
357.990853574100 ,
361.925137828170 ,
364.888985543530 ,
368.378008996760 ,
372.347812670187 ,
375.137166657224 ,
378.795071545264 ,
382.350428855712 ,
385.410327268222 ,
389.127576112145 ,
392.339640685607 ,
395.125514210474 ,
399.050450362965 ,
402.623398635805 ,
405.426657147253 ,
408.980327684375 ,
411.827925341617 ,
415.174560006299 ,
419.067520428735 ,
421.865974032207 ,
424.758560170521 ,
428.551415894984 ,
431.913082088719 ,
434.819632331765 ,
438.185870753979 ,
440.783000497539 ,
444.329537098422 ,
448.186010530809 ,
450.788883977938 ,
453.766489168972 ,
457.350758179614 ,
460.591364049534 ,
463.694936525159 ,
467.017429472901 ,
469.628877422600 ,
473.380084399420 ,
477.207717239915 ,
479.829396571935 ,
483.092006723393 ,
486.620254158866 ,
489.852238963273 ,
493.344390305095 ,
496.760146973661 ,
499.512525430942 ,
]
srp3=[
 13.9197580409484,
 17.6298368494607,
 15.6097119218449,
 19.9046813483254,
 21.5461011963805,
 22.8067035750513,
 26.1100809438055,
 25.3996071654195,
 30.9129577045352,
 32.3973738782197,
 36.2801577446668,
 38.4770506040813,
 39.5365243017877,
 44.6007521807867,
 46.4264486441249,
 51.7813452839044,
 52.9304548887123,
 55.5049682254844,
 59.1361693743122,
 61.7976771876575,
 67.3764709490098,
 68.1565767280193,
 71.5560870448082,
 73.8555695111739,
 77.4398137693056,
 82.4531361601170,
 83.4653248203574,
 87.0652702776825,
 88.6152119563647,
 92.8350103470340,
 97.0494483795311,
 98.5461688761367,
102.044604893790 ,
103.388573819431 ,
107.856463723209 ,
111.418709466080 ,
113.392504331916 ,
116.773272221266 ,
118.212315955273 ,
122.701912002341 ,
125.880403571822 ,
128.220046141062 ,
131.614901944590 ,
133.285791727145 ,
137.715637383645 ,
140.738687885541 ,
143.367789887570 ,
146.947160665995 ,
148.914185916186 ,
153.301991574634 ,
156.310238648683 ,
159.227761922881 ,
163.128127017090 ,
165.408992438863 ,
169.780670758463 ,
172.834731855201 ,
176.014618249939 ,
180.220576361400 ,
182.694290878195 ,
186.975439401558 ,
190.068360968982 ,
193.430932062509 ,
197.886931181133 ,
200.441764496805 ,
204.554205336924 ,
207.742705073575 ,
211.193109228813 ,
215.823109005223 ,
218.407541929961 ,
222.262421502076 ,
225.650859384889 ,
229.089303646243 ,
233.740528100369 ,
236.370844297374 ,
239.862852139277 ,
243.581201680759 ,
246.942652900559 ,
251.370631379428 ,
254.166380733324 ,
257.212350383855 ,
261.341391940995 ,
264.678715088134 ,
268.589028952664 ,
271.744958411610 ,
274.390085675165 ,
278.817671103443 ,
282.342024482158 ,
285.523286642608 ,
289.104239167335 ,
291.641985254968 ,
295.947496289351 ,
299.884643701662 ,
302.491352910812 ,
306.218731925019 ,
309.186769919727 ,
312.880932382944 ,
317.103101832704 ,
319.804496392904 ,
323.173876988391 ,
326.945690638918 ,
330.078317323432 ,
333.938353784767 ,
337.410148025166 ,
340.258466732026 ,
344.416591163116 ,
347.735919240886 ,
350.813766292886 ,
354.790359433065 ,
357.717833443332 ,
361.227972764798 ,
365.195020422293 ,
368.234158771922 ,
371.650634639871 ,
375.280260090261 ,
378.049792297101 ,
381.812755141061 ,
385.693834306947 ,
388.402376991829 ,
392.004202288797 ,
395.370788103527 ,
398.485196282383 ,
402.175424901035 ,
405.333844528212 ,
407.959449756010 ,
411.840535152374 ,
415.552656318534 ,
418.177896005244 ,
421.586468415747 ,
424.721434927449 ,
427.971157150828 ,
431.604325771526 ,
434.643989760037 ,
437.189517967298 ,
441.047132438126 ,
444.730475459362 ,
447.248694673772 ,
450.520305800849 ,
453.646915449695 ,
456.889339487358 ,
460.420621457747 ,
463.502038329261 ,
466.085530639692 ,
469.980796148448 ,
473.681280877147 ,
476.251744061298 ,
479.616225305649 ,
482.950828987851 ,
486.223600705399 ,
489.842614472554 ,
493.159260191847 ,
495.877976775582 ,
499.857091426702 ,
]


plt.plot(t,np.array(srp0),c='red')
#plt.plot(t,np.array(srp1),c='yellow')
plt.plot(t,np.array(srp2),c='green')
plt.plot(t,np.array(srp3),c='blue')
plt.plot(t,np.array(data),c='black')
plt.legend(('NoRP','SRP','ObsData'))
plt.xlabel('MJD')
plt.ylabel('$\\omega$/deg')
#plt.show()
#plt.savefig('omega.png',dpi=800)

print(np.array(srp2)[-1]-np.array(data)[-1])
print(np.array(srp3)[-1]-np.array(data)[-1])
print(abs(np.array(srp3)-np.array(srp2)).max())
plt.clf()
plt.plot(np.array(srp3)-np.array(srp2))
plt.show()
