from sage.all import FractionField
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import free_module_element as vector
from curve import WeierstrassCurve
from curve import PointWeierstrass

class  Lollipop956451(WeierstrassCurve):

    def __init__(self, p, a, b, r, cofactor,L):
        super().__init__(p, a,b, r, cofactor)
        self.a = a
        self.b = b
        self.D = -160807944 
        self.L = L
        self.cofactor = cofactor
        self.r = r
        M = Matrix([[-L,1], [r,0]])
        self.N = M.LLL()
        self.N_inv = self.N**-1
        self.Fpxz = self.Fp['x', 'z']
        self.FpXZ = FractionField(self.Fpxz)
        x, z = self.FpXZ.gens()
        
        a0_map = 36*x**7 + 52093178300019920831088149019314051232888168213855215466901398150856759165825142081300844754729266384812814186895830464818765409742963023209668068685551346824229758381291518549836095942317531113567498632709984830589066459540185798361769089894601095269793749238355573968396396248472966113*x**6*z + 101625270205915329150208244036819233114776677616667682143562350825890865178202618867968193383272169824102770255134215560272828592132856206249566232250400100147908327548050934015830847570278549715833061146739117054229084428040626154683019547224255401224433239151559275298375055242429392406*x**5*z**2 + 258434612088532982459344832459371716164045912907869304007571927894976984909716855974908565025897366164837992656504502017564296085382746721794012450276240627660364415626173564869675144871496806576080838085173300155092842544035299481590678456737686750811831784236306123852146610036252044571*x**4*z**3 + 169182265156315503932301417218794557353750714049107274750690822309494126528698399826952011883907685014506964714875663726886573181578889404419749469102754094812412448720439521490181430173244984277129958796043727359600252348021921086168670570254415080322602459633957659060856967392923559758*x**3*z**4 + 32645377557643449306671069693716494883272121329732367387366063547712140615654411978740340013400118010265778581978703766236105507766918402208177262324724438409523308667540034919444947946129938978990299003835351560181949670206743556207688093271021949880364576744650139720101999831902515510*x**2*z**5 + 252039864364971734406985414701800034827439835304182162328305103805921810303977505148164333856551626170723584501150114861683813479431717379988209083981417018171478469207255850300091172828650018553220539692191184094208385456992651976113387766319305806867898873394430586344415060994084601245*x*z**6 + 247308289805541814661115560309824463340238179875245956570532676049791316336961356169769559580359974376402551636054505296199468195763532031676306153443238851840893084361583172051694064585116432650517752211089147110223010401507471329295641141304632698693325812752135693833445626018245107073*z**7
        b0_map = 216*x**6 + 312559069800119524986528894115884307397329009283131292801408388905140554994950852487805068528375598308876885121374982788912592458457778139258008412113308080945378550287749111299016575653905186681404991796259908983534398757241114790170614539367606571618762495430133443810378377490837796678*x**5*z + 162125072814838451233029298951232469136864583214201686433151727172795925126055444333240669353645004663840271569490148096602470796759183106986042121724216713378942763898232097425152900006150475747544147258818009812119643054999184096996263783742729129692366008632420133444715140056900511954*x**4*z**2 + 279879028411507332634510735373739535464986379401760009218378678364781594708221238156200665469187007733747360706408015092833315993618887696104537307043391883070509553600558296005543637182151980139615014154928763814526567140060848431341158255450711540727440206262337594944078369992458871919*x**3*z**3 + 64540403583062258976072200252673540547841874808402209597757777804591842477837962770725903010027620629402183629109335413781216274427483652762762988613667671788104133803690571491851652001703935387850099335693672842751236608420222904049552546359559441219091347407474492466325601475438491930*x**2*z**4 + 315867002703168198541123467858808287884118081039117008022615737496189805220007592429815446224788134444717991016555616222726499006037302683153164771043727307869617917806911797663147322908119772997731099189944411457194313281792663631761292299610244065678367107630776450956517502032771362635*x*z**5 + 84750481658368488859932655998829222173793819073655313937935537243733989808543008188572076720586408324189517939578482567617779804017655060141227159121526999798881200178117305355737511723826500873401637902267653446298845370199015548236993320225640664985604263368070622103232971249064790076*z**6
        c0_map = x**6 + 11085104051803073590443870335240333701902504078788795317563944793589177799498182091056471337924385211085647418551935087587153456656940013095752642349724274403939241043717824589219553463582610407574003031072627603385125702860199044312701239504686841344319595387284411351523079044933046011*x**5*z + 115837262241353735792751913787449394622409735216007150107364382128551755152043536646937526282294396750080710375459675192591442649528548421282865929290900038930933257533631160572994837193178325832023752785040988830320854873587671463389472006134291756853912094146974988682195965930384367793*x**4*z**2 + 57793859054048863341336011619138897223580985365073783856937373584444582058701718819519450053740362523779398969751702848410742437637224478839123776373763402741757983355580071794603843751117306540020672265707907417658837814894153800124119804789808529187900649813138242850800857037748854663*x**3*z**3 + 322834518595137995953557433582027113937436436158189065789020628626594010055857372260487170958793738298829182301977100305842642509451594574025932569891237298350140731981271685974045499809827291031378133617000195235985608085332846788684174269059403202305022277116009024800039188768482916523*x**2*z**4 + 95588887505548544476440179723946625224311654176354146482422715551709889817658949434515859567482075925283206585584583334539224952219199356335327691322408710546308563062402866694020184360505847953180902438674401664714993291009556923789052031979438510994937259918467177181244534259216655119*x*z**5 + 307776330919120610757256276241369893466823808171064426370245272277935707585628373869777692080335589614229926092262526934277698875423224413765732766062387817621052919319343478396850793195473188418285308107093104410859606229340401953778985301702453975729834101966409165888307169415416524268*z**6

        a1_map = 36*x**7 + 273273281080736870653053436810408693914091736013522939915116148363116378228026843408875157492598038879811907050391905004747309875885634021327443962066886997849675859452134634891377606663629268743557880944897632051712918647901125609630990483069737008218226618552716118144210249821187934308*x**6*z + 104776769605450814481285983876507433627895013789236527531355271010428816902225577340539818839318902476605391214457602716073992266099560803651302682325313288009797927885185742808481045724728125390059882768244405322661515528838445125592113843080224711393897427629517451799121024123953717643*x**5*z**2 + 125946996076625134548737294177810713355419844980128575354801907503427709646708062259531588593062892113010937811389768553882633061916938361360381694992703931121757191283601723640290403419713247374839964973956612695011662001060736422874422029077956697503791860513771863409810048564921117533*x**4*z**3 + 335201731215499922379775993058967817804098557841407518385885080977611027778090517633731093050246692287512031045097307715385983592992221161111009411048358931726822868507126096644614910183455237960635826185337277278761040654359182414562534777040151475275440608373489810763296638203560437786*x**3*z**4 + 214395254717938006432706061630853711963076395930907029968065442858648295475516814431842475527124285557172603493988010682381025467612459391266622133656321141425079565479418035092217859163329229422851232433207588376631061686659435703786539043052098912590310797957179193281968126415952860380*x**2*z**5 + 69689421607828878268985248835962878658797569583646317862482604564181083755465771009949859307812000274535571724432768265543133574175160709116290029541316581462784301620162616635049626523770122424123380127110137427478001534720096221022186395152803160507462818070150149010041194801525559045*x*z**6 + 135334249385933808738419572863332632114214117675485219726999563479489490757357613884393606516110526565698429788557627134326776835307483115091656204398950303260963212396422059371563119051949525570176299928588498827175008260062538053458502074137397709992806709156401417430951718177521984239*z**7
        b1_map = 216*x**6 + 251757416224858310218755888665100315342142501590971975629094432505283702903723407666322451313393828421789470778456099275208823135686294335014955548783231860227723479942607142699994324995149838224960843725767356745175676513698834470204040769321921276806512972500763770119523701450660844716*x**5*z + 112927789715328148191212621489124984635153952305603865522441189618362086654189833539799705664845011553885644003124239325111680122858000314469904664797214315134641276553168273199137664639109117800512603497829903535605760216278998649863291007007189832652338926184103558835618195458442080946*x**4*z**2 + 137713240372089457353006776923646873720060290178277978968812762513049812690462732388667628497431410463639386792138610847707378730004053158445507356761544357478155359727806599451448136654240331733881247114124723504452044828101298308142731010462151906457940686099046061742312627629715150642*x**3*z**3 + 266190474099022952160150716618134138563061044476155818462494272822332176084964307998555895273918401640799456074633884998633101169472549844101512264355730744882800750720178794989955295686866408317601347824669953792544101841165974657155669519799598184845842994825454891033359595214350959228*x**2*z**4 + 25874584147990417296800484724586420421112314441613007993172653350309442341275945391856097370946438624260691903808119842626162263341110138412761085903122257097765026893863923558638548017934601084344717337693424377923722893900453448252796173856611788832358127588109996235523946804084093070*x*z**5 + 40298366705078848866560795372128456092173691846303411235633206436731620658930295015954188926271603255424385752040050302734409058957551103678216370109409913654942526945490680453454746640922885117928299499529586696642665198521110538225764692513949538719585575940208189666291594608265733415*z**6
        c1_map = x**6 + 123247780329428489215326323150068451065036541096778331430553431149537779933926105604712795511643923706973828378391052701793678284299115650004904684247809927179463080271187128745449765322875183095908811299326470520750511799139580088121574913065534425714632412083314303732817323395849672714*x**5*z + 69183139080552903293448583103917959253224745924385421231028297284311668712547995384510502280525094741694998595582573614445211556091743219417411095960020802563588875129718742429777133226172729465520803659544637690598873602400870306348694256427118181304572928257146759220580643720388811043*x**4*z**2 + 30308768211786096092394121217775345792589246116800200168564613425234387259683962058311406310065710741722320469981743527910465072184758561644946226819551143376389816471335439151684294135604344064487565820411801510049804072873469727932140486336957369263095393120847184795575214949973031345*x**3*z**3 + 286602558769953442942987695024234268189121740420356281522173354270880956339546368127417046321905419897598869188243563409090770844043855547400336918468726965943620587773301576843436750130493802298038343735883027743957774226773072520377154411225315029531672640614647327679465706160796403568*x**2*z**4 + 167826039331114358100718022215503321296270908102236805575844059632038881425494981092599378874516506987654774621583417336837286593888091241890148363528673865489124947268879883961643493407373161858930835131734488926992241803931324346496972420835219973218739232470087471562769516193176771073*x*z**5 + 220583958489697287871987179510025461091190005563308845435931282270134732827082180112496702350695842250282715606539046926652405673682427360913195271645661344427611645507504385675815488767090116591828476005959060449103540831480946454687723618435101301105074096104215336863419997875177342424*z**6

        a2_map = 36*x**7 + 160709482922905657596027866236216919620843588765293794244379512354106041235171911743993215419520728766804386107875485434977620911700989132497383732365241357469396485197954673983077409129896089174524087826855249070420271681732761871110567863582234376285937659109761410068402690713864727378*x**6*z + 281614629406100270982544946338658174671472686016059758949012399792623688673856686749029115450742557806587375817439143680711376102068814861533018305607206195083759716790204233738151859598561601620859529603884498011142104756408931087551557843340171793875778480565214775028217345893502971817*x**5*z**2 + 290757630225730524656732300080764843180520895191639959960452823016567383312641344464923372773122629794191708231564997482228308333985817367297749840191934701173538233474298441259165386386070143054851822482598386036263476021002375556734693418554072308392315356015377259842730306241875157968*x**4*z**3 + 23796325461642658986744309828699825470914289573639674401159876040463888841285313300596480845846705736369731601733031995501135599346902402572398219309176606209827650163008546960225455639260112775556357265246790505111453561238799869777290841784127732377460978991217265403110770581503487327*x**3*z**4 + 98643360086267892249701804652792926536764416834910477979835987315034455777809225507872736917844124397199943221539989510651771300627946463241866159542625117509943735114359521598075992939140338868063809709541529490623494266208877153084175291209189059094189747217493767080132301743113085339*x**2*z**5 + 113037619392772139952201998078298178759277673075561015416378871229168036172746948847043826101089137463234392084148632781127105785170527913695834357076019305364013372621840329149688149853656739091508376697241321784442703231776454383305970520077181783259081023057154204653697046620119160054*x*z**6 + 32458062314534341639627331840763941526197476752738952335164782980229448008028780668206177575301574409632900407245754666227061096781335768144266104426551910807808712411831364716174695785774990263506451739253761198976899758789042505998856860327423203489700605262844535669742392129644192313*z**7
        b2_map = 216*x**6 + 270315762407652488726384831318625593653857575346679933535475845287928964178812644070495045696027170172285330885305247233228207410392179898509448282382403081381213072802627710574328797286063647928951305989322276639970712403542611632872456116945155871464202585250801991037547245544954983702*x**5*z + 33555617998312942129758229926726435976492367207552820646726797266961533714425011130252312143819303241444882402628876771474013741408540897423260158311815148077518419660257431237851801190805363466514044114980184543208044032260486453014559900433696256126887334408642870215091399438167939273*x**4*z**2 + 152636937447283692213262378806592974945716837299384123743496206934130592031251838547697483563396624999169920189349523092336584693816115351968867610094483207765434344025283772299941134592935252237984222534799609987037090644238047644719722089673913287632604692230415616261672076340991358544*x**3*z**3 + 52398802628757849870623713369141808919522584487582478547031394658374324217428604625804715421273785046255918850997792600053604388168467864982271146323528093089881710542266486779919284126658473340969888279169791790401231131282297992500682791427397792348082190846629654059038786702748730572*x**2*z**4 + 141363272591780114871400264916356180054176560023601110665681756334008121516124191360513012042696664018644282792351740490979851209336970664029469069606658552531258929749256676327680923613654115262828380640189281653052712084991135314765210615705797606446541618963877773792427854151244182494*x*z**5 + 166534643694648153429940204510962059204553082065756350764998150466790482713330538235483576805984065668275910810217444232299748681080856809005236129776930055371195113582558925472827956340589988488568301640614558533657498091330836186836772629054779434575511256831416930681989873628498858321*z**6
        c2_map = x**6 + 332158577225810845112286891386491878578647523942547276029666677849392496005080332234246816982727032501444476223916272134383761886903745065904719545364416873218533770492907231680371766292117088060661356787433665315938496454395835415709909332247291415071114859417827538594643721368439783250*x**5*z + 289115265925859081386253235857577898163527849292375263986245415079949644535825015778852617628687565002655078958455959127798383624932321583775249941754769794162782437578140722219214912416308362121923536907950054002344556856452064585361615094717104561525484273814070555171455150464897247270*x**4*z**2 + 1262744787273165565319141833489531990331431370988183537701814092796261126081403169519437530532656377175437617770913016770894748804795518898573537435554796349590640986204261913007401115755696035257153357412562519819425462611541197902578778000967292641043894703220642748083792811086172542*x**3*z**3 + 119008143653832960269441297935985866794712474085306196218085706942673854648752987699394426398221552141881756438213888374525652865485784396160399352161995841751863470474533647223078531510897545604181139727420290024066692207115759584123876386409615436132049782978564056621961920972247749260*x**2*z**4 + 131426498723692698274277966377779329983833503414395186919025441374736004217768996242825566613041411228791469549174630391455222221150931964408383710316814124903745901139924553336546967931971907027501653787100889334091585911741382594027670385509324780480056597517328173000407234347704726306*x*z**5 + 101454985793787133726597844935240423993175497127231411771280069358745381684171895778216520459984704337578526541875112531399125187419604538976512826783530732931195720848286165291975612354580205048954319851481267517202534962030086361141437064198385423499475441007714452378753187698280589577*z**6
        
        a3_map = 36*x**7 + 72051954641070995332123494877613400437064939891721327430104194295131435030127293147064534382858646173012119351352364243444831021710793824341966892589109709921780456962694010102973114569124207605779625882297512866436894351980481797677507831875319179738922669313624931501482923367595286720*x**6*z + 197194319690574746285617381297133074746714601541847573495761748149381466177554910237326256478824401861425958904957814876515188176828724578710091763629321113788421538520646169726202037265890049884374298873322297375662132592885804592225962675210852386883299855430556035739953601350418429162*x**5*z**2 + 27315173387406020102089002449098244194914135313382944379738136291988809584257802345653571801860870387937276587789455765721055192080253900398795003720217224408067137033562629457026824686310657196410776250849091841399335378428001113608579165120038658985656971238962051375615785685233488406*x**4*z**3 + 330216963503818405854620227084529226083244412924072438962129073303667516560411056749885077924876092232546752299725022084096114998269467528150805649031452717694533930759605275713147889300654307836488870620674705126572758818232392842666144019213319203983432991429033130214454029555458082828*x**3*z**4 + 207786834891298999638339124537847893646175824646259609153069001784910094510928312072059427883931190394002957706020930552734544894404581668960655581299085950686416981697710361048028421711329894242589699003353866760457820937988046637903205210761997841181495815711637965883205973906238984700*x**2*z**5 + 24073143896910356653662174162412543013146181671089999278503860552854110211801426330913082029434474938602512578656894020232589493411896675676548328002417039734754501541334315664758571887787329837908103097028475394949390492476643327518728508639967208876764558444108790454672820302779069061*x*z**6 + 70705144135424704242676138800749014545770057333418704174884977260057716609974873353927928410274108767898982885117400864286052888152709103915484752985534442438940440621558845504091077649203393646712503481060104262073654354037869322148486260755993563852670676653456677190314074532765767706*z**7
        b3_map = 216*x**6 + 85341160281535243567849786216342440586787660727786548615224551352434968564654345685655082886603275823802223227140352772350227100357885497814374299630135727813099822583613893955770858668088802075581144807880468307345907268455910989169571458977789885307824331177866354322463090836455030037*x**5*z + 191135937920435070661956047506535903685662971609765232098670569882366297951990627765493912767886275444421341207871685130481832401138476674550506146751865176598399328552569174042241453789296937456342456280177591470556027638368147864862961332992914770944852414858077284681819267357025551387*x**4*z**2 + 210650597461186013923957173065883219376416099299704639155703620139203735057755232857417891274222697962729627662987114154199019380956288909774782586317860861972694568274228814888430955746225725953361723164465072257829625402172288462930775251903776377493084967274395227788961208140839561448*x**3*z**3 + 268401409243047154739092068057767168789291295822126977787327861114318197638781433303627485485966448381498792876226528165371375203015792500712719253371270503980301054639521225801032884182846476584249392345046710471079275918598469885400008106271220145993380598360984767233777837720453453863*x**2*z**4 + 123227115432877331759428010662927856302283783747567795779957744643169538569368946706316652325022703491159863118281486902413465064751505348029417426281953602956042472358615343563735559319109061508891465855951779111855801391215367548204367017662986548975027150076489603626547119941495527217*x*z**5 + 20678061402340664052229949328811736997770196329692764089525368836301752158426713139112505748726593079505246970048717612205965797942907047581278284479396512246782553776054796401229332611477367217298275297378750813576962483575681390321901012445939649084686354181894224041609780185562014118*z**6
        c3_map = x**6 + 117658299039437770456411602540824137357341352315617175528191987980982642622873340319662500425262273909562723164528843236201942760572036755644196987873427224737021541313147111390216307209138931397637193103143133876715066679808451093345089283865689486034651746826672881881630453216583210281*x**5*z + 114822247100888182406100774951729886747950522243967841417277607553078266851056778989664853253951593183788369456233753051078427006062886086694012540507538017981650357239446572588951115447152506693233835895974194471323651857481558257376596843449910719760922628984447183149388725318904527953*x**4*z**2 + 71917347354061818452056243794002304184623362035221540613870730757091980760524121182353329208426057785631100206189804744924738822273450773197160511884632279171479602775616744577379116551254885680677328665167347689248092194266700083796197327999595076445115298248964389114591403512472789641*x**3*z**3 + 329097996962139179980512172671425296273207867941533241439645889108619422143157766719870726975870984510797033114199007947565321435950001178068841212138618767561436153719102303052698384363960771154887705246838238963644711476871281137523718784262061528124373956388474100227981980484343358094*x**2*z**4 + 211436656227068744196273849557907208959254734671515852596204894270167673541933304200314322466985889085204358250494762252142634408722224664221555909906691891373383843597069580625095952623392177723046074982255122877940150673488213091823870521440129962714060016779017249176704629132728294646*x*z**5 + 329171517312602726918307103811030292493720718202231342504100071591035063032033344414726662609650391460930611211626724685840460040296977298858845997231495838158352274341786121705499719083807294408395458761659652830915044693679414080249251422646847627141140818744256067937920860862446745268*z**6
        
        a4_map = 78364164096*x**7 + 212754038764259230563310879050775869828969748318899353580799673684499802576250*x**6*z + 1395753593892944420376111654753141437022053235732109664942190513232673002665941*x**5*z**2 + 140744998288905480129594275235941748908322965268043296245400183979525245446664*x**4*z**3 + 1758013572475244011248874196368531615224954973624963267816585389001962108638816*x**3*z**4 + 660908732607138698406007542114583911273286252701443619819188841307157079466744*x**2*z**5 + 897564107490539783202813083049574233824324844214142170724095378095943384576538*x*z**6 + 1617181084958655874666188185346609704706511258980425164104290088667012395776992*z**7
        b4_map = 470184984576*x**6 + 1276524232585555383379865274304655218973818489913396121484798042106998815457500*x**5*z + 50280436130337387002440463352056069516521612309087191897279153028715920593753*x**4*z**2 + 1210058530201633688471755675064356569199513476201401908772777501051083932610365*x**3*z**3 + 1847907225308104973665379831053057191346170247191545755315744884851300589029704*x**2*z**4 + 388583100936910331125940990325457366115722078515087518166794332664142446404981*x*z**5 + 1830466435720113281286605528729977321229117884727385054516126448437671436338230*z**6
        c4_map = 2176782336*x**6 + 2048894856808388543112749115899850999900188306607993260994151491035244035694198*x**5*z + 455138368609339222282470856334730840858410267095612063306431594105002346869154*x**4*z**2 + 2075178815782899770587040729446765843301704872700780258421817931305870036905014*x**3*z**3 + 1904666006208540449822628498443936083239940345165714136315856089286329868640114*x**2*z**4 + 605505512525506039631276325696569569873109797877767782065638926230105103929764*x*z**5 + 1289019751980764068309417578261907228073621040354328226615903611050475735353556*z**6

        a5_map = 78364164096*x**7 + 918173365797421856982116374111134049022032106155629355182815080890769595449773*x**6*z + 1216277079351257670828668244955233675738397670599540190431275538131421480763008*x**5*z**2 + 2043904946128269296382132078103921334256194369830450024384869485743210001766814*x**4*z**3 + 1333833876684118301943629517887616254539072912896290763352661924324324312361922*x**3*z**4 + 38948850958505616778281492532610963011161774519715701882169566252382498882485*x**2*z**5 + 463030069114208815305817753390224280913301177052811648373836758752489907433350*x*z**6 + 2101530624268586931820635965791601616263698249292720226126070110932953751867961*z**7
        b5_map = 470184984576*x**6 + 1182718970882311828863541933293636051392321502253255069790499073369071367850704*x**5*z + 186789573026254427446743084768186189856317864702386360761090711676572003252500*x**4*z**2 + 820896646479199322923404760694906048872900648748280603642363935774597070495670*x**3*z**3 + 154096910148287948053445935637125062680364398994682197471398710422990017293988*x**2*z**4 + 2104541610928218050078371014561350556125368836340879927713855884263793264747164*x*z**5 + 1695764959416863056759158210535648933875315586621232656861941847471477310542990*z**6
        c5_map = 2176782336*x**6 + 1287348506021409017883007156764705572160852195008363902747220136295389020732030*x**5*z + 197588504813365157608665662405594019938899037760530107886258771657465121701605*x**4*z**2 + 2155689545164463403890172392049941236190092887982086736441679040121787027141963*x**3*z**3 + 407157455533946611278451029044236990428300331045196476517679841219965402700637*x**2*z**4 + 243201542432039435503693759490279184423983255117588136201779965723864731177317*x*z**5 + 192776563220037312948486207527521874632307899779784006738468618799063013717278*z**6

        a6_map = 78364164096*x**7 + 748723229029389554935339082109009080748457909497180776480983428543164589405846*x**6*z + 1391303938891551317232766995482978103355038723341132038605776604367092346654093*x**5*z**2 + 976456343641866673403460671255362609540693587256432563326693887872311625460274*x**4*z**3 + 1653961684829618698769715808731669093644298109201115468744470112666256714686811*x**3*z**4 + 1114653992257323466695402677531427401793623699266356277920090040191534323782640*x**2*z**5 + 2108607027386354826390837549863034441969749759739166544581539577047386461917759*x*z**6 + 1193614370663587251353382929683500784151470249566405174380508866261819192271619*z**7
        b6_map = 470184984576*x**6 + 166018150274118016582878181280886241750876322302563597579509159283441331587142*x**5*z + 1561609859329543545934114307480054966208701712177558959245381307544654633889526*x**4*z**2 + 1313807235922346076316039538736698115272167332182137072197512927704942796830192*x**3*z**3 + 1287612522338618938999960327168476648941438766313747654332719787622926078515289*x**2*z**4 + 505669846774700238617693208172928455055535949038459540967592424909229433170513*x*z**5 + 1205753046611077754895215355409024013073550600263723851381673594665083158626232*z**6
        c6_map = 2176782336*x**6 + 1102378173448593427005492941235208979594647170045051953562180726008974456473257*x**5*z + 764021305697784908827519174356463827659929162080053031776743068101707471927225*x**4*z**2 + 414873907231075961042811921023098852965714131361356813704037231832506997775188*x**3*z**3 + 622620192902363846175099846868260227695216996132954756473277628138129629904574*x**2*z**4 + 1480832065775769817787466844764020114600680598048823797448169609316071198490178*x*z**5 + 2051805517718069552242545345618923148734147689857553637814277725577494746473698*z**6

        self.a_maps = [a0_map, a1_map, a2_map, a3_map, a4_map, a5_map, a6_map]
        self.b_maps = [b0_map, b1_map, b2_map, b3_map, b4_map, b5_map, b6_map]
        self.c_maps = [c0_map, c1_map, c2_map, c3_map, c4_map, c5_map, c6_map]

        self.iso_x = 1591776315227191114928250266108102815068276791771315371416973227795476042238011
        self.iso_y = 362723301092365577775585079158435711361056655579945183829463270975417886340321
        
    def random_point(self):
        P = super().random_point()
        return Lollipop956451Point(P.X, P.Y,P.Z, self)
       
class Lollipop956451Point(PointWeierstrass):
      
        def __init__(self, x, y, z, curve):
            super().__init__(x, y, z, curve)

        def fast_scalar_mul(self, n):
            psiP = self.psi()
            beta = vector([n,0]) * self.curve.N_inv
            b = vector([int(beta[0]), int(beta[1])]) * self.curve.N
            k1 = n-b[0]
            k2 = -b[1]
            return self.multi_scalar_mul(k1, psiP, k2)

        def psi(self):
            x,y,z = self.X, self.Y, self.Z

            for i in range(7):
                new_x = self.curve.a_maps[i](x,z) 
                y = y* self.curve.b_maps[i](x,z) 
                z = z* self.curve.c_maps[i](x,z) 
                x = new_x

            return Lollipop956451Point(self.curve.iso_x*x,  self.curve.iso_y*y, z, self.curve)
