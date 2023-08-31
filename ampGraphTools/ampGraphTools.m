(* ::Package:: *)

(* Mathematica Package *)


BeginPackage["ampGraphTools`"]
$AmpGTVersion = .63;  
(* .63 Fixed some somePriveledgedLegs boundaries and Jacobi op clarity *)
(* Exported symbols added here with SymbolName::usage *) 


nGonGraph::usage = "nGonGraph[n] gives an n-gon 1-loop graph";
halfLadderGraph::usage = "halfLadderGraph[legs] gives a legs multiplicity half-ladder graph.";
multiLoopGraph::usage = "multiLoopGraph[m,l] gives an $m$-point $l$-loop graph.";


extractType::usage = "expr ~ extractType ~ head  returns a unique list of expressions of type head from expression."

getUniqDotsFromExpr::usage = "getUniqDotsFromExpr[expr] returns a list of all dot products in the expression."

DEFAULTLABELSET::usage="TBD";
DotPower::usage="TBD";
FermionArrowSign::usage="TBD";
LEGS::usage="TBD";
LOOPS::usage="TBD";
Lorentz::usage="TBD";
Lsq::usage="TBD";
SPECIALNUM::usage="TBD";
Say::usage="TBD";
aRule::usage="TBD";
addLeg::usage="TBD";
addLegUO::usage="TBD";
addMastersToWeb::usage="TBD";
allFunction::usage="TBD";
applyMomentum::usage="TBD";
atreeRule ::usage="TBD";
bcjExpress::usage="TBD";
bcjReorder::usage="TBD";
bipartiteGraph::usage="TBD";
blobForm::usage="TBD";
blowOutExt::usage="TBD";
buildADressingContainer::usage="TBD";
buildAllMomenta::usage="TBD";
buildAllMomentaD::usage="TBD";
buildGraphDenom::usage="TBD";
buildLoops::usage="TBD";
buildMomenta::usage="TBD";
buildRHSNumberFromMom::usage="TBD";
buildRules::usage="TBD";
canonNecklOrderRules ::usage="TBD";
canonVacuumGraph::usage="TBD";
clarifyTrees::usage="TBD";
collapseTrees::usage="TBD";
collectionAdd::usage="TBD";
colorForm::usage="TBD";
colorToStuff::usage="TBD";
complexSpinorProd::usage="TBD";
concatenateNecklaces::usage="TBD";
consistentGraphToTrees::usage="TBD";
corruptCutList::usage="TBD";
corruptGraph::usage="TBD";
corruptGraphNoMemo::usage="TBD";
crappyGraph::usage="TBD";
curGuys::usage="TBD";
cutDisplayRule::usage="TBD";
cutGraph3D::usage="TBD";
cutGraphData::usage="TBD";
cutSig::usage="TBD";
cutSumFormat::usage="TBD";
daLabelList::usage="TBD";
definitelyGood::usage="TBD";
doAColorOrderedCut::usage="TBD";
doADoubleCopyDressedCut::usage="TBD";
doContactCollapse::usage="TBD";
doCubicGravDressedCut::usage="TBD";
doGravDressedCutWithContacts::usage="TBD";
doGravDressedCutWithContactsList::usage="TBD";
doMomentaPlot::usage="TBD";
doOrderedPlot::usage="TBD";
doOrderedPlotRand::usage="TBD";
doublePol::usage="TBD";
dressGraph::usage="TBD";
dressTreeCO::usage="TBD";
dressTreeUO::usage="TBD";
dressTreeUOWithContacts::usage="TBD";
dressTreeUOWithContactsCNTR::usage="TBD";
dressedGraphPlotForm::usage="TBD";
dressedMomentaPlot::usage="TBD";
evaluateKltCut::usage="TBD";
evaluateKltCutNoLsq::usage="TBD";
evaluateLpRule::usage="TBD";
evaluateLsqRule::usage="TBD";
evaluateSRule::usage="TBD";
exSetRightOfALeftOfB::usage="TBD";
extLegNum::usage="TBD";
extLegs::usage="TBD";
externalLegsRightOfALeftOfB::usage="TBD";
fancyFormatOff ::usage="TBD";
fancyFormatOn ::usage="TBD";
fancyFormatQ ::usage="TBD";
fancySameTest::usage="TBD";
feynRule::usage="TBD";
findAllColl::usage="TBD";
findAllCorruptEdgeIso::usage="TBD";
findAllMhV::usage="TBD";
findGoodSet::usage="TBD";
flattenMomenta::usage="TBD";
fromCompressedDressing::usage="TBD";
fromCompressedGraph::usage="TBD";
functionKeys ::usage="TBD";
gammamu ::usage="TBD";
generateMassiveDecay::usage="TBD";
generateMasslessDecay::usage="TBD";
generateNullMomenta::usage="TBD";
getAConnectingEdge::usage="TBD";
getAllGravNumerators::usage="TBD";
getAllNumerators::usage="TBD";
getAvgGravNumerator::usage="TBD";
getCDGraphContribution::usage="TBD";
getCOGraphContribution::usage="TBD";
getColorlyDotRules::usage="TBD";
getColorlyRules::usage="TBD";
getConnectingEdges::usage="TBD";
getCutPrecursor::usage="TBD";
getDotRules::usage="TBD";
getDressing::usage="TBD";
getDressingLabelForm::usage="TBD";
getEveryIsomorphism::usage="TBD";
getEveryIsomorphismBasisProject::usage="TBD";
getExtLegs::usage="TBD";
getExtLegsFromTrees::usage="TBD";
getGoodUniqLegs::usage="TBD";
getGraphDenom::usage="TBD";
getGravGraphContribution::usage="TBD";
getGravGraphContributionAll::usage="TBD";
getGravNumerator::usage="TBD";
getHardCutId::usage="TBD";
getIndepLegs::usage="TBD";
getIndepRules::usage="TBD";
getIntLegs::usage="TBD";
getIntLegsFromTrees::usage="TBD";
getKRules::usage="TBD";
getKRules2::usage="TBD";
getKRulesNum::usage="TBD";
getMHVRules::usage="TBD";
getMHVRules2::usage="TBD";
getMyCycles::usage="TBD";
getMyCyclesOld::usage="TBD";
getMyMomHash::usage="TBD";
getMyUniqLegs::usage="TBD";
getNColorTree::usage="TBD";
getNUOTree::usage="TBD";
getNUOTreeWithContacts::usage="TBD";
getNonConnectingEdges::usage="TBD";
getNumerator::usage="TBD";
getOffShellDots::usage="TBD";
getOffShellDotsBetter::usage="TBD";
getRandomMassVector::usage="TBD";
getRandomNullVector::usage="TBD";
getSig::usage="TBD";
getTheCut::usage="TBD";
getTheCutCompare::usage="TBD";
getTheCutOld::usage="TBD";
getTheProcessedCut::usage="TBD";
getUniqDotsFromExpr::usage="TBD";
getVertices::usage="TBD";
goodMHV::usage="TBD";
graphDenomProd::usage="TBD";
graphExclusion::usage="TBD";
graphExclusionF::usage="TBD";
graphHashCode::usage="TBD";
graphPlotForm::usage="TBD";
graphSet::usage="TBD";
graphToWeb::usage="TBD";
gravPhase::usage="TBD";
importCompressedGraphsAndDressings::usage="TBD";
inSetRightOfALeftOfB::usage="TBD";
indepLegLabels::usage="TBD";
internalOneLoopBubbleQ::usage="TBD";
internalOneLoopTadpoleQ::usage="TBD";
internalOneLoopTriangleQ::usage="TBD";
invSoftLimitCO::usage="TBD";
invSoftLimitUO::usage="TBD";
isIsomorphic::usage="TBD";
isomorphicEdgeRule::usage="TBD";
isomorphicEdgeRulesAll::usage="TBD";
isomorphicEdgeRulesAllBad::usage="TBD";
isomorphicRulesAll::usage="TBD";
isomorphicRulesAllWithSig::usage="TBD";
isomorphicVertexRule::usage="TBD";
isomorphicVertexRulesAll::usage="TBD";
jacobiGraphOnLeg::usage="TBD";
kkExpress::usage="TBD";
kltCutWithPhase::usage="TBD";
mapCutWithInsToLs::usage="TBD";
mapToAbsGeneric::usage="TBD";
mapToAbsGenericTrees::usage="TBD";
matchGraphOrAdd::usage="TBD";
maxPowerCount::usage="TBD";
mergeEdge::usage="TBD";
mergeLegsTreeLevel::usage="TBD";
mergeLegsTreeLevelList::usage="TBD";
mergeLegsTreeLevelListWhoo::usage="TBD";
mergeLegsTreeLevelWhoo::usage="TBD";
mhvTreeRule ::usage="TBD";
mhvTreeRules ::usage="TBD";
minStripPriveledge::usage="TBD";
monsterPairing::usage="TBD";
monsterPairingOld::usage="TBD";
newColorSimplify::usage="TBD";
newEqns::usage="TBD";
newGraph::usage="TBD";
newGraphSig::usage="TBD";
newGraphSigRule::usage="TBD";
newGraphUO::usage="TBD";
newKLT::usage="TBD";
newRules::usage="TBD";
nextBatch::usage="TBD";
nmhvLevel::usage="TBD";
noCrappyGuysInJacEqn::usage="TBD";
nonZero::usage="TBD";
nonZeroUnion::usage="TBD";
numNEWGRAPHS::usage="TBD";
opTillClosure::usage="TBD";
opTillClosureOnCollection::usage="TBD";
outSideRules::usage="TBD";
planarGraphs::usage="TBD";
planarQ::usage="TBD";
plotBlob::usage="TBD";
plotCutGraph::usage="TBD";
plotSmallBlob::usage="TBD";
polGaugeAnsatz::usage="TBD";
polGaugeAnsatzMaxPowers::usage="TBD";
polMinus::usage="TBD";
polPlus::usage="TBD";
polRules::usage="TBD";
priveledgeExtLegs::usage="TBD";
priveledgeLegs::usage="TBD";
priveledgeSomeLegs::usage="TBD";
processCut::usage="TBD";
processCutList::usage="TBD";
refineEdges::usage="TBD";
reformatEqn::usage="TBD";
reformatEqnRule::usage="TBD";
refreshGRAPHRULES ::usage="TBD";
refreshHLP ::usage="TBD";
resetMathematicaGraph ::usage="TBD";
restEqns::usage="TBD";
rightOf::usage="TBD";
rotateList::usage="TBD";
scoreRulesForLegs::usage="TBD";
scramble::usage="TBD";
selectCutList::usage="TBD";
setLeftOf::usage="TBD";
setRightOf::usage="TBD";
sewGraphs::usage="TBD";
showCollection::usage="TBD";
singlePol::usage="TBD";
someAddlLabels::usage="TBD";
someLabels::usage="TBD";
spa::usage="TBD";
spaEval::usage="TBD";
spb::usage="TBD";
spbEval::usage="TBD";
specialRules::usage="TBD";
stagePlot::usage="TBD";
stringList::usage="TBD";
stripPriveledge::usage="TBD";
stripVacuum::usage="TBD";
stripVacuumRules::usage="TBD";
stupidRep::usage="TBD";
stupidSameTest::usage="TBD";
symCount::usage="TBD";
tHat::usage="TBD";
tadpoleQ::usage="TBD";
theGuys::usage="TBD";
theRulez::usage="TBD";
toGraph::usage="TBD";
toVacuum::usage="TBD";
treesToLoops::usage="TBD";
treesToLoopsExtFlag::usage="TBD";
treesToLoopsExtFlagNoTri::usage="TBD";
treesToLoopsExtFlagNoTriUO::usage="TBD";
treesToLoopsUO::usage="TBD";
treesToLoopsUOWithContacts::usage="TBD";
treesToLoopsUOWithContactsSlow::usage="TBD";
twoExternalVertices::usage="TBD";
twoFersEqual::usage="TBD";
twoFersMap::usage="TBD";
uHat::usage="TBD";
dot::usage="TBD";
unionSortExprs::usage="TBD";
uniqDots::usage="TBD";
uniqLegs::usage="TBD";
uniqPol::usage="TBD";
zeGraph::usage="TBD";
zeReals::usage="TBD";
zeroRules::usage="TBD";




StylePrint["Welcome to ampGraphTools, version "<>ToString[$AmpGTVersion]<>", 
a work in progress, but fairly simple implementation of
ideas in http://arxiv.org/abs/arXiv:1506.00974 and refs 
therein. - dr.jjmc@gmail.com."];

(* Implementation of the package *)


(* ::Subsubsection::Closed:: *)
(*Some raw graphs *)


simpleBoxGraph = {Atree[{k[ 1],l[ 1],-l[ 4]}],
Atree[{k[ 2],l[ 2],-l[ 1]}],
Atree[{k[ 3],l[ 3],-l[ 2]}],
Atree[{k[ 4],l[ 4],-l[ 3]}]}//toGraph;

nGonGraph[nn_] :=
    (Join[consistentGraphToTrees[simpleBoxGraph]/.-l[ 4]:>-l[ nn],
    Table[Atree[{k[ i],l[ i],-l[ i-1]}],{i,5,nn}]]/.l[ a_]:>l[ a])/.l[ a_]:>l[ a+nn]//toGraph;

lRuleForGraphBuilding[mm_,ll_] :=
    Module[ {nn = mm+2 (ll-1)},
        Thread[Rule[
        k[ #]&/@Flatten[Table[ {i+2, nn- Floor[ll-1/2]+i}, {i,1,ll-1}]],
        Flatten[Table[{l[ i+(2*(ll-1)+mm)*2],-l[ i+(2*(ll-1)+mm)*2]},{i,1,ll-1}]]]]
    ];

halfLadderGraph[legs__] :=If[ Length[legs]>3,
    toGraph[Flatten[
    	{Atree[{legs[[1]],legs[[2]],in[1]}],
    	Table[
    		Atree[{-in[i],legs[[i+2]],in[i+1]}],
         	{i,1,Length[legs]-4}],
   	    Atree[{-in[Length[legs]-3],legs[[-2]],legs[[-1]]}
   	    ]}]],
    toGraph[{ Atree[legs] }]
]
    
multiLoopGraph[mm_,ll_] :=
    Module[ {tmpGraph = If[ ll>0,
                            nGonGraph[ 2*(ll-1)+mm]/.lRuleForGraphBuilding[mm,ll],
                            halfLadderGraph[k[ #]&/@Range[mm]]
                        ], allKs},
        allKs = tmpGraph/. k[ a_]:>Sow[k[ a]]//Reap//Last//Flatten//Union;
        tmpGraph/.MapIndexed[#->k[ #2[[1]]]&,allKs]
    ];


(* ::Subsubsection:: *)
(*Expression tools*)


(*
extractType[data_,head_]:=data /. head[a___]:>Sow[head[a]]//Reap//Last//Flatten;
*)
extractType[data_, head_] := Union[Cases[data, head[___], {0,Infinity}]]


getUniqDotsFromExpr[expr_]:=(expr~extractType~dot)


functionKeys = DownValues[#][[All, 1, 1, 1]] &;  
(* gives argument to a function so you can use it as a hash *)


Say[a_]:=StylePrint[ Grid[{
Flatten[{DateString[]," -- ", Flatten[{a}] }]}]]


polRules[LEGS_]:={dot[k[ a_],\[Epsilon][ k[a_]]]:>0,
    dot[k[ 1],\[Epsilon][k[ LEGS]]]->Total[-dot[k[ #],\[Epsilon][ k[LEGS]]]&/@Range[2,LEGS-1]]};


scramble[list_] := (* Put list in bucket, remove list from bucket *)
    Module[ {hat = list,bucket = {},curN = 0},
        While[Length[hat]>0,
        curN = Random[Integer,{1,Length[hat]}];
        bucket = Join[bucket,hat[[{curN}]]];
        hat = Join[hat[[Range[curN-1]]],
        hat[[Range[curN+1,Length[hat]]]]]];
        bucket
    ]


(* ::Subsubsection::Closed:: *)
(*Access graph meta data *)


getMyCyclesOld[gr_] :=
    Module[ {g = System`Graph[gPF = graphPlotForm[gr]/.{a_->b_,c_}:>UndirectedEdge[a,b]],fc},
        fc = FindFundamentalCycles[g];
        Table[{Intersection@@(First[Complement[fc[[i]],
        Sequence@@fc[[Range[Length[fc]]/.i:>{}//Flatten]]]]/.neckl[a__]:>a /.-a_:>a
        )//First,Length[fc[[i]]]},{i,1,Length[fc]}]
    ]
 


getMyCycles[graph_]:=Module[{uniqLegs=getMyUniqLegs[graph],
extLegs=getExtLegs[graph],cycles=FindFundamentalCycles[mathematicaGraph[graph]],lLegs,cnt},
lLegs=Complement[uniqLegs,extLegs];
MapIndexed[(cnt[#2[[1]]]=Length[#])&,cycles];
Table[  {leg, Sort[Map[#[[1]]&,Position[cycles,leg]]//Tally ,!OrderedQ[Last/@{#1,#2}]&]//First//First//cnt,Sort[Map[#[[1]]&,Position[cycles,leg]]//Tally ,!OrderedQ[Last/@{#1,#2}]&]//First//First},{leg,
lLegs}]]


getExtLegs[graph_]:=(
	getExtLegsFromTrees[
		consistentGraphToTrees[ graph ]
	])


(* ::Subsubsection:: *)
(*Gauge Feynman Rules *)


feynRule[p_, q_, r_] := 
 g structureF[{c[p /. -a_ :> a], c[q /. -a_ :> a], 
    c[r /. -a_ :> a]}] (  
   component[p - r , q /. -a_ :> a] g[p /. -a_ :> a, r /. -a_ :> a] +
     component[r - q, p /. -a_ :> a] g[q /. -a_ :> a, r /. -a_ :> a] + 
    component[q - p, r /. -a_ :> a] g[p /. -a_ :> a, q /. -a_ :> a])
feynRule[aa_, bb_, cc_, dd_] := 
 Block[{\[Mu] = aa, \[Nu] = bb, \[Sigma] = cc, \[Rho] = 
    dd}, -I g^2  ( 
    structureF[c /@ {aa, bb, e}] structureF[c /@ {cc, dd, e}]*
      ( g[\[Mu], \[Sigma]] g[\[Nu], \[Rho]] - 
        g[\[Mu], \[Rho]] g[\[Nu], \[Sigma]]) + 
     structureF[c /@ {aa, cc, e}] structureF[c /@ {bb, dd, e}]*
      ( g[\[Mu], \[Nu]] g[\[Rho], \[Sigma]] - 
        g[\[Mu], \[Rho]] g[\[Nu], \[Sigma]]) + 
     structureF[c /@ {aa, dd, e}] structureF[c /@ {cc, bb, e}]*
      ( g[\[Mu], \[Sigma]] g[\[Nu], \[Rho]] - 
        g[\[Mu], \[Nu]] g[\[Rho], \[Sigma]]))]

theRulez[stuff_] := 
 stuff /. component[a_, b_] :> component[a][b] /. 
    g[a_, b_] :> component[\[Eta]][{a, b}] //. 
        {   component[a_, b_] :> component[a][ b] ,
    component[-(b_)][c_] :> -component[b][c] ,      
    component[\[Eta]][{a_, b_}]*component[\[Eta]][{b_, c_}] :> 
     component[\[Eta]][{a, c}], 
     component[\[Eta]][{b_, a_}]*component[\[Eta]][{b_, c_}] :> 
     component[\[Eta]][{a, c}] ,
     component[\[Eta]][{a_, b_}]*component[\[Eta]][{c_, b_}] :> 
     component[\[Eta]][{a, c}] ,   
    component[a_][b_]*component[\[Eta]][{c_, b_}] :> component[a][c],
     component[a_][b_]*component[c_][b_] :> dot[a, c] } /. \[Epsilon][
    k[ a_]] :> \[Epsilon][ a]


(* ::Subsubsection:: *)
(*Formatting*)


fancyFormatOn := (FFON = True;

Format[k[a_]] := Subscript[\[ScriptK], a];

Format[l[a_]] := Subscript[\[ScriptL], a];

Format[flat[a_]] := 
 Superscript[RawBoxes[RowBox[{"(", MakeBoxes[a], ")"}]], 
  Style[\[Flat], {Blue, Bold}]];

Format[spb[a_, b_]] := 
  RawBoxes[RowBox[{"[", MakeBoxes[a], ",", MakeBoxes[b], "]"}]];
  
Format[spa[a_, b_]] := 
  RawBoxes[RowBox[{"\[LeftAngleBracket]", MakeBoxes[a], ",", MakeBoxes[b], 
        "\[RightAngleBracket]"}]];

Format[AtreeF[a__]] := 
 DisplayForm[
  RowBox[{SuperscriptBox[SubscriptBox["\[ScriptCapitalA]", Length[a]], 
     "tree"], "(", Sequence @@ Riffle[a, ","], ")"}]];

Format[Atree[a__]] := 
 DisplayForm[
  RowBox[{SuperscriptBox[SubscriptBox["A", Length[a]], "tree"], "(", 
    Sequence @@ Riffle[a, ","], ")"}]];

Format[AtreeF[a__, h__]] := 
 DisplayForm[
  RowBox[{SuperscriptBox[SubscriptBox["\[ScriptCapitalA]", Length[a]], 
     "tree"], "(", Sequence @@ Riffle[Thread[SuperscriptBox[a, h]], ","], 
    ")"}]];

Format[Atree[a__, h__]] := 
 DisplayForm[
  RowBox[{SuperscriptBox[SubscriptBox["A", Length[a]], "tree"], "(", 
    Sequence @@ Riffle[Thread[SuperscriptBox[a, h]], ","], ")"}]];

  Format[numerator[a_, b_]] := 
   n @@ (Join[
       b /. -Total[Subscript[k, #] & /@ Range[LEGS - 1]] :> 
         Subscript[k, LEGS], 
       Complement[Subscript[k, #] & /@ Range[LEGS], 
        b /. -Total[Subscript[k, #] & /@ Range[LEGS - 1]] :> 
          Subscript[k, LEGS]]] /. Subscript[k, aa_] :> aa);
  Format[perm[a__]] := 
   DisplayForm[RowBox[{"(", a /. k[b_] :> b, ")"} // Flatten]];
  Format[dot[a_, b_]] := 
   DisplayForm[RowBox[{"(", a, "\[CenterDot]", b, ")"}]];
  Format[\[Epsilon][k[a_]]] := 
   DisplayForm[RowBox[{SubscriptBox["\[CurlyEpsilon]", a]}]];
  Format[aa[i_, j_]] := Subscript[\[ScriptA], j];
  Format[\[Tau][\[Epsilon][k[a_]], \[Epsilon][k[b_]]]] := 
   DisplayForm[
    RowBox[{SubscriptBox["\[CurlyEpsilon]", 
        RowBox[{a, "\[CenterDot]", b}]]} // Flatten]];
  SetAttributes[\[Tau], Orderless];
  Format[\[Tau][a_, b_]] := 
   DisplayForm[
    RowBox[{"(", If[Head[a]  ~SameQ~  Plus, {"(", a, ")"}, a], 
       "\[CenterDot]", If[Head[b]  ~SameQ~  Plus, {"(", b, ")"}, b], ")"} //
       Flatten]];
  Format[uLsq[a__]] := 
   DisplayForm[SuperscriptBox[RowBox[{"(", Plus @@ a, ")"}], "2"]];
  Format[k[a_]] := Subscript[k, a];
  Format[l[a_]] := Subscript[l, a];
  )
fancyFormatOn;


fancyFormatOff := (FFON = False;
  Format[spa[a_,b_]]=.;
  Format[spb[a_,b_]]=.;
  Format[Atree[a__]]=.;
  Format[AtreeF[a__]]=.;
  Format[numerator[a_, b_]] =.;
  Format[perm[a__]] =.;
  Format[dot[a_, b_]] =.;
  Format[\[Epsilon][k[a_]]] =.;
  Format[aa[i_, j_]] =.;
  Format[\[Tau][a_, b_]] = .;
  Format[\[Tau][\[Epsilon][k[a_]], \[Epsilon][k[b_]]]] =.;
  Format[uLsq[a__]] =.;
  Format[k[a_]] =.;
  Format[l[a_]] =.;)
fancyFormatOff;
fancyFormatQ := FFON


(* ::Subsubsection:: *)
(*Properties of dot products etc*)


Clear[Lorentz,Lsq];

Lorentz[a__, b__] :=
    First[a]*First[b] - Rest[a] . Rest[b]

Lsq[a_] :=
    If[ Head[a]  ~SameQ~  List,
        Lorentz[a, a],
        HoldForm[Lorentz[a, a]]
    ]


ulp[x___]:=dot[x];

SetAttributes[dot, Orderless]

dot[-(a_), b_] :=
    -dot[a, b]

dot[a_, (b_) + (c_)] :=
    dot[a, b] + dot[a, c]

dot[a_?NumberQ b_, c_] :=
    a dot[b,c]

uLsqCleaningRule :=
    uLsq[a__]:>uLsq[Flatten[a/.Plus:>List]]


(* ::Subsubsection:: *)
(*Scaffolding graph operations, dressing, etc.*)


Clear[leftOf];
leftOf[el_,neckl[list__]] :=
    Module[ {loc},
        If[ (loc = Position[list,el,1]) ~SameQ~ {},
            Throw[
                Error[{"Element Not In Necklace", 
                    el, 
                    neckl[list]}]
            ],
            If[ Length[loc = Flatten[loc]]>1,
                Throw[Error[
                    {"Ambiguous element in Necklace", 
                        el,
                        neckl[list]}]
                        ],
                loc = loc[[1]];
                If[ loc ~SameQ~ 1,
                    list[[-1]],
                    list[[loc-1]]
                ]
            ]
        ]
    ];

leftOf[el_,list__] :=
    Module[ {elCount,elCanon = el/.-a_:>a,nextVal,connedge,StylePrint},
        StylePrint[el];
        StylePrint[list];
        If[ NumberQ[el]&&el<0,
            el = -el
        ];
        If[ Head[list[[1]]] ~UnsameQ~ neckl,
            Throw[Error["Not a Necklace in leftOf",el,list]]
        ];
        elCount = Count[list,elCanon,5];
        If[ elCount === 0,
            elCanon = -elCanon,
            elCount = Count[list,elCanon,5];
        ];
        If[ elCount === 0,
            Throw[Error["Not Present in leftOf[...,list]",list,el,elCanon]];
        ];
        nextVal = If[ elCount === 1,
                      leftOf[el,Select[list,MemberQ[#[[1]],el]&,1][[1]]],
                      leftOf[-el,Select[list,MemberQ[#[[1]],-el]&,1][[1]]]
                  ];
        connedge = getConnectingEdges[list];
        If[ MemberQ[connedge,nextVal] || MemberQ[connedge,-nextVal],
            leftOf[nextVal,list],
            nextVal
        ]
    ];


Clear[rightOf];
rightOf[el_,neckl[list__]] :=
    Module[ {loc},
        If[ (loc = Position[list,el,1]) ~SameQ~ {},
            Throw[Error[{"Element Not In Necklace", el, neckl[list]}]],
            If[ Length[loc = Flatten[loc]]>1,
                Throw[Error[{"Ambiguous element in Necklace", el, neckl[list]}]],
                loc = loc[[1]];
                If[ loc ~SameQ~ Length[list],
                    list[[1]],
                    list[[loc+1]]
                ]
            ]
        ]
    ];

rightOf[el_,list__] :=
    Module[ {elCount,elCanon = el/.-a_:>a,nextVal,connedge},
        If[ NumberQ[el]&&el<0,
            el = -el
        ];
        If[ Head[list[[1]]] ~UnsameQ~ neckl,
            Throw["Not a Necklace in rightOf"]
        ];
        elCount = Count[list,elCanon,5];
        If[ elCount === 0,
            elCanon = -elCanon,
            elCount = Count[list,elCanon,5];
        ];
        If[ elCount === 0,
            Throw[Error["Not Present in rightOf[...,list]",list,el,elCanon]];
        ];
        nextVal = If[ elCount === 1,
                      rightOf[el,Select[list,MemberQ[#[[1]],el]&,1][[1]]],
                      rightOf[-el,Select[list,MemberQ[#[[1]],-el]&,1][[1]]]
                  ];
        connedge = getConnectingEdges[list];
        If[ MemberQ[connedge,nextVal] || MemberQ[connedge,-nextVal],
            rightOf[nextVal,list],
            nextVal
        ]
    ];

externalLegsRightOfALeftOfB[a_,b_,list__] :=
    Module[ {tmp,StylePrint},
        StylePrint[a];
        StylePrint[b];
        StylePrint[list];
        StylePrint[leftOf[b,list]];
        NestWhileList[(tmp = rightOf[#,list];
                       StylePrint[{#,list,tmp}];
                       tmp)&,a,(# ~UnsameQ~ leftOf[b,list])&]
    ]

setRightOf[el_,neckl[list__]] :=
    Rest[NestWhileList[rightOf[#,neckl[list]]&,el,(# ~UnsameQ~ leftOf[el,neckl[list]])&]]

setLeftOf[el_,neckl[list__]] :=
    Rest[NestWhileList[leftOf[#,neckl[list]]&,el,(# ~UnsameQ~ rightOf[el,neckl[list]])&]]

inSetRightOfALeftOfB[elA_,elB_,neckl[list__]] :=
    If[ elA ~SameQ~ elB,
        NestWhileList[rightOf[#,neckl[list]]&,
        elA,(# ~UnsameQ~ leftOf[elA,neckl[list]])&],
        NestWhileList[rightOf[#,neckl[list]]&,
        elA,(# ~UnsameQ~ elB)&]
    ]

exSetRightOfALeftOfB[elA_,elB_,neckl[list__]] :=
    If[ rightOf[elA,neckl[list]] ~SameQ~ elB,
        {},
        NestWhileList[rightOf[#,neckl[list]]&,
        rightOf[elA,neckl[list]],(# ~UnsameQ~ leftOf[elB,neckl[list]])&]
    ]

getConnectingEdges[list__] :=
    Module[ {l = Flatten[list/.neckl[a___]:>List[a]]},
        Select[l,Count[l,#,\[Infinity]] ~SameQ~ 2&]
    ]

getAConnectingEdge[list__] :=
    Module[ {l = Flatten[list/.neckl[a___]:>List[a],2]},
        Select[l,Count[l,#,\[Infinity]] ~SameQ~ 2&,1]
    ]

(* getNonConnectingEdges[list__] :=
    Module[ {l = Flatten[list/.neckl[a___]:>List[a],2]},
        Select[l,Count[l,#,\[Infinity]]+Count[l,-#,\[Infinity]] \[Equal] 1&]
    ] *)

getNonConnectingEdges[list__] :=
    Module[ {a,l = Flatten[list/.neckl[a___]:>List[a],2]},
        Select[l,(a = #/.-b_:>b;
                  (Count[l/.-b_:>b,a,\[Infinity]]+Count[l/.-b_:>b,-a,\[Infinity]]) ~SameQ~ 1)&]
    ]

concatenateNecklaces[neckl[listA__],neckl[listB__],ell_] :=
    Module[ {locA,locB,ret,el},
        el = If[ (Head[ell] ~SameQ~ Times && ell[[1]] ~SameQ~ -1)||(NumberQ[ell]&&ell<0),
                 -ell,
                 ell
             ];
        If[ listA ~SameQ~ listB,
            ret = neckl[listA/.{-el:>merged[-el,el],el:>merged[el,-el]}],
            locA = Flatten[Position[listA,el]];
            locB = Flatten[Position[listB,el]];
            If[ Length[locA]>1,
                locA = Flatten[Position[listA,-el]]
            ];
            If[ Length[locB]>1,
                locB = Flatten[Position[listB,-el]]
            ];
            If[ Length[locA] ~UnsameQ~ 1,
                Throw[Error[{"Bad element in A",el,neckl[listA]}]]
            ];
            If[ Length[locB] ~UnsameQ~ 1,
                Throw[Error[{"Bad element in B",el,neckl[listB]}]]
            ];
            locA = locA[[1]];
            locB = locB[[1]];
            ret = listA;
            ret[[locA]] = Append[Prepend[
            setRightOf[listB[[locB]],neckl[listB] ],
            merged[listA[[locA]],listB[[locB]]]
            ],merged[listB[[locB]],listA[[locA]]]];
            neckl[Flatten[ret,1]]
        ]
    ]

concatenateNecklaces[neckl[listA__],neckl[listB__],ell_] :=
    Module[ {locA,locB,ret,el,StylePrint},
        el = If[ (Head[ell] ~SameQ~ Times && ell[[1]] ~SameQ~ -1)||(NumberQ[ell]&&ell<0),
                 -ell,
                 ell
             ];
        If[ listA ~SameQ~ listB,
            neckl[listA/.{el:>merged[el,-el],-el:>merged[-el,el]}],
            locA = Flatten[Position[listA,el]];
            locB = Flatten[Position[listB,el]];
            If[ Length[locA]>1,
                locA = Flatten[Position[listA,-el]]
            ];
            If[ Length[locB]>1,
                locB = Flatten[Position[listB,-el]]
            ];
            If[ Length[locA] ~UnsameQ~ 1,
                Throw[Error[{"Bad element in A",el,neckl[listA]}]]
            ];
            If[ Length[locB] ~UnsameQ~ 1,
                Throw[Error[{"Bad element in B",el,neckl[listB]}]]
            ];
            locA = locA[[1]];
            locB = locB[[1]];
            StylePrint[{locA,locB}];
            ret = listA;
            ret[[locA]] = Append[Prepend[
            setRightOf[listB[[locB]],neckl[listB] ],
            merged[listA[[locA]],listB[[locB]]]
            ],merged[listB[[locB]],listA[[locA]]]];
            neckl[Flatten[ret,1]]
        ]
    ]

mergeEdge[necklaces__,nextEdge_] :=
    Module[ {newForm = necklaces,concatForm,locs,StylePrint},
        locs = Map[First,Position[newForm,nextEdge]];
        If[ Length[locs] ~UnsameQ~ 2,
            Throw[Error["tried to concatenate a bad edge",
            newForm,nextEdge]]
        ];
        StylePrint[{"Necklaces", necklaces}];
        StylePrint[{newForm,nextEdge,locs}];
        concatForm = concatenateNecklaces[
        newForm[[locs[[1]]]],
        newForm[[locs[[2]]]],
        nextEdge];
        newForm[[locs[[1]]]] = newForm[[locs[[2]]]] = {};
        newForm = Append[Flatten[newForm],concatForm]
    ]

concatenateNecklaces[necklaces__] :=
    Module[ {newForm = necklaces,nextEdge},
        While[(nextEdge = getAConnectingEdge[newForm]) ~UnsameQ~ {},
        newForm = mergeEdge[newForm,nextEdge[[1]]];
];
        If[ Length[newForm]>1,
            Throw[Error["Distinct Graphs present in concatenate",necklaces,newForm]],
            newForm[[1]]
        ]
    ]

Clear[newGraph];
newGraph[vertexFormGraph[neck_],old_,new_,{leftBound_,rightBound_}] :=
    Module[ {newA, newB,numConnecting,replEdge,necklaces = neck,StylePrint},
        numConnecting = Length[getConnectingEdges[necklaces]];
        replEdge = in[numConnecting+1];
        necklaces = Replace[necklaces,old->replEdge,3];
        newA = Append[necklaces, neckl[{-replEdge,new,old}]];
        newB = Append[necklaces, neckl[{-replEdge,old,new}]];
        ewA = newA;
        ewB = newB;
        StylePrint[{leftBound,rightBound,newA}];
        StylePrint[{leftBound,rightBound,newB}];
        vertexFormGraph/@Select[{newA,newB},MemberQ[ 
        externalLegsRightOfALeftOfB[leftBound,rightBound,#],new]&
        ]
    ]

Clear[newGraphUO];
newGraphUO[vertexFormGraph[neck_],old_,new_,{leftBound_,rightBound_}] :=
    Module[ {newA, newB,numConnecting,replEdge,necklaces = neck},
        numConnecting = Length[getConnectingEdges[necklaces]];
        replEdge = in[numConnecting+1];
        necklaces = Replace[necklaces,old->replEdge,3];
        newA = Append[necklaces, neckl[{-replEdge,new,old}]];
        newB = Append[necklaces, neckl[{-replEdge,old,new}]];
        ewA = newA;
        ewB = newB;
        vertexFormGraph/@Select[{newA,newB},MemberQ[ 
        externalLegsRightOfALeftOfB[leftBound,rightBound,#],new]&,1
        ]
    ]

Clear[newGraphUO];
newGraphUO[vertexFormGraph[neck_],old_,new_,{leftBound_,rightBound_}] :=
    Module[ {newA,numConnecting,replEdge,necklaces = neck},
        numConnecting = Length[getConnectingEdges[necklaces]];
        replEdge = in[numConnecting+1];
        necklaces = Replace[necklaces,old->replEdge,3];
        newA = Append[necklaces, neckl[{-replEdge,new,old}]];
        vertexFormGraph/@{newA}
    ]

Protect[vertexFormGraph]

addLeg[vertexFormGraph[necklaces__],new_,{i_,j_}] :=
    Module[ {StylePrint,edges = Union[inSetRightOfALeftOfB[i,j,concatenateNecklaces[necklaces]]/.merged[a_,b_]:>merged@@Sort[{a,b}]]}(* only consider each edge once *),
        StylePrint[edges];
        Map[newGraph[vertexFormGraph[necklaces],#/.merged[a_,b_]:>a/.-a_:>a,new,{i,j}]&,
        edges]
    ]

Clear[refineEdges];
refineEdges[edges__,necklaces__] :=
    Module[ {gpf = graphPlotForm[vertexFormGraph[necklaces]],revHash,revSame},
        revSame[a_,b_] :=
            revHash[a/.merged[c_,d_]:>d] ~SameQ~ revHash[b/.merged[c_,d_]:>d];
        Map[(revHash[#[[2]]] = Sort[List@@#[[1]]])&,gpf];
        Union[edges,SameTest->revSame]
    ]

addLegUO[vertexFormGraph[necklaces__],new_,{i_,j_}] :=
    Module[ {
    edges = Union[
    inSetRightOfALeftOfB[i,j,concatenateNecklaces[necklaces]]/.merged[a_,b_]:>merged@@Sort[{a,b}]
    ]
    }(* only consider each edge once *),
        edges = refineEdges[edges,necklaces];
        Map[newGraphUO[vertexFormGraph[necklaces],#/.merged[a_,b_]:>a/.-a_:>a,new,{i,j}]&,
        edges]
    ]

Clear[invSoftLimitCO,invSoftLimitUO]

invSoftLimitCO[graphSet_] :=
    invSoftLimitCO[graphSet] = Map[
    Module[ {StylePrint,g = #,numEx},
        numEx = Length[concatenateNecklaces[g[[1]]][[1]]/.merged[a___]:>{}//Flatten];
        StylePrint[numEx];
        addLeg[g,numEx+1,{numEx,1}]
    ]&,graphSet]//Flatten

invSoftLimitCO[graphSet_,exLabel_] :=
    invSoftLimitCO[graphSet,exLabel] = Map[
    Module[ {StylePrint,g = #,numEx},
        numEx = Length[concatenateNecklaces[g[[1]]][[1]]/.merged[a___]:>{}//Flatten];
        StylePrint[numEx];
        addLeg[g,exLabel[numEx+1],{exLabel[numEx],exLabel[1]}]
    ]&,graphSet]//Flatten

invSoftLimitUO[graphSet_,exLabel_] :=
    invSoftLimitUO[graphSet,exLabel] = Map[
    Module[ {g = #,numEx},
        numEx = Length[concatenateNecklaces[g[[1]]][[1]]/.merged[a___]:>{}//Flatten];
        addLegUO[g,exLabel[numEx+1],{exLabel[1],exLabel[1]}]
    ]&,graphSet]//Flatten

Clear[getNColorTree];
getNColorTree[n_] :=
    Nest[invSoftLimitCO[#]&,{ vertexFormGraph[{neckl[{1,2,3}]}]},n-3];

getNColorTree[n_,ex_] :=
    Nest[invSoftLimitCO[#,ex]&,{vertexFormGraph[{neckl[{ex[1],ex[2],ex[3]}]}]},n-3];

getNUOTree[n_,ex_] :=
    Nest[invSoftLimitUO[#,ex]&,{vertexFormGraph[{neckl[{ex[1],ex[2],ex[3]}]}]},n-3];

dressTreeCO[Atree[list__]] :=
    Module[ {graphs = getNColorTree[Length[list],ex]},
        graphs/.Thread[
        Rule[(ex/@Range[Length[list]]),list]
        ]
    ];

dressTreeUO[Atree[list__]] :=
    Module[ {graphs = getNUOTree[Length[list],ex]},
        graphs/.Thread[
        Rule[(ex/@Range[Length[list]]),list]
        ]
    ];

sewGraphs[vertexFormGraph[a__],vertexFormGraph[b__]] :=
    Module[ {na = Length[getConnectingEdges[a]],
    bConEdges = getConnectingEdges[b]},
        vertexFormGraph[
        Join[a, b/.
        Map[
        #->in[#[[1]]+na]&,
        bConEdges]
        ]
        ]
    ]

dressGraph[g_] :=
    g /.graphAlgRules[g]

refreshGRAPHRULES :=
    (Clear[graphAlgRules];
     graphAlgRules[g_] :=
         graphAlgRules[g,{}];
     graphAlgRules[g_,hold_] :=
         graphAlgRules[g,hold] = Module[ {rules = (Rule@@#&/@(List@@ Reduce[g/. 
                         vertexFormGraph[necklaces_]:> necklaces/.neckl[a__]:>Plus@@a==0, 
                         Complement[Union[Flatten[Reap[g/.in[a_]:>Sow[in[a]]][[2]]]],hold], 
                         Backsubstitution->True]))},
                                     Select[rules,Head[#[[1]]] ~SameQ~ in&]/. Flatten[Map[{#,-#[[1]]->-#[[2]]}&,(Reverse/@
                                     Select[rules,Head[#[[1]]] ~UnsameQ~ in&])]]
                                 ];)
refreshGRAPHRULES;

graphDenomProd[graph_] :=
    1/getGraphDenom[stripPriveledge[graph]]

graphDenomProd[graph_,hold_] :=
    (Times@@Map[uLsq[Flatten[{# /.Plus:>List}]]&,Select[#[[2]]&/@graphPlotForm[graph],Head[#] ~SameQ~ in&]/.graphAlgRules[graph,hold]])^(-1)

getGraphDenom[graph_] :=
    (Times@@Map[uLsq[Flatten[{# /.Plus:>List}]]&,
    Select[#[[2]]&/@graphPlotForm[graph],Head[#] ~SameQ~ in&]/.graphAlgRules[graph]])



treesToLoops[treeListUS_] :=
    Module[ {treeGraphs,negRules,bbb, newForms, next,treeList},
        negRules = MapIndexed[-#1->bbb[#2[[1]]]&, Reap[treeListUS/.-a_:>Sow[a]]//Last//Flatten//Union];
        treeList = treeListUS/.negRules;
        treeGraphs = dressTreeCO /@ treeList;
        newForms = First[treeGraphs];
        treeGraphs = Rest[treeGraphs];
        While[Length[treeGraphs] > 0, next = First[treeGraphs];
                                      treeGraphs = Rest[treeGraphs];
                                      newForms = Flatten[Outer[sewGraphs[#1, #2] & , newForms, next]]; ];
        newForms /. Reverse/@negRules
    ]


treesToLoopsUO[treeListUS_] :=
    Module[ {treeGraphs,negRules,bbb, newForms, next,treeList},
        negRules = MapIndexed[-#1->bbb[#2[[1]]]&, Reap[treeListUS/.-a_:>Sow[a]]//Last//Flatten//Union];
        treeList = treeListUS/.negRules;
        treeGraphs = dressTreeUO/@treeList;
        newForms = First[treeGraphs];
        treeGraphs = Rest[treeGraphs];
        While[Length[treeGraphs]>0,next = First[treeGraphs];
                                   treeGraphs = Rest[treeGraphs];
                                   newForms = Flatten[Outer[sewGraphs[#1,#2]&,newForms,next]];];
        newForms /. Reverse/@negRules
    ]



dressTreeUOWithContactsCNTR[Atree[list__]] :=
    Module[ {StylePrint,boo,graphs = getNUOTreeWithContacts[Length[list]],gooz},
        gooz =  graphs /. Thread[k /@ Range[Length[list]] -> list]/.
          in[a_]:>in[a+CNTR];
        CNTR = If[ (boo = (gooz/.in[a_]:>Sow[a]//Reap//Last//Flatten//Max)+1)>0,
                   boo,
                   CNTR
               ];
        StylePrint[{CNTR}];
        gooz
    ]


treesToLoopsUOWithContacts[treeListUS_] :=
    Module[ {treeGraphs,negRules,bbb, newForms, next,treeList},
        negRules = MapIndexed[-#1->bbb[#2[[1]]]&, Reap[treeListUS/.-a_:>
           Sow[a]]//Last//Flatten//Union];
        treeList = treeListUS/.negRules;
        Block[ {CNTR = 0},
            treeGraphs = dressTreeUOWithContactsCNTR /@treeList;
        ];
        newForms = First[treeGraphs];
        treeGraphs = Rest[treeGraphs];
        While[Length[treeGraphs] > 0, next = First[treeGraphs];
                                      treeGraphs = Rest[treeGraphs];
                                      newForms = Flatten[Outer[sewGraphs[#1, #2] & , 
                                                   newForms, next]]; ];
        newForms/. Reverse/@negRules
    ]
                


treesToLoopsExtFlag[treeList_] :=
    Module[ {treeGraphs = dressTreeCO/@treeList,newForms,next},
        newForms = First[treeGraphs];
        treeGraphs = Rest[treeGraphs];
        While[Length[treeGraphs]>0,
        next = First[treeGraphs];
        treeGraphs = Rest[treeGraphs];
        newForms = Flatten[Outer[sewGraphs[#1,#2]&,Select[newForms,
        And@@Map[Count[#,k,\[Infinity]]<2&,#[[1]]]&],Select[next,
        And@@Map[Count[#,k,\[Infinity]]<2&,#[[1]]]&]
        ]
        ];
        ];
        newForms
    ]

treesToLoopsExtFlagNoTri[treeList_] :=
    Module[ {treeGraphs = dressTreeCO/@treeList,newForms,next},
        newForms = First[treeGraphs];
        treeGraphs = Rest[treeGraphs];
        While[Length[treeGraphs]>0,
        next = First[treeGraphs];
        treeGraphs = Rest[treeGraphs];
        newForms = Flatten[Outer[sewGraphs[#1,#2]&,Select[newForms,
        And@@Map[Count[#,k,\[Infinity]]<2&,#[[1]]]&],Select[next,
        And@@Map[Count[#,k,\[Infinity]]<2&,#[[1]]]&]
        ]
        ];
        ];
        Select[newForms,(!internalOneLoopBubbleQ[#])&&(!internalOneLoopTriangleQ[#])&]
    ]

treesToLoopsExtFlagNoTriUO[treeList_] :=
    Module[ {treeGraphs = dressTreeUO/@treeList,newForms,next},
        newForms = First[treeGraphs];
        treeGraphs = Rest[treeGraphs];
        While[Length[treeGraphs]>0,
        next = First[treeGraphs];
        treeGraphs = Rest[treeGraphs];
        newForms = Flatten[Outer[sewGraphs[#1,#2]&,Select[newForms,
        And@@Map[Count[#,k,\[Infinity]]<2&,#[[1]]]&],Select[next,
        And@@Map[Count[#,k,\[Infinity]]<2&,#[[1]]]&]
        ]
        ];
        ];
        Select[newForms,(!internalOneLoopBubbleQ[#])&&(!internalOneLoopTriangleQ[#])&]
    ]


buildLoops[treeList__] :=
    treesToLoopsExtFlagNoTri[treeList]

graphPlotForm[vertexFormGraph[necklaces_]] :=
    Sort[
    Module[ {tmp},
        Flatten[{Table[
        Map[ 
        If[ Head[#] ~UnsameQ~ Times&&(
        tmp = Position[necklaces,-#]) ~UnsameQ~ {},
            Ray[neckl1,necklaces[[tmp[[1,1]]]],#],
            {}
        ]&,neckl1[[1]] ],{neckl1,necklaces}],
        Table[
        Ray[
        necklaces[[Position[necklaces,edge][[1,1]]]],
        neckl[{-edge}],edge],{edge,
        getNonConnectingEdges[necklaces]}]
        }]/.Ray[a_,b_,c_]:>{Rule[a,b],c}
    ],OrderedQ[{#2[[2]],#1[[2]]}]&]

internalOneLoopTriangleQ[vertexFormGraph[necklaces__]] :=
    Module[ {intEdges = getConnectingEdges[necklaces]},
        Or@@Map[internalOneLoopBubbleQ[vertexFormGraph[ mergeEdge[necklaces,#]/.merged[a___]:>{}/.neckl[a__]:>neckl[Flatten[a]]]]&,intEdges]
    ]

tadpoleQ[g_] :=
    Module[ {
    gpf = graphPlotForm[g]},
        Length[gpf] ~UnsameQ~ Length[Select[gpf,#[[1]][[1]] ~UnsameQ~ #[[1]][[2]]&]]
    ]

internalOneLoopBubbleQ[g_] :=
    Module[ {
    gpf = graphPlotForm[g]},
        Length[gpf] ~UnsameQ~ Length[Union[Map[#[[1]]/.Rule[a_,b_]:>Sort[{a,b}]&,gpf]]]
    ]

dressedGraphPlotForm[g_] :=
    Sort[graphPlotForm[g]/.graphAlgRules[g],OrderedQ[{#2[[2]],#1[[2]]}]&]

dressedMomentaPlot[g_,opts___] :=
    GraphPlot[dressedGraphPlotForm[g],DirectedEdges->True,MultiedgeStyle->.24,opts]

doMomentaPlot[g_,opts___] :=
    GraphPlot[graphPlotForm[g],DirectedEdges->True,MultiedgeStyle->.24,Method-> {"SpringElectricalEmbedding"},
    opts]

blobForm[treeList__] :=
    vertexFormGraph[treeList/.Atree[a__]:>neckl[a]]

plotCutGraph[graph_,opts___] :=
    GraphPlot[graphPlotForm[graph],MultiedgeStyle->.24,Method-> {"SpringElectricalEmbedding","RepulsiveForcePower"->-2,"InferentialDistance"->.5},ImageSize->Large,opts]

Clear[stupidRep];
stupidRep[graph_] :=
    stupidRep[graph] = (graph /.k[a_]:>Map[k[a,#]&,Range[a]]/.neckl[a__]:>neckl[Flatten[a]])

stupidSameTest[a_,b_] :=
    stupidSameTest[a,b] = isIsomorphic[stupidRep[a],stupidRep[b]]

plotBlob[graph_] :=
    (plotCutGraph[graph,BaseStyle->{PointSize[.06],Thickness[Large]},ImageSize->Medium]/.
    Tooltip[Point[a_],neckl[{-k[#]}]]:>{}/@Range[20])

plotSmallBlob[graph_] :=
    (plotCutGraph[graph,BaseStyle->{PointSize[.06],Thickness[Large]},ImageSize->Small]/.
    Tooltip[Point[a_],neckl[{-k[#]}]]:>{}/@Range[20])

resetMathematicaGraph :=
    (Clear[mathematicaGraph,mmaGraphRule];
     mathematicaGraph[graph_] :=
         Module[ {gPF = graphPlotForm[graph],necklaces,rules},
             necklaces = Union[Flatten[Reap[gPF/.neckl[a__]:>Sow[neckl[a]]][[2]]]];
             rules = Thread[Rule[necklaces,Range[Length[necklaces]]]];
             mmaGraphRule[graph] = Reverse/@rules;
             Graph[Map[#[[1]]&,gPF]/.Rule :> UndirectedEdge]
         ]);
resetMathematicaGraph  

graphHashCode[graph_] :=
    graphHashCode[graph] = CanonicalGraph[mathematicaGraph[graph]]

isIsomorphic[graphA_, graphB_] :=
    IsomorphicGraphQ[mathematicaGraph[graphA], mathematicaGraph[graphB]]

isIsomorphic[graphA_, graphB_] :=
    graphHashCode[graphA]  ~SameQ~  graphHashCode[graphB]

isomorphicVertexRule[graphA_,graphB_] :=
    Module[ {rule = FindGraphIsomorphism[mathematicaGraph[graphA],
    mathematicaGraph[graphB],1],vC = VertexCount[mathematicaGraph[graphA]]},
        Normal[First[rule]]
    ]
    
isomorphicVertexRulesAll[graphA_,graphB_] :=
    Module[ {rule = FindGraphIsomorphism[mathematicaGraph[graphA],
    mathematicaGraph[graphB],All],vC = VertexCount[mathematicaGraph[graphA]]}, 
        Normal[#]&/@rule
    ]

isomorphicEdgeRule[graphA_,graphB_] :=
    Module[ {
    gPFA = graphPlotForm[graphA],
    gPFB = graphPlotForm[graphB]},
        ((Reverse[(Rule@@#)]&/@Join[ gPFA, Map[ {Reverse[#[[1]]],-#[[2]]}&,gPFA]])/.(isomorphicVertexRule[graphA,graphB])/.
        (((Rule@@#)&/@Join[ gPFB, Map[ {Reverse[#[[1]]],-#[[2]]}&,gPFB]])))
    ]

isomorphicEdgeRulesAllBad[graphA_,graphB_] :=
    Module[ {
    gPFA = graphPlotForm[graphA],
    gPFB = graphPlotForm[graphB],rules = isomorphicVertexRulesAll[graphA,graphB]},
        Table[((Reverse[(Rule@@#)]&/@Join[ gPFA, Map[ {Reverse[#[[1]]],-#[[2]]}&,gPFA]])/.rule/.
        (((Rule@@#)&/@Join[ gPFB, Map[ {Reverse[#[[1]]],-#[[2]]}&,gPFB]]))),{rule,rules}]
    ]


isomorphicEdgeRulesAll[graphA_,graphB_] :=
    Module[ {StylePrint,
    gPFA = graphPlotForm[graphA],
    gPFB = graphPlotForm[graphB],(* Assumes first rule is identity fully spelled out. *) 
    rules = isomorphicVertexRulesAll[graphA,graphB]},
        StylePrint[rules];
(* Sometimes return value only tells us the difference btetween the input and output, confusing the hell out of us later. *)
        rules = Map[   Join[#,First[rules]]&,rules];
        zibbles = Table[{(Reverse[(Rule@@#)]&/@Join[ gPFA, Map[ {Reverse[#[[1]]],-#[[2]]}&,gPFA]]),rule,
        (((Rule@@#)&/@Join[ gPFB, Map[ {Reverse[#[[1]]],-#[[2]]}&,gPFB]]))},{rule,rules}];
        Map[#[[1]]/.#[[2]]/.#[[3]]&,zibbles]
    ]


isomorphicRulesAll[graphA_,graphB_] :=
    isomorphicEdgeRulesAll[graphA,graphB]

rotateList[l_, ip_] :=
    Join[l[[Range[ip, Length[l]]]], l[[Range[ip - 1]]]]



isomorphicEdgeRulesAll[graphA_,graphB_] :=
    Module[ {StylePrint,
    gPFA = graphPlotForm[graphA],
    gPFB = graphPlotForm[graphB],(* Assumes first rule is identity fully spelled out. *) 
    rules = isomorphicVertexRulesAll[graphA,graphB]},
        StylePrint[rules];
(* Sometimes return value only tells us the difference btetween the input and output, confusing the hell out of us later. *)
        rules = Map[   Join[#,First[rules]]&,rules];
        zibbles = Table[{(Reverse[(Rule@@#)]&/@Join[ gPFA, Map[ {Reverse[#[[1]]],-#[[2]]}&,gPFA]]),rule,
        (((Rule@@#)&/@Join[ gPFB, Map[ {Reverse[#[[1]]],-#[[2]]}&,gPFB]]))},{rule,rules}];
        Map[#[[1]]/.#[[2]]/.#[[3]]&,zibbles]
    ]


getGravGraphContributionAll[graph_, graphList__, dressingStorage_] :=
    getAllGravNumerators[graph, graphList, dressingStorage]/getGraphDenom[stripPriveledge[graph]];
   


getAllGravNumerators[graph_, graphList__, dressingStorage_] :=
    Module[ {myIsoRule,
    myIsoRules,
    myGraph = Select[graphList, isIsomorphic[#1, graph] & ]},
        If[ Length[myGraph]  ~SameQ~  0,
            Return[NoDressing[graph]],
            myGraph = First[myGraph];
            myIsoRules = isomorphicEdgeRulesAll[myGraph, graph];
      (*   myIsoRulesBad= isomorphicEdgeRulesAllBad[myGraph, graph];
       clearingHouse[graph]={"good"\[Rule]myIsoRules,
"bad"\[Rule]myIsoRulesBad,
 "from"\[Rule]myGraph,
"to"\[Rule]graph}; *)
            1/Max[Length[myIsoRules],1] Sum[(
            dressingStorage[myGraph] /. myIsoRule //. graphAlgRules[graph]
            ), {myIsoRule, myIsoRules}]
        ]
    ];


doGravDressedCutWithContactsList[cut_, gravDressing_, gravGraphs_] :=
    Module[ {allCutGrphs, cutGraphs, cutValue, isGood },
        StylePrint[{"Doing cut: ", cut}];
        ((isGood[graphHashCode[#1]] = True) & ) /@ 
        gravGraphs;
        isGood[_] = False;
        StylePrint[
        stringList[{"CutGraphs: ",
        Timing[allCutGrphs = corruptGraph/@treesToLoopsUOWithContacts[cut]; ]}]];
        aCGZ = allCutGrphs;
        StylePrint[stringList[{"N=4Only: ", 
                Timing[cutGraphs = Select[allCutGrphs, isGood[graphHashCode[#1]] & ]; ]}]];
        StylePrint[stringList[{"Expression: ", 
           Timing[cutValue = (getGravGraphContributionAll[mee = #1, Select[gravGraphs, 
                   isIsomorphic[mee, #1] & , 1], gravDressing] & ) /@ cutGraphs; ]}]];
        cutValue
    ]


doCubicGravDressedCut[cut_, gravDressing_, gravGraphs_] :=
    Module[ {allCutGrphs, cutGraphs, cutValue, isGood,StylePrint},
        StylePrint[{"Doing cut: ", cut}];
        ((isGood[graphHashCode[#1]] = True) & ) /@ 
        gravGraphs;
        isGood[_] = False;
        StylePrint[
        stringList[{"CutGraphs: ",
        Timing[allCutGrphs = corruptGraph/@treesToLoopsUO[cut]; ]}]];
        aCGZ = allCutGrphs;
        StylePrint[stringList[{"N=4Only: ", 
                Timing[cutGraphs = Select[allCutGrphs, isGood[graphHashCode[#1]] & ]; ]}]];
        StylePrint[stringList[{"Expression: ", 
           Timing[cutValue = (getGravGraphContributionAll[mee = #1, Select[gravGraphs, 
                   isIsomorphic[mee, #1] & , 1], gravDressing] & ) /@ cutGraphs; ]}]];
        cutValue
    ]


canonNecklOrderRules = neckl[a__]:>neckl[ 
    Select[
       Map[rotateList[a,#]&, Range[Length[a]] ],
       (#[[1]] ~SameQ~ Sort[a][[1]] )&][[1]]];

getSig[graph_] :=
    (graph/.neckl:>Signature)/.vertexFormGraph[a__]:>Times@@a

newGraphSig[graphA_,graphB_] :=
    getSig[graphA]/getSig[graphB]

newGraphSigRule[graphA_,graphB_,isoRule_] :=
    newGraphSig[graphA/.isoRule,graphB]

getDressing[graph_,graphList__,dressingStorage_] :=
    Module[ {
    myIsoRule,
    myGraph = Select[graphList,isIsomorphic[#,graph]&]},
        If[ Length[myGraph] ~SameQ~ 0,
            Return[NoDressing[graph]],
            myGraph = First[myGraph];
            myIsoRule = isomorphicEdgeRule[myGraph,graph];
            dressed[{graph,myGraph,dressingStorage[myGraph] /. myIsoRule/.graphAlgRules[graph],newGraphSig[graph,myGraph/.isomorphicVertexRule[myGraph,graph]]}
            ]
        ]
    ]

getDressing[graph_, graphList__, dressingStorage_] :=
    Module[ {myIsoRule, myGraph = Select[graphList, 
            isIsomorphic[#1, graph] & ]},
        If[ Length[myGraph]  ~SameQ~  0,
            Return[NoDressing[graph]],
            myGraph = First[myGraph];
            myIsoRule = isomorphicEdgeRule[myGraph, graph];
            dressed[{graph, myGraph, dressingStorage[myGraph] /. 
                    myIsoRule /. graphAlgRules[graph], newGraphSig[graph, 
                  myGraph /. myIsoRule]}]
        ]
    ]

getNumerator[graph_,graphList__,dressingStorage_] :=
    (getDressing[graph,graphList,dressingStorage]/.
    dressed[{a_,b_,c_,d_}]:>c d /. graphAlgRules[graph])

getAllNumerators[graph_, graphList__, dressingStorage_] :=
    Module[ {myIsoRule, myIsoRules, myGraph = Select[graphList, 
        isIsomorphic[#1, graph] & ]},
        If[ Length[myGraph]  ~SameQ~  0,
            Return[NoDressing[graph]],
            myGraph = First[myGraph];
            myIsoRules = isomorphicEdgeRulesAll[myGraph, graph];
            Table[(dressingStorage[myGraph] /. 
                 myIsoRule /. graphAlgRules[graph])*(newGraphSig[graph, 
                myGraph /. myIsoRule]), {myIsoRule, myIsoRules}]
        ]
    ]


getCOGraphContribution[graph_, graphList__, dressingStorage_] :=
    getNumerator[graph,graphList,dressingStorage]/getGraphDenom[stripPriveledge[graph]]
   
getCDGraphContribution[graph_, graphList__, dressingStorage_] :=
    colorFactor[graph] getNumerator[graph,graphList,dressingStorage]/getGraphDenom[stripPriveledge[graph]]


getDressingLabelForm[graph_,graphList__,dressingStorage_] :=
    Module[ {
    myIsoRule,myLabel = Select[graphList,isIsomorphic[#["vertexForm"],graph]&,1],
    myGraph},
        If[ Length[myLabel] ~SameQ~ 0,
            Return[NoDressing[graph]],
            myLabel = First[myLabel];
            myGraph = myLabel["vertexForm"];
            myIsoRule = isomorphicEdgeRule[myGraph,graph];
            dressed[{graph,myLabel,dressingStorage[myGraph] /. myIsoRule/.graphAlgRules[graph],newGraphSig[graph,myGraph/.isomorphicVertexRule[myGraph,graph]]}
            ]
        ]
    ]

Clear[getMyMomHash];
getMyMomHash[myTree_] :=
    getMyMomHash[myTree] = Module[ {myRule,legs,myMomHash},
                               myRule = Rule@@#&/@List@@Reduce[
                                 Join[myTree/.Atree[a__]:>Plus@@a == 0,
                                 {leg[k, 1]+leg[k, 2]+leg[k, 3]+leg[k, 4] == 0}]];
                               legs = Flatten[myTree/.Atree[a__]:>a];
                               Clear[myMomHash];
                               (myMomHash[#/.myRule] = #;
                                myMomHash[-(#/.myRule)] = #;)&/@(moosh = Sort[
                               Union[Flatten[Outer[#1+#2+#3+#4&,legs,legs,legs,legs]]],
                               OrderedQ[{#1,#2}/.-a_:>a/.leg[p, a_]:>.12*a/.leg[l, a_]:>.12*a/.leg[k, a_]:>100*(5-a)]&
                               ]);
                               (myMomHash[#/.myRule] = #;
                                myMomHash[-(#/.myRule)] = #;)&/@(moosh = Sort[
                               Union[Flatten[Outer[#1+#2+#3&,legs,legs,legs]]],
                               OrderedQ[{#1,#2}/.-a_:>a/.leg[p, a_]:>.12*a/.leg[l, a_]:>.12*a/.leg[k, a_]:>100*(5-a)]&
                               ]);
                               (myMomHash[#/.myRule] = #;
                                myMomHash[-(#/.myRule)] = #;)&/@Sort[Union[Flatten[Outer[#1+#2&,legs,legs]]],OrderedQ[{#1,#2}/.-a_:>a/.leg[p, a_]:>.12*a/.leg[l, a_]:>.12*a/.leg[k, a_]:>100*(5-a)]&];
                               (myMomHash[#/.myRule] = #;
                                myMomHash[-(#/.myRule)] = #;)&/@Sort[legs,OrderedQ[{#1,#2}/.-a_:>a/.leg[p, a_]:>.12*a/.leg[l, a_]:>.12*a/.leg[k, a_]:>100*(5-a)]&];
                               myMomHash
                           ]

clarifyTrees[expr_,myTree_] :=
    Module[ {myRule,legs,myMomHash},
        myRule = Rule@@#&/@List@@Reduce[Join[myTree/.Atree[a__]:>Plus@@a == 0,{leg[k, 1]+leg[k, 2]+leg[k, 3]+leg[k, 4]==0}]];
        legs = Flatten[myTree/.Atree[a__]:>a];
        Clear[myMomHash];
        myMomHash = getMyMomHash[myTree];
        expr/.uLsq[b__]:>If[ (Head[tmp = myMomHash[Plus@@b/.myRule]] ~UnsameQ~ myMomHash)&&Length[tmp]<=Length[b],
                             uLsq[Flatten[{tmp/.Plus:>List}]],
                             uLsq[b]
                         ]
    ]

fancySameTest[treeA_,treeB_] :=
    isIsomorphic[toGraph[treeA],toGraph[treeB]]


getIndepRules[treeList_] :=
    ToRules[getIndepReq[treeList,{}]];

getIndepRules[treeList_,req_] :=
    ToRules[getIndepReq[treeList,req]];



getGoodUniqLegs[treeList_] :=
    (Module[ {Print,Subscript},
         Subscript[a_,b_] :=
             a[b];
         (Print[Complement[Union[
         Flatten[treeList/.Atree[a__]:>a/.-a_:>a]],{Subscript[k, 1],Subscript[k, 2],Subscript[k, 3],Subscript[k, 4]}]];
          Rule@@#&/@List@@Reduce[  treeList/.Atree[a__]:>Plus@@a==0,
          Complement[Union[Flatten[treeList/.Atree[a__]:>a/.-a_:>a]],{Subscript[k, 1],Subscript[k, 2],Subscript[k, 3]}]])
     ]);
    
    getIndepReq[treeList_,req_] :=
        ( Module[ {me,legs,loops,doo,eqns2,Subscript,
        myTrees = treeList,
        eqns,
        free = Flatten[{List@@Expand[req/.uLsq[{a_,b_}]:>dot[a,b]/.-1:>1]}],ks},
              Subscript[a_,b_] :=
                  a[b];
              legs = Union[Flatten[myTrees/.Atree[a__]:>a /.-a_:>a]];
              loops = Select[legs,#[[1]] ~SameQ~ l&];
              ks = Select[legs,#[[1]] ~SameQ~ k&];
              doo = Union[Flatten[Outer[dot[#1,#2]&,legs,legs]/.dot[a_,a_]:>{}]];
              doo = Join[Complement[doo,free],
              Complement[{dot[Subscript[k, 2],Subscript[k, 3]],
              dot[Subscript[k, 2],Subscript[k, 4]],
              dot[Subscript[k, 3],Subscript[k, 4]],
              dot[Subscript[k, 1],Subscript[k, 2]]},free]];
              eqns = myTrees/.Atree[a__]:>(Plus@@a);
              eqns2 = Flatten[Map[(me = #;
                                   Map[(dot[#,me]==0)&,eqns])&,legs]]/.dot[a_,a_]:>0;
              Reduce[eqns2,doo,Backsubstitution->True]
          ]);

consistentGraphToTrees[vertexFormGraph[necklaces__]] :=
    necklaces/.neckl[a___]:>Atree[a]


corruptGraphNoMemo[graph_] := If[Head[graphHashCode[graph]] === Graph, graph, 
     Module[{goodSet, cG, StylePrint}, cG = minStripPriveledge[If[Head[graphHashCode[graph]] === Graph, graph, 
          Module[{intLegs = getIntLegs[graph], n = 0}, goodSet = {}; While[goodSet === {} && n < 24, 
             n++; StylePrint[StringJoin["Trying ", ToString[n]]]; goodSet = Select[Subsets[intLegs, {n}], 
                Head[graphHashCode[priveledgeSomeLegs[graph, Flatten[{#1}]]]] === Graph & , 1]; StylePrint[{"Got ", goodSet}]; ]; 
            If[n >= 24, n = 0; intLegs = Sort[getIntLegs[graph]]; goodSet = {}; While[goodSet === {} && n < 24, n++; 
                StylePrint[StringJoin["Trying ", ToString[n]]]; goodSet = Select[Subsets[intLegs, {n}], 
                  Head[graphHashCode[priveledgeSomeLegs[graph, Flatten[{#1}]]]] === Graph & , 1]; StylePrint[{"Got ", goodSet}]; ]; ]; 
            If[n >= 24, Throw["Some crappy corrupt graphs here.! n>24!!!"]]; priveledgeSomeLegs[graph, goodSet[[1]]]]]]; 
       If[Head[graphHashCode[cG]] =!= Graph, Print[{"Had to box oldschool!  minStripPriveledge insufficient for crazy graph.", 
           stripPriveledge[graph]}]; priveledgeSomeLegs[graph, goodSet[[1]]], cG]]]



Clear[corruptGraph];
corruptGraph[graph_] :=
    corruptGraph[graph] = If[ Head[graphHashCode[graph]] ~SameQ~ Graph,
                              graph,
                              Module[ {goodSet,cG,StylePrint},
                                  cG = minStripPriveledge[If[ Head[graphHashCode[graph]] ~SameQ~ Graph,
                                                              graph,
                                                              Module[ {intLegs = getIntLegs[graph],n = 0},
                                                                  goodSet = {};
                                                                  While[goodSet ~SameQ~ {}&&n<24,
                                                                     n++;
                                                                     StylePrint["Trying "<>ToString[n]];
                                                                     goodSet = Select[Subsets[intLegs,{n}],
                                                                        Head[graphHashCode[priveledgeSomeLegs[graph,
                                                                        Flatten[{#}]]]] ~SameQ~ Graph&,1];
                                                                     StylePrint[{"Got ", goodSet}];  
                                                                   ];
                                                                  If[ n>=24,
                                                                      n = 0;
                                                                      intLegs = Sort[getIntLegs[graph]];
                                                                      goodSet = {};
                                                                      While[goodSet ~SameQ~ {}&&n<24,
                                                                       n++;
                                                                       StylePrint["Trying "<>ToString[n]];
                                                                       goodSet = Select[Subsets[intLegs,{n}],
                                                                          Head[graphHashCode[priveledgeSomeLegs[graph,
                                                                          Flatten[{#}]]]] ~SameQ~ Graph&,1];
                                                                       StylePrint[{"Got ", goodSet}];  
                                                                      ];
                                                                  ];
                                                                  If[ n>=24,
                                                                      Throw["Some crappy corrupt graphs here.! n>24!!!"]
                                                                  ];
                                                                  priveledgeSomeLegs[graph,goodSet[[1]]]
                                                              ]
                                                          ]
                                      ];
                                  If[ Head[graphHashCode[cG]] ~UnsameQ~ Graph,
                                      Print[{"Had to box oldschool!  minStripPriveledge insufficient for crazy graph.",
                                      stripPriveledge[graph]}];
                                      priveledgeSomeLegs[graph,goodSet[[1]]],
                                      cG
                                  ]
                              ]
                          ]


(* ::Subsubsection:: *)
(*More plot code*)


doOrderedPlot[expr_, myExtLegs_, options___] := 
   Module[{labels, edges,styles,gpf, lgpf,labelHash},
gpf=graphPlotForm[expr];
edges=First/@gpf /. Rule:>DirectedEdge;
labels=Last/@gpf;
styles = Map[If[Head[#]===l,{"DashedLine",Red},Blue]&,labels];
labelHash=Association[ Thread[Rule[edges,Head/@(labels/.-a_:>a)]]];
lgpf=MapThread[Style[Labeled[#,#2],#3]&,{edges,labels/.highlight[a_]:>a,styles}]; GraphPlot[lgpf, DirectedEdges->True,Method->{"SpringElectricalEmbedding"},     MultiedgeStyle -> .50 , VertexCoordinateRules -> 
      MapIndexed[neckl[-{#1}] -> N[{-Cos[(#2[[1]] - 3/2)*2*(Pi/Length[myExtLegs])], 
           Sin[(#2[[1]] - 3/2)*2*(Pi/Length[myExtLegs])]}] & , myExtLegs] ,  EdgeRenderingFunction -> 
(Flatten[{If[ labelHash[#3/.DirectedEdge[a_,b_,c_]:>DirectedEdge[a,b]] === l, {Dashed, Red, 
           Arrowheads[{{0.05, 0.75}}], Thick, Arrow[#1] (*, Text[#3, Mean[#1], 
            Background -> Opacity[0.6, White]]*) }, If[ labelHash[#3/.DirectedEdge[a_,b_,c_]:>DirectedEdge[a,b]] ===highlight, {RGBColor[0.074189364461738, 0.5290608072022583, 
             0.0075379568169680325], Thick, Arrow[#1]} ,{Blue, Arrowheads[1/18], 
            Thick, Arrow[#1] (*,  Text[#3[[-1,1,1]], Mean[#1] ,Background -> Opacity[0.6, 
               White]]*) }]]}] & ), options]]

doOrderedPlotRand[expr_, myExtLegs_, options___] := 
  Module[{firstRules = 
     MapIndexed[
      neckl[-{#}] -> 
        N[{-Cos[(#2[[1]] - 3/2) 2 \[Pi]/Length[myExtLegs]], 
          Sin[(#2[[1]] - 3/2) 2 \[Pi]/Length[myExtLegs]]}] &, myExtLegs], 
    fisk,labels, edges,styles,gpf, lgpf,labelHash},
   fisk = graphPlotForm[expr] /. firstRules /. neckl[a___] :> Sow[neckl[a]] //
         Reap // Last // Flatten // Union;
   firstRules = 
    Join[firstRules, 
     MapIndexed[#1 -> 
        N[3/4 Mean[# /. neckl[a__] :> (neckl[{-#}] & /@ a) /. firstRules /. 
            neckl[___] :> RandomReal[{-2, -2}/3, 2]]
         ] &, fisk]];
   gpf=graphPlotForm[expr];
   edges=First/@gpf /. Rule:>DirectedEdge;
   labels=Last/@gpf;
   styles = Map[If[Head[#]===l,{"DashedLine",Red},Blue]&,labels];
   labelHash=Association[ Thread[Rule[edges,Head/@(labels/.-a_:>a)]]];
   lgpf=MapThread[Style[Labeled[#,#2],#3]&,{edges,labels/.highlight[a_]:>a,styles}]; 
   GraphPlot[lgpf, VertexCoordinateRules -> firstRules,
    MultiedgeStyle -> 50, 
    EdgeRenderingFunction -> (Flatten[{If[ labelHash[#3/.DirectedEdge[a_,b_,c_]:>DirectedEdge[a,b]] === l, {Dashed, Red, 
           Arrowheads[{{0.05, 0.75}}], Thick, Arrow[#1] (*, Text[#3, Mean[#1], 
            Background -> Opacity[0.6, White]]*) }, If[ labelHash[#3/.DirectedEdge[a_,b_,c_]:>DirectedEdge[a,b]] ===highlight, {RGBColor[0.074189364461738, 0.5290608072022583, 
             0.0075379568169680325], Thick, Arrow[#1]} ,{Blue, Arrowheads[1/18], 
            Thick, Arrow[#1] (*,  Text[#3[[-1,1,1]], Mean[#1] ,Background -> Opacity[0.6, 
               White]]*) }]]}] & ), 
    options]];

doOrderedPlot[expr_] := 
 Module[{myExtLegs = 
    mergeLegsTreeLevelList[consistentGraphToTrees[expr], 
       expr /. in[a_] :> Sow[in[a]] // Reap // Last // Flatten // Union] // 
      First // First},
  myExtLegs = 
   rotateList[myExtLegs, 
    Position[myExtLegs, First[Sort[myExtLegs]]] // Flatten // First]; 
  doOrderedPlot[expr, myExtLegs]]

cutDisplayRule[expr_] := 
 expr /. in[a_] :> "" /. 
  Graphics[a___] :> 
   Graphics @@ (List[a] /. {-l[b_] :> Style[-l[b], {Blue, Bold, Large}] , 
       l[b_] :> Style[l[b], {Blue, Bold, Large}],
       k[b_] :> Style[b, {Purple, Bold, Large}]})


cutDisplayRule[expr_] := 
 expr /. in[a_] :> "" /. 
  Graphics[a___] :> 
   Graphics @@ (List[a] /. {-l[b_] :> Style[-l[b], {Blue, Bold, Large}], 
        l[b_] :> Style[l[b], {Blue, Bold, Large}], 
        k[b_] :> Style[b, {Purple, Bold, Large}]})



cutGraph3D[graph_] := 
 Graph3D[HighlightGraph[
   mathematicaGraph[
    graph], (First /@ 
      Select[graphPlotForm[graph], (#[[-1]] // Head)  ~SameQ~  l &]) /. 
    Rule :> UndirectedEdge]]
 
stagePlot[sGraph_, order_] := 
 Manipulate[Module[{myLegs, StylePrint, dasRed, myTrees},
   myTrees = consistentGraphToTrees[sGraph][[ 1 ;; Round[x] ]];
   StylePrint[x];
   myLegs = myTrees /. Atree :> List /. -a_ :> a // Flatten // Union;
   StylePrint[myLegs];
          (dasRed[#] = True) & /@ myLegs;
   StylePrint[{#, dasRed[#]}] & /@ myLegs;
   Show[GraphicsRow[{doOrderedPlot[sGraph, k /@ Range[4], 
       EdgeRenderingFunction -> (({If[True  ~SameQ~  dasRed[#3], {Red, Arrow[#]}, 
             Arrow[#]], 
            Inset[#3, Mean[#1], Automatic, Automatic, #[[1]] - #[[2]], 
             Background -> White]}) &), PlotLabel -> x], 
      MatrixForm[myTrees /. Atree :> neckl]}]]], 
  {{x, 1}, 1, 
          Length[sGraph[[1]]], 1}] 

canonVertexList[
  vertexFormGraph[
   graph_]] := (# /. 
     neckl[a__] :> 
      neckl[rotateList[a, 
        Position[a, Sort[a] // First] // Flatten // First]]) & /@ graph


getOffShellDotsBetter[graph_] := 
    Module[{trees = consistentGraphToTrees[graph], ext, 
        int = getIntLegs[graph], legs}, 
      ext = getExtLegsFromTrees[trees]; 
       legs = Join[int, ext]; ToRules[
         (Reduce[ (daRInp = Flatten[#1]), 
                  (daList = 
                  Sort[Union[Flatten[Outer[dot[#1, #2] & , 
                           Join[ext[[{-1}]], int], Join[ext, int]]]], 
                     OrderedQ[(Abs[#1] /. {k[a_] :> a, l[b_] :> 
                                    2^b, dot[a_, b_] :> a*b} & ) /@ 
                           {#1, #2}] & ] // Reverse ), 
                   Backsubstitution -> True] & )[Join[(dot[#1, #1] == 0 & ) /@ ext, 
             Append[trees, Atree[ext]] /. Atree[a__] :> 
                 (dot[#1, Plus @@ a] == 0 & ) /@ legs]]]]


polGaugeAnsatzMaxPowers[gr_,param_,powerCycleDrop_] :=
    Module[{eLegs,LEGS,LOOPS,uniqLegs,
uniqDots,uniqPol,singlePol,doublePol,
maxPowerCount,curGuys,\[Epsilon]Used,nextBatch,lMaps,StylePrint},
lMz=lMaps=Map[Rule[#[[1]],#[[2]]-powerCycleDrop]&,getMyCycles[gr]];
StylePrint[lMaps];
      eLegs=getExtLegs[gr];
	LEGS=Length[eLegs];
LOOPS=Length[getMyCycles[gr]];
uniqLegs=Select[getMyUniqLegs[gr],(#/.Append[lMaps,a_:>2])>=1&];
StylePrint[uniqLegs];
uniqDots=Outer[dot[#1,#2]&, uniqLegs,uniqLegs]/.getIndepRules[{Atree[k/@Range[LEGS]]}]/.dot[k[a_],k[a_]]:>{}/.dot[a_,b_]:>Sow[dot[a,b]]//Reap//Last//Flatten//Union;
uniqDots=
Select[uniqDots,(me=#; And@@((0 ~SameQ~ (D[me/.dot[a_,b_]:> a*b,{#[[1]],1+#[[2]]}]))&/@lMaps))&];
StylePrint[uniqDots];
uniqPol=Map[\[Epsilon][#]&,eLegs];
singlePol=Outer[dot[#1,#2]&,uniqLegs,uniqPol]/.dot[First[uniqLegs],Last[uniqPol]]:>{}/.
dot[k[a_],\[Epsilon][k[a_]]]:>{}//Flatten//Union;
doublePol=
Outer[dot[#1,#2]&,uniqPol,uniqPol]/.dot[a_,a_]:>{}//Flatten//Union;
maxPowerCount=2*(LEGS+LOOPS-1)-2;
curGuys=Flatten[{singlePol,doublePol}];
cG=curGuys=If[Length[lMaps]>0,Select[curGuys,(me=#; And@@((0 ~SameQ~  
  (D[me/.dot[a_,b_]:> a*b,{#[[1]],1+#[[2]]}]))&/@lMaps))&],curGuys];
While[maxPowerCount>0,
maxPowerCount-=2;
curGuys=Table[\[Epsilon]Used=expr/.\[Epsilon][a_]:>Sow[\[Epsilon][a]]//Reap//Last//Flatten//Union;
lUsed=((Rule@@#)&/@(expr/.l[a_]:>Sow[l[a]]//Reap//Last//Flatten//Tally));
lUR=Map[#[[1]]->(#[[2]]-(#[[1]]/.lUsed/.l[a_]:>0))&,lMaps];
nextBatch=If[(Length[uniqPol]-Length[\[Epsilon]Used]) ~SameQ~ 0,
uniqDots,
If[(Length[uniqPol]-Length[\[Epsilon]Used])-1<=maxPowerCount,
(* i.e. I can afford one mistake :-) *)
{singlePol,doublePol},
doublePol]];
nextBatch=(nextBatch/.(dot[_,#]:>{}&/@\[Epsilon]Used))//Flatten//Union;
nextBatch=Select[nextBatch, ((#/.dot[a_,b_]:>a/.lUR/.a_[b_]:>10)>0&&(#/.dot[a_,b_]:>b/.lUR/.a_[b_]:>10)>0)&];

If[nextBatch ~SameQ~ {},{},
Map[expr*#&,nextBatch]],
{expr,curGuys}]//Flatten//Union;
];
MapIndexed[a[param,#2[[1]]]*#&,curGuys]//Total
]


polGaugeAnsatz[gr_,param_] :=
    Module[{eLegs,LEGS,LOOPS,uniqLegs,
uniqDots,uniqPol,singlePol,doublePol,
maxPowerCount,curGuys,\[Epsilon]Used,nextBatch},
      eLegs=getExtLegs[gr];
	LEGS=Length[eLegs];
LOOPS=Length[getMyCycles[gr]];
uniqLegs=getMyUniqLegs[gr];
uniqDots=Map[#[[2]]/.dot[a_,b_]:>Sow[dot[a,b]]&,getOffShellDots[gr]]  //Reap//Last//Flatten//Union;
uniqPol=Map[\[Epsilon][#]&,eLegs];
singlePol=Outer[dot[#1,#2]&,uniqLegs,uniqPol]/.dot[First[uniqLegs],Last[uniqPol]]:>{}//Flatten//Union;
doublePol=
Outer[dot[#1,#2]&,uniqPol,uniqPol]/.dot[a_,a_]:>{}//Flatten//Union;
maxPowerCount=2*(LEGS+LOOPS-1)-2;
curGuys=Flatten[{singlePol,doublePol}];
While[maxPowerCount>0,
maxPowerCount-=2;
curGuys=Table[\[Epsilon]Used=expr/.\[Epsilon][a_]:>Sow[\[Epsilon][a]]//Reap//Last//Flatten//Union;
nextBatch=If[(Length[uniqPol]-Length[\[Epsilon]Used]) ~SameQ~ 0,
uniqDots,
If[(Length[uniqPol]-Length[\[Epsilon]Used])-1<=maxPowerCount,
(* i.e. I can afford one mistake :-) *)
{singlePol,doublePol},
doublePol]];
nextBatch=(nextBatch/.(dot[_,#]:>{}&/@\[Epsilon]Used))//Flatten//Union;If[nextBatch ~SameQ~ {},{},
Map[expr*#&,nextBatch]],
{expr,curGuys}]//Flatten//Union;
];
MapIndexed[a[param,#2[[1]]]*#&,curGuys]//Total
]


(* ::Subsubsection::Closed:: *)
(*Color rep code*)


colorForm[vertexFormGraph[gr_]] :=
    colorForm[consistentGraphToTrees[vertexFormGraph[gr]]]

colorForm[trees__] :=
    newColorSimplify[Times @@ trees /. Atree[a__] :> structureF[a] /. 
    l[a_] :> b[a] /. -(a_) :> a /. k[i_] :> a[i]]

tr[a_,a_] :=
    Nc^2-1
tr[] :=
    Nc
tr[x_] :=
    0
tr[x___,a_,y___,a_,z___] :=
    tr[x,z] tr[y]-tr[x,y,z]/Nc
tr /:tr[w___,a_]^2:=tr[w,w]-(tr[w] tr[w])/Nc
tr/:tr[w___,a_,x___] tr[y___,a_,z___]:=tr[x,w,z,y]-(tr[x,w] tr[z,y])/Nc

ClearAll[newColorSimplify]

newColorSimplify[expr_] :=
    FixedPoint[ ExpandAll[#]&,expr/.structureF[{a_,b_,c_}]:>tr[a,b,c]-tr[b,a,c]]/.tr[x___,a[2],y___]:>tr[a[2],y,x]/.
    tr[x___,a[1],y___]:>tr[a[1],y,x] /.{
    tr[a[4],a[3]]->tr[a[3],a[4]],tr[a[4],a[2]]->tr[a[2],a[4]],tr[a[3],a[2]]->tr[a[2],a[3]]}/.tr[a___]:>ctr[a]


(* ::Subsubsection:: *)
(*Momenta and spinor code*)


spa[x__] :=
    Signature[{x}]*spa @@ Sort[{x}] /;  !OrderedQ[{x}] && FreeQ[{x}, Pattern] 
spa[jj_, jj_] :=
    0


spb[x__] :=
    Signature[{x}]*spb @@ Sort[{x}] /;  !OrderedQ[{x}] && FreeQ[{x}, Pattern]
spb[jj_, jj_] :=
    0


(* Complex spinor code *)
iV[p_] :=
    Module[ {q = Sqrt[p[[1]]+p[[4]]]},
        {q,(p[[2]]+I p[[3]])/q,0,0}
    ]
iU[p_] :=
    Module[ {q = Sqrt[p[[1]]+p[[4]]]},
        {0,0,q,(p[[2]]-I p[[3]])/q}
    ]

Ui[p_] :=
    Module[ {a = iU[p]},
        {0,0,-a[[4]],a[[3]]}
    ]

Vi[p_] :=
    Module[ {b = iV[p]},
        {b[[2]],-b[[1]],0,0}
    ]

up[p_] :=
    Module[ {q = Sqrt[p[[1]]+p[[4]]]},
        {q,(p[[2]]+I p[[3]])/q,0,0}
    ]

um[p_] :=
    Module[ {q = Sqrt[p[[1]]+p[[4]]]},
        {0,0,q,(p[[2]]-I p[[3]])/q}
    ]

upb[p_] :=
    Module[ {a = iU[p]},
        {0,0,-a[[4]],a[[3]]}
    ]

umb[p_] :=
    Module[ {b = iV[p]},
        {b[[2]],-b[[1]],0,0}
    ]

gammamu = {{
{0,0,0,1},
{0,0,-1,0},
{0,-1,0,0},
{1,0,0,0}},
-1{
{0,0,-1,0},
{0,0,0,1},
{1,0,0,0},
{0,-1,0,0}
},
-1{{0,0,I,0},
{0,0,0,I},
{I,0,0,0},
{0,I,0,0}},
-1{{0,0,0,1},
{0,0,1,0},
{0,-1,0,0},
{-1,0,0,0}}};

spaEval[a_,b_] :=
    umb[a] . up[b]
spbEval[a_,b_] :=
    upb[a] . um[b]
polMinus[k_,q_] :=
    - (Map[upb[q] . # . up[k]&,gammamu]/(Sqrt[2]  spbEval[q,k]))
polPlus[k_,q_] :=
    +  (Map[umb[q] . # . um[k]&,gammamu]/(Sqrt[2] spaEval[q,k]))
complexSpinorProd[q_] :=
    {spa[a_, b_] :> umb[q[a]] . up[q[b]], spb[a_, b_] :> upb[q[a]] . um[q[b]]}

(* Cached Momenta Tools  -- needs updating to the land of mostly s*)

refreshHLP :=
    (Clear[hLP,hLsq];
     hLP[a_,b_,myMomenta_] :=
         hLP[a,b,myMomenta] = Lorentz[myMomenta[a],myMomenta[b]];
     hLsq[a__,myMomenta_] :=
         hLsq[a,myMomenta] = Lsq[Plus@@(myMomenta/@a)];)
evaluateLpRule[myMomenta_] :=
    {dot[a_,b_]:> hLP[a,b,myMomenta]}
evaluateLsqRule[myMomenta_] :=
    {uLsq[a___]:>hLsq[Flatten[{a/.Plus:>List}],myMomenta]}
evaluateSRule[myMomenta_] :=
    {s[a___]:>hLsq[Flatten[{a}/.Plus:>List],myMomenta]}
applyMomentum[mom_][expr_] :=
    expr/.evaluateLpRule[mom]/.evaluateLsqRule[mom]/.complexSpinorProd[mom]
refreshHLP

getRandomNullVector[prec_,D_] :=
    Module[ {v = Table[ Random[Real, {-1, 1}, prec],{i,1,D-1}]},
        Prepend[ -1^(RandomInteger[{0,1}])*v,Sqrt[v . v]]
    ]

generateNullMomenta[l_, prec_,D_] :=
    Module[ {length = Length[l], Kvecs, K, \[Rho], l1, l2, momenta, StylePrint},
        StylePrint[length];
        Kvecs = (getRandomNullVector[prec,D] & ) /@ Range[length - 2];
        StylePrint[N[Kvecs]];
        K = Plus @@ Kvecs;
        \[Rho] = getRandomNullVector[prec,D];
        StylePrint[N[\[Rho]]];
        l1 = -((Lsq[K]*\[Rho])/(2*Lorentz[\[Rho], K]));
        l2 = -K - l1;
        allVals = Join[Kvecs, {l1, l2}];
        ((momenta[l[[#1]]] = allVals[[#1]];
          momenta[-l[[#1]]] = -allVals[[#1]]; ) & ) /@ 
        Range[length];
        momenta
    ]

generateMasslessDecay[massiveMom_, l__, prec_,D_] :=
    Module[ {length = Length[l], Kvecs, K, \[Rho], l1, l2, StylePrint, momenta, allVals},
        StylePrint[l];
        Kvecs = If[ length ~SameQ~ 2,
                    {0*massiveMom},
                    (getRandomNullVector[prec,D] & ) /@ Range[length - 2]
                ];
        StylePrint[Kvecs];
        K = Plus @@ Kvecs + massiveMom;
        StylePrint[K];
        \[Rho] = getRandomNullVector[prec,D];
        StylePrint[Lsq[K]];
        l1 = (-(Lsq[K]/(2*Lorentz[\[Rho], K])))*\[Rho];
        l2 = -K - l1;
        allVals = If[ length  ~SameQ~  2,
                      {l1, l2},
                      Join[Kvecs, {l1, l2}]
                  ];
        ((momenta[l[[#1]]] = allVals[[#1]];
          momenta[-l[[#1]]] = -allVals[[#1]]; ) & ) /@ Range[length];
        momenta
    ]

buildMomenta[mtrees__,D_] :=
    Module[ {myLegs, firstCut, known, notknown, mom, StylePrint, myDiscover, myTree, myKnown, myNextMom, trees = mtrees, Print, bob},
        StylePrint[trees];
        myLegs = (#1[[1]] & ) /@ Rest[trees];
        firstCut = First[trees][[1]];
        mom = generateNullMomenta[firstCut, PREC,D];
        known = Join[firstCut, -firstCut];
        notknown = Complement[Union[Flatten[myLegs] /. -(a_) :> a], known];
        ((myTree = #1;
          StylePrint[myTree];
          myKnown = Select[myTree, Count[known, #1, Infinity] + Count[known, -#1, Infinity] > 0 & ];
          StylePrint[{"Known!", myKnown}];
          myDiscover = Complement[myTree, myKnown];
          StylePrint[{"Discover!", myDiscover}];
          If[ Length[myKnown] ~UnsameQ~ Length[myTree] && Length[myDiscover] > 1,
              myNextMom = generateMasslessDecay[Plus @@ (mom[#1] & ) /@ myKnown, myDiscover, 850,D];
              StylePrint[N[myNextMom /@ myDiscover]];
              ((mom[#1] = myNextMom[#1];
                mom[-#1] = -myNextMom[#1]; ) & ) /@ Complement[myTree, myKnown];
              known = Join[known, myDiscover, -myDiscover];,
              Print[Style[{"Warning!", myTree, myKnown, Lsq[bob = Plus @@ (mom[#1] & ) /@ myKnown]}]];
              If[ Chop[Lsq[bob]]  ~UnsameQ~  0 && Length[myDiscover] ~SameQ~ 1,
                  Print[Style["Warning Dude"]];
                  Throw["Bad Sort!"];
              ];
              If[ Chop[Lsq[bob]] ~SameQ~ 0 && Length[myDiscover] ~SameQ~ 1,
                  mom[myDiscover[[1]]] = -bob;
                  mom[-myDiscover[[1]]] = bob;
                  known = Join[known, myDiscover, -myDiscover];
                  StylePrint[{"Now known", known}];
              ];
          ]; ) & ) /@ myLegs;
        mom
    ]

buildAllMomentaD[mtrees_,D_] :=
    Module[ {rmom, good, trees = mtrees},
        good = Select[Select[(Catch[buildMomenta[#1,D]] & ) /@ Permutations[trees],  !StringQ[#1] & ], 
           (rmom = #1;
            Chop[trees /. Atree[a__] :> Plus @@ rmom /@ a] == (trees /. Atree[a__] :> Table[0,{D}])) & ];
        If[ Length[good]  ~SameQ~  0,
            Throw["bad Momenta"]
        ];
        good[[1]]
    ]



generateNullMomenta[l_, prec_] :=
    Module[ {length = Length[l], Kvecs, K, \[Rho], l1, l2, 
      momenta,StylePrint},
        StylePrint[length];
        Kvecs = (getRandomNullVector[prec] & ) /@ Range[length - 2];
        StylePrint[N[Kvecs]];
        K = Plus @@ Kvecs;
        \[Rho] = getRandomNullVector[prec];
        StylePrint[N[\[Rho]]];
        l1 = -((Lsq[K]*\[Rho])/(2*Lorentz[\[Rho], K]));
        l2 = -K - l1;
        allVals = Join[Kvecs, {l1, l2}];
        ((momenta[l[[#1]]] = allVals[[#1]];
          momenta[-l[[#1]]] = 
          -allVals[[#1]]; ) & ) /@ Range[length];
        momenta
    ]

buildAllMomenta[mtrees_,n__] :=
    Module[ {rmom, good, trees = mtrees},
        good = 
        Select[
        Select[
        (Catch[buildMomenta[#1]] & ) /@ Permutations[trees],  (!StringQ[#1]) & ,1], 
           (rmom = #1;
            Chop[trees /. Atree[a__] :> Plus @@ rmom /@ a] == (trees /. Atree[a__] :> {0, 0, 0, 0})) & ];
        If[ Length[good]  ~SameQ~  0,
            Throw["bad Momenta"]
        ];
        good[[1]]
    ]

buildMomenta[mtrees__] :=
    Module[ {myLegs, firstCut, known, notknown, mom, StylePrint, myDiscover, myTree, myKnown, myNextMom, trees = mtrees, Print, bob},
        StylePrint[trees];
        myLegs = (#1[[1]] & ) /@ Rest[trees];
        firstCut = First[trees][[1]];
        mom = generateNullMomenta[firstCut, PREC];
        known = Join[firstCut, -firstCut];
        notknown = Complement[Union[Flatten[myLegs] /. -(a_) :> a], known];
        ((myTree = #1;
          StylePrint[myTree];
          myKnown = Select[myTree, Count[known, #1, Infinity] + Count[known, -#1, Infinity] > 0 & ];
          StylePrint[{"Known!", myKnown}];
          myDiscover = Complement[myTree, myKnown];
          StylePrint[{"Discover!", myDiscover}];
          If[ Length[myKnown] ~UnsameQ~ Length[myTree] && Length[myDiscover] > 1,
              myNextMom = generateMasslessDecay[Plus @@ (mom[#1] & ) /@ myKnown, myDiscover, 850];
              StylePrint[N[myNextMom /@ myDiscover]];
              ((mom[#1] = myNextMom[#1];
                mom[-#1] = -myNextMom[#1]; ) & ) /@ Complement[myTree, myKnown];
              known = Join[known, myDiscover, -myDiscover];,
              Print[Style[{"Warning!", myTree, myKnown, Lsq[bob = Plus @@ (mom[#1] & ) /@ myKnown]}]];
              If[ Chop[Lsq[bob]] ~UnsameQ~ 0 && Length[myDiscover] ~SameQ~ 1,
                  Print[Style["Warning Dude"]];
                  Throw["Bad Sort!"];
              ];
              If[ Chop[Lsq[bob]] ~SameQ~ 0 && Length[myDiscover] ~SameQ~ 1,
                  mom[myDiscover[[1]]] = -bob;
                  mom[-myDiscover[[1]]] = bob;
                  known = Join[known, myDiscover, -myDiscover];
                  StylePrint[{"Now known", known}];
              ];
          ]; ) & ) /@ myLegs;
        mom
    ]

generateMassiveDecay[mSq_,l__,prec_] :=
    Module[ {length = Length[l], Kvecs, K, \[Rho],b, l1, l2, momenta,allVals},
        StylePrint[length];
        Kvecs = (getRandomMassVector[prec] & ) /@ Range[length - 2];
        StylePrint[N[Kvecs]];
        K = Plus @@ Kvecs;
        \[Rho] = If[ mSq ~SameQ~ 0,
                     getRandomNullVector[prec],
                     b = getRandomMassVector[prec];
                     b*Sqrt[mSq]/Sqrt[Lsq[b]]
                 ];
        StylePrint[{Lsq[\[Rho]],Sqrt[Lsq[\[Rho]]]}];
        l1 = \[Rho]-K/Random[Real,{1/100,99/100},prec];
        l2 = \[Rho]-l1-K;
        StylePrint[Lsq[l1+l2+K]];
        allVals = Join[Kvecs, { l1,l2}];
        ((momenta[l[[#1]]] = allVals[[#1]];
          momenta[-l[[#1]]] = -allVals[[#1]]; ) & ) /@ Range[length];
        momenta
    ]

generateMasslessDecay[massiveMom_,l__,prec_] :=
    Module[ {length = Length[l], Kvecs, K, \[Rho], l1, l2,StylePrint, momenta,allVals},
        StylePrint[l];
        Kvecs = If[ length~SameQ~ 2,
                    {0*massiveMom},
                    (getRandomNullVector[prec] & ) /@ Range[length - 2]
                ];
        StylePrint[Kvecs];
        K = Plus @@ Kvecs+massiveMom;
        StylePrint[K];
        \[Rho] = getRandomNullVector[prec];
        StylePrint[Lsq[K]];
        l1 = (-(Lsq[K]/(2*Lorentz[\[Rho], K])))*\[Rho];
        l2 = -K - l1;
        allVals = If[ length ~SameQ~ 2,
                      {l1,l2},
                      Join[Kvecs, {l1, l2}]
                  ];
        ((momenta[l[[#1]]] = allVals[[#1]];
          momenta[-l[[#1]]] = -allVals[[#1]]; ) & ) /@ Range[length];
        momenta
    ]

getRandomNullVector[prec_] :=
    Module[ {x = Random[Real,{-1,1},prec],
    y = Random[Real,{-1,1},prec],
    z = Random[Real,{-1,1},prec]},
        {Sqrt[{x,y,z} . {x,y,z}],x,y,z}
    ];

getRandomMassVector[prec_] :=
    (#/#[[1]])&[getRandomNullVector[prec]+getRandomNullVector[prec]]

buildAllMomenta[mtrees_,n__] :=
    Module[ {rmom, good, trees = mtrees},
        good = 
        Select[
        Select[
        ({#,Catch[buildMomenta[#1]]} & ) /@ Permutations[trees],  (!StringQ[#1[[2]]]) & ], 
           (rmom = #1[[2]];
            Chop[trees /. Atree[a__] :> Plus @@ rmom /@ a] == (trees /. Atree[a__] :> {0, 0, 0, 0})) &,1 ];
        If[ Length[good]  ~SameQ~  0,
            Throw["bad Momenta"]
        ];
        Table[buildMomenta[good[[1,1]]],{i,1,n}]
    ]


(* ::Subsubsection:: *)
(*4D SUSY cuts *)


mergeLegsTreeLevel[trees_, a_] :=
    Module[ {bt = Select[trees, Count[#, a /. -b_ :> b, \[Infinity]] > 0 &], gt = a /. -b_ :> b},
        If[ Length[bt]  ~UnsameQ~  2,
            Throw[{"Bad collapse ", trees, a, bt}]
        ];
        rules = {bt[[2]] -> {},
          bt[[1]] -> (bt[[1]] /. (-gt :> gt) /. gt :> Rest[rotateList[bt[[2]][[1]], Flatten[Position[bt[[2]][[1]], gt]][[1]]]])};
        Flatten[trees /. rules /. Atree[b__] :> Atree[Flatten[b]]]
    ]

Clear[mergeLegsTreeLevelList];
mergeLegsTreeLevelList[trees_, a__] :=
    Module[ {foo = trees,rest = a},
        While[Length[rest]>0,foo = mergeLegsTreeLevel[foo,First[rest]];
                             rest = Rest[rest];
        ];
        foo
    ]

getIntLegs[graph_] :=
    (#[[1]]&/@Select[Tally[ graph[[1]]/.neckl[a__]:>a/.-a_:>a//Flatten],#[[2]] ~SameQ~ 2&])

getIntLegsFromTrees[trees_] :=
    (#[[1]]&/@Select[Tally[ trees/.Atree[a__]:>a/.-a_:>a//Flatten],#[[2]] ~SameQ~ 2&])

getExtLegsFromTrees[trees_] :=
    (#1[[1]] & ) /@ 
    Select[Tally[Flatten[trees /. Atree[a__] :> a /. 
    -(a_) :> a]], #1[[2]]  ~SameQ~  1 & ]


getKRules[trees_] :=
    ToRules[Reduce[trees /. Atree[a__] :> 0 == Plus @@ a, 
    Append[Sort[getIntLegsFromTrees[trees]], Last[Sort[getExtLegsFromTrees[trees]]]]]]



(* ::Subsubsection:: *)
(*KK & BCJ basis*)


kkExpress[Atree[legs_],first_, second_] :=
    Module[ {\[Alpha],\[Beta],pems,legsA,legsU, one, nnn,sign = 1,pos,StylePrint},
        pos[a_,l_] :=
            Position[l,a][[1,1]];
        one = second;
        nnn = first;
        If[ Intersection[legs,{first,second}] ~UnsameQ~ Sort[{first,second}],
            Throw[notInKK[{first,second,Atree[legs]}]];
        ];
        If[ (pos[first,legs]-pos[second,legs]) ~SameQ~ -1,
            Return[Atree[rotateList[legs,pos[first,legs]]]],
            If[ (pos[first,legs]-pos[second,legs]) ~SameQ~ 1,
                Return[(-1)^Length[legs] Atree[rotateList[Reverse[legs],pos[first,Reverse[legs]]]]],
                If[ Position[legs,nnn][[1,1]]<Position[legs,one][[1,1]],
                    legsU = Reverse[legs];
                    sign = (-1)^Length[legs];
                    StylePrint[{"rotatedLegs",legsU}];,
                    legsU = legs;
                ];
                legsU = rotateList[legsU,Position[legsU,one][[1,1]]];
                StylePrint[{"rotated: legs U",legsU}];
                StylePrint[{"conditional eq: ",legsU[[-1]],nnn}];
                If[ legsU[[-1]] ~SameQ~ nnn,
                    Return[ sign * Atree[rotateList[legsU,Length[legsU]]]],
                    StylePrint[{legsU,pos[one,legsU], pos[nnn,legsU]}];
                    \[Alpha] = legsU[[ Range[pos[one,legsU]+1, pos[nnn,legsU]-1]]];
                    \[Beta] = Reverse[legsU[[ Range[pos[nnn,legsU]+1, Length[legsU]]]]];
                    StylePrint[\[Alpha]];
                    StylePrint[\[Beta]];
                    StylePrint[sign];
                    sign*Plus@@Table[(-1)^Length[\[Beta]] Atree[Join[{nnn},{one},perm]],
                    {perm,Select[Permutations[Join[\[Alpha],\[Beta]]],(\[Alpha] ~SameQ~  Flatten[#/.Thread[\[Beta]->Map[{}&,\[Beta]]]] )&&
                    (\[Beta] ~SameQ~  Flatten[#/.Thread[\[Alpha]->Map[{}&,\[Alpha]]]] )&]}
                    ]
                ]
            ]
        ]
    ]



bcj\[ScriptCapitalG][i_,j_] :=
    If[ (i<j)||(j ~SameQ~ 1)||(j ~SameQ~ 3),
        uLsq[{i,j}],
        0
    ]


bcj\[ScriptCapitalF][3,\[Sigma]_,1,k_,m_,n_] :=
    Module[ {\[Rho] = Join[{3},\[Sigma],{1}],
    t},
        t[kk_] :=
            If[ kk ~SameQ~ (1+m),
                0,
                If[ kk ~SameQ~ 3,
                    t[5],
                    Position[\[Rho],kk][[1,1]]
                ]
            ];
        If[ t[k-1]<t[k],
            Sum[bcj\[ScriptCapitalG][k,\[Rho][[l]]],{l,t[k],n-1}],
            -Sum[bcj\[ScriptCapitalG][k,\[Rho][[l]]],{l,1,t[k]}]
        ]+
        If[ (t[k-1]<t[k])&&(t[k+1]>t[k]),
            uLsq[Prepend[Range[4,k],2]],
            0
        ]+
        If[ (t[k-1]>t[k])&&(t[k+1]<t[k]),
            -uLsq[Prepend[Range[4,k],2]],
            0
        ]
    ]


bcjReorder[Atree[{first_,second_,\[Alpha]___,third_,\[Beta]___}],first_,second_,third_] :=
    Module[ {rule = Reverse/@Join[{first->1,second->2,third->3},Thread[{\[Alpha]}->(Range[Length[{\[Alpha]}]]+3)],
    Thread[{\[Beta]}->Range[Length[{\[Beta]}]]+Length[{\[Alpha]}]+3]],
    m = Length[{\[Alpha]}]+3,n = Length[{\[Beta]}]+Length[{\[Alpha]}]+3,
    \[Alpha]\[Alpha] = (Range[Length[{\[Alpha]}]]+3),
    \[Beta]\[Beta] = Range[Length[{\[Beta]}]]+Length[{\[Alpha]}]+3},
        Plus@@Table[ Atree[Join[{1,2,3},\[Sigma] ]/.rule] *
        Product[bcj\[ScriptCapitalF][3,\[Sigma],1,k,m,n]/uLsq[Prepend[Range[4,k],2]],{k,4,m}],
        {\[Sigma],Select[Permutations[Join[\[Alpha]\[Alpha],\[Beta]\[Beta]]],
        \[Beta]\[Beta] ~SameQ~  Flatten[#/.Thread[\[Alpha]\[Alpha]->Map[{}&,\[Alpha]\[Alpha]]]]&]}]/.uLsq[a__]:>uLsq[a/.rule]
    ]


bcjExpress[Atree[legs__],first_,second_,third_] :=
    Module[ {expr = kkExpress[Atree[legs],first,second]},
        expr/.Atree[l_]:>bcjReorder[Atree[l],first,second,third]
    ]


bcjExpress[Atree[legs__,spins__],first_,second_,third_] :=
    bcjExpress[Atree[legs],first,second,third]/.Atree[a__]:>Atree[a,a/.Thread[Rule[legs,spins]]]


(* ::Subsubsection:: *)
(*4d SUSY Cut*)


collapseTrees[trees_] :=
    Module[ {propTarget,legs,subsets,happyHash,happySame,StylePrint,goodCandidates,candidates,uniqCandidates},
        StylePrint[trees];
        goodCandidates = {};
        propTarget = Length[trees[[1]]]-1;
        StylePrint[propTarget];
        While[propTarget>=1,
                candidates = Flatten[Table[
        legs = getIntLegsFromTrees[tree];
        StylePrint[{tree, legs}];
        subsets = Subsets[legs,{propTarget}];
        StylePrint[{propTarget,subsets,tree}];
        Table[mergeLegsTreeLevelList[tree,subset],{subset,subsets}],
        {tree,trees}],1];
                StylePrint[{"candidates", candidates}];
                happyHash[tree_] :=
                    happyHash[tree] = Sort[tree/.Atree[a_]:>Atree[a/.getKRules[tree]]/.Atree[a_]:>rotateList[a,Ordering[a,1][[1]]]];
                happySame[a_,b_] :=
                    (happyHash[a] ~SameQ~ happyHash[b]);
                uniqCandidates = Union[candidates,SameTest->happySame];
                StylePrint[{"uniq Cand",uniqCandidates}];
                goodCandidates = Append[goodCandidates,uniqCandidates];
                StylePrint[{"good Cand",goodCandidates}];
                propTarget--;
        ];
        goodCandidates
    ]

      
getCutPrecursor[cut_] :=
    Module[ {StylePrint,candStuff =
    collapseTrees[consistentGraphToTrees/@treesToLoops[{#}]]&/@cut,sewCut,cutCandidates},
        candStuff = MapIndexed[#1/. in[a_]:>in[#2[[1]]*100+a]&,candStuff];
        cutCandidates = (Outer[sewCut,##,2]&@@candStuff) /.
        sewCut[a__]:>sewCut[Flatten[{a}]]//Flatten;
        StylePrint[Length[cutCandidates]];
        cutCandidates/.sewCut[a__]:>a
    ]     

nmhvLevel[list_] :=
    -2-Plus@@(list/. 1:>0)

getMHVRules[treeList_,hashRule_] :=
    Module[ {candidates = Select[Union[
    Flatten[treeList/.Atree[a__]:>a/.hashRule/.-a_:>a]],!NumberQ[#]&],
    maxDigits,goodRules,rules,lC},
        zeTL = treeList;
        zeHR = hashRule;
(*maxDigits=FromDigits[Table[1,{Length[candidates]}],2]; *)
        maxDigits = 2^(lC = Length[candidates])-1;
        rules =
        (Thread[candidates->(2IntegerDigits[#,2,lC]-1)])&/@Range[0,maxDigits];
        (* goodRules=Select[rules,(me=#;(Map[nmhvLevel[#[[1]]/.me]&,treeList]//Union) ~SameQ~ {0})&] ;*)
        tL = treeList/.hashRule;
        goodRules = Map[Join[hashRule,#]&,Select[rules,(And@@Map[nmhvLevel[#[[1]]] ~SameQ~ 0&,tL/.#])&] ]
    ]

goodMHV[treeList__] :=
    nmhvLevel[First[treeList][[1]]] ~SameQ~ 0 &&(Length[treeList] ~SameQ~ 1||
    goodMHV[Rest[treeList]])

goodMHV[treeList__,rules_] :=
    treeList ~SameQ~ {}||(nmhvLevel[treeList[[1,1]]/.rules] ~SameQ~ 0 &&
    goodMHV[Rest[treeList],rules])

getMHVRules2[treeList_,hashRule_] :=
    Module[ {candidates = Select[Union[Flatten[treeList/.Atree[a__]:>a/.hashRule/.-a_:>a]],!NumberQ[#]&],maxDigits,goodRules,rules,lC},
        zeTL = treeList;
        zeHR = hashRule;
(*maxDigits=FromDigits[Table[1,{Length[candidates]}],2]; *)
        maxDigits = 2^(lC = Length[candidates])-1;
        rules =
        ((Thread[candidates->(2IntegerDigits[#,2,lC]-1)])&/@Range[0,maxDigits]);
        (* goodRules=Select[rules,(me=#;(Map[nmhvLevel[#[[1]]/.me]&,treeList]//Union) ~SameQ~ {0})&] ;*)
        tL = treeList/.hashRule;
        goodRules = Map[Join[hashRule,#]&,Select[rules,goodMHV[tL,#]&] ]
    ]


selectCutList[cutList_,extRule_] :=
    Module[ {ccntr = 0},
        Select[(ccntr++;
        (* If[Mod[ccntr,1000] ~SameQ~ 0,
        StylePrint[{ccntr,Length[cutList],Date[]}]];*)
                {#,getMHVRules[#,extRule]})&/@(cutList),#[[2]] ~UnsameQ~ {}&]
    ]

cutSig[cut_] :=
    Signature[(Flatten[(cut /. Atree[a_,b_]:>(a*(b/.(1->{}))))/.-a_:>a])] (* Take only negative helicity guys from the cut, and that's your signature *)

getKRules2[trees_,guys_] :=
    (Rule @@ #1 & ) /@ 
    Flatten[{Reduce[trees /. Atree[a__] :> 0 == Plus @@ a,guys ] /. 
    And :> List}]

processCut[{cut__,ruleSets__},SUSY_] :=
    Module[ {
    flatGuys = Select[cut/.Atree[a__]:>a/.-a_:>a//Flatten//Union,Head[#] ~SameQ~ in&],num,kRules,
    denom = Times@@(cut/.Atree[a__]:>{Times@@Thread[foo[a,rotateList[a,2]]]/.foo[i_,j_]:>spa[i,j]})},
        kRules = getKRules2[cut,flatGuys];
        num = (Sum[cutSig[cut /. Atree[a__]:>Atree[a,a/.rule]]*Times@@ cut/.Atree[a__]:>spa@@Select[a,(#/.rule) ~SameQ~ -1&],{rule,ruleSets}]^SUSY)*
        If[ SUSY<4,
            (
            Sum[(cutSig[cut /. Atree[a__]:>Atree[a,a/.rule]]*Times@@ cut/.Atree[a__]:>spa@@Select[a,(#/.rule) ~SameQ~ -1&])^(4-SUSY),{rule,ruleSets}]),
            1
        ];
        (num/denom/.in[a_]:>flat[in[a]/.kRules])*1/(Times@@(uLsq[Flatten[{#/.kRules}/.Plus:>List]]&/@flatGuys))
    ]

processCut[{cut__,ruleSets__}] :=
    Module[ {
    flatGuys = Select[cut/.Atree[a__]:>a/.-a_:>a//Flatten//Union,Head[#] ~SameQ~ in&],num,kRules,
    denom = Times@@(cut/.Atree[a__]:>{Times@@Thread[foo[a,rotateList[a,2]]]/.foo[i_,j_]:>spa[i,j]})},
    StylePrint[{"Using default globabl SUSY", SUSY}];
        kRules = getKRules2[cut,flatGuys];
        num = (Sum[cutSig[cut /. Atree[a__]:>Atree[a,a/.rule]]*Times@@ cut/.Atree[a__]:>spa@@Select[a,(#/.rule) ~SameQ~ -1&],{rule,ruleSets}]^SUSY)*
        If[ SUSY<4,
            (
            Sum[(cutSig[cut /. Atree[a__]:>Atree[a,a/.rule]]*Times@@ cut/.Atree[a__]:>spa@@Select[a,(#/.rule) ~SameQ~ -1&])^(4-SUSY),{rule,ruleSets}]),
            1
        ];
        (num/denom/.in[a_]:>flat[in[a]/.kRules])*1/(Times@@(uLsq[Flatten[{#/.kRules}/.Plus:>List]]&/@flatGuys))
    ]

(* Resolve sign ambiguity of spinor products based on helicity.
   cf. Blackhat and  Freedman arxiv 0808.1720 eq. 4.6
*)

FermionArrowSign[expr_] :=
    Module[ { result},
        result = expr;
        result = result//. {spa[a1_,a2_] :> -spa[-a1,a2] /; a1 < 0, 
                            spa[a1_,a2_] :> -spa[a1,-a2] /; a2 < 0,
                            spb[a1_,a2_] :>  spb[-a1,a2] /; a1 < 0,
                            spb[a1_,a2_] :>  spb[a1,-a2] /; a2 <0};
        result = result //.{spa[-a1_, a2_] :> -spa[a1,a2],
                            spa[a1_, -a2_] :> -spa[a1,a2],
                            spb[-a1_, a2_] :>  spb[a1,a2],
                            spb[a1_, -a2_] :>  spb[a1,a2]};
        Return[result]
    ]       

processCutList[goodCuts_] :=
    Flatten[(processCut/@goodCuts)]


buildGraphDenom[graph_] :=
    Times@@(uLsq[Flatten[{#1/.Plus:>List}]]&)/@(Select[(#1[[2]]&)/@graphPlotForm[graph],Head[#1] ~SameQ~ in&]/.graphAlgRules[graph])^(-1)

buildRHSNumberFromMom[trees_, parents_,mom_] :=
    (CCCCjj = 0;
     Clear[recBeast];
     Plus @@ Flatten[(CCCCjj++;
                      ( getDressing[mEE = #1, recBeast[mEE] = Select[parents,graphHashCode[#] ~SameQ~ graphHashCode[mEE]&], 
                            DEFAULTDRESS]) /. dressed[{a_, b_, c_, d_}] :> 
                           buildGraphDenom[a]*c*d /. NoDressing[a___] :> {} /. 
                        uLsq[a__] :> uLsq[Flatten[a /. Plus :> List]] /.evaluateLsqRule[mom]/.evaluateLpRule[mom]& ) /@ 
        (RHSPARENTS = treesToLoops[trees])])

getTheCut[cut_] :=
    Module[ {},
        Clear[mom];
        {mom} = buildAllMomenta[cut,1];
        \[Xi] = getRandomNullVector[64,4];
        mom[-flat[a_]] :=
            -mom[flat[a]]; (*
        mom[flat[a_]]:=(mom/@a)- \[Xi] Lorentz[(mom/@a),(mom/@a)]/(2 Lorentz[(mom/@a),\[Xi]]); *)
        mom[flat[a_]] :=
            If[ Head[a] ~SameQ~ Plus,
                (mom/@a)- \[Xi] Lorentz[(mom/@a),(mom/@a)]/(2 Lorentz[(mom/@a),\[Xi]]),
                (mom[a])- \[Xi] Lorentz[(mom[a]),(mom[a])]/(2 Lorentz[(mom[a]),\[Xi]])
            ];
        cutGraphSoln = buildRHSNumberFromMom[cut,all,mom];
        candStuff = (collapseTrees[consistentGraphToTrees/@treesToLoops[{#1}]]&)/@cut;
        cutPrecursorList = getCutPrecursor[cut];
        StylePrint[Timing[goodCuts = selectCutList[cutPrecursorList,extRule];]];
        processedCut = Plus@@processCutList[goodCuts];
        refreshHLP;
        StylePrint[cutVal = (processedCut//FermionArrowSign)/.complexSpinorProd[mom]/.evaluateLsqRule[mom]];
        StylePrint[ansatzVal = cutGraphSoln * I uLsq[{leg[k, 1],leg[k, 2]}] uLsq[{leg[k, 2],leg[k, 3]}] Atree[{leg[k, 1],leg[k, 2],leg[k, 3],leg[k, 4]},{leg[k, 1],leg[k, 2],leg[k, 3],leg[k, 4]}/.extRule]//.treeRules//.complexSpinorProd[mom]/.evaluateSRule[mom]/.evaluateLsqRule[mom]];
        If[ NumberQ[cutVal],
            {cutVal,ansatzVal},
            cutVal = (Map[(Numerator[#]/. complexSpinorProd[mom]/. evaluateLsqRule[mom])/Denominator[#]&,FermionArrowSign[processedCut]]//Rationalize)/. complexSpinorProd[mom]/. evaluateLsqRule[mom];
            StylePrint[{cutVal,ansatzVal}];
            {cutVal,ansatzVal}
        ]
    ]

getTheCutOld[cut_] :=
    Module[ {},
        Clear[mom];
        {mom} = buildAllMomenta[cut,1];
        \[Xi] = getRandomNullVector[64,4];
        mom[-flat[a_]] :=
            -mom[flat[a]];
        mom[flat[a_]] :=
            (mom/@a)- \[Xi] Lorentz[(mom/@a),(mom/@a)]/(2 Lorentz[(mom/@a),\[Xi]]);
        cutGraphSoln = buildRHSNumberFromMom[cut,all,mom];
        candStuff = (collapseTrees[consistentGraphToTrees/@treesToLoops[{#1}]]&)/@cut;
        cutPrecursorList = getCutPrecursor[cut];
        StylePrint[Timing[goodCuts = selectCutList[cutPrecursorList,extRule];]];
        processedCut = Plus@@processCutList[goodCuts];
        refreshHLP;
        StylePrint[cutVal = (processedCut//FermionArrowSign)/.complexSpinorProd[mom]/.evaluateLsqRule[mom]];
        StylePrint[ansatzVal = cutGraphSoln * I uLsq[{leg[k, 1],leg[k, 2]}] uLsq[{leg[k, 2],leg[k, 3]}] Atree[{leg[k, 1],leg[k, 2],leg[k, 3],leg[k, 4]},{leg[k, 1],leg[k, 2],leg[k, 3],leg[k, 4]}/.extRule]//.treeRules//.complexSpinorProd[mom]/.evaluateSRule[mom]/.evaluateLsqRule[mom]];
        If[ NumberQ[cutVal],
            {cutVal,ansatzVal},
            cutVal = (Map[(Numerator[#]/. complexSpinorProd[mom]/. evaluateLsqRule[mom])/Denominator[#]&,FermionArrowSign[processedCut]]//Rationalize)/. complexSpinorProd[mom]/. evaluateLsqRule[mom];
            {cutVal,ansatzVal}
        ]
    ]

Clear[findAllMhV];
findAllMhV[tree_] :=
    findAllMhV[tree] = Module[ {nextSteps,runHash,legs = Select[Flatten[tree/.Atree[a__]:>a/.-a_:>a]//Union,#[[1]] ~UnsameQ~ k&]},
                           If[ Max[Length[#[[1]]]&/@tree]<=5 &&Min[Length[#[[1]]]&/@tree]>=4,
                               Return[validMHV[tree]],
                               Map[(runHash[#] = Plus@@((-1+Length[tree[[ (#[[1]][[1]][[1]]) ]][[1]] ])&/@Position[tree,#]))&,legs];
                               nextSteps = Select[Map[Catch[mergeLegsTreeLevel[tree,#]]&,
                                                                                   Select[legs,5>=runHash[#]>=4&]],!StringQ[#]&];
                               If[ nextSteps ~UnsameQ~ {},
                                   Return[Sort[Union[Flatten[Map[findAllMhV,nextSteps]]],OrderedQ[{Length[#1[[1]]],Length[#2[[1]]]}]&]];,
                                   {}
                               ]
                           ]
                       ]


Clear[findAllColl];
findAllColl[tree_,max_] :=
    findAllColl[tree,max] = Module[ {nextSteps,runHash,legs = Select[Flatten[tree/.Atree[a__]:>a/.-a_:>a]//Union,#[[1]] ~UnsameQ~ k&]},
                                If[ Max[Length[#[[1]]]&/@tree]<=max &&Min[Length[#[[1]]]&/@tree]>=4,
                                    Return[validCol[tree]],
                                    Map[(runHash[#] = Plus@@((-1+Length[tree[[ (#[[1]][[1]][[1]]) ]][[1]] ])&/@Position[tree,#]))&,legs];
                                    nextSteps = Select[Map[Catch[mergeLegsTreeLevel[tree,#]]&,
                                                                                        Select[legs,max>=runHash[#]>=4&]],!StringQ[#]&];
                                    If[ nextSteps ~UnsameQ~ {},
                                        Return[Sort[Union[Flatten[Map[findAllColl[#,max]&,nextSteps]]],OrderedQ[{Length[#1[[1]]],Length[#2[[1]]]}]&]];,
                                        {}
                                    ]
                                ]
                            ]


getTheProcessedCut[cut_] :=
    Module[ {},
        candStuff = (
        collapseTrees[consistentGraphToTrees/@treesToLoops[{#1}]]&)/@cut;
        cutPrecursorList = getCutPrecursor[cut];
        StylePrint[Timing[goodCuts = selectCutList[cutPrecursorList,extRule];]];
        processedCut = Plus@@processCutList[goodCuts]
    ]

getTheCutCompare[cut_,cutGraphSoln_] :=
    Module[ {},
        Clear[mom];
        {mom} = buildAllMomenta[cut,1];
        \[Xi] = getRandomNullVector[64,4];
        mom[-flat[a_]] :=
            -mom[flat[a]]; (*
        mom[flat[a_]]:=(mom/@a)- \[Xi] Lorentz[(mom/@a),(mom/@a)]/(2 Lorentz[(mom/@a),\[Xi]]); *)
        mom[flat[a_]] :=
            If[ Head[a] ~SameQ~ Plus,
                (mom/@a)- \[Xi] Lorentz[(mom/@a),(mom/@a)]/(2 Lorentz[(mom/@a),\[Xi]]),
                (mom[a])- \[Xi] Lorentz[(mom[a]),(mom[a])]/(2 Lorentz[(mom[a]),\[Xi]])
            ];
        candStuff = (collapseTrees[consistentGraphToTrees/@treesToLoops[{#1}]]&)/@cut;
        cutPrecursorList = getCutPrecursor[cut];
        StylePrint[Timing[goodCuts = selectCutList[cutPrecursorList,extRule];]];
        processedCut = Plus@@processCutList[goodCuts];
        refreshHLP;
        StylePrint[cutVal = (processedCut//FermionArrowSign)/( uLsq[{leg[k, 1],leg[k, 2]}] uLsq[{leg[k, 2],leg[k, 3]}]  Atree[{leg[k, 1],leg[k, 2],leg[k, 3],leg[k, 4]},{leg[k, 1],leg[k, 2],leg[k, 3],leg[k, 4]}/.extRule])/.complexSpinorProd[mom]/.evaluateLsqRule[mom]];
        StylePrint[ansatzVal = cutGraphSoln * I (1/(uLsq[{leg[k, 1],leg[k, 2]}]uLsq[{leg[k, 1],leg[k, 3]}]uLsq[{leg[k, 2],leg[k, 3]}]))^2//.treeRules//.complexSpinorProd[mom]/.evaluateSRule[mom]/.evaluateLsqRule[mom]/.evaluateLpRule[mom]];
        If[ NumberQ[cutVal],
            {cutVal,ansatzVal},
            cutVal = (Map[(Numerator[#]/. complexSpinorProd[mom]/. evaluateLsqRule[mom])/Denominator[#]&,FermionArrowSign[processedCut]]//Rationalize)/. complexSpinorProd[mom]/. evaluateLsqRule[mom];
            StylePrint[{cutVal,ansatzVal}];
            {cutVal,ansatzVal}
        ]
    ]



(* ::Subsubsection:: *)
(*Jacobi graph op code...*)


(* tHat[graph_, leg_] := 
 Module[{trees = consistentGraphToTrees[graph], a, b, rest, neg, 
        pos, counter, as, at, au, k1, k2, k3, k4, ls, myEqnFs, 
   myNormCrap, zat, i, StylePrint,
        allDaFudges, myTree, rules, fancyRules, Print, allFunction}, 
      StylePrint[{"graph,leg", graph, leg}]; 
       {a, b} = Select[trees, Count[#1, leg, Infinity] > 0 & ]; 
       rest = Complement[trees, {a, b}]; {pos, neg} = 
   If[Count[a, -leg, Infinity] > 0, 
           {b, a}, {a, b}]; counter = 0; While[pos[[1]][[-1]]  ~UnsameQ~  leg, 
         pos = Atree[rotateList[pos[[1]], 2]]]; counter = 0; 
  While[neg[[1]][[1]]  ~UnsameQ~  -leg, 
         neg = Atree[rotateList[neg[[1]], 2]]]; {k1, k2, a} = 
   pos[[1]]; 
       {a, k3, k4} = neg[[1]]; 
  myTree = Join[rest, {Atree[{k1, k2, k3, k4}]}]; 
       StylePrint[myTree]; rules = mapToAbsGenericTrees[myTree]; 
       StylePrint[{"rules", "InputForm"[rules]}]; 
       toGraph[
   Join[rest, {Atree[{k4, k1, at}], Atree[{-at, k2, k3}]}] /. 
    at :> leg]]

uHat[graph_, leg_] :=
    Module[ {trees = consistentGraphToTrees[graph], a, b, rest, neg, 
    pos, counter, as, at, au, k1, k2, k3, k4, ls, myEqnFs, myNormCrap, zat, i, 
    allDaFudges, myTree, rules, fancyRules, Print,StylePrint },
        Clear[fancyRules, SS, TT, UU, 
          nGraphs, theNUMS, myTree, rules];
        StylePrint[{"graph,leg", graph, leg}];
        {a, b} = Select[trees, Count[#1, leg, Infinity] > 0 & ];
        StylePrint[{"ab",a,b}];
        rest = Complement[trees, {a, b}];
        {pos, neg} = If[ Count[a, -leg, Infinity] > 0,
                         {b, a},
                         {a, b}
                     ];
        StylePrint[{"+-",pos,neg}];
        counter = 0;
        While[pos[[1]][[-1]]  ~UnsameQ~  leg, 
        pos = Atree[rotateList[pos[[1]], 2]]];
        counter = 0;
        While[neg[[1]][[1]]  ~UnsameQ~  -leg, 
        neg = Atree[rotateList[neg[[1]], 2]]];
        StylePrint[{"+-",pos,neg}];
        counter = 0;
        {k1, k2, a} = pos[[1]];
        {a, k3, k4} = neg[[1]];
        myTree = Join[rest, {Atree[{k1, k2, k3, k4}]}];
        StylePrint[myTree];
        rules = mapToAbsGenericTrees[myTree];
        StylePrint[{"rules", "InputForm"[rules]}];
        toGraph[Join[rest, {Atree[{k3, k1, au}], Atree[{-au, k4, k2}]}] /. au :> leg]
    ]
*)

jacobiTripleGraph[graph_, leg_] := 
 Module[{ a, b, rest, neg, 
        pos, counter, as, at, au, k1, k2, k3, k4, ls, myEqnFs, 
   myNormCrap, zat, i, StylePrint,
        allDaFudges, myTree, rules, fancyRules, Print, allFunction}, 
      StylePrint[{"graph,leg", graph, leg}]; 
       {neg, pos} =graph[[Sequence@@#[[1;;2]]]]&/@Sort[
Position[graph,leg],
OrderedQ[Length/@{#2,#1}]&]; 
       rest = Complement[graph[[1]], {neg,pos}]; 
 While[pos[[1]][[-1]]  ~UnsameQ~  leg, 
         pos = neckl[rotateList[pos[[1]], 2]]]; counter = 0; 
  While[neg[[1]][[1]]  ~UnsameQ~  -leg, 
         neg = neckl[rotateList[neg[[1]], 2]]]; {k1, k2, a} = 
   pos[[1]]; 
       {a, k3, k4} = neg[[1]]; (* return shat, that, uhat *)
      { graph,vertexFormGraph[
   Join[rest, {neckl[{k4, k1, leg}], neckl[{-leg, k2, k3}]}] ],
vertexFormGraph[
   Join[rest, {neckl[{k3, k1, leg}], neckl[{-leg, k4, k2}]}] ]}]
   
tHat[graph_,leg_]:=jacobiTripleGraph[graph,leg][[2]]
uHat[graph_,leg_]:=jacobiTripleGraph[graph,leg][[2]]
   


(* ::Text:: *)
(**)


(* ::Subsubsection:: *)
(*Tree operations*)


(* so 4pt oriented. ?? getDotRules[trees_] :=
    Module[ {legs = 
       Flatten[trees /. Atree[a__] :> a /. -a_ :> a] // Union, vars, 
      expr, Subscript},
        Subscript[a_, b_] :=
            a[b];
        vars = Union[
          Flatten[Outer[dot[#1, #2] &, 
            Complement[
             legs, {Subscript[k, 1], Subscript[k, 2], Subscript[k, 3], 
              Subscript[k, 4]}], legs]]];
        expr = Flatten[{dot[Subscript[k, 1], Subscript[k, 1]] == 0, 
           dot[Subscript[k, 2], Subscript[k, 2]] == 0, 
           dot[Subscript[k, 3], Subscript[k, 3]] == 0, 
           dot[Subscript[k, 4], Subscript[k, 4]] == 0,
           trees /. 
            Atree[a__] :> 
             Map[(0 == dot[#, Plus @@ a] /. 
                  dot[Subscript[k, z_], Subscript[k, z_]] :> 0 /. 
                 dot[0, z_] :> 0) &, legs]}];
        Rule @@ # & /@ (List @@ Reduce[expr, vars])
    ]  *)
    
getDotRules[trees_] := Module[{legs = Union[Flatten[trees /. Atree[a__] :> a /. -(a_) :> a]], vars, expr,extLegs, Subscript}, 
    Subscript[a_, b_] := a[b]; extLegs=Select[legs,Head[#]===k&];vars = Union[Flatten[Outer[dot[#1, #2] & , Flatten[legs/.k[a_]:>{}], legs]]]; 
     expr = Flatten[{Map[dot[#,#]==0&,extLegs], 
        trees /. Atree[a__] :> (0 == dot[#1, Plus @@ a] /. 
        dot[Subscript[k, z_], Subscript[k, z_]] :> 0 /. 
             dot[0, z_] :> 0 & ) /@ legs}]; (Rule @@ #1 & ) /@ List @@ 
             Reduce[expr, vars]]

getIndepLegs[trees_] :=
    Module[ {StylePrint, 
    rules = getKRules[trees], stuff, 
    intLegs = getIntLegsFromTrees[trees]},
        StylePrint[{"indpe legs krules", rules}];
        stuff = Union[Complement[intLegs, 
              (#1[[1]] & ) /@ rules], 
            Select[Flatten[(#1[[2]] /. Plus :> List /. 
                       -(a_) :> a & ) /@ rules], #1  ~UnsameQ~  0 & ]];
        If[ Length[stuff] < NUMINDEPLEGS,
            Throw[{trees, stuff}]
        ];
        stuff
    ]

Clear[getMyUniqLegs];
getMyUniqLegs[graph_] :=
    getMyUniqLegs[graph] = getIndepLegs[consistentGraphToTrees[graph]]

getColorlyRules[num_, listA_, allFunction_] :=
    Join[getKRules[consistentGraphToTrees[allFunction[num]]] /. 
       Thread[Rule[getMyUniqLegs[ allFunction[num]], listA]], 
      Thread[Rule[getMyUniqLegs[ allFunction[num]], listA]]] // Union


getColorlyDotRules[num_, list_, allFunction_] :=
    (Rule @@ #1 & ) /@ 
    List @@ 
    Reduce[(Equal @@ #1 & ) /@ getDotRules[consistentGraphToTrees[
    allFunction[num]]] /. 
    getColorlyRules[num, list, allFunction]]

getKRulesNum[tree_, list_] :=
    ToRules[Reduce[tree /. Atree[a__] :> 0 == Plus @@ a, 
           Complement[Flatten[tree /. Atree[a__] :> a /. -(a_) :> a], 
       list]]]

colorToStuff[num_, listA_, allFunction_] := 
   daDress[num] /. getKRulesNum[consistentGraphToTrees[allFunction[num]], 
       getMyUniqLegs[allFunction[num]]] /. Thread[getMyUniqLegs[allFunction[num]] -> 
       listA] /. getColorlyDotRules[num, listA, allFunction];


someLabels={a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q,
              r, s, t, u, v, w, x, y, 
             z};


someAddlLabels=Outer[ToExpression[StringJoin[ToString[#1],ToString[#2]]]&,
someLabels,
someLabels]//Flatten//Union;


daLabelList=Flatten[{someLabels,someAddlLabels}];


indepLegLabels[trees_] :=
    Module[ {internalIndep, externalIndep, indep, Subscript},
        indep = getIndepLegs[trees];
        Subscript[a_, b_] :=
            a[b];
        internalIndep = 
         Subscript[
            l, #] & /@ (daLabelList[[Range[
              Length[Intersection[getIntLegsFromTrees[trees], indep]]]]]);
        externalIndep =
         intIndep = 
          Subscript[k, #] & /@ (Range[
             Length[Intersection[getExtLegsFromTrees[trees], indep]]]);
        Join[externalIndep, internalIndep]
    ];


mapToAbsGeneric[graph_] :=
    Module[ {trees = consistentGraphToTrees[graph]},
        mapToAbsGenericTrees[trees]
    ]


mapToAbsGenericTrees[trees_] :=
    Module[ {},
        Join[(getKRules[
            trees] /.
            (Thread[
             Rule[getIndepLegs[trees], indepLegLabels[trees]]])), (Thread[
           Rule[getIndepLegs[trees], indepLegLabels[trees]]])]
    ]


DEFAULTLABELSET[LEGS_,LOOPS_] :=
    Module[ {Subscript},
        Subscript[a_, b_] :=
            a[b];
        indepLegLabels[{Atree[Join[{
             Subscript[k, 1], Subscript[k, 2]}, 
            Subscript[l, #] & /@ Range[LOOPS + 1]]], 
          Atree[Join[-Subscript[l, #] & /@ Range[LOOPS + 1], 
            Subscript[k, #] & /@ Range[3, LEGS]]]}]
    ]



internalOneLoopTadpoleQ[g_] :=
    Module[ {gpf = 
       consistentGraphToTrees[
        g]},
        ! (And @@ ((Length[# /. Atree[a__] :> a /. -a_ :> a]  ~SameQ~  
        Length[Union[# /. Atree[a__] :> a /. -a_ :> a]]) & /@ gpf))
    ]

matchGraphOrAdd[nG_, rule_, metaHolder_, badMaxSYMCheck_] :=
    If[ badMaxSYMCheck[nG],
        0,
        Module[ {all = metaHolder["graphSet"], 
          toWorry = metaHolder["worryList"], num, StylePrint},
            num = Select[Range[Length[all]], isIsomorphic[corruptGraph[all[[#]]],corruptGraph[ nG]] &, 1];
            StylePrint[{"Here we are now", num, InputForm[nG], InputForm[all]
              }];
            If[ num  ~SameQ~  {},
                num = Length[all] + 1;
                metaHolder["graphSet"] = Append[all, nG /. rule];
                metaHolder["worryList"] = Append[toWorry, num];
                num,
                num[[1]]
            ]
        ]
    ]



scoreRulesForLegs[rules_] :=
    Module[ {Subscript},
        Subscript[a_, b_] :=
            a[b];
        Total[Map[If[ #  ~SameQ~  (# /. rules),
                      0,
                      100
                  ] &, 
          Subscript[k, #] & /@ Range[LEGS]]]
    ];

isomorphicRulesAllWithSig[graphA_, graphB_] :=
    Module[ { rules = isomorphicEdgeRulesAll[corruptGraph[graphA],corruptGraph[ graphB]], nl},
        Table[{rule, newGraphSig[graphA /. rule, graphB]}, {rule, rules}]
    ]

getMatchingGraphsOne[graphA_, graphFunc_, 
  num_] :=
    (If[ Flatten[{num}][[1]]  ~UnsameQ~  0,
         Module[ {extLabelRules, indepLegs, StylePrint},
             StylePrint["About to extLabelRules"];
             extLabelRules = 
              isomorphicRulesAllWithSig[graphFunc[num], graphA];
             StylePrint["just ext labeled rules"];
             StylePrint[{"ext Label rules", extLabelRules}];
             extLabelRules = 
              Sort[extLabelRules, OrderedQ[scoreRulesForLegs /@ {#1[[1]], #2[[1]]}] &];
             preUnionA = ((StylePrint[{num, #1[[1]], #1[[2]]}];
                          {getMyUniqLegs[graphFunc[num]] /. #1[[1]], #1[[2]]}) &) /@ 
               extLabelRules[[{1}]];
              preUnion=preUnionA/.graphAlgRules[graphA]/.
         red[a_,1,b_]:>-k[a]/.red[a_,0,b_]:>l[a];
             StylePrint[{"preUnion", preUnion}];
             extLabelRules = 
              Union[preUnion];
             (color[{num}][#1[[1]]]*#1[[2]] &) /@ 
             extLabelRules
         ],
         {0}
     ])

toGraph[trees_] :=
    toGraph[trees] = vertexFormGraph[trees /. Atree[a__] :> neckl[a]]

fromCompressedGraph[graph_]:=Module[{extLegNum=Select[graph,Length[#]===1&]//Length},
If[extLegNum===0, Select[Tally[Abs[Flatten[graph]]],#[[2]]===1&]//Length];
tree[ a/.{n_:>-l[-n]/;
NumberQ[n]&&n<0,
n_:>k[n]/;NumberQ[n]&&0<n<=extLegNum,
n_:>l[n]/;NumberQ[n]&&n>extLegNum}]//toGraph]

fromCompressedDressing[dressing_]:=Throw["Not Implemented Yet"];

symCount[graph_]:=Length[findAllCorruptEdgeIso[graph]]

zeReals[trees_] := (StylePrint["zeReals is deprecated, use
toGraph"];
    toGraph[trees])
    


fromCompressedDressing[dressing_,extLegNum_]:=dressing /.t[a_,b_]:>2dot@@(k/@{a,b})/.k[a_]:>If[a>extLegNum,l[a],k[a]]


fromCompressedGraph[graph_]:=Module[{extLegNum=Select[graph,Length[#]===1&]//Length},
If[extLegNum===0,extLegNum= Select[Tally[Abs[Flatten[graph]]],#[[2]]===1&]//Length];
Atree/@graph/. Atree[a__]:>Atree[ a/.{n_:>-l[-n]/;
NumberQ[n]&&n<0,
n_:>k[n]/;NumberQ[n]&&0<n<=extLegNum,
n_:>l[n]/;NumberQ[n]&&n>extLegNum}]//toGraph]


importCompressedGraphsAndDressings[graphNameFunction_,dressingNameFunction_]:=Module[{graphSet,dressingFunction,list=Range[1,functionKeys[graphNameFunction]//Length],extLegNum},
Print[list//Length," graphs"];
graphSet=corruptGraph@*fromCompressedGraph@*graphNameFunction/@list;
extLegNum=getExtLegs[graphSet[[1]]]//Length;
(dressingFunction[graphSet[[#]]]=fromCompressedDressing[dressingNameFunction[#],extLegNum])&/@list;
{graphSet,dressingFunction}]


planarQ[graph_]:=planarQ[graph]=PlanarGraphQ[
	mathematicaGraph[
		vertexFormGraph[Append[graph[[1]],
		neckl[-getExtLegs[graph]]]]]];




jacobiGraphOnLeg[graph_, leg_, metaHolder_, excludeGraphFunc_] :=
    Module[ {trees = consistentGraphToTrees[graph], a, b, rest, neg, pos, 
      counter, as, at, au, k1, k2, 
           k3, k4, ls, myEqnFs, myNormCrap, zat, i, allDaFudges, myTree, 
      rules, fancyRules, Print, allFunction, StylePrint},
        allFunction[num_] :=
            metaHolder["graphSet"][[num]];
        Clear[fancyRules, SS, TT, UU, nGraphs, theNUMS, myTree, rules];
        StylePrint[{"graph,leg", graph, leg}];
        {a, b} = Select[trees, Count[#1, leg, Infinity] > 0 & ];
        rest = Complement[trees, {a, b}];
        {pos, neg} = If[ Count[a, -leg, Infinity] > 0,
                         {b, a},
                         {a, b}
                     ];
        counter = 0;
        While[pos[[1]][[-1]]  ~UnsameQ~  leg, 
        pos = Atree[rotateList[pos[[1]], 2]]];
        counter = 0;
        While[neg[[1]][[1]]  ~UnsameQ~  -leg, neg = Atree[rotateList[neg[[1]], 2]]];
        {k1, k2, a} = pos[[1]];
        {a, k3, k4} = neg[[1]];
        myTree = Join[rest, {Atree[{k1, k2, k3, k4}]}];
        StylePrint[myTree];
        rules = mapToAbsGenericTrees[myTree];
        StylePrint[{"rules", rules // InputForm}];
        fancyRules = {as -> -(k1 + k2), at -> -(k1 + k4), 
        au -> -(k1 + k3)};
        SS = Join[rest, {Atree[{k1, k2, as}], Atree[{-as, k3, k4}]}];
        TT = Join[rest, {Atree[{k1, k4, at}], Atree[{-at, k3, k2}]}];
        UU = Join[rest, {Atree[{k1, k3, au}], Atree[{-au, k2, k4}]}];
        nGraphs = toGraph /@ {SS, TT, UU};
        StylePrint["Just built toGraph"];
        theNUMS = (matchGraphOrAdd[#1, {as -> leg, at -> leg, 
         au -> leg}, metaHolder, excludeGraphFunc] & ) /@ nGraphs;
        StylePrint[theNUMS];
        StylePrint["Just did theNUMS"];
        theSTUFF = 
         Table[getMatchingGraphsOne[nGraphs[[i]], allFunction, 
           theNUMS[[i]] /. {} -> 0], {i, 1, 3}];
        StylePrint[
        retVal = Flatten[(Outer[#1 == #2 + #3 & , ##1] & ) @@ theSTUFF]];
        StylePrint[retVal /. fancyRules /. rules];
        retVal /. fancyRules /. rules
    ]

Clear[twoFersMap]; 
twoFersMap[eqn_] :=
    twoFersMap[
      eqn] = ({(eqn[[2]] - eqn[[1]] // ExpandAll) /. 
             Plus :> List /. -a_ :> a /. (i_ color[{a_}][b__] :> 
             color[{a}][b]) /. color[{a_}][b__] :> a} // Flatten // Sort)

twoFersEqual[a_, b_] :=
    (twoFersMap[a]  ~SameQ~  twoFersMap[b])

unionSortExprs[jacEqns_] :=
    #[[1]] & /@ 
    Tally[Sort[Complement[(jacEqns // Flatten), {True}], 
    OrderedQ[(StringLength /@ 
    ToString /@ ({#1, #2} /. 
    color[{a_}][{b1_, b2_, b3_, b4___}] :> 
    color[{a}][{b1, b2, b3}])) +
    (1/4) (StringLength /@ ToString /@ ({#1, #2}))] &], 
    twoFersEqual];

specialRules[SPECIALNUM_] :=
    Thread[
    Rule[SPECIALNUM, Table[{}, {Length[SPECIALNUM]}]]];

buildRules[rules_, eqns_, SN_] :=
    Module[ {newRules = rules, 
      newEqns = Flatten[#] & /@ (eqns /. rules /. specialRules[SN]),
      cntr = Length[rules], rle, nextEqns, num},
        While[(nextEqns = Select[newEqns, Length[#]  ~SameQ~  1 &, 1])  ~UnsameQ~  {},
         nextEqns = nextEqns[[1]];
         num = nextEqns[[1]];
         rle = num -> {};
         newEqns = Select[
           Map[Flatten[# /. rle] &, newEqns]
           , #  ~UnsameQ~  {} &];
         cntr++;];
        cntr
    ];

mergeLegsTreeLevelWhoo[trees_, a_] :=
    Module[ {bt = 
    Select[trees, Count[#1, a /. -(b_) :> b, Infinity] > 0 & ], 
        gt = a /. -(b_) :> b},
        If[ Length[bt]  ~UnsameQ~  2,
            If[ Length[bt]  ~SameQ~  1,
                trees /. {-a -> {}, a -> {}} /. Atree[b_] :> Atree[Flatten[b]],
                Throw[{"Bad collapse ", trees, a, bt}]
            ],
            rules = {bt[[2]] -> {}, bt[[1]] -> (bt[[1]] /. -gt :> gt /. 
                     
            gt :> Rest[
              rotateList[bt[[2]][[1]], 
               Flatten[Position[bt[[2]][[1]], gt]][[
                              1]]]])};
            Flatten[trees /. rules /. Atree[b__] :> Atree[Flatten[b]]]
        ]
    ]

mergeLegsTreeLevelListWhoo[trees_, a__] :=
    Module[ {foo = trees, rest = a},
        While[Length[rest] > 0, 
        foo = mergeLegsTreeLevelWhoo[foo, First[rest]];
        rest = Rest[rest]; ];
        foo
    ]

crappyGraph[graph_] :=
    Module[ {tree = consistentGraphToTrees[graph], intLegs, extLegs, 
      otherLegs, otherOther},
        intLegs = 
         tree /. Subscript[k, a_] :> {} /. -a_ :> a /. v[a_, b_] :> {} /. 
            Atree[a__] :> a // Flatten // Union;
        otherOther = 
         Reap[tree /. v[a_, b_] :> Sow[v[a, b]]][[2]] // Flatten // Union;
        Or @@ Table[
          otherLegs = Complement[intLegs, {ll}];
          res = mergeLegsTreeLevelListWhoo[tree /. -ll :> pp, 
            Flatten[{otherLegs, otherOther}]];
          Length[res] > 1, {ll, intLegs}]
    ]

twoExternalVertices[graph_] :=
    Max[Flatten[
       consistentGraphToTrees[graph] /. -a_ :> a /. 
            Subscript[l, a_] :> 0 /. Subscript[k, a_] :> 1 /. 
          v[a_, 0] :> 1 /. v[a_, b_] :> 0 /. Atree[a__] :> Plus @@ a]] > 1

noCrappyGuysInJacEqn[jacEqn_] :=
    Module[ {guys = 
       Reap[jacEqn /. color[{a_}][b__] :> Sow[a]][[2]] // Flatten // 
     Union},
        Intersection[guys, definitelyGood]  ~SameQ~  guys
    ]

reformatEqnRule[eqn_, num_,LEGS_,LOOPS_] :=
    Module[ {StylePrint},
        RuleDelayed[
           ToExpression[
            StringReplace[ ToString[#[[1]]], 
             Table[ToString[leg] -> 
             	    StringJoin[ToString[leg], "_"],
                   {leg, DEFAULTLABELSET[LEGS,LOOPS] /. {k[b_] :> 
                  ToExpression[StringJoin[ToString /@ {k, b}]]
                 , 
                 l[b_] :> 
                  ToExpression[StringJoin[ToString /@ {l, b}]]}}]]],
           Evaluate[#[[2]]]] &[
         reformatEqn[eqn, num,LEGS,LOOPS] /. 
           k[b_] :> ToExpression[StringJoin[ToString /@ {k, b}]] /. 
          l[b_] :> ToExpression[StringJoin[ToString /@ {l, b}]]]
    ]


reformatEqnRule[eqn_, num_] :=
    Module[ {},
    	StylePrint[{"USING DEFAULT: LEGS/LOOPS",LEGS,LOOPS}];
        RuleDelayed[
           ToExpression[
            StringReplace[ ToString[#[[1]]], 
             Table[ToString[leg] -> 
             	    StringJoin[ToString[leg], "_"],
                   {leg, DEFAULTLABELSET[LEGS,LOOPS] /. {k[b_] :> 
                  ToExpression[StringJoin[ToString /@ {k, b}]]
                 , 
                 l[b_] :> 
                  ToExpression[StringJoin[ToString /@ {l, b}]]}}]]],
           Evaluate[#[[2]]]] &[
         reformatEqn[eqn, num, LEGS,LOOPS] /. 
           k[b_] :> ToExpression[StringJoin[ToString /@ {k, b}]] /. 
          l[b_] :> ToExpression[StringJoin[ToString /@ {l, b}]]]
    ]


reformatEqn[eqn_, num_,LEGS_,LOOPS_] :=
     Module[ {StylePrint,colors = 
           
    Reap[(List @@ eqn) /. color[a_][b_] :> Sow[color[a][b]]][[2]] // 
            Flatten // Union,
         rules, aa, eqns, col, 
         Subscript},
          Subscript[a_, b_] := a[b];
          StylePrint[colors];
          StylePrint[num];
          colors = 
            Select[colors, # /. c_ color[a_][b_] :> color[a][b] /. 
                  color[{a_}][b_] :> a  ~SameQ~  num &, 1];
          StylePrint[{"Got the guy", colors}];
          If[ Length[colors] > 1,
               
   Print[Style[{"Wacky greater colors in reformatEQN", eqn, num}]]
           ];
          col = colors[[1]];
          StylePrint[{"Now a col[[1]]", col[[1]]}];
          
  StylePrint[
   eqns = Thread[col[[1]] == (aa /@ Range[Length[col[[1]]]]) ]];
          rules = (Rule @@ #) & /@ 
              List @@ Reduce[
                  Thread[col[[1]] == (aa /@ Range[Length[col[[1]]]]) ],
                  Flatten[{Map[k[#] &, Range[LEGS - 1]], { 
                      Subscript[l, a], Subscript[l, b], 
                      Subscript[l, c], Subscript[l, d], 
                      Subscript[l, e]
                      }
                     }]];
          StylePrint[{"rules", rules}];
          StylePrint[{eqn, eqn /. rules}];
          StylePrint[{"Reformat Goods",
              theReformatGoods = 
                Reduce[eqn /. rules /. 
                    
       Thread[Rule[(aa /@ Range[Length[col[[1]]]]), 
         Flatten[{Map[k[#] &, Range[LEGS - 1]], { 
                          Subscript[l, a], Subscript[l, b], 
                          Subscript[l, c], Subscript[l, d], 
                          Subscript[l, e]
                          }
                         }][[Range[Length[col[[1]]]]]]]], 
                  col /. rules /. 
                    
       Thread[Rule[(aa /@ Range[Length[col[[1]]]]), 
         Flatten[{Map[k[#] &, Range[LEGS - 1]], { 
                          Subscript[l, a], Subscript[l, b], 
                          Subscript[l, c], Subscript[l, d], 
                          Subscript[l, e]
                          }
                         }][[Range[Length[col[[1]]]]]]]]]}];
          theReformatGoods
      ]


Clear[SPECIALNUM];
iterate[{oldRules_, oldEqns_}] := (StylePrint["DEFAULT LEGS/LOOPS GLOBVAL USED IN iterate"];iterate[{oldRules,oldEqns},LEGS,LOOPS])

iterate[{oldRules_, oldEqns_},LEGS_,LOOPS_] :=
    If[ ! ListQ[SPECIALNUM],
        Throw["Must define SPECIALNUM for iterate"],
        Module[ {nextRules = oldRules, cand, rulez, candidateSearch, nZ, nzU,
           newNumz = {},StylePrint},
            candidateSearch = 
             Select[oldEqns, (Length[Flatten[twoFersMap[#]]] > 2 &&
                 
                 Length[Flatten[twoFersMap[# /. specialRules[SPECIALNUM]]]]  ~SameQ~ 
                   1) &];
            StylePrint["CANDIDATE SEARCH"];
            StylePrint[candidateSearch // Short];
            If[ candidateSearch  ~SameQ~  {},
                {oldRules, oldEqns},
                Do[StylePrint[{"sh eqn", eqn // Short}];
                   StylePrint[{"sh eqn, ruled", 
                     eqn //. nextRules // Expand // Short}];
                   cand = If[ (eqn //. nextRules // Expand)  ~SameQ~  True,
                              True,
                              StylePrint[{"mapped", 
                                 twoFersMap[eqn //. nextRules // Expand]} // Short];
                              Flatten[
                                twoFersMap[eqn //. nextRules // Expand] /. 
                                 specialRules[SPECIALNUM]][[1]]
                          ];
                   StylePrint[{"cand", cand}];
                   If[ cand  ~SameQ~  True || cand  ~SameQ~  True[[1]] || cand  ~SameQ~  {}[[1]],
                       StylePrint["Already got one!"],
                       rulez = reformatEqnRule[eqn /. nextRules, cand,LEGS,LOOPS];
                       StylePrint[rulez // Short];
                       If[ Head[rulez]  ~SameQ~  RuleDelayed,
                           nextRules = Append[nextRules, rulez],
                           Throw[{"bogus dude!!", eqn}]
                       ];
                       StylePrint["updated rulez"];
                   ],
                 {eqn, candidateSearch}];
                StylePrint["Made it out of do loop!"];
                nZ = List @@ (And @@ (oldEqns /. nextRules));
                StylePrint["Made it past nZ!"];
                StylePrint[nZ // Short];
                nzU = 
                 Select[unionSortExprs[nZ], Length[Union[twoFersMap[#]]] > 1 &];
                {nextRules, nzU}
            ]
        ]
    ]



findGoodSet[graph_] :=
    If[ Head[graphHashCode[graph]] ~SameQ~ Graph,
        {},
        Module[ {StylePrint,intLegs = getIntLegs[graph],n = 0,goodSet = {}},
            While[goodSet ~SameQ~ {}&&n<10,
               n++;
               StylePrint["Trying "<>ToString[n]];
               goodSet = Map[{#,priveledgeSomeLegs[graph,Flatten[{#}]],graphHashCode[priveledgeSomeLegs[graph,Flatten[{#}]]]}&, Select[Subsets[intLegs,{n}],Head[graphHashCode[priveledgeSomeLegs[graph,Flatten[{#}]]]] ~SameQ~ Graph&]];
               StylePrint[{"Got ", goodSet}];  
             ];
            goodSet
        ]
    ]


monsterPairingOld[poR_,rule_] :=
    Module[ {StylePrint,strippedRule =
    Select[rule,Head[First[#]] ~UnsameQ~ Times&&#[[1]] ~UnsameQ~ #[[2]]&],expr},
        StylePrint["Monster pairing"];
        StylePrint[{"Pairs of Redundancy", poR}];
        StylePrint[{"Stripped rule",strippedRule}];
        safetyCheckA = strippedRule/.
          Rule[c___]:>Rule@@((blam = {Rule[c],"In",poR/.-a_:>a,
          "the thing",#/.-a_:>a,
          "and it's", Position[poR/.-a_:>a,# /.-a_:>a]};
                              If[ blam[[-1]] ~UnsameQ~ {},
                                  StylePrint[{blam}]
                              ];
                              Position[poR/.-a_:>a,# /.-a_:>a])&/@
          List[c]);
        StylePrint[{"SaftyCheckA", safetyCheckA}];
        safetyCheck = safetyCheckA/.Rule[{},{}]->{}/. Rule[a_,{}]:>BLAMM/.Rule[{},a_]:>BLAMM//Flatten;
        StylePrint[{"SaftyCheck", safetyCheck}];
        If[ Position[safetyCheck,BLAMM] ~UnsameQ~ {},
            {},
            If[ safetyCheck ~SameQ~ {},
                rule,
                expr = Position[poR/.-a_->a,#[[1]]]&/@strippedRule;
                StylePrint[{"expr",expr}];
                StylePrint[{"strippedRule",strippedRule}];
                If[ Flatten[expr] ~SameQ~ {},
                    rule,
                    somePairz = First[Flatten[expr[[#]]]]&/@(
                      goodLy = Select[Range[Length[expr]],expr[[#]] ~UnsameQ~ {}&]);
                    StylePrint[{"Some pairz",somePairz}];
                    StylePrint[{"GoodLy",goodLy}];
                    If[ Length[somePairz//Union] ~SameQ~ 1,
                        StylePrint[{"1] sr@goodly",strippedRule[[#]]&/@goodLy}];
                        sstrippedRule = strippedRule;
                        ppoR = poR;
                        Union[Flatten[{#,#/.Rule[a_,b_]:>-a->-b}]]&[
                         Join[rule,
                             Thread[poR[[somePairz[[1]]]]->
                                strippedRule[[goodLy[[1]],2]]/
                                poR[[somePairz[[1]]]][[1]] *(poR[[somePairz[[1]],1]]/
                                  strippedRule[[goodLy[[1]],1]]) poR[[somePairz[[1]]]]]]],
                        If[ Length[somePairz]>2,
                            If[ OddQ[Length[somePairz]],
                                Throw["Problem -- odd pairing monsterPairing!"],
                                Print[Style["Even higher mult pairing!"]];
                                Print[Style["strippedRule!"]];
                                Throw["FILL IN HIGHER MULT PAIRING!"];
                                Do[vishLookup[ poR[[i,1]] /.-a_:>a] = i,{i,1,Length[poR]}];
                                Thread[poR[[sP[[1]]]]->
                                       strippedRule[[goodLy[[1]],2]]/
                                       poR[[sP[[2]]]][[1]] poR[[somePairz[[2]]]]]
                            ],
                            If[ ((Length/@somePairz)//Union//Length) ~UnsameQ~ 1,
                                Throw["Bad play dude!"];,
                                StylePrint[{"2] sr@goodly",strippedRule[[#]]&/@goodLy}];
                                sstrippedRule = strippedRule;
                                ppoR = poR;
                                Union[Flatten[{#,#/.Rule[a_,b_]:>-a->-b}]]&[
                                 Join[rule,
                                     Thread[poR[[somePairz[[1]]]]->
                                        strippedRule[[goodLy[[1]],2]]/
                                        poR[[somePairz[[2]]]][[1]] *(ppoR[[somePairz[[1]],1]]/sstrippedRule[[goodLy[[1]],1]]) poR[[somePairz[[2]]]]],
                                        Thread[poR[[somePairz[[2]]]]->
                                        strippedRule[[goodLy[[2]],2]]/
                                        poR[[somePairz[[1]]]][[1]] (ppoR[[somePairz[[2]],1]]/sstrippedRule[[goodLy[[2]],1]]) poR[[somePairz[[1]]]]]]]
                            ]
                        ]
                    ]
                ]
            ]
        ]  //Sort
    ]


monsterPairing[poR_,rule_] :=
    Module[ {StylePrint,
    strippedRule = Select[rule,Head[First[#]] ~UnsameQ~ Times&&#[[1]] ~UnsameQ~ #[[2]]&],
    safetyCheck,safetyCheckA,porSwitchingLooksLike,
    blowUpStrippedRule,relStrippedRule,nullIt,doob,daab},
(* stripped Rule is the goodIso rule taking poR someplace else,
either to it's own negative or to a friend somewhere else *)
        StylePrint["Monster pairing"];
        StylePrint[{"Pairs of Redundancy", poR}];
        StylePrint[{"Stripped rule",strippedRule}];
        (* SafetyCheckA is a mapping from |poR|\[Rule]pos{|poR|,|poR|}*)
        safetyCheckA = strippedRule/.
          Rule[c___]:>Rule@@((blam = {Rule[c],"In",poR/.-a_:>a,
          "the thing",#/.-a_:>a,
          "and it's", Position[poR/.-a_:>a,# /.-a_:>a]};
                              If[ blam[[-1]] ~UnsameQ~ {},
                                  StylePrint[{blam}]
                              ];
                              Position[poR/.-a_:>a,# /.-a_:>a])&/@
          List[c]);
        StylePrint[{"SaftyCheckA", safetyCheckA}];
        safetyCheck = safetyCheckA/.Rule[{},{}]->{}/. 
           Rule[a_,{}]:>BLAMM/.Rule[{},a_]:>BLAMM//Flatten;
        StylePrint[{"SaftyCheck", safetyCheck}];
        porSwitchingLooksLike = safetyCheck/.Rule[a__,b__]:>Rule[
            First[a//Flatten],First[b//Flatten]];
        If[ Position[safetyCheck,BLAMM] ~UnsameQ~ {},
            {},
            If[ safetyCheck ~SameQ~ {},
                StylePrint["Empty safety"];
                blowUpStrippedRule = Map[Thread[Rule[#,#]]&,pairsOfRedundancy];
                Flatten[
                 Join[blowUpStrippedRule,
                      blowUpStrippedRule/.Rule[a_,b_]:>(-a -> -b),
                      rule]]//Union,
                (*Relevant Stripped Rule *)
                relStrippedRule = Select[strippedRule,
                   Position[poR/.-a_:>a,#[[1]]/.-a_:>a] ~UnsameQ~ {}&];
                StylePrint[{"relStrippedRule",relStrippedRule}];
                nullIt = Map[#->1&,Flatten[relStrippedRule/.(-a_ :> a)/.
                   Rule:>List]]//Flatten//Union;
             (* note sign correction from Rule[a_,b_]\[RuleDelayed]b/a,
              and sign correction from poR (doob/daab) *)
                blowUpStrippedRule = relStrippedRule/.Rule[a_,b_]:>
                Thread[
                 Rule[
                   (*from*) (doob = poR[[Flatten[
                    Position[
                     poR,
                     a/. -g_ :> g]
                     ][[1]]]]),
                    (*to*) ( (b/a) /.nullIt ) (daab = poR[[ Flatten[
                   Position[poR, b /. -g_ :> g]
                   ][[1]]]]) *((First[doob]/First[daab]) /.nullIt)
                   ]
                   ]//Flatten//Union;
                StylePrint[{"blowUpStrippedRule", blowUpStrippedRule}];
                Flatten[
                Join[blowUpStrippedRule,
                     blowUpStrippedRule/.Rule[a_,b_]:>(-a -> -b),
                     rule]]//Union
            ]
        ]  //Sort
    ]


findAllCorruptEdgeIso[graph_] :=
    Module[ {StylePrint,styleprint,comb,myGoodset, myFirstClust, myRules},
        myStrippedGuy = stripPriveledge[graph];
        mmyGoodset = myGoodset = findGoodSet[myStrippedGuy];
        If[ Flatten[myGoodset]  ~SameQ~  {},
            isomorphicEdgeRulesAll[corruptGraph[graph], corruptGraph[graph]],
            makingTheWorldSafe = First/@mmyGoodset;
            redGuyz = Flatten[ makingTheWorldSafe]//Union;
            redGuyz = -redGuyz;
            While[Min[(Length/@(pairsOfRedundancy = Map[#[[1]]&,#]&/@GatherBy[zzzz =
            Table[{myGuy,styleprint[{"rrrr",rrrr = First/@Select[
            styleprint[{"qqqq",myGuy,qqqq = Join[Position[
            consistentGraphToTrees[myStrippedGuy]
              ,myGuy], Position[
            consistentGraphToTrees[myStrippedGuy]
              ,-myGuy]]}];
            qqqq,Length[#] ~SameQ~ 3&]} ];
                         rrrr},{myGuy,redGuyz}],#[[2]]&])//Union)]<2,
            redGuyz = redGuyz/.Map[#->-#&,Flatten[Select[pairsOfRedundancy,Length[#] ~SameQ~ 1&,1]]];
];
            strippedGraph = myStrippedGuy /.Flatten[Map[ Sort[{-#->{},#->{}}]&,Rest[#]]&/@pairsOfRedundancy]/.neckl[a__]:>neckl[Flatten[a]];
            goodIso = isomorphicEdgeRulesAll[strippedGraph,strippedGraph];
            Clear[mmmyLength];
            Do[mmmyLength[First[pair]] = mmmyLength[-First[pair]] = Length[pair],{pair,pairsOfRedundancy}];
            mmmyLength[_] :=
                0;
            goodIso = Table[Join[rules,rules/.Rule[a_,b_]:>Rule[-a,-b]]//Flatten//Union,{rules,goodIso}];
            preGIso = goodIso;
            notBadIso = Select[preGIso, 
              And@@(#/.Rule[a_,b_]:>Equal[mmmyLength[a],mmmyLength[b]])&];
            moR = monsterPairing[pairsOfRedundancy,#]&/@notBadIso;
            StylePrint[moR];
            goodIso = Select[
            moR,# ~UnsameQ~ {}&];
            stupidIsos = Outer[comb[Sort[{##}//Flatten]]&,##,1]&[Sequence@@Table[Table[Thread[pair->perm],
            {perm,Permutations[pair]}],{pair,pairsOfRedundancy}]]//Flatten//Union;
            stupidIsos = stupidIsos/.comb[rules__]:>comb[Join[rules,rules/.Rule[a_,b_]:>Rule[-a,-b]]//Flatten//Union];
            Map[( Join[#,Map[#/.Rule[a_,b_]:>Rule[-a,-b]&,#]]//Flatten//Union)&,
               Map[Union[Flatten[{# (* Don't get to even if yoy wanna
                 ,Reverse/@# *) }]]&,((Sort/@Flatten[stupidIsos/.comb[rules__]:>Flatten[
            {bbbb = (#/.Rule[a_,b_]:>Rule[a,b/.rules]),
            you = #;
            Select[rules,Position[First/@you,#[[1]]] ~SameQ~ {}&&
            Position[First/@you,-#[[1]]] ~SameQ~ {}&]}] &/@goodIso,1])//Union) ]  ]
        ]
    ]


getEveryIsomorphismBasisProject[graph_,dressing_] :=
    Module[ {indepRules,myIsoRule,myIsoRules},
        myIsoRules = findAllCorruptEdgeIso[graph]/.
         red[a_,1,b_]:>-k[a]/.red[a_,0,b_]:>l[a];
        indepRules = getIndepRules[consistentGraphToTrees[graph]];
        Table[dressing/.indepRules/.myIsoRule//.graphAlgRules[graph]/.
         getKRules[graph//consistentGraphToTrees]/.dot[k[a_],k[a_]]:>0/.
         getIndepRules[{Atree[getExtLegsFromTrees[consistentGraphToTrees[graph]]]}],
           {myIsoRule,myIsoRules}]
    ]


getEveryIsomorphism[graph_,dressing_] :=
    Module[ {indepRules,myIsoRule,myIsoRules},
        myIsoRules = findAllCorruptEdgeIso[graph]/.
         red[a_,1,b_]:>-k[a]/.red[a_,0,b_]:>l[a];
        Table[dressing/.myIsoRule//.graphAlgRules[graph]/.
         getKRules[graph//consistentGraphToTrees]/.dot[k[a_],k[a_]]:>0/.
         getIndepRules[{Atree[getExtLegsFromTrees[consistentGraphToTrees[graph]]]}],
           {myIsoRule,myIsoRules}]
    ]


getGravNumerator[graph_, graphList__, dressingStorage_] :=
    getGravDressing[graph, graphList, dressingStorage] /. 
      dressed[{a_, b_, c_, d_}] :> c /. graphAlgRules[graph]


getAvgGravNumerator[graph_,graphList_,dressingStorage_] :=
    getAllGravNumerators[graph,graphList,dressingStorage]/. graphAlgRules[graph]


getOffShellDots[graph_] :=
    Module[ {trees = consistentGraphToTrees[graph],ext,int = getIntLegs[graph],legs},
        ext = getExtLegsFromTrees[trees];
        legs = Join[int,ext];
        Join[Map[dot[#,#]==0&,ext], Append[trees,Atree[ext]] /.Atree[a__]:>
           Map[dot[#,Plus@@a]==0&,legs]]//Reduce[Flatten[#],
               (daList = Sort[Outer[dot[#1,#2]&,Join[ext[[{-1}]],int],
                  Join[ext,int]]//Flatten//Union,OrderedQ[(Abs[#]/.{k[a_]:>a ,
        l[b_]:>2^b,dot[a_,b_]:>a b})&/@{#1,#2}]&]),Backsubstitution->True]&//ToRules
    ]


getGravGraphContribution[graph_, graphList__, dressingStorage_] :=
    getNumerator[graph, graphList, dressingStorage]^2/
         getGraphDenom[stripPriveledge[graph]]

getDblCopyGravGraphContribution[graph_, graphList__, 
  dressingStorage_] :=
    getNumerator[graph, graphList, dressingStorage]^2/
          getGraphDenom[stripPriveledge[graph]]

getGravGraphContribution[graph_, graphList__, dressingStorage_] :=
    getGravNumerator[graph, graphList, dressingStorage]/
          getGraphDenom[stripPriveledge[graph]]

gravPhase[x_] :=
    I

kltCutWithPhase[cut_] :=
    Module[ {expr, nn, cc, cutList, thisKLTCut, Print},
        expr =
        Table[Module[ {n = Length[cutExpr /. Atree[a__] :> a], first, 
        nm1label, nlabel, StylePrint, cntr = 0, labT, lab, perms, mat,negRules,bbb},
                  first = cutExpr[[1, 1]];
                  nm1label = cutExpr[[1, -2]];
                  nlabel = cutExpr[[1, -1]];
                  StylePrint[{cutExpr, n, first, nm1label, nlabel}];
                  perms = Permutations[cutExpr[[1, 2 ;; -3]]];
                  labT = 
                   atreeTwiddle[{first, #, nlabel, nm1label} // Flatten] & /@ perms;
                  lab = 
                   labT /. (atreeTwiddle[a__] :> 
                      atree[{a[[1 ;; -3]], a[[-1]], a[[-2]]} // Flatten]);
                  labT = labT /. atreeTwiddle[a__] :> atreeTwiddle[Atree[a]];
                  lab = lab /. atree[a__] :> atree[Atree[a]];
                  Print[Style[{labT, lab}]];
                  negRules = MapIndexed[-#1->bbb[#2[[1]]]&, 
                     Reap[perms/.-a_:>Sow[a]]//Last//Flatten//Union];
                  StylePrint[{negRules,perms/.negRules}];
                  mat = Outer[newKLT[#1, #2, first] &, perms/.negRules, perms/.negRules, 1]/.
                     (Reverse/@negRules);
                  Print[Style[mat]];
                  {labT, lab, mat}
              ],
        {cutExpr, cut}];
        hmm = expr;
        zeTwiddle = (Outer[Times, ##] &[Sequence @@ (#[[1]] & /@ expr)] // 
      Flatten) //. {atree[a_] atree[b_] :> atree[Flatten[{a, b}]],
           atreeTwiddle[a_] atreeTwiddle[b_] :> 
            atreeTwiddle[Flatten[{a, b}]]};
        zeNT = (Outer[Times, ##] &[Sequence @@ (#[[2]] & /@ expr)] // 
      Flatten) //. {atree[a_] atree[b_] :> atree[Flatten[{a, b}]],
           atreeTwiddle[a_] atreeTwiddle[b_] :> 
            atreeTwiddle[Flatten[{a, b}]]};
        zeMat = KroneckerProduct[##] &[Sequence @@ (#[[3]] & /@ expr)];
        thisKLTCut["leftMovers"] = 
         zeTwiddle /. atreeTwiddle :> someCut /. 
          someCut[a__] :> 
           someCut[a]*Times @@ (gravPhase /@ (a /. Atree[b__] :> Length[b]));
        thisKLTCut["rightMovers"] = 
         zeNT /. atree :> someCut /. 
          someCut[a__] :> 
           someCut[a]*Times @@ (gravPhase /@ (a /. Atree[b__] :> Length[b]));
        thisKLTCut["kltMatrix"] = zeMat;
        thisKLTCut
    ]

stringList[list__] :=
    StringJoin[ToString[list]]

evaluateKltCut[cut_, ymDressing_, ymGraphs_] :=
    Module[ {rawKLT, left, right, matrix, leftVal, rightVal, 
      neededFunctions, basisDots, myFunction, StylePrint},
        StylePrint["Getting KLT Matrix and color-ordered Cuts"];
        rawKLT = kltCutWithPhase[cut];
        left = rawKLT["leftMovers"];
        right = rawKLT["rightMovers"];
        matrix = rawKLT["kltMatrix"];
        StylePrint[{"I have to do ", Length[left], " + " , Length[right], 
           " cuts."} // stringList];
        StylePrint[{"**************** Lets do the left movers ", 
           DateString[]} // stringList];
        leftVal = 
         left /. someCut[cutt_] :> 
           doAColorOrderedCut[cutt, ymDressing, ymGraphs];
        StylePrint[{"**************** Lets do the right movers ", 
           DateString[]} // stringList];
        rightVal = 
         right /. 
          someCut[cutt_] :> doAColorOrderedCut[cutt, ymDressing, ymGraphs];
        neededFunctions = 
         Map[# -> (# /. uLsq[a__] :> dot[Plus @@ a, Plus @@ a] /. 
               dot[a_, a_] :> 0 /. 
              getIndepRules[
               cut]) &, {leftVal, rightVal, matrix} /. 
                dot[a_, b_] :> Sow[dot[a, b]] /. 
               uLsq[aa_] :> Sow[uLsq[aa]] // Reap // Last // Flatten // 
     Union];
        basisDots = 
        Map[#[[2]] /. dot[a_, b_] :> Sow[dot[a, b]] &, neededFunctions] // 
       Reap // Last // Flatten // Union;
        StylePrint[{"**************** Lets build the package ", 
           DateString[]} // stringList];
        myFunction["leftVal"] = leftVal;
        myFunction["rightVal"] = rightVal;
        myFunction["matrix"] = matrix;
        myFunction["neededFunctions"] = neededFunctions;
        myFunction["evaluate"][rules_] :=
            Block[ {dot, uLsq, myRules},
                myRules = (#[[1]] -> (#[[2]] //. rules)) & /@ neededFunctions;
                myRules /. Rule :> Set;
                leftVal . matrix . rightVal
            ];
        myFunction["dots"] = basisDots;
        myFunction
    ]

doAColorOrderedCut[cut_, ymDressing_, ymGraphs_] :=
    Module[ {allCutGrphs, cutGraphs, cutValue, isGood, StylePrint},
        StylePrint[{"Doing cut: ", cut}];
        (isGood[graphHashCode[#]] = True) & /@ ymGraphs;
        isGood[_] = False;
        StylePrint[{"CutGraphs: ", 
           Timing[allCutGrphs = treesToLoops[cut];]} // stringList];
        StylePrint[{"N=4Only: ", 
           Timing[cutGraphs = 
              Select[allCutGrphs, isGood[graphHashCode[#]] &];]} // 
    stringList];
        StylePrint[{"Expression: ", 
           Timing[cutValue = 
              Total[(getCOGraphContribution[mee = #, 
                   Select[ymGraphs, isIsomorphic[mee, #] &, 1], 
                   ymDressing] & /@ cutGraphs)];]} // stringList];
        cutValue
    ]

doADoubleCopyDressedCut[cut_, ymDressing_, ymGraphs_] :=
    Module[ {allCutGrphs, cutGraphs, cutValue, isGood, StylePrint},
        StylePrint[{"Doing cut: ", cut}];
        (isGood[graphHashCode[#]] = True) & /@ ymGraphs;
        isGood[_] = False;
        StylePrint[{"CutGraphs: ", 
           Timing[allCutGrphs = treesToLoopsUO[cut];]} // stringList];
        StylePrint[{"N=4Only: ", 
           Timing[cutGraphs = 
              Select[allCutGrphs, isGood[graphHashCode[#]] &];]} // 
    stringList];
        StylePrint[{"Expression: ", 
           Timing[cutValue = 
              Total[(getDblCopyGravGraphContribution[mee = #, 
                   Select[ymGraphs, isIsomorphic[mee, #] &, 1], 
                   ymDressing] & /@ cutGraphs)];]} // stringList];
        cutValue
    ]


doGravDressedCutWithContacts[cut_, gravDressing_, gravGraphs_] :=
    Module[ {allCutGrphs, cutGraphs, cutValue, isGood},
        StylePrint[{"Doing cut: ", cut}];
        (isGood[graphHashCode[#]] = True) & /@ gravGraphs;
        isGood[_] = False;
        StylePrint[{"CutGraphs: ", 
           Timing[allCutGrphs = treesToLoopsUOWithContacts[cut];]} // 
    stringList];
        StylePrint[{"N=4Only: ", 
           Timing[CGGG = cutGraphs = 
              Select[ corruptGraph/@ allCutGrphs, isGood[graphHashCode[#]] &];]} // 
    stringList];
        StylePrint[{"Expression: ", 
           Timing[cutValue = 
              Total[(getGravGraphContribution[mee = #, 
                   Select[gravGraphs, isIsomorphic[mee, #] &, 1], 
                   gravDressing] & /@ cutGraphs)];]} // stringList];
        cutValue/.red[a_,1,b_]:>-k[a]/.red[a_,0,b_]:>l[a]
    ]


  getGravDressing[graph_, graphList__, dressingStorage_] :=
      Module[ {myIsoRule, myGraph = Select[graphList, isIsomorphic[#1, graph] & ]},
          If[ Length[myGraph]  ~SameQ~  0,
              Return[NoDressing[graph]],
              myGraph = First[myGraph];
              myIsoRule = isomorphicEdgeRule[myGraph, graph];
              dressed[{graph, myGraph, dressingStorage[myGraph] /. myIsoRule /. 
                 graphAlgRules[graph],dsntMatter}]
          ]
      ]


evaluateKltCutNoLsq[cut_, ymDressing_, ymGraphs_] :=
    Module[ {rawKLT, left, right, matrix, leftVal, rightVal, 
      neededFunctions, basisDots, myFunction, foop, lC, mC, rC, varNames,
       data,StylePrint},
        StylePrint[a_] :=
            Print[Style["              <<KLT>> "<>ToString[a]]];
        StylePrint["Getting KLT Matrix and color-ordered Cuts"];
        rawKLT = kltCutWithPhase[cut];
        left = rawKLT["leftMovers"];
        right = rawKLT["rightMovers"];
        matrix = rawKLT["kltMatrix"];
        StylePrint[{"I have to do ", Length[left], " + " , Length[right], 
           " cuts."} // stringList];
        StylePrint[{"Lets do the left movers ", 
           DateString[]} // stringList];
        leftVal = 
         left /. someCut[cutt_] :> 
           doAColorOrderedCut[cutt, ymDressing, ymGraphs];
        StylePrint[{"Lets do the right movers ", 
           DateString[]} // stringList];
        rightVal = 
         right /. 
          someCut[cutt_] :> doAColorOrderedCut[cutt, ymDressing, ymGraphs];
        StylePrint[{"Lets build the package ", 
           DateString[]} // stringList];
        {leftVal,rightVal,matrix} = Block[ {uLsq,dot},
                                        uLsq[a__] :=
                                            dot[Plus@@a,Plus@@a];
                                        getIndepRules[cut]/.Rule:>Set;
                                        {leftVal,rightVal,matrix}
                                    ];
        neededFunctions = Map[# -> (#  /. 
              dot[a_, a_] :> 0 /. getIndepRules[cut]) &, 
         Reap[{ getIndepRules[cut]} /. dot[a___] :> Sow[dot[a]] /. 
               uLsq[a__] :> Sow[uLsq[a]];] // Last // Flatten // Union];
        basisDots = 
         Map[#[[2]] /. dot[a_, b_] :> Sow[dot[a, b]] &, neededFunctions] // 
       Reap // Last // Flatten // Union;
        StylePrint[{"Start building package ", 
           DateString[]} // stringList];
        myFunction["dots"] = basisDots;
        myFunction["leftVal"] = leftVal;
        myFunction["rightVal"] = rightVal;
        myFunction["matrix"] = matrix;
        myFunction["neededFunctions"] = neededFunctions;
        myFunction["evaluate"][rules_] :=
            Block[ {dot, uLsq, myRules},
                myRules = (#[[1]] -> (#[[2]] /. rules)) & /@ neededFunctions;
                dot[a_,a_] = 0;
                myRules /. Rule :> Set;
                leftVal . matrix . rightVal
            ];
        varNames = 
         ToExpression["XX" <> ToString[#]] & /@ Range[Length[basisDots]];
        myFunction
    ]



evaluateKltCut[cut_, ymDressing_, ymGraphs_] :=
    Module[ {rawKLT, left, right, matrix, leftVal, rightVal, 
      neededFunctions, basisDots, myFunction, foop, lC, mC, rC, varNames,
       data,StylePrint},
        StylePrint["Getting KLT Matrix and color-ordered Cuts"];
        rawKLT = kltCutWithPhase[cut];
        left = rawKLT["leftMovers"];
        right = rawKLT["rightMovers"];
        matrix = rawKLT["kltMatrix"];
        StylePrint[{"I have to do ", Length[left], " + " , Length[right], 
           " cuts."} // stringList];
        StylePrint[{"**************** Lets do the left movers ", 
           DateString[]} // stringList];
        leftVal = 
         left /. someCut[cutt_] :> 
           doAColorOrderedCut[cutt, ymDressing, ymGraphs];
        StylePrint[{"**************** Lets do the right movers ", 
           DateString[]} // stringList];
        rightVal = 
         right /. 
          someCut[cutt_] :> doAColorOrderedCut[cutt, ymDressing, ymGraphs];
        StylePrint[{"**************** Lets build the package ", 
           DateString[]} // stringList];
        neededFunctions = 
         Map[# -> (# /. uLsq[a__] :> dot[Plus @@ a, Plus @@ a] /. 
               dot[a_, a_] :> 0 /. getIndepRules[cut]) &, 
          Reap[{leftVal, rightVal, matrix, getIndepRules[cut]} /. dot[a___] :> Sow[dot[a]] /. 
                uLsq[a__] :> Sow[uLsq[a]];] // Last // Flatten // Union];
        basisDots = 
         Map[#[[2]] /. dot[a_, b_] :> Sow[dot[a, b]] &, neededFunctions] // 
       Reap // Last // Flatten // Union;
        StylePrint[{"**************** Start building package ", 
           DateString[]} // stringList];
        myFunction["dots"] = basisDots;
        myFunction["leftVal"] = leftVal;
        myFunction["rightVal"] = rightVal;
        myFunction["matrix"] = matrix;
        myFunction["neededFunctions"] = neededFunctions;
        myFunction["evaluate"][rules_] :=
            Block[ {dot, uLsq, myRules},
                myRules = (#[[1]] -> (#[[2]] /. rules)) & /@ neededFunctions;
                myRules /. Rule :> Set;
                leftVal . matrix . rightVal
            ];
        varNames = 
         ToExpression["XX" <> ToString[#]] & /@ Range[Length[basisDots]];
        myFunction
    ]





Clear[newKLT];
newKLT[\[Rho]_, \[Tau]_, ki_] :=
    If[ \[Rho]  ~SameQ~  {},
        1,
        Module[ {j = \[Rho][[-1]], l, \[Beta], \[Gamma]},
            l = Flatten[Position[\[Tau], j]][[1]];
            \[Beta] = \[Tau][[1 ;; l - 1]];
            \[Gamma] = \[Tau][[l + 1 ;; -1]];
            2 dot[Plus @@ \[Beta] + ki, j] newKLT[\[Rho][[1 ;; -2]], 
              Join[\[Beta], \[Gamma]], ki]
        ]
    ]

treesToLoopsUOWithContactsSlow[cut_] :=
    Module[ {cubicCut = treesToLoopsUO[cut], grz, StylePrint},
        If[ VERBOSETIMING  ~SameQ~  True,
            StylePrint[a_] :=
                Print[Style[a]]
        ];
        Join[stripPriveledge /@ 
           Union[ Timing[
              grz = priveledgeLegs /@ (Table[
                    Map[toGraph[mergeLegsTreeLevelList[curGraph, #]] &, 
                     Select[Subsets[
                       curGraph /. in[a_] :> Sow[in[a]] // Reap // Last // 
                   Flatten // Union], Length[#] > 0 &]], {curGraph, 
                     consistentGraphToTrees /@ cubicCut}] // Flatten);] // 
       StylePrint;
                  StylePrint[{"About to union ", Length[grz], 
                    "somewhat redundant graphs including contacts"}];
                  Module[ {j = 0},
                      Do[ graphHashCode[grz[[i]]];
                          If[ IntegerQ[i/Length[grz] *10] || IntegerQ[i/50],
                              StylePrint[{(i/Length[grz])*100 // N, 
                                "percent has been graphHashCoded", i, 
                                DateString[]}]
                          ];, {i, 1, Length[grz]}];
                  ];
                  grz, SameTest -> isIsomorphic], cubicCut] // Flatten
    ]

priveledgeLegs[graph_] :=
    Module[ {allKs, allLs, newLs, newKs, 
      tree = consistentGraphToTrees[graph], StylePrint},
        allKs = tree /. k[a_] :> Sow[k[a]] // Reap // Last // Flatten // 
    Union;
        allLs = 
         tree /. l[a_] :> Sow[l[a]] // Reap // Last // Flatten // Union;
        StylePrint[{allKs, allLs}];
        newLs = 
         Table[  Table[
            Atree[{red[num, 0, i], -red[num, 0, i + 1]}], {i, 1, num}] /. 
           red[num, 0, num + 1] -> l[num], {num, allLs /. l[a_] :> a}];
        StylePrint[allLs];
        newKs = 
         Table[  Table[
            Atree[{red[num, 1, i], -red[num, 1, i + 1]}], {i, 1, num}] /. 
           red[num, 1, num + 1] -> -k[num], {num, allKs /. k[a_] :> a}];
        Join[newLs, newKs, 
           tree /. Thread[-allLs -> (allLs /. l[a_] :> -red[a, 0, 1])] /. 
            Thread[allKs -> (allKs /. k[a_] :> -red[a, 1, 1])]] // Flatten //
    toGraph
    ]

blowOutExt[graph_] := 
  Module[{tree = consistentGraphToTrees[graph], ks, myNum},
   Say["Depends oddly on LOOPS, probably replacable by \
priveledgeExtLegs"];
   ks = Reap[tree /. k[ a_] :> Sow[k[ a]] ][[2]] // Flatten // Union;
   Join[tree, Flatten[Table[myNum = myK /. k[a_] :> a;
          
        Table[Atree[ {in[
            Prime[LOOPS + myNum]*2^j], -in[
             Prime[LOOPS + myNum]*2^(j + 1)]}],
         {j, 1, myNum}], {myK, ks}]]] /. 
     Flatten[Table[
       myNum = myK /. 
         k[a_] :> a; {myK -> -in[Prime[LOOPS + myNum]*2],
        -in[Prime[LOOPS + myNum]*2^(myNum + 1)] -> myK}, {myK, ks}]] //
     toGraph];
    
priveledgeExtLegs[graph_] :=
    Module[ {allKs, allLs, newLs, newKs, 
      tree = consistentGraphToTrees[graph], StylePrint},
        allKs = tree /. k[a_] :> Sow[k[a]] // Reap // Last // Flatten // 
    Union;
        allLs = {};
        newLs = 
         Table[  Table[
            Atree[{red[num, 0, i], -red[num, 0, i + 1]}], {i, 1, num}] /. 
           red[num, 0, num + 1] -> l[num], {num, allLs /. l[a_] :> a}];
        StylePrint[allLs];
        newKs = 
         Table[  Table[
            Atree[{red[num, 1, i], -red[num, 1, i + 1]}], {i, 1, num}] /. 
           red[num, 1, num + 1] -> -k[num], {num, allKs /. k[a_] :> a}];
        Join[newLs, newKs, 
           tree /. Thread[-allLs -> (allLs /. l[a_] :> -red[a, 0, 1])] /. 
            Thread[allKs -> (allKs /. k[a_] :> -red[a, 1, 1])]] // Flatten //
    toGraph
    ]

stripPriveledge[graph_] :=
    consistentGraphToTrees[graph] /. 
         Atree[{a_, b_}] :> {} /. -red[a_, 1, 1] :> 
         k[a] /. -red[a_, 0, 1] :> -l[a] // Flatten // toGraph

mapCutWithInsToLs[cut_] :=
    Module[ {maxL = 
       Max[cut /. l[a_] :> Sow[a] // Reap // Last // Flatten // Union] + 
        1},
        cut /. in[a_] :> l[a + maxL]
    ]

dressTreeUOWithContacts[Atree[list__]] :=
    Module[ {graphs = getNUOTreeWithContacts[Length[list]]},
        graphs /. Thread[k /@ Range[Length[list]] -> list]
    ]

doContactCollapse[cubicCut_] :=
    Module[ {grz, StylePrint},
        If[ VERBOSETIMING  ~SameQ~  True,
            StylePrint[a_] :=
                Print[Style[a]]
        ];
        Join[stripPriveledge /@ 
           Union[ Timing[
              grz = priveledgeLegs /@ (Table[
                    Map[toGraph[mergeLegsTreeLevelList[curGraph, #]] &, 
                     Select[Subsets[
                       curGraph /. in[a_] :> Sow[in[a]] // Reap // Last // 
                   Flatten // Union], Length[#] > 0 &]], {curGraph, 
                     consistentGraphToTrees /@ cubicCut}] // Flatten);] // 
       StylePrint;
                  StylePrint[{"About to union ", Length[grz], 
                    "somewhat redundant graphs including contacts"}];
                  Module[ {j = 0},
                      Do[ graphHashCode[grz[[i]]];
                          If[ IntegerQ[i/Length[grz] *10] || IntegerQ[i/50],
                              StylePrint[{(i/Length[grz])*100 // N, 
                                "percent has been graphHashCoded", i, 
                                DateString[]}]
                          ];, {i, 1, Length[grz]}];
                  ];
                  grz, SameTest -> isIsomorphic], cubicCut] // Flatten
    ]

Clear[getNUOTreeWithContacts]; 
getNUOTreeWithContacts[n_] :=
    getNUOTreeWithContacts[n] = Block[ {VERBOSETIMING = False, k},
                                    doContactCollapse[ Nest[invSoftLimitUO[#1, k] & , 
                                          {vertexFormGraph[{neckl[{k[1], k[2], k[3]}]}]}, n - 3]]
                                ]

                    
 (* Henrik expand -- refactor to match other code. *)                   
Clear[getNUOTreeWithContacts]; 
getNUOTreeWithContacts[n_] :=
    getNUOTreeWithContacts[n] = (Module[ {graf = Join[Table[{-i},{i,1,n}],{Range[n]}],StylePrint,
    g1,gg,ng,t,q,tt,r,rr,k,kk,GExp,lm,ll},
                                     GExp[g_,i_] :=
                                         (g1 = g;
                                          gg = g[[i]];
                                          ng = Length[gg];
                                          t = {{q[1,2,3]}};
                                          Do[t = Flatten[
                                          Table[ll = Union[Abs[Flatten[l/.q->List]]];
                                                lm = Max[{ll,ng}]+1;
                                                Join[
                                                Table[Append[l/.{k->-lm},q[k,j,lm]],{k,ll}],
                                                Table[l/.k:>Append[k,j],{k,l}]],
                                          {l,t}],1],
                                          {j,4,ng}];
                                          Table[g1[[i]] = (t[[j]]/.Join[
                                          {List->Sequence,q->List},Thread[Range[ng]->gg],
                                          Thread[Range[ng+1,2 ng-3]->Range[r+1,r+ng-3]],
                                          Thread[-Range[ng+1,2 ng-3]->-Range[r+1,rr = r+ng-3]]]);
                                                g1,
                                          {j,Length[t]}]);
                                     tt = {graf};
                                     r = Max[graf];
                                     Do[StylePrint[{Length[tt],r,i}];
                                        tt = If[ Length[graf[[i]]]>3,
                                                 Flatten[Table[GExp[tt[[j]],i],{j,Length[tt]}],r = rr;
                                                                                               1],
                                                 tt
                                             ];,
                                         {i,Length[graf],1,-1}];
                                     Map[vertexFormGraph[(neckl/@Select[#,Length[#]>1&])]&,tt]
                                 ]/.
        neckl[a__]:>neckl[k/@a]/.k[a_]:>If[ a<0,
                                            -k[-a],
                                            k[a]
                                        ]/.k[a_]:>If[ a>n,
                                                      in[a-n],
                                                      k[a]
                                                  ])


(* ::Subsubsection:: *)
(*Graph Web from Jacobi*)


graphExclusionF[a_]:=internalOneLoopTriangleQ[a]||
internalOneLoopBubbleQ[a]||
internalOneLoopTadpoleQ[a];


graphToWeb[aGraph_,graphExclusionF_, zcutLeg___]:=Module[{web,cutGraphData,
bad,jacEqns,jacEqn,jacEqns2,jacEqns2b,jacEqns2a,numNEWGRAPHS,val,zeNum,zeGraph,zenum,gnL,allFunction,graphExclusion,nextEqn,StylePrint},
bad=jacEqns=jacEqn=jacEqns2=jacEqns2b=jacEqns2a={};
LEGS=getExtLegs[aGraph]//Length;
LOOPS=Length[getMyCycles[aGraph]];
web["LEGS"]=LEGS;
web["LOOPS"]=LOOPS;
StylePrint[web];
cutGraphData["graphSet"]={aGraph} ;
numNEWGRAPHS:=Length[cutGraphData["graphSet"]];
cutGraphData["worryList"]=Range[numNEWGRAPHS];
jacEqns={};
allFunction[i_]:=cutGraphData["graphSet"][[i]];
graphExclusion[a_]:=graphExclusion[a]=graphExclusionF[a];
Timing[While[Length[cutGraphData["worryList"]]>0,zenum=First[cutGraphData["worryList"]];
StylePrint["Doing "<>ToString[zenum]];
cutGraphData["worryList"]=Rest[cutGraphData["worryList"]];
zeGraph=allFunction[zenum];
gnL=Select[getIntLegs[zeGraph],Head[#] ~UnsameQ~ zcutLeg&];
jacEqns=Join[jacEqns,Table[StylePrint[{zenum,leg,val=jacobiGraphOnLeg[zeGraph,leg,cutGraphData,graphExclusion]}];
val,{leg,gnL}]];
Say[{zenum,Length[cutGraphData["worryList"]],Length[jacEqns]}];
web["graphSet"]=cutGraphData["graphSet"];
];
];
web["graphSet"]=cutGraphData["graphSet"];
web["jacEqns"]=Sort[Complement[(jacEqns//Flatten),{True}],OrderedQ[(StringLength/@ToString/@({#1,#2}/.color[{a_}][{b1_,b2_,b3_,b4___}]:>color[{a}][{b1,b2,b3}]))+(1/4) (StringLength/@ToString/@({#1,#2}))]&];
web
]


addMastersToWeb[web_]:=Module[{StylePrint,jacEqns,jacEqns2=web["jacEqns"],jacEqns2a, jacEqns2b,specialRules,buildRules,zeroRules,jacEqns3,bipartiteGraph,bi1,bish,bi2,bish2,ubi2,candSingleMasters,allFunction,
theGuys,others,jEq3,bad={},planarGraphs,newRules={},restEqns,myNumz,aRule,known,nonZero,nonZeroUnion, nr,nzU,allDa,numNEWGRAPHS,outSideRules,
all=web["graphSet"],nextEqn},jacEqns=jacEqns2;numNEWGRAPHS=Length[all];
allFunction[i_]:=all[[i]];
StylePrint[{"Length of jacEqns",Length[jacEqns2]}];
jacEqns2a=Tally[jacEqns2,twoFersEqual];
jacEqns2b=#[[1]]&/@jacEqns2a;
jacEqns2=jacEqns2b;
StylePrint[{"jacEqns2b",Length[jacEqns2b]}];
specialRules[SPECIALNUM_]:=Thread[Rule[SPECIALNUM,Table[{},{Length[SPECIALNUM]}]]];
buildRules[rules_,eqns_,SN_]:=Module[{newEqns=Flatten[#]&/@(eqns/.rules/.specialRules[SN]),cntr=Length[rules],rle,nextEqns,num},While[(nextEqns=Select[newEqns,Length[#] ~SameQ~ 1&,1]) ~UnsameQ~ {},nextEqns=nextEqns[[1]];
num=nextEqns[[1]];
rle=num->{};
newEqns=Select[Map[Flatten[#/.rle]&,newEqns],# ~UnsameQ~ {}&];
cntr++;];
cntr];
zeroRules=reformatEqnRule[color[{#}][DEFAULTLABELSET[LEGS,LOOPS]]==0,#,LEGS,LOOPS]&/@bad;
jacEqns3=(jacEqns2/.zeroRules)//Flatten//Union;
jacEqns3=Union[Select[jacEqns3,# ~UnsameQ~ True&],SameTest->twoFersEqual];
bipartiteGraph=Map[twoFersMap,jacEqns2]//Union;
bi1=Select[bipartiteGraph,(Length[#]>1&&Length[Union[#]]>1)&];
bish=Map[#->{}&,Select[bipartiteGraph,(Length[#] ~SameQ~ 1)&]//Flatten];
bi2=Select[Flatten[#/.bish]&/@bipartiteGraph,(Length[#]>1&&Length[Union[#]]>1)&];
bish2=Map[#->{}&,Select[Flatten[#/.bish]&/@bipartiteGraph,(Length[#] ~SameQ~ 1)&]//Flatten];
ubi2=Union[Flatten[bi2]];
StylePrint[ubi2];
SPECIALNUM=If[Length[all] ~SameQ~ 1,
{1},
candSingleMasters=Select[ubi2,buildRules[bish,bi2,{#}] ~SameQ~ (numNEWGRAPHS-1)&];
If[candSingleMasters ~UnsameQ~ {}, Say[" *** YeS Cand Single Masters***"];
myFoundMaster=
   First[ Sort[candSingleMasters, OrderedQ[ If[planarQ[ web["graphSet"][[#]] ],
      -100*1/#,100*1/#]&/@{#1,#2}]&]];
Say[{"I found this",myFoundMaster}];
{myFoundMaster},
   Say[" *** No Cand Single Masters, trying dbl planar***"];
   planarGraphs=Select[ubi2,planarQ[all[[#]]]&];
   pGC=Flatten[Table[Table[{planarGraphs[[i]],
   planarGraphs[[j]]},
   {j,i+1,Length[planarGraphs]}],
   {i,1,Length[planarGraphs]}],1];
   candidatePlanarPairs=Select[pGC,( buildRules[bish,bi2,#] ~SameQ~ (numNEWGRAPHS-Length[#]) )&];
   If[Length[candidatePlanarPairs]>0, 
      Say[" *** Yes dbl planar***"];
      First[candidatePlanarPairs],
      Say[" *** No dbl planar, try 1 np+1pl ***"];
      pGC=Flatten[Table[Table[{planarGraphs[[i]],
         ubi2[[j]]},{j,1,Length[ubi2]}],{i,1,Length[planarGraphs]}],1];
      candidatePlanarPairs=Select[pGC,
         (buildRules[bish,bi2,#] ~SameQ~ (numNEWGRAPHS-Length[#]))&];
      If[Length[candidatePlanarPairs]>0,
          Say[" *** Yes 1 np+1pl ***"];
           First[candidatePlanarPairs],
          Say["***** go for 3 np!!"];
          pGT=Flatten[Table[Table[Table[{planarGraphs[[i]],planarGraphs[[j]],planarGraphs[[k]]},{k,j+1,Length[planarGraphs]}],
               {j,i+1,Length[planarGraphs]}],{i,1,Length[planarGraphs]}],2];
          candidatePlanarTrips=Select[pGT, ( buildRules[bish,bi2,#] ~SameQ~ 
          (numNEWGRAPHS-Length[#])) &];
          If[Length[candidatePlanarTrips]>0,
             candidatePlanarTrips//First,
             Throw["Fault--nomastersfound"]
             ]
      ]
    ]
   ]
];
StylePrint[{"Candidate Single Masters!",candSingleMasters,"THIS SHOULD BE 1"}];
StylePrint[{"SPECIALNUM RIGHT HERE",SPECIALNUM}];
theGuys=Select[Range[Length[all]],twoExternalVertices[allFunction[#]]&];
others=Complement[Range[numNEWGRAPHS],theGuys];
StylePrint[{"pre crappyGraph defGood",definitelyGood}];
definitelyGood=Select[others,!crappyGraph[allFunction[#]]&];
StylePrint[{"post crappyGraph defGood",definitelyGood}];
jEq3=Select[jacEqns2,noCrappyGuysInJacEqn[#]&];
Print[{"jEq3a",jEq3//Length}];
jEq3=Flatten[jacEqns];
Print[{"jEq3b",jEq3//Length}];
bad=((twoFersMap[#]&/@Select[jEq3,(Length[twoFersMap[#]]) ~SameQ~ 1&]))//Flatten//Union;
Print[{"definitelyGood",definitelyGood}];
Print[{"bad",bad}];
Print[{"jEq3",jEq3//Length}];
LEGS=web["LEGS"];
LOOPS=web["LOOPS"];
zeroRules=reformatEqnRule[color[{#}][DEFAULTLABELSET[LEGS,LOOPS]]==0,#,LEGS,LOOPS]&/@bad;
Print[{"DEFAULTLABELSET",DEFAULTLABELSET[LEGS,LOOPS]}];
Print[{"zeroRules",zeroRules//Length}];
jacEqns3=(jEq3/.zeroRules)//Flatten//Union;
jacEqns3=Union[Select[jacEqns3,# ~UnsameQ~ True&],SameTest->twoFersEqual];
StylePrint[{"nR first",newRules//Length}];
planarGraphs=Select[Range[numNEWGRAPHS],planarQ[all[[#]]]&];
        StylePrint[planarGraphs];
planarGraphs={};

Module[{c=0},known=Join[SPECIALNUM];
newRules=zeroRules;
restEqns=unionSortExprs[Select[jacEqns3,(twoFersMap[#]//Length) ~SameQ~ 2&&(twoFersMap[#]//Union//Length) ~SameQ~ 2&]];
While[Length[restEqns=unionSortExprs[Select[restEqns//.newRules,(# ~UnsameQ~ True&&(twoFersMap[#]//Length) ~SameQ~ 2&&(twoFersMap[#]//Union//Length) ~SameQ~ 2)&]]]>0&&restEqns ~UnsameQ~ {True},
StylePrint[{"indapreloop",newRules//Length,Length[restEqns]}];
StylePrint[{"First eqn",nextEqn=First[restEqns]}]; 
(myNumz=Reap[(nextEqn/.color[{a_}][b__]:>Sow[a])][[2]]/.Map[#->{}&,known ] /.Map[#->{}&, planarGraphs]//Flatten);
 Say[{"myNumz",myNumz}];
StylePrint[{"Da my known",known}];
StylePrint[{"Da my rules",newRules}];
StylePrint[{"Da my planar",planarGraphs}]; 

If[myNumz ~UnsameQ~ {},
 StylePrint[{"whoopdee",c++;myNumz . c,c}];
myNumz=First[Flatten[myNumz/.Map[#->{}&,planarGraphs]]];
StylePrint[{"Decided",myNumz}];
aRule=reformatEqnRule[nextEqn,myNumz,LEGS,LOOPS];
StylePrint[{"aRule",aRule}];
newRules=Append[newRules,aRule]/.True:>{}//Flatten//Union;
known=Join[SPECIALNUM,Map[#[[1]]&,newRules]/.color[{a_}][b__]:>Sow[a]//Reap//Last//Flatten//Union]//Flatten//Union];
Say[{"nR",newRules//Length,restEqns//Length}];]
];
nonZero=List@@(And@@(jacEqns3//.newRules));
nonZeroUnion=Select[unionSortExprs[nonZero],Length[Union[twoFersMap[#]]]>1&];
nr=newRules//Flatten//Sort;
nzU=nonZeroUnion//Flatten//Sort;
StylePrint[{"nR Second",nr}];
StylePrint[{"nzU Second",nzU}];

{nr,nzU}=FixedPoint[Map[Union[Flatten[#]]&,iterate[#,LEGS,LOOPS]]&,{nr,nzU}];
StylePrint[{"nR Third",nr}];
StylePrint[{"nzU Third",nzU}];
Say[{"Did we get it all?",Length[all]-Length[nr],Length[SPECIALNUM]}];
outSideRules=nr;
web["jacSoln"]=nr;
web["masterIds"]=SPECIALNUM;
web
]


(* ::Subsubsection:: *)
(*Delayed Graph Container Code*)


priveledgeSomeLegs[graph_,legs_] :=
    Module[ {allKs, allLs, newLs, newKs, 
    tree = consistentGraphToTrees[graph],numLz,StylePrint,dz},
        allLs = legs;
        Do[numLz[num] = If[Length[dz=Flatten[Position[allLs,l[num]]]]>0,Last[dz],0],{num,allLs/.l[a_]:>a}];
        
        StylePrint[{#,numLz[#]}&/@(allLs/.l[a_]:>a)];
        StylePrint[{allKs, allLs}];
        newLs = Table[ Table[
         Atree[{red[num, 0, i], -red[num, 0, i + 1]}], {i, 1,numLz[num]}] /. 
           red[num, 0, numLz[num] + 1] -> l[num], {num, allLs /. l[a_] :> a}];
        StylePrint[allLs];
        toGraph[Flatten[Join[newLs, 
           tree /. Thread[-allLs -> (allLs /. l[a_] :> -red[a, 0, 1])] ]]]
    ]


minStripPriveledge[graph_] :=
    toGraph[Flatten[consistentGraphToTrees[graph] /. Atree[{red[a___],-red[b___]}] :> {} /. 
    red[a_, 0, b_] :> red[a, 0, 1]]]


Clear[corruptCutList];
corruptCutList[cutList_] :=
    corruptCutList[cutList] = Module[ {n = 0,
    StylePrint,ret},
                                  Print[Style[{"Working out a corrupt list.", 
                                  DateString[]}]];
                                  ret = consistentGraphToTrees/@(corruptGraph[toGraph[#]]&/@
                                     cutList);
                                  Print[Style[{"Worked out a corrupt list.", 
                                  DateString[]}]];
                                  ret
                              ]


getHardCutId[graph_,cutList_] :=
    Module[ {expr = Select[Range[Length[cutList]],isIsomorphic[toGraph[corruptCutList[cutList][[#]]],corruptGraph[graph]]&,1]},
        If[ expr ~SameQ~ {},
            Throw[badCutThrowObject[{"cutlist Incomplete!!",graph//stripPriveledge//consistentGraphToTrees,
                     graph}]],
            If[ Length[expr]>1,
                Throw[{"too many cutlist matches",graph}],
                expr[[1]]
            ]
        ]
    ]


buildADressingContainer[graphs_,graphDressings_] :=
    Module[ {myDressingObject},
        myDressingObject["graphs"] = graphs;
        myDressingObject["delayedContainers"] = {};
        If[ Length[graphs] ~UnsameQ~ Length[graphDressings],
            Throw["buildDressing length mismatch"];
        ];
        Do[myDressingObject["dressings"][graphs[[i]]] = graphDressings[[i]],{i,1,Length[graphs]}];
        myDressingObject["loadDelayedDressings"][configRule_] :=
            (myDressingObject["delayedContainers"] =
            Append[myDressingObject["delayedContainers"],
            Module[ {newContainer, contactDir,
            filenames,goodList,strReplaceRule,theFGPath,
            filename,filenum,str,myGraph,testString,test,myGraphs,contList,myCuts,numZ = {},impTable = {}},
                myGraphs = {};
                myCuts = "canonicalGraphLabels"/.configRule;
                newContainer["configRule"] = configRule;
                newContainer["inclusionQ"][graph_] :=
                    (ToExpression["isIncluded"/. configRule])[graph];
                newContainer["labelGenerator"][graph_] :=
                    (ToExpression["labelGenerator"/. configRule])[graph];
                newContainer["canonicalGraphLabels"] :=
                    (Evaluate["canonicalGraphLabels"/. configRule]);
                contactDir = "graphListPath"/.configRule;
                If[ !FileExistsQ[contactDir],
                    CreateDirectory[contactDir]
                ];
                theFGPath = contactDir<>"*"<>("graphPattern"/.configRule);
                StylePrint[{"theFGPath",theFGPath,DateString[]}];
                impTable = {};
                If[ ("FastLoad"/.configRule) ~UnsameQ~ False,
                    StylePrint["Going speedy here!"];
                    iimptTable = impTable = Import[someStr =
                     "!find "<>contactDir<>"  -type f -name '*"<>("graphPattern"/.configRule)<>"' -exec grep -c True {} + | grep -v :0 | grep -v "<>("dressPattern"/.configRule),
                       "Table"];
                    StylePrint[someStr];
                    StylePrint[{"Got some import",Length[impTable],"graphPattern"/.configRule}];
                    If[ Length[impTable]>1,
                        numZ = ToExpression[StringReplace[#,("graphPattern"/.configRule)->""]]&/@ Map[StringSplit[#,{"/",":"}][[1,-2]]&,impTable ];
                        StylePrint[{"Got some numZ!", Length[numZ], numZ//First,"graphPattern"/.configRule,DateString[]}];
                        If[ Length[numZ]>0,
                            gNumZ = Select[numZ,#<=Length[myCuts]&];
                            If[ (gNumZ//Sort) ~UnsameQ~ (numZ//Sort),
                                StylePrint[{"Cuts no longer part of package.",
                                 Complement[numZ,gNumZ]}];
                                numZ = gNumZ;
                            ];
                            myyListOfNonZero[ "graphPattern"/.configRule] = numZ;
                            myGraphs = corruptGraph[toGraph[myCuts[[#]]]]&/@numZ;
                            StylePrint[{"Maybe got some graphs!", myGraphs//First,DateString[],"graphPattern"/.configRule}];
                            (newContainer["dressingPath"][corruptGraph[toGraph[myCuts[[#]]]]] = contactDir<>ToString[#]<>(
                            "dressPattern"/.configRule))&/@numZ;
                            StylePrint[{"Maybe got some load paths!", newContainer["dressingPath"][myGraphs//First],DateString[]}];
                        ]
                    ]
                ];
                StylePrint[{"Now looking at: ", Length[impTable],Length[impTable] ~SameQ~ 1}];
                If[ ("FastLoad"/.configRule) ~SameQ~ False || Length[impTable] ~SameQ~ 1,
                    StylePrint["Going slow -- fast is false or length 1"];
                    filenames = FileNames[theFGPath];
                    StylePrint[{"Filenames found:", Length[filenames]}];
                    Do[ 
                       ClearAll[myDelayedObject];
                       Get[filename];
                       myGraph = corruptGraph[myDelayedObject["graphName"/.configRule]];
                         (* Either no contList, or graph matches canonical form at a string level. *)
                       If[ ("verifyLoadSwitch"/.configRule) ~UnsameQ~ "verifyLoadSwitch",
                           test = myDelayedObject["verifyLoadSwitch"/.configRule];,
                           test = True
                       ];
                       If[ test,
                           If[ ListQ[contList = ("canonicalGraphLabels"/.configRule)]&& 
                           !isIsomorphic[
                           priveledgeLegs[stripPriveledge[myGraph]],
                              priveledgeLegs[stripPriveledge[ 
                                  toGraph[contList[[ newContainer["labelGenerator"][myGraph] ]]
                                         ]
                                       ]
                              ]
                              ],
                               If[ ("deleteNonCanon"/.configRule) ~SameQ~ True,
                                   DeleteFile[filename];
                                   strReplaceRule = ("graphPattern"->"dressPattern")/.configRule;
                                   DeleteFile[filename/.strReplaceRule];,
                                   StylePrint[{"Not Cannon.",filename}];
                               ],
                               myGraphs = Append[myGraphs,myGraph];
                               strReplaceRule = ("graphPattern"->"dressPattern")/.configRule;
                               newContainer["dressingPath"][myGraph] = StringReplace[filename,strReplaceRule];
                           ];
                       ];,
                        {filename,filenames}
                     ];
                ];
                StylePrint[{"Done loading", "labelName"/.configRule,DateString[]}];
                newContainer["myGraphs"] = myGraphs;
                StylePrint[{"Graphs loaded:", Length[myGraphs]}];
                newContainer
            ]
            
            ];
             myDressingObject["graphs"] = Join[myDressingObject["graphs"],Last[myDressingObject["delayedContainers"]]["myGraphs"]];
            );
        myDressingObject["dressings"][myGraph_] :=
            Module[ {myContainer,myDress,myFNz,cg1,myFNzU},
                myContainer = Select[myDressingObject["delayedContainers"],
                  (#["inclusionQ"][myGraph] ~SameQ~ True&&MemberQ[#["myGraphs"],corruptGraph[myGraph]])&,1];
                If[ myContainer ~UnsameQ~ {},
                    myContainer = First[myContainer];
                    ClearAll[myDelayedDressing];
                    StylePrint[{"Loading for the very first time!!",
                      myContainer["dressingPath"][corruptGraph[myGraph]]}];
                    Clear[completeGuy,currentDressHolderVar,myDelayedDressing];
                    completeGuy = False;
                    myFNz = FileNames[myContainer["dressingPath"][corruptGraph[myGraph]]<>"*"];
                    If[ Length[myFNz] ~SameQ~ 1,
                        Get[myContainer["dressingPath"][corruptGraph[myGraph]]];,
                        Get[myContainer["dressingPath"][corruptGraph[myGraph]]];
                        currentDressHolderVar = (myDelayedDressing["dressingName"/.myContainer["configRule"]])[corruptGraph[myGraph]];
                        cg1 = True&& completeGuy;
                        completeGuy = False;
                        myFNzU = Complement[myFNz, FileNames[myContainer["dressingPath"][corruptGraph[myGraph]]]];
                        Map[(Get[#];
                             currentDressHolderVar = (myDelayedDressing["dressingName"/.myContainer["configRule"]])[corruptGraph[myGraph]];
                             cg1 = cg1&&completeGuy;
                             completeGuy = False)&,myFNzU];
                        completeGuy = cg1;
                    ];
                    If[ completeGuy ~UnsameQ~ True,
                        If[ ("RequireConfirmedCompleteDressing"/.myContainer["configRule"]) ~SameQ~ True,
                            Throw[{"Dressing not confirmed with completeGuy=True", 
                                 myContainer["dressingPath"][corruptGraph[myGraph]]}];,
                            StylePrint[{"WARNING: Dressing not confirmed with completeGuy=True.  
       If you want me to halt,  set RequireConfirmedCompleteDressing -> True in config rules.  ", 
                            myContainer["dressingPath"][corruptGraph[myGraph]]}];
                        ];
                    ];
                    myDress = (myDelayedDressing["dressingName"/.myContainer["configRule"]])[corruptGraph[myGraph]];
                    myDressingObject["dressings"][corruptGraph[myGraph]] = If[ Head[myDress] ~SameQ~ (myDelayedDressing["dressingName"/.myContainer["configRule"]]),
                                                                               Throw[{"!!!!!!!!!WARNING!!!!!!! Need to redo graph dressing in:",
                                                                               myContainer["dressingPath"][corruptGraph[myGraph]], "and assign to graph",
                                                                                corruptGraph[myGraph], "probably because corruptGraph has been updated and this dressing has not." }];
                                                                               0,
                                                                               myDress
                                                                           ],
                    0
                ]
            ];
        myDressingObject["memberQ"][graph_] :=
            MemberQ[myDressingObject["graphs"],corruptGraph[graph]]||
            Select[myDressingObject["graphs"],isIsomorphic[corruptGraph[graph],#]&,1] ~UnsameQ~ {};
        myDressingObject["addDressing"][myGraphA_,myDressing_] :=
            (Module[ {myGraph = corruptGraph[myGraphA]},
                 If[ MemberQ[myDressingObject["graphs"],myGraph],
                     Throw["Graph already defined"]
                 ];
                 If[ (Select[myDressingObject["graphs"],isIsomorphic[myGraph,#]&] ~UnsameQ~ {} ),
                     Throw["< GRAPH CONFLICT >:  Graph really already defined"]
                 ];
                 myDressingObject["graphs"] = Append[myDressingObject["graphs"],myGraph];
                 myDressingObject["dressings"][myGraph] = myDressing/.red[a_,1,b_]:>-k[a]/.red[a_,0,b_]:>l[a];
             ]
            );
        myDressingObject["modifyDressing"][a_,b_] :=
            (
            If[ !MemberQ[myDressingObject["graphs"],corruptGraph[a]],
                Throw["Graph not defined for modify"]
            ];
            myDressingObject["dressings"][corruptGraph[a]] = b/.red[aa_,1,bb_]:>-k[aa]/.red[aa_,0,bb_]:>l[aa];
            );
        myDressingObject["findALabel"][myGraph_] :=
            (
            Module[ {myContainer,myId,configRule},
                myContainer = Select[myDressingObject["delayedContainers"],(#["inclusionQ"][corruptGraph[myGraph]] ~SameQ~ True)&,1];
                If[ myContainer ~UnsameQ~ {},
                    myContainer = First[myContainer];
                    configRule = myContainer["configRule"];
                    myId = (ToExpression["labelGenerator"/. configRule])[corruptGraph[myGraph]],
                    Throw[badCutThrowObject[{"no label possible!",consistentGraphToTrees[stripPriveledge[myGraph]],myGraph}]];
                ]
            ]
            );
        myDressingObject["findCompleteLabel"][myGraph_] :=
            (
            Module[ {myContainer,myId,configRule},
                myContainer = Select[myDressingObject["delayedContainers"],
                      (#["inclusionQ"][corruptGraph[myGraph]] ~SameQ~ True)&,1];
                If[ myContainer ~UnsameQ~ {},
                    myContainer = First[myContainer];
                    configRule = myContainer["configRule"];
                    { "labelName"/.configRule, myId = (ToExpression["labelGenerator"/. configRule])[corruptGraph[myGraph]]},
                    Throw[badCutThrowObject[{"no label possible!",consistentGraphToTrees[stripPriveledge[myGraph]],myGraph}]];
                ]
            ]
            );
        myDressingObject["getProfilePath"][myGraph_] :=
            (
            Module[ {myContainer,myId,configRule,path,dir,profilePath},
                myContainer = Select[myDressingObject["delayedContainers"],(#["inclusionQ"][corruptGraph[myGraph]] ~SameQ~ True)&,1];
                If[ myContainer ~UnsameQ~ {},
                    myContainer = First[myContainer];
                    configRule = myContainer["configRule"];
                    dir = ("graphListPath"/.configRule);
                    StylePrint[{"Container Label","labelName"/.configRule}];
                    myId = ToString[(ToExpression["labelGenerator"/. configRule])[myGraph]];
                    profilePath = dir<>"/"<>myId<>("profilePattern"/.configRule/.{"profilePattern"->".PROFILE"}),
                    Throw["no label possible!"];
                ]
            ]
            );
        myDressingObject["getCubicKLTDIFFPath"][myGraph_] :=
            (
            Module[ {myContainer,myId,configRule,path,dir,profilePath},
                myContainer = Select[myDressingObject["delayedContainers"],(#["inclusionQ"][corruptGraph[myGraph]] ~SameQ~ True)&,1];
                If[ myContainer ~UnsameQ~ {},
                    myContainer = First[myContainer];
                    configRule = myContainer["configRule"];
                    dir = ("graphListPath"/.configRule);
                    StylePrint[{"Container Label","labelName"/.configRule}];
                    myId = ToString[(ToExpression["labelGenerator"/. configRule])[myGraph]];
                    profilePath = dir<>"/"<>myId<>("cubicKLTPattern"/.configRule/.{"cubicKLTPattern"->".CUBICKLTDIFF"}),
                    Throw["no label possible!"];
                ]
            ]
            );
        myDressingObject["myContainer"][myGraph_] :=
            (
            Module[ {myContainer,myId,configRule},
                myContainer = Select[myDressingObject["delayedContainers"],(#["inclusionQ"][corruptGraph[myGraph]] ~SameQ~ True)&,1];
                If[ myContainer ~UnsameQ~ {},
                    myContainer,
                    Null
                ]
            ]
            );
        myDressingObject["myCutPath"][myGraphA_] :=
            Module[ {dir,
            DEATHALTCOUNTER,altID,myGraph = corruptGraph[myGraphA],
            myContainer,graphFile,graphPath,dressPath,tmpId,dressFile,loadFlag,contList,theDress,myId,
            configRule,verifyString,grapthPath},
                myContainer = Select[myDressingObject["delayedContainers"],((#["inclusionQ"][myGraph] ~SameQ~ True))&,1];
                If[ myContainer ~UnsameQ~ {},
                    myContainer = First[myContainer];
                    configRule = myContainer["configRule"];
                    dir = ("graphListPath"/.configRule);
                    myId = ToString[(ToExpression["labelGenerator"/. configRule])[myGraph]];
                    graphPath = dir<>"/"<>myId;
                ];
                graphPath
            ];
        myDressingObject["saveDelayedDressing"][myGraphA_] :=
            Module[ {dir,
            DEATHALTCOUNTER,altID,myGraph = corruptGraph[myGraphA],
            myContainer,graphFile,graphPath,dressPath,tmpId,dressFile,loadFlag,contList,theDress,myId,
            configRule,verifyString,grapthPath},
                myContainer = Select[myDressingObject["delayedContainers"],((#["inclusionQ"][myGraph] ~SameQ~ True))&,1];
                If[ NOSAVE ~SameQ~ False,
                    If[ myContainer ~UnsameQ~ {},
                        myContainer = First[myContainer];
                        configRule = myContainer["configRule"];
                        If[ ListQ[contList = ("canonicalGraphLabels"/.configRule)]&& 
                        !isIsomorphic[
                              priveledgeLegs[stripPriveledge[myGraph]],
                                 priveledgeLegs[stripPriveledge[ 
                                     toGraph[contList[[ myContainer["labelGenerator"][myGraph] ]]
                                            ]
                                          ]
                                 ]],
                            StylePrint[{"This is label",(myContainer["labelGenerator"][myGraph])}];
                            StylePrint[ stripPriveledge[myGraph]//InputForm];
                            StylePrint[ 
                                         stripPriveledge[
                                              toGraph[
                                                 contList[[ (myContainer["labelGenerator"][myGraph]) ]]
                                              ]
                                         ]
                         //InputForm];
                            StylePrint["Not cannonical storage, will not save!"];,
                            dir = ("graphListPath"/.configRule);
                            StylePrint[{"Container Label","labelName"/.configRule}];
                            myId = ToString[(ToExpression["labelGenerator"/. configRule])[myGraph]];
                            graphPath = dir<>"/"<>myId<>("graphPattern"/.configRule);
                            StylePrint[{"graph path",graphPath}];
                            dressPath = dir<>"/"<>myId<>("dressPattern"/.configRule);
                            StylePrint[{"dress path",dressPath}];
                            If[ FileExistsQ[graphPath]||FileExistsQ[dressPath],
                                StylePrint[{"WARNING --save file already exists",
                                            myId, 
                                            graphPath,
                                            dressPath}];
                                DEATHALTCOUNTER = 0;
                                While[DEATHALTCOUNTER<=3&(
                                      FileExistsQ[graphPath]||
                                      FileExistsQ[dressPath]
                                      ),
                                  DEATHALTCOUNTER++;
                                  altID = "_alt_"<>ToString[RandomPrime[{1,10^9}]];
                                  graphPath = graphPath<>altID;
                                  dressPath = dressPath<>altID;
                                 ];
                                If[ DEATHALTCOUNTER>=3,
                                    Throw[{"Couldn't create alt files",graphPath,dressPath}]
                                ];
                                StylePrint[{"altPath", myId, graphPath,dressPath}];
                            ];
                            verifyString = "verifyLoadSwitch"/.configRule;
                            theDress = myDressingObject["dressings"][myGraph];
                            loadFlag = If[ theDress ~SameQ~ 0 || 
                                        Head[theDress] ~SameQ~ 
                                     myContainer["dressings"],
                                           theDress = 0;
                                           False,
                                           True
                                       ];
                            ClearAll[myDelayedObject];
                            myDelayedObject["verifyLoadSwitch"/.configRule] = loadFlag;
                            myDelayedObject["graphName"/.configRule] = myGraph;
                            If[ Head[GRAPHPROPERTIES[myGraph]] ~UnsameQ~ GRAPHPROPERTIES,
                                myDelayedObject["graphProperties"] = GRAPHPROPERTIES[myGraph]
                            ];
                            If[ ToExpression[("cleanFile"/.configRule)] ~SameQ~ True && 
                            FileExistsQ[graphPath],
                                DeleteFile[graphPath];
                            ];
                            If[ ToExpression[("cleanFile"/.configRule)] ~SameQ~ True && 
                            FileExistsQ[dressPath],
                                DeleteFile[dressPath];
                            ];
                            tmpIDGraph = graphPath<>(myId<>"_tmp_"<>
                               ToString[RandomPrime[{1,10^9}]]);
                            Save[tmpIDGraph, myDelayedObject];
                            RenameFile[tmpIDGraph, graphPath, OverwriteTarget->True];
                            ClearAll[myDelayedDressing];
                            myDelayedDressing["dressingName"/. 
                                myContainer["configRule"]][myGraph] = theDress;
                            tmpIDDress = dressPath<>(myId<>"_tmp_"<>
                              ToString[RandomPrime[{1,10^9}]]);
                            Save[tmpIDDress, myDelayedDressing];
                            Clear[completeGuy];
                            completeGuy = True;
                            Save[tmpIDDress, completeGuy];
                            RenameFile[tmpIDDress, dressPath, OverwriteTarget->True];,
                            Throw["No container suitable for saving delayed dressing!"]
                        ];
                    ]
                ]
            ];
        myDressingObject["saveAdjunctDelayedDressing"][myGraphA_, adjunctDressing_,adjunctLabel_] :=
            Module[ {dir,DEATHALTCOUNTER,altID,myGraph = corruptGraph[myGraphA],
            myContainer,graphFile,graphPath,dressPath,tmpId,dressFile,loadFlag,contList,theDress,myId,configRule,verifyString,grapthPath},
                myContainer = Select[myDressingObject["delayedContainers"],((#["inclusionQ"][myGraph] ~SameQ~ True))&,1];
                If[ NOSAVE ~SameQ~ False,
                    If[ myContainer ~UnsameQ~ {},
                        myContainer = First[myContainer];
                        configRule = myContainer["configRule"];
                        If[ ListQ[contList = ("canonicalGraphLabels"/.configRule)]&& 
                        !isIsomorphic[
                              priveledgeLegs[stripPriveledge[myGraph]],
                                 priveledgeLegs[stripPriveledge[ 
                                     toGraph[contList[[ myContainer["labelGenerator"][myGraph] ]]
                                            ]
                                          ]
                                 ]],
                            StylePrint[{"This is label",(myContainer["labelGenerator"][myGraph])}];
                            StylePrint[ stripPriveledge[myGraph]//InputForm];
                            StylePrint[ 
                                         stripPriveledge[
                                              toGraph[
                                                 contList[[ (myContainer["labelGenerator"][myGraph]) ]]
                                              ]
                                         ]
                         //InputForm];
                            StylePrint["Not cannonical storage, will not save!"];,
                            dir = ("graphListPath"/.configRule);
                            StylePrint[{"Container Label","labelName"/.configRule}];
                            myId = ToString[(ToExpression["labelGenerator"/. configRule])[myGraph]];
                            dressPath = dir<>"/"<>myId<>("dressPattern"/.configRule)<>adjunctLabel;
                            StylePrint[{"dress path",dressPath}];
                            verifyString = "verifyLoadSwitch"/.configRule;
                            ClearAll[myDelayedDressing,currentDressHolderVar];
                            myDelayedDressing["dressingName"/.myContainer["configRule"]][myGraph] = (
                               currentDressHolderVar+adjunctDressing);
                            tmpIDDress = StringReplace[dressPath,myId->myId<>"_tmp_"<>ToString[RandomPrime[{1,10^9}]]];
                            Save[tmpIDDress, myDelayedDressing];
                            Clear[completeGuy];
                            completeGuy = True;
                            Save[tmpIDDress, completeGuy];
                            RenameFile[tmpIDDress, dressPath, OverwriteTarget->True];,
                            Throw["No container suitable for saving delayed dressing!"]
                        ];
                    ]
                ]
            ];
        myDressingObject
    ];


DotPower[expr_, a_, n_] :=
    If[ n > 1,
        DotPower[expr . a, a, n - 1],
        expr . a
    ];


(* ::Subsubsection:: *)
(*Vacuum Code*)


toVacuum[graph_] :=
    graph/.neckl[a__]:>neckl[a/.k[___]:>{}//Flatten]


toVacuum[graph_,extLabels__] :=
    graph/.neckl[a__]:>neckl[a/.Map[#->{}&,extLabels]//Flatten]


stripVacuum[graph_] :=
    Module[ {origGraph = graph,myRule,StylePrint,nextRules = graph/.neckl[{a_,b_}]:>Sow[neckl[{a,b}]]//Reap//Last//Flatten//Union},
        StylePrint[nextRules];
        While[nextRules ~UnsameQ~ {},
        myRule = First[nextRules]/.neckl[{a_,b_}]:>Sort[{-a->b,a->-b}];
        StylePrint[myRule];
        origGraph = origGraph//.myRule  /.neckl[{a_,a_}]:>{}/.neckl[{a_,-a_}]:>{}/.neckl[{-a_,a_}]:>{}/.vertexFormGraph[a__]:>vertexFormGraph[Flatten[a]];
        StylePrint[origGraph];
        nextRules = origGraph/.neckl[{a_,b_}]:>Sow[neckl[{a,b}]]//Reap//Last//Flatten//Union;
            ];
        origGraph
    ]



canonVacuumGraph[graph_] :=
    Module[ {myExtra,myV,v},
        myExtra = Block[ {cntr = 0},
                      (graph//toVacuum//stripVacuum//corruptGraph//graphHashCode//EdgeRules)/.Rule[a_,b_]:>{Rule[a,b],l[++cntr]}
                  ]/.Rule[a_,b_]:>Rule[v[a],v[b]];
        Map[ #/.{Rule[a_,b_],ll_}:>(myV[a] = Flatten[{myV[a],ll}/.myV[___]:>{}];
                                    myV[b] = Flatten[{myV[b],-ll}/.myV[___]:>{}];)&,myExtra];
        vertexFormGraph[Map[ neckl[myV[#]]&,myExtra/.v[a_]:>Sow[v[a]]//Reap//Last//Flatten//Union]]
    ]


stripVacuumRules[graph_] :=
    Module[ {origGraph = graph, myRule, cacheRules = {},StylePrint, 
    nextRules = Union[Flatten[Last[Reap[graph /. neckl[{a_, b_}] :> 
    Sow[neckl[{a, b}]]]]]]},
        StylePrint[nextRules];
        While[nextRules  ~UnsameQ~  {}, myRule = First[nextRules] /. neckl[{a_, b_}] :> 
             Sort[{-a -> b, a -> -b}];
                                cacheRules = Flatten[{cacheRules,myRule} ];
                                StylePrint[myRule];
                                origGraph = origGraph //. myRule /. neckl[{a_, a_}] :> {} /. 
                                    neckl[{a_, -(a_)}] :> {} /. neckl[{-(a_), a_}] :> {} /. 
                                  vertexFormGraph[a__] :> vertexFormGraph[Flatten[a]];
                                StylePrint[origGraph];
                                nextRules = Union[Flatten[Last[Reap[origGraph /. neckl[{a_, b_}] :> 
                                       Sow[neckl[{a, b}]]]]]]; ];
        cacheRules
    ]


collectionAdd[collection_, item_, canonForm_] :=
 
 Module[{pos, newCol, maxItem},
  maxItem = collection["maxItem"] /. collection[___] :> 0;
  collection[canonForm[item]] = 
   collection[canonForm[item]] /. collection[___] :> (++maxItem);
  collection[  collection[canonForm[item]]] = 
   collection[  collection[canonForm[item]]] /. ( collection[___] :> item);
  collection["maxItem"] = maxItem;
  collection["allItems"];
  collection[canonForm[item]]
  ]

opTillClosure[initial_, op_] := 
 Module[{collection, oldMaxId, myCurrentID, uHatList},
  Clear[collection];
  myCurrentID = collectionAdd[collection, initial, graphHashCode];
  uHatList = {};
  uHatList = FixedPoint[  (oldMaxId = collection["maxItem"];
      uHatList = Flatten[{#,
           
           Table[StylePrint[{leg, 
              myRow = (myCurrentID <-> 
                 collectionAdd[collection, 
                  op[collection[myCurrentID] // stripPriveledge, leg] // 
                   priveledgeLegs, graphHashCode])}];
            
            myRow, {leg, 
             getIntLegs[collection[myCurrentID] // stripPriveledge]}]}] /. 
         UndirectedEdge[a_, b_] :> UndirectedEdge @@ Sort[{a, b}] // Union;
      If[oldMaxId < collection["maxItem"] || myCurrentID < oldMaxId,
       myCurrentID++];
      uHatList) &, uHatList];
  {uHatList, collection}]

opTillClosureOnCollection[collection_, {op1_, op2_}] := 
   Module[{oldMaxId, myCurrentID, graphList}, myCurrentID = 1; graphList = {}; 
     graphList = FixedPoint[(oldMaxId = collection["maxItem"]; 
         graphList = Union[Flatten[{#1, Table[DYNSTAT = {leg, myRow = myCurrentID <-> 
                    collectionAdd[collection, priveledgeLegs[op1[stripPriveledge[
                        collection[myCurrentID]], leg]], graphHashCode]}; 
                {myRow, myCurrentID <-> collectionAdd[collection, priveledgeLegs[
                    op2[stripPriveledge[collection[myCurrentID]], leg]], 
                   graphHashCode]}, {leg, getIntLegs[stripPriveledge[collection[
                   myCurrentID]]]}]}] /. UndirectedEdge[a_, b_] :> UndirectedEdge @@ 
              Sort[{a, b}]]; If[oldMaxId < collection["maxItem"] || 
           myCurrentID < oldMaxId, myCurrentID++]; graphList) & , graphList]; 
     {graphList/.TwoWayRule[a_,b_]:>TwoWayRule@@Sort[{a,b}]//Union, collection}]

getVertices[vertexFormGraph[necklaceList__]] := necklaceList

flattenMomenta[myMom4_, \[Xi]_] := (myMom4[-flat[a_]] := -myMom4[flat[a]]; 
  myMom4[flat[a_]] := If[Head[a]  ~SameQ~  Plus,
    (myMom4 /@ 
       a) - \[Xi] Lorentz[(myMom4 /@ a), (myMom4 /@ 
          a)]/(2 Lorentz[(myMom4 /@ a), \[Xi]]),
    (myMom4[
       a]) - \[Xi] Lorentz[(myMom4[a]), (myMom4[
          a])]/(2 Lorentz[(myMom4[a]), \[Xi]])];)

atreeRule = {Atree[lab_List, h_List] :> Module[{s = Plus @@ h, l = Length[h]},    
     If[(Abs[s] >= l - 2 && l > 3) || (Abs[s] ~SameQ~ l), 0, 
       If[ (Abs[s] ~SameQ~ Abs[l - 4]), 
                  MHVtree[lab, h], 
                 Print[Style[{Abs[s], l - 2, l, h}]]; Atree[lab, h]]]]}; 

mhvTreeRule = {MHVtree[lab_, hel_] :> Module[{negLegs, posLegs}, 
             
     negLegs = 
      Transpose[Select[Transpose[{lab, hel}], #1[[2]] < 0 & ]][[1]]; 
              
     posLegs = 
      Transpose[Select[Transpose[{lab, hel}], #1[[2]] > 0 & ]][[1]]; 
              If[Length[negLegs] ~SameQ~ 2, (I*spa[negLegs[[1]], negLegs[[2]]]^4)/
                  (Times @@ 
          Table[spa[lab[[i]], lab[[i + 1]]], {i, 1, Length[lab] - 1}]*
         spa[Last[lab], lab[[1]]]), 
                -(I*spb[posLegs[[2]], posLegs[[1]]]^4)/(Times @@ 
          Table[spb[lab[[i + 1]], lab[[i]]], 
                      {i, 1, Length[lab] - 1}]*spb[Last[lab], lab[[1]]]), 
      Null]]}; 

mhvTreeRules = Join[atreeRule, mhvTreeRule]; 

showCollection[list_] := Manipulate[
  Show[list[[Round[x]]], ImageSize -> Large, 
   PlotLabel -> Style[Round[x], {Large, Bold}]], {x, 1, Length[list]}]

cutSumFormat[cut__] := 
 TraditionalForm[
  Style[DisplayForm[
    RowBox[{UnderscriptBox["\[Sum]", RowBox[{"s", " \[Element] ", "States"}]],
        Sequence @@ cut} /. Module[{myRule = Map[# -> SuperscriptBox[#, "s"] &,
          First /@ 
           Select[ Flatten[cut /. -a_ :> a /. Atree :> List] // 
             Tally, #[[2]]  ~SameQ~  2 &]], StylePrint}, StylePrint[myRule]; 
       myRule]]], Large]]


(* ::Subsubsection::Closed:: *)
(**)


(*Private Methods & Close*)
Begin["`Private`"]


End[]

EndPackage[]


