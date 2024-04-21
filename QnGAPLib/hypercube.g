#######################################################################
##
##  hypercube.g  (Version 0.530)  HYPERCUBE Library   Anton Bondarenko
##
##
##  Copyright 2001 Anton Bondarenko      ( 21.03.01 )
##
##
##  03.05.01 - Added draft of perfect matching enumerator in 4-cubes
##  05.05.01 - Added draft of perfect colorings of 4-cubes
##
##
##
##
##
######################################################################

if not QUIET and BANNER  then
   Print("\nLoading  HYPERCUBE 0.530,");
   Print("\nby A.V.Bondarenko (a_bond@cs.urao.edu)\n\n"); 
fi;

RequirePackage("grape");

######################################################################
Hypercube := function(n)
######################################################################
#
# Returns the n-cube, whose vertices are the subsets
# of {1,...,n},  with x joined to y iff  Intersection(x,y)
# has size  e-1 and x has size e and y has size e-1 or vv.
#
#####################################################################
local rel,H;
if not IsInt(n) then 
   Error("usage: Hypercube( <Int> )");
fi;
if n<1 then
   Error("must have 1 <= <n>");
fi;
rel := function(x,y)
   return   AbsInt(Length(x)-Length(y))=1 and 
   Length(Intersection(x,y))=Minimum(Length(x),Length(y));
end;
H:=Graph(SymmetricGroup(2),Combinations([1..n]),OnSets,rel,true);
H.isSimple:=true;
return H;
end;


######################################################################
LoadMatchingsLibrary := function(n)
######################################################################
#
# Loads sets of matchings of Hn (n = 0,2,3,4)
#
# n=0 loads all matchings
#
# ATT! Currently released only for n = 4 !!!
#
######################################################################
if not IsInt(n) then
   Error("usage: LoadMatchingsLibrary(< Int n >) where n= 1, 2, 3 or 0");
fi;
if not n in [0,2,3,4] then
   Error("n must be equal to 0 or 2 or 3 or 4");
fi;
if n=4 then
   Read(Filename(DirectoriesPackageLibrary("hypercube"),"hn_4.lib"));

fi;
if n=3 then
   Read(Filename(DirectoriesPackageLibrary("hypercube"),"hn_3.lib"));
fi;
if n=2 then
   Read(Filename(DirectoriesPackageLibrary("hypercube"),"hn_2.lib"));
fi;
if n=0 then
   Read(Filename(DirectoriesPackageLibrary("hypercube"),"hn_4.lib"));
   Read(Filename(DirectoriesPackageLibrary("hypercube"),"hn_3.lib"));
   Read(Filename(DirectoriesPackageLibrary("hypercube"),"hn_2.lib"));
fi;
end;


######################################################################
AutGroupHypercubeSubgraph:=function(gamma, n)
######################################################################
#
# Returns automorphism group of gamma as subgroup of Aut(Hn)
#
#
######################################################################
return
  (Intersection(AutGroupGraph(Hypercube(n)),AutGroupGraph(gamma)));
end;

######################################################################
AutGroupH4Colorings:=function(colorings, m_list, aut_h)
######################################################################
#
# Returns automorphism group of H4 coloring presented by
# signature, aut group and matchings rep list
#
#
######################################################################
local aut_c, col_bs, i;
aut_c:=[];
col_bs:=H4ColoringBySignature(colorings,m_list,aut_h);
for i in [1..4] do
  Add(aut_c, AutGroupHypercubeSubgraph(GraphByVL(col_bs[i]),4));
od;
return Intersection(Intersection(aut_c[1],aut_c[2]),
                    Intersection(aut_c[3],aut_c[4]));
end;


######################################################################
AllH4MatchingsAutRep := function()
######################################################################
#
# Finds all non-isomorpfic perfect matchings reprezentatives in H4
#
# 
######################################################################
local Matchings, Pretender, n, EdgesList, SizeOfMatching, 
AutH4List, NumberOfVertexes, NumberOfEdges, p, j, i;
Matchings := [];
Pretender := [];
SizeOfMatching := 8;
NumberOfVertexes := 16;
NumberOfEdges := 32;
AutH4List := List(AutGroupGraph(Hypercube(4)));
n := [1..SizeOfMatching];
EdgesList := List(UndirectedEdges(Hypercube(4)));
p := SizeOfMatching;
repeat
  Pretender := [];
  for j in [1..SizeOfMatching] do
    Add(Pretender, EdgesList[n[j]]);     
  od;
  if IsMatchingVL(Pretender) then
    if not ExistInSVL(Pretender, Matchings, AutH4List) then
      Add(Matchings, Pretender);
    fi;
  fi;
  if n[SizeOfMatching] = NumberOfEdges then p := p - 1;
    else p := SizeOfMatching;
  fi;
  if p >= 1 then
    for i in [p..SizeOfMatching] do
      n[SizeOfMatching-i+p] := n[p] + (SizeOfMatching-i+p) - p + 1;
    od;
  fi;
until p < 1;
return Matchings;
end;


######################################################################
AllH4MatchingsAutRepLib := function()
######################################################################
#
# Reads all non-isomorpfic perfect matchings reprezentatives in H4
#
# 
######################################################################
local matchings;
matchings:=ReadAsFunction(Filename(
  DirectoriesPackageLibrary("hypercube"),"h4_matchings.g"));
return matchings();
end;


######################################################################
AllH4MatchingsPairsAutRepLib := function()
######################################################################
#
# Reads all non-isomorpfic perfect matchings pairs 
# reprezentatives in H4
#
# 
######################################################################
local pairs;
pairs:=ReadAsFunction(Filename(
  DirectoriesPackageLibrary("hypercube"),"h4_pairs.g"));
return pairs();
end;


######################################################################
IsMatchingVL := function(gammaVL)
######################################################################
#
# Check if gammaVL is a list of edge of matching
#
#
######################################################################
local VertexesList, SSVertexesList, i;
VertexesList := [];
SSVertexesList := [];
for i in gammaVL do
  Append(VertexesList, [i[1],i[2]]);
od;
SSVertexesList := SSortedList(VertexesList);
if Size(VertexesList) = Size(SSVertexesList) then return true;
  else return false; fi;
end;

######################################################################
ExistInMVL := function(matching, m_list, aut_g_list)
######################################################################
#
# Check if m_list already contains matching, isomorphical to   
# given one under group listed in aut_g_list.
#
#
######################################################################
local p, pm, s, sm, mr, i, j, k;
for k in [1..Size(aut_g_list)] do
  mr:=OnTuplesTuples(matching, aut_g_list[k]);
  p:=[];
  for i in [1..Size(matching)] do
    p[mr[i][1]]:=mr[i][2];
    p[mr[i][2]]:=mr[i][1];
  od;
  pm:=PermList(p);
  for i in [1..Size(m_list)] do
    s:=[];
    for j in [1..Size(m_list[i])] do
      s[m_list[i][j][1]]:=m_list[i][j][2];
      s[m_list[i][j][2]]:=m_list[i][j][1];
    od;
    sm:=PermList(s);
    if pm = sm then return true; fi;
  od;
od;
return false;
end;


AllH4ProperColoringsAutRep := function()
#############################################################################
#
# Finds all non-isomorphic H4 proper colorings
#
#
#############################################################################
local aut_h, aut_p, aut_m, ag, agp, agp1, agp2, matchings, pairs, colorings, 
NOM, SOC, SOHAG, p, n, s, i, j, k, l, m1, m2, m3, m4, ms1, ms2, p1, p2, 
c1, c2, c3, c4, a1, a2, a12, colors, dn, df, ds, pretender, p_list, b, tmp;
#
# data definitions and initialization
#
matchings:=AllH4MatchingsAutRepLib();
NOM:=Size(matchings);
pairs:=AllH4MatchingsPairsAutRepLib();
colorings:=[];
aut_h:=List(AutGroupGraph(Hypercube(4)));
aut_m:=[];
for i in [1..Size(matchings)] do
  Add(aut_m, List(AutGroupHypercubeSubgraph(GraphByVL(matchings[i]), 4)));
od;
aut_p:=[];
for i in [1..NOM] do Add(aut_p, []); od;
for i in [1..NOM] do for j in [1..NOM] do Add(aut_p[i], []); od; od;
for i in [1..NOM] do for j in [i..NOM] do for k in [1..Size(matchings[i][j])] do
  Add(aut_p[i][j], List(AutGroupHypercubeSubgraph(GraphByVL(
  Concatenation(matchings[i],OnTuplesTuples(matchings[j], aut_h[k]))),4))); od; od; od;
SOC:=4;
SOHAG:=Size(aut_h);
#
# for class 1-1-1-1:
#
Add(colorings, []);
n := [1..SOC];
p := SOC;
#
# for each four different matchings:
#
repeat
#
# choose smallest matchings complect in class
#
  if ( Size(pairs[n[1]][n[2]]) * Size(pairs[n[3]][n[4]]) ) <
     ( Size(pairs[n[1]][n[3]]) * Size(pairs[n[2]][n[4]]) ) then
    if ( Size(pairs[n[1]][n[2]]) * Size(pairs[n[3]][n[4]]) ) <
       ( Size(pairs[n[1]][n[4]]) * Size(pairs[n[2]][n[3]]) ) then
      c1:=n[1]; c2:=n[2]; c3:=n[3]; c4:=n[4];
    else
      c1:=n[1]; c2:=n[4]; c3:=n[2]; c4:=n[3]; fi;
  else
    if ( Size(pairs[n[1]][n[3]]) * Size(pairs[n[2]][n[4]]) ) <
       ( Size(pairs[n[1]][n[4]]) * Size(pairs[n[2]][n[3]]) ) then
      c1:=n[1]; c2:=n[3]; c3:=n[4]; c4:=n[4];
    else
      c1:=n[1]; c2:=n[4]; c3:=n[2]; c4:=n[3]; fi;
  fi;
#
# search for complement pairs of choosen matchings set
#
  for a1 in [1..Size(pairs[c1][c2])] do
    p1:=Concatenation(matchings[c1], OnTuplesTuples(matchings[c2], aut_h[pairs[c1][c2][a1]]));
    for a2 in [1..Size(pairs[c3][c4])] do
      p2:=Concatenation(matchings[c3], OnTuplesTuples(matchings[c4], aut_h[pairs[c3][c4][a2]]));
      l:=0;
      a12:=0;
      while (a12 < SOHAG and l = 0) do
        a12:=a12 + 1;
        colors:=Concatenation(p1, OnTuplesTuples(p2, aut_h[a12]));
        if IsSimpleGraphVL(colors) then
          l:=a12;
        fi;
      od;
      if l <> 0 then
#
# add finded coloring
#
        Add(colorings[1], [a12, pairs[c1][c2][a1], c1, c2, pairs[c3][c4][a2], c3, c4]);
#
# generate other designs of choosen matchings set
#
        ag:=List(Difference(
              Difference(
                AutGroupHypercubeSubgraph(GraphByVL(OnTuplesTuples(p2, aut_h[a12])), 4),
                AutGroupHypercubeSubgraph(
                  GraphByVL(OnTuplesTuples(matchings[c3], aut_h[a12])), 4)),
              Intersection(
                AutGroupHypercubeSubgraph(GraphByVL(matchings[c1]), 4),
                AutGroupHypercubeSubgraph(GraphByVL(
                  OnTuplesTuples(matchings[c2], aut_h[pairs[c1][c2][a1]])), 4))));
        for i in [1..Size(ag)] do
          Add(colorings[1], [ag[i], pairs[c1][c2][a1], c1, c2, pairs[c3][c4][a2], c3, c4]);
        od;
      fi;
    od;
  od;
#
# generate next complect
#
  if n[SOC] = NOM then p := p - 1;
    else p := SOC;
  fi;
  if p >= 1 then
    for i in [p..SOC] do
      n[SOC - i + p] := n[p] + (SOC - i + p) - p + 1;
    od;
  fi;
until p < 1;
#
# for class 2-1-1
#
Add(colorings, []);
SOC:=3;
n := [1..SOC];
p := SOC;
#
# for each three matchings (first applied twice):
#
repeat
for dn in [1..3] do
  if dn = 1 then df:=2; ds:=3; fi;
  if dn = 2 then df:=1; ds:=3; fi;
  if dn = 3 then df:=1; ds:=2; fi;
  c1:=n[dn]; c2:=n[dn]; c3:=n[df]; c4:=n[ds];
#
# search for complement pairs of choosen matchings set
#
  for a1 in [1..Size(pairs[c1][c2])] do
    p1:=Concatenation(matchings[c1], OnTuplesTuples(matchings[c2], aut_h[pairs[c1][c2][a1]]));
    for a2 in [1..Size(pairs[c3][c4])] do
      p2:=Concatenation(matchings[c3], OnTuplesTuples(matchings[c4], aut_h[pairs[c3][c4][a2]]));
      l:=0;
      a12:=0;
      while (a12 < SOHAG and l = 0) do
        a12:=a12 + 1;
        colors:=Concatenation(p1, OnTuplesTuples(p2, aut_h[a12]));
        if IsSimpleGraphVL(colors) then
          l:=a12;
        fi;
      od;
      if l <> 0 then
#
# add finded coloring
#
        Add(colorings[2], [a12, pairs[c1][c2][a1], c1, c2, pairs[c3][c4][a2], c3, c4]);
#
# generate other designs of choosen matchings set
#
        agp:=List(Difference(
              Difference(
                AutGroupHypercubeSubgraph(GraphByVL(OnTuplesTuples(p2, aut_h[a12])), 4),
                AutGroupHypercubeSubgraph(
                  GraphByVL(OnTuplesTuples(matchings[c3], aut_h[a12])), 4)),
              Intersection(
                AutGroupHypercubeSubgraph(GraphByVL(matchings[c1]), 4),
                AutGroupHypercubeSubgraph(GraphByVL(
                  OnTuplesTuples(matchings[c2], aut_h[pairs[c1][c2][a1]])), 4))));
        ag:=[];
#
# exclude also isomorphism m1 <-> Image(m2)
#
        for i in [1..Size(agp)] do
          ms1:=[];
          ms2:=[];
          for j in [1..Size(matchings[c1])] do Add(ms1, Set(matchings[c1][j])); od;
          for j in [1..Size(m2)] do Add(ms1, Set(OnTuplesTuples(matchings[c2], agp[i])[j])); od;
          if Set(ms1) <>  Set(ms2) then
            Add(ag, agp[i]);
          fi;
        od;
        for i in [1..Size(ag)] do
          Add(colorings[2], [Position(aut_h,ag[i]), pairs[c1][c2][a1], c1, c2, 
                                                    pairs[c3][c4][a2], c3, c4]);
        od;
      fi;
    od;
  od;
od;
#
# generate next complect
#
  if n[SOC] = NOM then p := p - 1;
    else p := SOC;
  fi;
  if p >= 1 then
    for i in [p..SOC] do
      n[SOC - i + p] := n[p] + (SOC - i + p) - p + 1;
    od;
  fi;
until p < 1;
#
# for class 2-2
#
Add(colorings, []);
SOC:=2;
n := [1..SOC];
p := SOC;
#
# for each two matchings (each applied twice):
#
repeat
  c1:=n[1]; c2:=n[1]; c3:=n[2]; c4:=n[2];
#
# search for complement pairs of choosen matchings set
#
  tmp:=[];
  for a1 in [1..Size(pairs[c1][c2])] do
    p1:=Concatenation(matchings[c1], OnTuplesTuples(matchings[c2], aut_h[pairs[c1][c2][a1]]));
    for a2 in [1..Size(pairs[c3][c4])] do
      p2:=Concatenation(matchings[c3], OnTuplesTuples(matchings[c4], aut_h[pairs[c3][c4][a2]]));
      l:=0;
      a12:=0;
      while (a12 < SOHAG and l = 0) do
        a12:=a12 + 1;
        colors:=Concatenation(p1, OnTuplesTuples(p2, aut_h[a12]));
        if IsSimpleGraphVL(colors) then
          l:=a12;
        fi;
      od;
      if l <> 0 then
#
# add finded coloring
#
        Add(tmp, [a12, pairs[c1][c2][a1], c1, c2, pairs[c3][c4][a2], c3, c4]);
#
# generate other designs of choosen matchings set
#
        agp1:=List(Difference(
              Difference(
                AutGroupHypercubeSubgraph(GraphByVL(OnTuplesTuples(p2, aut_h[a12])), 4),
                AutGroupHypercubeSubgraph(
                  GraphByVL(OnTuplesTuples(matchings[c3], aut_h[a12])), 4)),
              Intersection(
                AutGroupHypercubeSubgraph(GraphByVL(matchings[c1]), 4),
                AutGroupHypercubeSubgraph(GraphByVL(
                  OnTuplesTuples(matchings[c2], aut_h[pairs[c1][c2][a1]])), 4))));
        agp2:=[];
#
# exclude also isomorphism m3 <-> Image(m4)
#
        for i in [1..Size(agp1)] do
          ms1:=[];
          ms2:=[];
          for j in [1..Size(matchings[c3])] do 
            Add(ms1, Set(OnTuplesTuples(matchings[c3], aut_h[a12])[j]));
          od;
          for j in [1..Size(matchings[c4])] do 
            Add(ms1, Set(OnTuplesTuples(
              OnTuplesTuples(OnTuplesTuples(matchings[c4], aut_h[pairs[c3][c4][a2]]),aut_h[a12])
                                            , agp1[i])[j])); 
          od;
          if Set(ms1) <>  Set(ms2) then
            Add(agp2, agp1[i]);
          fi;
        od;
        ag:=[];
#
# exclude also isomorphism m1 <-> Image(m2)
#
        for i in [1..Size(agp1)] do
          ms1:=[];
          ms2:=[];
          for j in [1..Size(matchings[c1])] do Add(ms1, Set(matchings[c1][j])); od;
          for j in [1..Size(matchings[c2])] do Add(ms1, Set(OnTuplesTuples(matchings[c2], agp1[i])[j])); od;
          if Set(ms1) <>  Set(ms2) then
            Add(ag, agp2[i]);
          fi;
        od;
#
# exclude also isomorphisms between brothers
#
        p_list:=[tmp[Size(tmp)]];
        for i in [1..Size(ag)] do
          pretender:=[Position(aut_h,aut_h[a12]*ag[i]), pairs[c1][c2][a1], c1, c2, 
                                                               pairs[c3][c4][a2], c3, c4];
          b:=false;
          for k in [1..Size(p_list)] do
            ms1:=[];
            ms2:=[];
            for j in [1..Size(matchings[pretender[6]])] do 
              Add(ms1, Set(OnTuplesTuples(matchings[pretender[6]], aut_h[pretender[1]])[j]));
            od;
            for j in [1..Size(matchings[p_list[k][6]])] do 
              Add(ms2, Set(OnTuplesTuples(matchings[p_list[k][6]], aut_h[p_list[k][1]])[j]));
            od;
            if Set(ms1) = Set(ms2) then
              b:=true;
              break;
            fi;
          od;
          if not b then
            Add(p_list, pretender);
            Add(tmp, pretender);
          fi;
        od;
      fi;
    od;
  od;
  Append(colorings[3], H4ColoringsSetAutRep(tmp, matchings, aut_h));
#
# generate next complect
#
  if n[SOC] = NOM then p := p - 1;
    else p := SOC;
  fi;
  if p >= 1 then
    for i in [p..SOC] do
      n[SOC - i + p] := n[p] + (SOC - i + p) - p + 1;
    od;
  fi;
until p < 1;
#
# for class 3-1
#
Add(colorings, []);
SOC:=2;
n := [1..SOC];
p := SOC;
#
# for each two matchings (first applied three times):
#
repeat
tmp:=[];
for df in [1..2] do
  if df = 1 then ds:=2; fi;
  if df = 2 then ds:=1; fi;
  c1:=n[df]; c2:=n[df]; c3:=n[df]; c4:=n[ds];
#
# search for complement pairs of choosen matchings set
#
  for a1 in [1..Size(pairs[c1][c2])] do
    p1:=Concatenation(matchings[c1], OnTuplesTuples(matchings[c2], aut_h[pairs[c1][c2][a1]]));
    for a2 in [1..Size(pairs[c3][c4])] do
      p2:=Concatenation(matchings[c3], OnTuplesTuples(matchings[c4], aut_h[pairs[c3][c4][a2]]));
      l:=0;
      a12:=0;
      while (a12 < SOHAG and l = 0) do
        a12:=a12 + 1;
        colors:=Concatenation(p1, OnTuplesTuples(p2, aut_h[a12]));
        if IsSimpleGraphVL(colors) then
          l:=a12;
        fi;
      od;
      if l <> 0 then
#
# add finded coloring
#
        Add(tmp, [a12, pairs[c1][c2][a1], c1, c2, pairs[c3][c4][a2], c3, c4]);
#
# generate other designs of choosen matchings set
#
        agp1:=List(Difference(
              Difference(
                AutGroupHypercubeSubgraph(GraphByVL(OnTuplesTuples(p2, aut_h[a12])), 4),
                AutGroupHypercubeSubgraph(
                  GraphByVL(OnTuplesTuples(matchings[c3], aut_h[a12])), 4)),
              Intersection(
                AutGroupHypercubeSubgraph(GraphByVL(matchings[c1]), 4),
                AutGroupHypercubeSubgraph(GraphByVL(
                  OnTuplesTuples(matchings[c2], aut_h[pairs[c1][c2][a1]])), 4))));
        agp2:=[];
#
# exclude also isomorphism m3 <-> Image(m4)
#
        for i in [1..Size(agp1)] do
          ms1:=[];
          ms2:=[];
          for j in [1..Size(matchings[c3])] do 
            Add(ms1, Set(OnTuplesTuples(matchings[c3], aut_h[a12])[j]));
          od;
          for j in [1..Size(matchings[c4])] do 
            Add(ms1, Set(OnTuplesTuples(
              OnTuplesTuples(OnTuplesTuples(matchings[c4], aut_h[pairs[c3][c4][a2]]),aut_h[a12])
                                            , agp1[i])[j])); 
          od;
          if Set(ms1) <>  Set(ms2) then
            Add(agp2, agp1[i]);
          fi;
        od;
        ag:=[];
#
# exclude also isomorphism m1 <-> Image(m2)
#
        for i in [1..Size(agp1)] do
          ms1:=[];
          ms2:=[];
          for j in [1..Size(matchings[c1])] do Add(ms1, Set(matchings[c1][j])); od;
          for j in [1..Size(matchings[c2])] do Add(ms1, Set(OnTuplesTuples(matchings[c2], agp1[i])[j])); od;
          if Set(ms1) <>  Set(ms2) then
            Add(ag, agp2[i]);
          fi;
        od;
#
# exclude also isomorphisms between brothers
#
        p_list:=[tmp[Size(tmp)]];
        for i in [1..Size(ag)] do
          pretender:=[Position(aut_h,aut_h[a12]*ag[i]), pairs[c1][c2][a1], c1, c2, 
                                                               pairs[c3][c4][a2], c3, c4];
          b:=false;
          for k in [1..Size(p_list)] do
            ms1:=[];
            ms2:=[];
            for j in [1..Size(matchings[pretender[6]])] do 
              Add(ms1, Set(OnTuplesTuples(matchings[pretender[6]], aut_h[pretender[1]])[j]));
            od;
            for j in [1..Size(matchings[p_list[k][6]])] do 
              Add(ms2, Set(OnTuplesTuples(matchings[p_list[k][6]], aut_h[p_list[k][1]])[j]));
            od;
            if Set(ms1) = Set(ms2) then
              b:=true;
              break;
            fi;
          od;
          if not b then
            Add(p_list, pretender);
            Add(tmp, pretender);
          fi;
        od;
      fi;
    od;
  od;
od;
Append(colorings[4], H4ColoringsSetAutRep(tmp, matchings, aut_h));
#
# generate next complect
#
  if n[SOC] = NOM then p := p - 1;
    else p := SOC;
  fi;
  if p >= 1 then
    for i in [p..SOC] do
      n[SOC - i + p] := n[p] + (SOC - i + p) - p + 1;
    od;
  fi;
until p < 1;
#
# for class 4
#
Add(colorings, []);
for df in [1..NOM] do
  tmp:=[];
  c1:=df; c2:=df; c3:=df; c4:=df;
#
# search for complement pairs of choosen matchings set
#
  for a1 in [1..Size(pairs[c1][c2])] do
    p1:=Concatenation(matchings[c1], OnTuplesTuples(matchings[c2], aut_h[pairs[c1][c2][a1]]));
    for a2 in [a1..Size(pairs[c3][c4])] do
      p2:=Concatenation(matchings[c3], OnTuplesTuples(matchings[c4], aut_h[pairs[c3][c4][a2]]));
      l:=0;
      a12:=0;
      while (a12 < SOHAG and l = 0) do
        a12:=a12 + 1;
        colors:=Concatenation(p1, OnTuplesTuples(p2, aut_h[a12]));
        if IsSimpleGraphVL(colors) then
          l:=a12;
        fi;
      od;
      if l <> 0 then
#
# add finded coloring
#
        Add(tmp, [a12, pairs[c1][c2][a1], c1, c2, pairs[c3][c4][a2], c3, c4]);
#
# generate other designs of choosen matchings set
#
        agp1:=List(Difference(
              Difference(
                AutGroupHypercubeSubgraph(GraphByVL(OnTuplesTuples(p2, aut_h[a12])), 4),
                AutGroupHypercubeSubgraph(
                  GraphByVL(OnTuplesTuples(matchings[c3], aut_h[a12])), 4)),
              Intersection(
                AutGroupHypercubeSubgraph(GraphByVL(matchings[c1]), 4),
                AutGroupHypercubeSubgraph(GraphByVL(
                  OnTuplesTuples(matchings[c2], aut_h[pairs[c1][c2][a1]])), 4))));
        agp2:=[];
#
# exclude also isomorphism m3 <-> Image(m4)
#
        for i in [1..Size(agp1)] do
          ms1:=[];
          ms2:=[];
          for j in [1..Size(matchings[c3])] do 
            Add(ms1, Set(OnTuplesTuples(matchings[c3], aut_h[a12])[j]));
          od;
          for j in [1..Size(matchings[c4])] do 
            Add(ms1, Set(OnTuplesTuples(
              OnTuplesTuples(OnTuplesTuples(matchings[c4], aut_h[pairs[c3][c4][a2]]),aut_h[a12])
                                            , agp1[i])[j])); 
          od;
          if Set(ms1) <>  Set(ms2) then
            Add(agp2, agp1[i]);
          fi;
        od;
        ag:=[];
#
# exclude also isomorphism m1 <-> Image(m2)
#
        for i in [1..Size(agp1)] do
          ms1:=[];
          ms2:=[];
          for j in [1..Size(matchings[c1])] do Add(ms1, Set(matchings[c1][j])); od;
          for j in [1..Size(matchings[c2])] do Add(ms1, Set(OnTuplesTuples(matchings[c2], agp1[i])[j])); od;
          if Set(ms1) <>  Set(ms2) then
            Add(ag, agp2[i]);
          fi;
        od;
#
# exclude also isomorphisms between brothers
#
        p_list:=[tmp[Size(tmp)]];
        for i in [1..Size(ag)] do
          pretender:=[Position(aut_h,aut_h[a12]*ag[i]), pairs[c1][c2][a1], c1, c2, 
                                                               pairs[c3][c4][a2], c3, c4];
          b:=false;
          for k in [1..Size(p_list)] do
            ms1:=[];
            ms2:=[];
            for j in [1..Size(matchings[pretender[6]])] do 
              Add(ms1, Set(OnTuplesTuples(matchings[pretender[6]], aut_h[pretender[1]])[j]));
            od;
            for j in [1..Size(matchings[p_list[k][6]])] do 
              Add(ms2, Set(OnTuplesTuples(matchings[p_list[k][6]], aut_h[p_list[k][1]])[j]));
            od;
            if Set(ms1) = Set(ms2) then
              b:=true;
              break;
            fi;
          od;
          if not b then
            Add(p_list, pretender);
            Add(tmp, pretender);
          fi;
        od;
      fi;
    od;
  od;
  Append(colorings[5], H4ColoringsSetAutRep(tmp, matchings, aut_h));
#
# generate next complect
#
od;
#
# FINISHED
#
return colorings;
end;

######################################################################
AllH4ProperColoringsAutRepLib := function()
######################################################################
#
# Reads all non-isomorpfic perfect colorings 
# reprezentatives in H4
#
# 
######################################################################
local colorings;
colorings:=ReadAsFunction(Filename(
  DirectoriesPackageLibrary("hypercube"),"h4_colorings.g"));
return colorings();
end;




PrintColorings := function(colorings, matchings, aut_h)
#############################################################################
#
# Prints list of colorings
#
#
#############################################################################
local i, j, d, a;
d:=["ABCD","AABC","AABB","AAAB","AAAA"];
for i in [1..5] do
  Print("\n\n***************** Design Type - ", d[i], " *****************");
  for j in [1..Size(colorings[i])] do
    Print("\n\n***************** Coloring - ",colorings[i][j], "*****************");
#    Print("\nMatching # 1 - ", colorings[i][j][3]);
#    Print("\n", matchings[colorings[i][j][3]]);
#    Print("\nMatching # 2 - ", colorings[i][j][4]);
#    Print("\n", OnTuplesTuples(matchings[colorings[i][j][4]], aut_h[colorings[i][j][2]]));
#    Print("\nMatching # 3 - ", colorings[i][j][6]);
#    Print("\n", OnTuplesTuples(matchings[colorings[i][j][6]], aut_h[colorings[i][j][1]]));
#    Print("\nMatching # 4 - ", colorings[i][j][7]);
#    Print("\n", 
#      OnTuplesTuples(
#        OnTuplesTuples( matchings[colorings[i][j][7]], aut_h[colorings[i][j][5]] ), 
#                        aut_h[colorings[i][j][1]] 
#      )
#    );
#    Print("\nAut 1-2 # ", colorings[i][j][1]);
#    Print("\n", aut_h[colorings[i][j][1]]);
#    Print("\nAut 1 # ", colorings[i][j][2]);
#    Print("\n", aut_h[colorings[i][j][2]]);
#    Print("\nAut 2 # ", colorings[i][j][5]);
#    Print("\n", aut_h[colorings[i][j][5]]);
#    Print("\nAut Group of coloring :");
    a:=AutGroupH4Colorings(colorings[i][j],matchings,aut_h);
#    Print("\n", a);
    Print("\nSize of Aut Group : ", Size(a));
  od;
od;
Print("\n\n");
return true;
end;

PrintPairs := function(pairs, matchings, aut_h)
#############################################################################
#
# Prints list of colorings
#
#
#############################################################################
local i, j, d;
for i in [1..8] do
  for j in [i..8] do
    Print("\n\n***************** Pair type - ",i," - ",j, " *****************\n");
    for d in [1..Size(pairs[i][j])] do
      Print("\nPair - [", pairs[i][j][d], ", ",i,", ", j, " ]");
      Print("\nMatching 1 - ");
      Print("\n", matchings[i]);
      Print("\nMatching 2 - ");
      Print("\n", OnTuplesTuples(matchings[j], aut_h[pairs[i][j][d]]));
      Print("\nAut - ", pairs[i][j][d]);
      Print("\n", aut_h[pairs[i][j][d]]);
    od;
  od;
od;
return true;
end;


AllH4ComplementMatchingsPairsAutRep := function()
#############################################################################
#
# Finds all non-isomorphic pairs of H4 matchings listed in h_matchings
# Returns results as upper-triangular matrix of automorphisms
#
#
#############################################################################
local a_matrix, h_matchings, aut_g_list, aut_m_list, m1, m2, m3, p1, i, j, k,
NOM, NOHA, dummy, temp, count;
#
# let's define data and allocate ctructures
#
h_matchings:=AllH4MatchingsAutRepLib();
NOM:=Size(h_matchings);
count:=0;
dummy:=[];
a_matrix:=[];
for i in [1..NOM] do Add(a_matrix, []); od;
for i in [1..NOM] do for j in [1..NOM] do Add(a_matrix[i], []); od; od;
aut_g_list:=List(AutGroupGraph(Hypercube(4)));
NOHA:=Size(aut_g_list);
aut_m_list:=[];
for i in [1..NOM] do
  Add(aut_m_list, 
  List(AutGroupHypercubeSubgraph(GraphByVL(h_matchings[i]), 4))); 
od;
#
# function body begins here
#
for i in [1..NOM] do
  m1:=h_matchings[i];
  for j in [i..NOM] do
    m2:=h_matchings[j];

Print("\ni = ", i, " j = ", j);

#
# for each pair of matchings m1, m2 do
#
    temp:=[];
    for k in [1..NOHA] do
#
# check if aut_g_list[k] is not a m2 automorphism
#
      if not aut_g_list[k] in aut_m_list[j] then
        m3:=OnTuplesTuples(m2, aut_g_list[k]);
        p1:=Concatenation(m1, m3);
#
# check if matchings are complement
#
        if IsSimpleGraphVL(p1) then
#
# check if a_matrix already contains isomorph
#
          if not ExistInSVL(p1, temp, aut_g_list) then
            Add(temp, p1);
            Add(a_matrix[i][j], k); 
            count:=count+1;

Print("\nAdded block # ", count);

          fi;
        fi;
      fi;
    od;
  od;
od;
return a_matrix;
end;


AllH4Colorings := function()
#############################################################################
#
# Finds all non-isomorphic colorings
#
#
#############################################################################
local Matchings, AutH4List, Colorings, ok;

Matchings:=AllH4MatchingsAutRepLib();
AutH4List:=List(AutGroupGraph(Hypercube(4)));
Colorings:=AllH4ProperColoringsAutRepLib();
ok:=PrintColorings(Colorings, Matchings, AutH4List);
return Colorings;
end;


######################################################################
ExistInSVL := function(h_subgraph, hs_list, aut_g_list)
######################################################################
#
# Check if m_list already contains matching, isomorphical to   
# given one under group listed in aut_g_list.
#
#
######################################################################
local p, pm, s, sm, mr, i, j, k;
for k in [1..Size(aut_g_list)] do
  mr:=OnTuplesTuples(h_subgraph, aut_g_list[k]);
  p:=[];
  for i in [1..Size(mr)] do
    Add(p, Set(mr[i]));
  od;
  pm:=Set(p);
  for i in [1..Size(hs_list)] do
    s:=[];
    for j in [1..Size(hs_list[i])] do
      Add(s, Set(hs_list[i][j]));
    od;
    sm:=Set(s);
    if pm = sm then return true; fi;
  od;
od;
return false;
end;


######################################################################
GraphByVL := function(vl)
######################################################################
#
# Generates graph from vertex list
#
#
######################################################################
local nn, i;
nn:=NullGraph(Group(()),16);
for i in [1..Size(vl)] do
  AddEdgeOrbit(nn, [vl[i][1],vl[i][2]]);
  AddEdgeOrbit(nn, [vl[i][2],vl[i][1]]);
od;
nn.group:=AutGroupGraph(nn);
return nn;
end;


######################################################################
IsSimpleGraphVL := function(gammaVL)
######################################################################
#
# Check if gammaVL doesn't contain parallel edges
#
#
######################################################################
local EdgesList, SSEdgesList, i;
EdgesList := [];
SSEdgesList := [];
for i in [1..Size(gammaVL)] do
  Add(EdgesList, SortedList(gammaVL[i]));
od;
SSEdgesList := SSortedList(EdgesList);
if Size(EdgesList) = Size(SSEdgesList) then return true;
  else return false; fi;
end;


######################################################################
H4ColoringsSetAutRep := function(col_list, matchings, aut_h)
######################################################################
#
# Returns set of colorings representstive based on aut_h equivalence
#
#
######################################################################
local isomorph, i, j, col_set;
col_set:=[];
for i in [1..Size(col_list)] do
  isomorph:=false;
  for j in [1..Size(col_set)] do
    if IsomorphismsH4Colorings(col_list[i], col_set[j], 
                               matchings, aut_h) <> [] then
      isomorph:=true;
      break;
    fi;
  od;
  if not isomorph then
    Add(col_set, col_list[i]);
  fi;
od;
return col_set;
end;

######################################################################
IsomorphismsH4Colorings := function(col_1, col_2, matchings, aut_h)
######################################################################
#
# Finds all elements of aut_h presenting isomorphism col1 -> col2
#
#
######################################################################
local col_1_ord, col_2_ord, isomorphisms, i;
isomorphisms:=[];
col_1_ord:=OrderedH4Coloring(
             H4ColoringBySignature(col_1, matchings, aut_h));
col_2_ord:=OrderedH4Coloring(
             H4ColoringBySignature(col_1, matchings, aut_h));
for i in [1..Size(aut_h)] do
  if OrderedH4Coloring(H4ColoringImage(col_2_ord, aut_h[i])) =
     col_1_ord then
    Add(isomorphisms, i);
  fi;
od;
return isomorphisms;
end;


######################################################################
H4ColoringImage := function(col_bs, aut)
######################################################################
#
# Returns an image of col_bs under aut
#
#
######################################################################
local image, i;
image:=[];
for i in [1..4] do
  Add(image, OnTuplesTuples(col_bs[i], aut));
od;
return image;
end;


######################################################################
H4ColoringBySignature := function(col, matchings, aut_h);
######################################################################
#
# Returns colorings as list of matchings based on signature
#
#
######################################################################
return [
matchings[col[3]],
OnTuplesTuples(matchings[col[4]], aut_h[col[2]]),
OnTuplesTuples(matchings[col[6]], aut_h[col[1]]),
OnTuplesTuples(
  OnTuplesTuples(matchings[col[7]], aut_h[col[5]]), aut_h[col[1]])];
end;


######################################################################
OrderedH4Coloring := function(col_bs)
######################################################################
#
# Returns totally ordered representation of given coloring
#
#
######################################################################
local col_ord, mat_ord, i, j;
col_ord:=[];
for i in [1..4] do
  mat_ord:=[];
  for j in [1..8] do
    Add(mat_ord, SortedList(col_bs[i][j]));
  od;
  Add(col_ord, SortedList(mat_ord));
od;
return SortedList(col_ord);
end;






