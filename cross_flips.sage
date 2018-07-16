import random
import numpy
import itertools
from copy import deepcopy
from fractions import Fraction
from sage.probability.probability_distribution import GeneralDiscreteDistribution
from sage.combinat.posets.posets import FinitePoset
from sage.graphs.graph_coloring import first_coloring
import networkx as nx
from networkx.algorithms import isomorphism

####### some cross polytopes
cr3=SimplicialComplex([[0,1,2],[0,1,5],[0,4,2],[0,4,5],[3,1,2],[3,1,5],[3,4,2],[3,4,5]])
cr4=SimplicialComplex([[0,1,2,3],[0,1,2,7],[0,1,6,3],[0,1,6,7],[0,5,2,3],[0,5,2,7],[0,5,6,3],[0,5,6,7],[4,1,2,3],[4,1,2,7],[4,1,6,3],[4,5,2,3],[4,5,2,7],[4,5,6,3],[4,5,6,7],[4,1,6,7]])
cr5=SimplicialComplex([[0,1,2,3,4],[0,1,2,3,9],[0,1,2,8,4],[0,1,7,3,4],[0,6,2,3,4],[5,1,2,3,4],[0,1,2,8,9],[0,1,7,3,9],[0,6,2,3,9],[5,1,2,3,9],[0,1,7,8,4],[0,6,2,8,4],[5,1,2,8,4],[0,6,7,3,4],[5,1,7,3,4],[5,6,2,3,4],[0,1,7,8,9],[0,6,2,8,9],[5,1,2,8,9],[0,6,7,3,9],[5,1,7,3,9],[5,6,2,3,9],[0,6,7,8,4],[5,1,7,8,4],[5,6,2,8,4],[5,6,7,3,4],[0,6,7,8,9],[5,1,7,8,9],[5,6,2,8,9],[5,6,7,3,9],[5,6,7,8,4],[5,6,7,8,9]])
######WRITE A FUNCTION TO PRODUCE d-CROSS POLYTOPES ON [0,...,2d-1]

def cross_polytope(dim, zero_index=False):
	cp = []
	for perm in itertools.product('01',repeat=dim):
		index = 1
		if zero_index:
			index = 0
		face = []
		for i in perm:
			if i == '0':
				face.append(index)
			else:
				face.append(index+dim)
			index+=1
		cp.append(face)
	return SimplicialComplex(cp)

########################################################################################
########################################################################################
########################################################################################
##################### Aux functions. Probably some are useless.

###### INPUT: simp_comp = a complex 
###### OUTPUT: the list of facets of the boundary		
def boundary_complex(simp_comp):
	d=simp_comp.dimension()
	fac=list(set(simp_comp.facets()))
	r=list(simp_comp.n_faces(d-1))
	b=[]
	for i in range(len(r)):
		f=[]
		n=0
		for j in range(len(fac)):
			if r[i].is_face(fac[j]):
				n=n+1
		if n==1:
			b.append(r[i])	
	return b

	

def	is_in(Map,v):
	boo=False
	for i in range(len(Map)):
		if v==Map[i]:
			boo=True
	return boo					
		
def Vindex(V,w):
	for i in range(len(V)):
		if (V[i]==w):
			s=i
	return s	

########################### Relabeling vertices to 1,....,n ex: barycentric subdivision
def vertex_relabeling(S):
	ver=list(S.vertices())
	fac=list(S.facets())
	d=S.dimension()
	newfac=[]
	for i in range(len(fac)):
		newfac.append([])
	for j in range(len(fac)):
		for k in range(d+1):
			for l in range(len(ver)):
				if fac[j][k]==ver[l]:
					newfac[j].append(range(len(ver))[l])
	newS=SimplicialComplex(newfac)				
	return newS
	
	
########################### Relabeling vertices to 1,....,n
def random_vertex_relabeling(S):
	fix=set((1,2,3,19,20,35))
	ver=list(set(S.vertices())-set((1,2,3,19,20,35)))
	perm=list(numpy.random.permutation(ver))
	fac=list(S.facets())
	d=S.dimension()
	newfac=[]
	for i in range(len(fac)):
		newfac.append([])
	for j in range(len(fac)):
		for k in range(d+1):
			if fac[j][k] in fix:
				newfac[j].append(fac[j][k])
			else:	
				for l in range(len(ver)):
					if fac[j][k]==ver[l]:
						newfac[j].append(perm[l])
	newS=SimplicialComplex(newfac)				
	return newS	

def technical(IT,V):
	Map=[]
	fa=[]
	for t in range(len(IT)):
		fa.append(IT[t][0])
	S=SimplicialComplex([list(fa)])	
	VS=S.vertices()
	for k in range(len(V)):
		Map.append(-1)
	for i in range(len(IT)):
		for j in range(i+1,len(IT)):
			L=list(set(IT[j][0])-set(IT[i][0]))
			if len(L)==1:
				L1=list(set(IT[j][1])-set(IT[i][1]))
				L2=list(set(IT[i][0])-set(IT[j][0]))
				L3=list(set(IT[i][1])-set(IT[j][1]))
				Map[L1[0]]=L[0]
				Map[L3[0]]=L2[0]
	for s in range(len(Map)):
		if Map[s]==-1:
			candidates=[]
			for n in range(len(IT)):
				if is_in(list(IT[n][1]),s)==True:
					candidates.append(IT[n][0])
			if len(candidates)!=0:		
				intersection=set(candidates[0])
				for m in range(len(candidates)):
					intersection=intersection & set(candidates[m])	
				if len(list(intersection))==1:
					Map[s]=list(intersection)[0]
				else:
					for ll in range(len(list(intersection))):
						if is_in(Map,list(intersection)[len(list(intersection))-ll-1])==False:
							Map[s]=list(intersection)[len(list(intersection))-ll-1]
	for l in range(len(Map)):
		if Map[l]==-1:
			for ll in range(len(Map)):
				if is_in(Map,len(Map)-ll-1)==False:
					Map[l]=len(Map)-ll-1			
	return Map		
	
def powerset_n(s,n):
	return itertools.combinations(s,n)

###CHECK AND SUBS IN MATCH CROSS
def powerset_NEW(s):
	del s[-1]
	x=len(s)
	p=[]
	for i in range(1,1 << x):
		p.append([s[j] for j in range (x) if (i & (1 << j))])
	del p[len(p)-1]	
	return p

def powerset(s):
	x=len(s)
	p=[]
	for i in range(1,1 << x):
		p.append([s[j] for j in range (x) if (i & (1 << j))])
	del p[len(p)-1]	
	return p

def mirror_sum(h1,h2):
	g=[]
	l=len(h1)
	g.append(1)
	for i in range(1,len(h1)-1):
		g.append(h1[i]+h2[l-i-1])
	g.append(1)
	return g

def fvec_list(L):
	f=[]
	for i in range(len(L[0])):
		f.append(L[0][i].f_vector())
	return f, L[1]

###### the utility GM.subgraph_isomorphisms_iter() contains redundancies, due to the possible simmetry of move.flip_graph().
def remove_redundancies_list(t):
	s=[]
	for i in t:
		if i not in s:
			s.append(i)
	return s

#########################################################################################
#########################################################################################
#########################################################################################



#########################################################################################
#########################################################################################
#########################################################################################
#################### Functions to produce cross_flips

def cross(s):
	d=s.dimension()
	c=s
	falis=[]
	for w in range(d):
		falis.append(list(s.n_faces(d-w)))
	for i in range(d):
		fa=falis[i]
		for j in range(len(fa)):
			count=0
			for k in range(i+1):
				if (k in fa[j]):
					count=count+1
			if count == 0:		
				c.stellar_subdivision(fa[j],inplace=True)
	return c	
	
	
def all_cross(n):
	s=Simplex(n)
	p=powerset(s.faces())
	allc=[]
	fvec=[]
	for i in range(len(p)):
		allc.append(cross(SimplicialComplex(p[i])))
		fvec.append((cross(SimplicialComplex(p[i]))).h_vector())
	return allc, fvec



def no_redundant(L):
	A=L
	NR1=[]
	NR2=[]
	NR1.append(A[0][0])
	NR2.append(A[1][0])
	for i in range(1,len(A[0])):
		count=0
		for j in range(i):
			if A[1][j]==A[1][i]:
				count=count+1
		if count==0:
			NR1.append(A[0][i])
			NR2.append(A[1][i])
	return NR1,NR2
	
def match_cross(L):
	l=len(L[1][0])
	simp=[]
	m1=[]
	m2=[]
	n1=[]
	n2=[]
	for j in range(l):
		simp.append(binomial(l-1,j))	
	A=L
	p=[] 
	i=0
	while (len(A[1])>0): 
		m1.append(A[1][0])
		n1.append(A[0][0])
		for k in range(len(A[1])): 	
			if mirror_sum(A[1][0],A[1][k])==simp:
				m2.append(A[1][k])
				n2.append(A[0][k])
				p=k
		if p==0:	
			A[0].remove(A[0][0])
			A[1].remove(A[1][0])
		else:
			A[0].remove(A[0][0])
			A[1].remove(A[1][0])
			A[0].remove(A[0][p-1])
			A[1].remove(A[1][p-1])
	return n1,n2,m1,m2

#### Order cross-flips (separate those adding vertices from those eliminating vertices).

def ord_cross(c, trivial_move = false):
	trivial_index = -1;
	for i in xrange(len(c[0])):
		if c[0][i].f_vector()[1]>c[1][i].f_vector()[1]:
			c[0][i],c[1][i]=c[1][i],c[0][i]
			c[2][i],c[3][i]=c[3][i],c[2][i]
		else:
			if c[2][i] == c[3][i] and not trivial_move :
				trivial_index = i
	if trivial_index >= 0:
		for j in range(4):
			c[j].remove(c[j][trivial_index])
	return c	



###### A non redundant matched double list of flips. If trivial_move=True also includes the trivial flip	
###### INPUT: d a positive integer. If trivial_move = true the trivial flip is included.
###### OUTPUT: a double list c. c[0] contains the moves increasing f0 and c[1] those decreasing f0. For brevity lets call the move c[i][j] a move of type i,j.  
def	cross_flips(d, trivial_move = False):
	al=all_cross(d)
	n=no_redundant(al)
	mc=match_cross(n)
	o=ord_cross(mc,trivial_move)
	return o	
	
################################################################################
################################################################################
################################################################################	


################################################################################
################################################################################
################################################################################
################################################################################
################## Functions to test appliability of a move
####### INPUT: move = a cross flip, simp_comp = the source simplicial complex. If stop_first = true it stops when one applicable flip is found. 
####### OUTPUT: a list of all applicable flips of type move on simp_comp.
def is_applicable_graphs_VF2(move,simp_comp, stop_first=False, check_induced=True):
	G=simp_comp.flip_graph()
	nG = G.networkx_graph()
	g=move.flip_graph()
	ng = g.networkx_graph()
	GM = isomorphism.GraphMatcher(nG,ng)
	sub=[]	
	if stop_first==True:	
		for i in GM.subgraph_isomorphisms_iter():
			if is_induced(SimplicialComplex(set(i)),simp_comp)==True or check_induced ==False:
				sub.append(set(i))
				return sub
	else:
		for i in GM.subgraph_isomorphisms_iter():
			if is_induced(SimplicialComplex(set(i)),simp_comp)==True or check_induced ==False:
				sub.append(set(i))				
	return sub	

####### Same as above, but takes a graph as first argument
def is_applicable_graphs_VF2NOSTOP_for_update(H,move,c1):
	nH = H.networkx_graph()
	g=move.flip_graph()
	ng = g.networkx_graph()
	GM = isomorphism.GraphMatcher(nH,ng)
	sub=[]		
	for i in GM.subgraph_isomorphisms_iter():
		if is_induced(SimplicialComplex(set(i)),c1)==True:
			sub.append(set(i))
	return sub	



####### INPUT: cflips = list of cross flips, simp_comp = source simplicial complex, if stop_first = true applies the first found move.	
####### OUTPUT: A list app of two lists containing lists. Each app[i][j] contains all the applicable flips on simp_comp of type i,j. (See cross_flips).
def applicable_list(cflips,simp_comp,stop_first=False):
	app=[[],[]]
	for i in range(len(cflips[0])):
		isa0=list(is_applicable_graphs_VF2(cflips[0][i],simp_comp,stop_first))
		app[0].append(remove_redundancies_list(isa0))
	for j in range(len(cflips[1])):
		isa1=list(is_applicable_graphs_VF2(cflips[1][j],simp_comp,stop_first))
		app[1].append(remove_redundancies_list(isa1))	
	return app



###### INPUT: sub = a subcomplex of simp_comp, simp_comp = the source simplicial complex.
###### Output: true if sub is an induced subcomplex of simp_comp
def is_induced(sub,simp_comp):
	V1=sub.vertices()
	V2=simp_comp.vertices()
	is_ind=True
	for k in range(2,simp_comp.dimension()+2):
		p=powerset_n(V1,k)
		##p is an iterator
		for i in p:
			if (Simplex(i) in simp_comp.n_faces(k-1))==True:
				if (Simplex(i) not in sub.n_faces(k-1))==True:
					return 	False
	return is_ind

#################################################################################
#################################################################################
#################################################################################

#################################################################################
#################################################################################
#################################################################################
################### Functions to flip
###### INPUT: 2 complexes.
###### OUTPUT: mo with all the interior vertices relabeled, while those on the boundary are to be identified with some from simp_comp.
def vertex_rename(mo,simp_comp):
	ver=simp_comp.vertices()
	sub_comp=SimplicialComplex(mo)
	dd=sub_comp.dimension()
	fac=list(sub_comp.n_faces(dd))
	copy=fac
	vtc=list(set(sub_comp.vertices())-set((SimplicialComplex(boundary_complex(sub_comp))).vertices()))
	for i in range(len(vtc)): 
		for k in range(len(sub_comp.n_faces(dd))):
			fac[k]=[max(ver)+i+1 if x==vtc[i] else x for x in (copy)[k]]
		copy=fac
	return fac	

###### INPUT: two lists of facets.
######OUTPUT: The complex obtained by deleting the facets in mo1 and adding mo2	
def move(mo1,mo2,c):
	fac=set(c.facets())
	fmo1=set((SimplicialComplex(mo1)).facets())
	fmo2=set((SimplicialComplex(mo2)).facets())
	dele=fac-fmo1
	un=dele|fmo2
	return SimplicialComplex(list(un))


####TO SUBSTITUTE WITH SOMETHING EASIER	
###### INPUT: mo = an element of cross_flips that is applicable on c, c = the source simplicial complex.
###### OUTPUT: the complement of mo in the cross_polytope, 				
def right_complement(mo,c):
	V=c.vertices()
	D=SimplicialComplex(mo)
	G=D.flip_graph()
	VD=D.vertices()
	F=c.facets()
	isa=is_applicable_graphs_VF2(D, c, stop_first = True, check_induced = False)[0]
	G1=(SimplicialComplex(isa).flip_graph())
	M=list(G.is_isomorphic(G1, certificate = True))
	IT=M[1].items()
	Map=technical(IT,V)
	NF=[]
	for t in range(len(F)):
		NF.append([])
		for y in range(len(list(F[t]))):
			NF[t].append(Map[list(F[t])[y]])
	C=list(set(SimplicialComplex(NF).facets())-set(SimplicialComplex(mo).facets()))									
	return C					
		
####################################################################################
####################################################################################
####################################################################################

#####################################################################################
#####################################################################################
#####################################################################################
############ Main 
###### INPUT: i = integer 0,...,len(cross_flips[0]), s = integer 0,1, cflips = output of cross_flips, cross_poly = d-cross polytope, simp_comp = source simplicial complex.
###### OUTPUT: a simplicial complex obtained applying a flip of type s,i in a non specified "part" of simp_comp.
def mainVF2(i,s,cflips,cross_poly,simp_comp):
	isa=list(is_applicable_graphs_VF2(cflips[s][i],simp_comp))
	if len(isa)==0:
		return 0
	else:
	    prov3=right_complement(isa[0],cross_poly)
	    vr=vertex_rename(prov3,simp_comp)
	    mo=move(isa[0],vr,simp_comp)
	return mo

###### INPUT: s = integer 0,1, i = integer 0,...,len(cross_flips[0]), j = integer 0,..., len(app_list[s][i]), cross_poly = d-cross polytope, simp_comp = the source simplicial complex, app_list = the output of applicable_list.
###### OUTPUT: the simplicial complex obtained applying the flip app_list[s][i][j] to simp_comp.
def mainVF2_one_move(s,i,j,cross_poly,simp_comp,app_list):
	complement=right_complement(app_list[s][i][j],cross_poly)
	vr=vertex_rename(complement,simp_comp)
	new=move(app_list[s][i][j],vr,simp_comp)
	return new
	
###### INPUT: cflips = the output of cross_flips, sub = app[i][j][k] for some i,j,k, simp_comp = the source simp. complex, cross_poly = the d-cross polytope.
###### OUTPUT: the updated set of applicable moves on the new s.c. obtained applying sub to simp_comp.
###### Supposed to be faster than using applicable_list, especially for large complexes.
def update_moves(cflips,sub,app,simp_comp,cross_poly,only_down=False):
	complement=right_complement(sub,cross_poly)
	vr=vertex_rename(complement,simp_comp)
	vrs=[Simplex(i) for i in vr]
	new_complex=move(sub,vr,simp_comp)
	G=new_complex.flip_graph()
	new_app=[[],[]]
	facets_to_check=list(set([i for i in G.vertices() for j in vrs if G.distance(j,i)<=4]))
	H=G.subgraph(facets_to_check)
	if only_down==False:
		for i0 in range(len(app[0])):
			a=[]
			for j0 in range(len(app[0][i0])):
				if app[0][i0][j0]&sub==set() and is_induced(SimplicialComplex(app[0][i0][j0]),new_complex)==True:
					a.append(app[0][i0][j0])
			new_isa=is_applicable_graphs_VF2NOSTOP_for_update(H,cflips[0][i0],new_complex)
			for i in new_isa:
				if is_induced(SimplicialComplex(i),new_complex)==True:
					a.append(i)		
			new_app[0].append(remove_redundancies_list(a))	
	for i1 in range(len(app[1])):
		b=[]
		for j1 in range(len(app[1][i1])):
			if app[1][i1][j1]&sub==set() and is_induced(SimplicialComplex(app[1][i1][j1]),new_complex)==True:
				b.append(app[1][i1][j1])
		new_isa=is_applicable_graphs_VF2NOSTOP_for_update(H,cflips[1][i1],new_complex)
		for i in new_isa:
			if is_induced(SimplicialComplex(i),new_complex)==True:
				b.append(i)		
		new_app[1].append(remove_redundancies_list(b))	
	return new_app	






	
	


