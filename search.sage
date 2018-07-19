import random
import math

load("cross_flips.sage")

# use simulated annealing to reduce vertices
def search_cross_flip(S, initial_temp=1000, final_temp=1, cooling_rate = 0.003, accept_scale=200):

	# auxillary function to compute energy
	def accept(S,Snew,temp):
		if Snew == 0:
			return False
		prev_min = len(S.vertices())
		vert = len(Snew.vertices())
		if vert < prev_min:
			return True
		prob = math.exp( accept_scale*(prev_min - vert) / temp )
		print "Accept Probability: " + str(prob)
		if random.random() < prob :
			return True
		return False

	dim = S.dimension()+1
	cp = cross_polytope(dim,True)
	cf = cross_flips(dim)
	app = applicable_list(cf,S)
	mini = len(S.vertices())
	minval = S
	temp = initial_temp
	# begin simulated annealing
	while ( temp > final_temp ):
		print "Temperature: " + str(temp)
		# generate neighbors
		options_up   = [len(app[0][i]) for i in range(len(app[0]))]
		options_down = [len(app[1][i]) for i in range(len(app[1]))]
		opts = []
		for j in range(len(options_down)):
			if options_up[j] != 0:
				opts.append((0,j))
			if options_down[j] != 0:
				opts.append((1,j))
		choice = random.choice(opts)
		face = random.choice(range(len(app[choice[0]][choice[1]]))) #0
		Snew = mainVF2_one_move(choice[0],choice[1],face,cp,S,app)
		# check if neighbor is acceptable
		if accept(S, Snew, temp):
			print "Accepted: " + str(len(Snew.vertices()))
			if len(Snew.vertices()) < mini:
				mini = len(Snew.vertices())
				minval = Snew
				print "New minimum: " + str(mini)
			app = update_moves(cf,app[choice[0]][choice[1]][face],app,S,cp)
			S = Snew
		temp *= (1-cooling_rate)
	print "Final min: " + str(mini)
	return minval
