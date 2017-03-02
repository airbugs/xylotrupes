#power analysis GeoSSE
library('ape');
library('diversitree');

dynastes <- read.tree('xylosp.nwk');
geodata <- read.csv('geodata.csv');
geo <- geodata$state;
names(geo) <- geodata$species;

p <- starting.point.geosse(dynastes);
#Full model
lik.base <- make.geosse(dynastes, geo);
ml.base <- find.mle(lik.base, p);
p <- coef(ml.base);
#the previous part is all the same as Regular GeoSSE analysis
#
#then we specifically stimate only assymmatric lambda rates for different state
lik.var.sp <- constrain(lik.base, sAB~0, xA~xB, dA~dB);#only delet sA~sB so that sA and sB can differ
ml.var.sp <- find.mle(lik.var.sp, p[argnames(lik.var.sp)]);
#the let's see the parameter estimation
coef(ml.var.sp, TRUE);#this coef will be the parameter values for GeoSSE simulations
lik.con <- constrain(lik.base, sAB~0, sA~sB, xA~xB, dA~dB);
ml.con <- find.mle(lik.con, p[argnames(lik.con)]);
lik.free.sp <- constrain(lik.base, sAB~0, xA~xB, dA~dB);#delete sA~sB
ml.free.sp <- find.mle(lik.free.sp, p[argnames(lik.free.sp)]);
anova(ml.con, ml.free.sp);#the original dataset shows that freeing the lambda parameter does not help to statistically increase the model fit

#simulating GeoSSE with exact lambda par values and 36 tips.
#to investigate if the tree lambda values are different, how often can we identify this pattern with a tree of 36 tips
#thus, test type 2 error.
par <- coef(ml.var.sp, TRUE);
#a while loop since I do not know how many sims will it take to produce 1000 36-tips trees
count <- 1;
n.sims <- 0;#just to see how many sims does it take to get 1000 trees with exactly 36 tips
p.value.vector <- NULL;#create a null vector to store the p values from each sim
while(count < 11){
	sim.tre <- tree.geosse(par, max.t = 30, x0=0);#the root height of the real tree is 19.2438 (95%HPD = 13-27)
		if(length(sim.tre$tip.label) == 36){
			if(max(cophenetic(sim.tre))>26 & max(cophenetic(sim.tre))<54){#chose thos that have a tree depth similar to the real tree
			sim.states <- sim.tre$tip.state;
			p <- starting.point.geosse(sim.tre);
			lik.base <- make.geosse(sim.tre, sim.states);
			ml.base <- find.mle(lik.base, p);
			p <- coef(ml.base);
			lik.con <- constrain(lik.base, sAB~0, sA~sB, xA~xB, dA~dB);
			ml.con <- find.mle(lik.con, p[argnames(lik.con)]);
			lik.free.sp <- constrain(lik.base, sAB~0, xA~xB, dA~dB);#delete sA~sB
			ml.free.sp <- find.mle(lik.free.sp, p[argnames(lik.free.sp)]);
			temp.anova <- anova(ml.con, ml.free.sp);
			temp.p.value <- temp.anova[2,5];
			p.value.vector <- c(p.value.vector, temp.p.value);#save the p values
			count <- count + 1;
			#if you want to save the sim trees, run the following codes
			#sim.num <- paste(count, '.nex', sep = '');
			#write.nexus(sim.tre, file = sim.num);
			print(count);
			}
		}
	n.sims <- n.sims + 1;
	}
sum(p.value.vector<0.05);#to see how often can we infer sig different lambda